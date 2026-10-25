//! A proof that a spiral never stops growing.
//!
//! The analyser notices a spiral when a trajectory keeps coming back to one sign
//! pattern, larger each time. That alone proves nothing: a trajectory can grow
//! for a while and still settle into a cycle. This module turns the observation
//! into a proof.
//!
//! Take the last three returns. If they lie on a line `p + x·q` (they usually
//! do: a spiral re-enters its pattern through the same face, at a growing
//! distance), treat the position `x` as an unknown and run one lap of the event
//! map *symbolically*. Inside a pattern the state moves in a straight line, and
//! the time to the next sign change is `⌈|v_e| / |d_e|⌉`, which is an affine
//! function of `x` once `x` is restricted to a residue class modulo a suitable
//! `m`. So the whole lap is affine on each residue class, valid for every `x`
//! above some threshold, and ends back on the line at `x' = (α/m)·x + const`.
//!
//! A lap may also land on a *parallel* line (the offset can depend on a residue
//! class, e.g. `v₀ + v₂` alternating between 3 and 7), so the family is a finite
//! union of parallel lines, each split into residue classes.
//!
//! The certificate is then the induction: if every (line, residue class) state
//! reachable from the current one laps back onto a line of the family, lands in
//! a reachable state, and moves strictly outward (`x' > x`) for every `x` above
//! a common threshold `X`, and the current return is above `X`, then the
//! trajectory laps forever with strictly growing `x`. It never repeats a state, so it is unbounded and
//! is neither fixed, nor a cycle, nor a ray. `α > m` is a multiplicative spiral
//! (`x ↦ 2x + 13`); `α = m` with a positive offset is an additive one, where the
//! lap itself gets longer each time. (A slope-1 lap whose length does not grow
//! is a pure translation, which is a helix; it is left to the helix test.)
//!
//! Some spirals also carry a *steady* coordinate: one that keeps its sign for
//! the whole lap, so it never decides when the next event happens, but whose
//! value drifts in a way that puts the returns on a curve instead of a line.
//! Such a coordinate `z` is left out of the line. The lap needs it to stay
//! clear of zero, which holds when `|z| ≥ B(x)` for an affine bound `B` read off
//! the lap, and the induction carries that bound forward: the lap must add at
//! least `B(x') − B(x)` to `|z|`. A steady coordinate that is never pushed
//! towards zero needs no bound at all.

// Index loops mirror the formulas.
#![allow(clippy::needless_range_loop)]

use crate::exact::Reduced;

/// One proven lap, for display. Returns lie on `base + x₁·dir[0]` (a line,
/// `rank` 1) or `base + x₁·dir[0] + x₂·dir[1]` (a plane, `rank` 2), ignoring the
/// coordinates in the `steady` bitmask, and one lap maps `x ↦ (slope·x + off)/den`.
/// The induction covers `classes` states (a parallel line or plane and a residue
/// class of `x` modulo `modulus`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Lap<const N: usize> {
    pub rank: u8,
    pub base: [i64; N],
    pub dir: [[i64; N]; 2],
    pub steady: u32,
    pub slope: [[i64; 2]; 2],
    pub off: [i64; 2],
    pub den: i64,
    pub modulus: [u32; 2],
    pub classes: u32,
    /// Events in one lap from the current return.
    pub events: u32,
}

/// `a·y + b`, with `y` the unknown (`x = c + m·y`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
struct Aff {
    a: i128,
    b: i128,
}

impl Aff {
    fn at(&self, y: i128) -> i128 {
        self.a * y + self.b
    }
}

enum LapErr {
    /// The ceiling divisions need a finer modulus: multiply it by this.
    Refine(i128),
    Fail(&'static str),
}

struct LapOk<const N: usize> {
    /// Smallest `y` for which the symbolic lap is exact (`None`: every `y`).
    y_min: Option<i128>,
    /// The lap ends on the line through `line` at `x' = alpha·y + beta`.
    alpha: i128,
    beta: i128,
    line: [i128; N],
    events: u32,
    /// Some event time depends on `x`, so the lap gets longer as it moves out.
    grows_in_time: bool,
    /// Per steady coordinate: the bound `B` that `|z|` must meet at the start of
    /// the lap, the change `Δ` of `|z|` over the lap, and whether it is ever
    /// pushed towards zero (if not, `B` is irrelevant).
    bound: [Aff; 8],
    gain: [Aff; 8],
    pushed_in: [bool; 8],
}

/// Limits that keep failed attempts cheap (everything lives on the stack).
const MAX_MODULUS: i128 = 4096;
const MAX_LAP_EVENTS: u32 = 512;

fn gcd(mut a: i128, mut b: i128) -> i128 {
    a = a.abs();
    b = b.abs();
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn div_floor(a: i128, b: i128) -> i128 {
    let q = a / b;
    if (a % b != 0) && ((a < 0) != (b < 0)) {
        q - 1
    } else {
        q
    }
}

fn div_ceil(a: i128, b: i128) -> i128 {
    -div_floor(-a, b)
}

fn raise(y_min: &mut Option<i128>, y: i128) {
    *y_min = Some(y_min.map_or(y, |m| m.max(y)));
}

/// `hi ≥ lo` for every `y ≥ y_min`: raise `y_min` where needed, or fail if `hi`
/// falls behind for large `y`.
fn dominate(y_min: &mut Option<i128>, hi: Aff, lo: Aff) -> Result<(), LapErr> {
    if hi.a > lo.a {
        raise(y_min, div_ceil(lo.b - hi.b, hi.a - lo.a));
        Ok(())
    } else if hi.a == lo.a && hi.b >= lo.b {
        Ok(())
    } else {
        Err(LapErr::Fail("a bound falls behind as the spiral grows"))
    }
}

/// Run one lap symbolically from `p + (c + m·y)·q` until the pattern returns.
#[allow(clippy::too_many_arguments)]
fn lap<const N: usize>(
    w: &Reduced<N>,
    pattern: &[i8; N],
    steady: &[bool; N],
    p: &[i128; N],
    q: &[i128; N],
    m: i128,
    c: i128,
) -> Result<LapOk<N>, LapErr> {
    let mut v: [Aff; N] = std::array::from_fn(|e| Aff {
        a: m * q[e],
        b: p[e] + c * q[e],
    });
    // Steady coordinates: |z| + gain must stay ≥ 1, i.e. |z| ≥ 1 − gain at every event.
    let mut gain = [Aff::default(); 8];
    let mut bound = [Aff { a: 0, b: 1 }; 8];
    let mut pushed_in = [false; 8];
    let mut y_min = None;
    let mut grows_in_time = false;
    for events in 0..=MAX_LAP_EVENTS {
        // Signs, valid for every y above the threshold. Steady coordinates keep theirs.
        let mut s = *pattern;
        for e in 0..N {
            if steady[e] {
                continue;
            }
            let Aff { a, b } = v[e];
            if a == 0 {
                s[e] = b.signum() as i8;
            } else {
                s[e] = a.signum() as i8;
                // sign(a)·(a·y + b) > 0  ⇔  y > −sign(a)·b / |a|
                raise(&mut y_min, div_floor(-a.signum() * b, a.abs()) + 1);
            }
        }
        if events == 0 {
            if s != *pattern {
                return Err(LapErr::Fail("the family does not start in the pattern"));
            }
        } else if s == *pattern {
            // Back in the pattern: the state must be on a line parallel to q.
            let e0 = (0..N)
                .find(|&e| q[e] != 0)
                .ok_or(LapErr::Fail("no direction"))?;
            if v[e0].a % q[e0] != 0 {
                return Err(LapErr::Fail("the lap ends off the line"));
            }
            let alpha = v[e0].a / q[e0];
            if (0..N).any(|e| !steady[e] && v[e].a != alpha * q[e]) {
                return Err(LapErr::Fail("the lap ends off the line"));
            }
            let b: [i128; N] = std::array::from_fn(|e| if steady[e] { 0 } else { v[e].b });
            let (line, beta) = on_line(&b, q, e0);
            return Ok(LapOk {
                y_min,
                alpha,
                beta,
                line,
                events,
                grows_in_time,
                bound,
                gain,
                pushed_in,
            });
        }
        let d = w.drift(&s);
        if d.iter().all(|&x| x == 0) {
            return Err(LapErr::Fail("the lap reaches a fixed pattern"));
        }
        if (0..N).all(|e| d[e] == 0 || (s[e] != 0 && d[e].signum() as i8 == s[e])) {
            return Err(LapErr::Fail("the lap reaches an absorbing pattern"));
        }
        // Time to the next sign change: the smallest of the per-coordinate times.
        let mut times = [Aff::default(); N];
        let mut n = 0;
        for e in 0..N {
            let de = d[e] as i128;
            if steady[e] {
                continue;
            }
            if s[e] == 0 {
                if de != 0 {
                    times[n] = Aff { a: 0, b: 1 };
                    n += 1;
                }
            } else if de.signum() as i8 == -s[e] {
                let (a, b) = (s[e] as i128 * v[e].a, s[e] as i128 * v[e].b);
                let dd = de.abs();
                if a % dd != 0 {
                    return Err(LapErr::Refine(dd / gcd(a, dd)));
                }
                times[n] = Aff {
                    a: a / dd,
                    b: div_ceil(b, dd),
                };
                n += 1;
            }
        }
        let times = &times[..n];
        let k = *times
            .iter()
            .min_by_key(|t| (t.a, t.b))
            .ok_or(LapErr::Fail("only steady coordinates move"))?;
        for t in times {
            if t.a > k.a {
                dominate(&mut y_min, *t, k)?;
            }
        }
        if k.a < 0 {
            return Err(LapErr::Fail("an event time shrinks as the spiral grows"));
        } else if k.a == 0 {
            if k.b < 1 {
                return Err(LapErr::Fail("an event time is not positive"));
            }
        } else {
            raise(&mut y_min, div_ceil(1 - k.b, k.a));
            grows_in_time = true;
        }
        for e in 0..N {
            let de = d[e] as i128;
            if steady[e] {
                // |z| changes by s·d per step; it must stay ≥ 1 at the segment's end.
                let sd = s[e] as i128 * de;
                gain[e].a += k.a * sd;
                gain[e].b += k.b * sd;
                if sd < 0 {
                    pushed_in[e] = true;
                    let need = Aff {
                        a: -gain[e].a,
                        b: 1 - gain[e].b,
                    };
                    // Keep the larger of the old bound and this one, for large y.
                    if (need.a, need.b) > (bound[e].a, bound[e].b) {
                        if need.a > bound[e].a {
                            dominate(&mut y_min, need, bound[e])?;
                        }
                        bound[e] = need;
                    } else {
                        dominate(&mut y_min, bound[e], need)?;
                    }
                }
                continue;
            }
            v[e].a += k.a * de;
            v[e].b += k.b * de;
            if v[e].a.abs() > 1 << 100 || v[e].b.abs() > 1 << 100 {
                return Err(LapErr::Fail("numbers too large"));
            }
        }
    }
    Err(LapErr::Fail("the lap does not return to the pattern"))
}

/// Write `b = line + x·q` with `line[e0]` in `[0, |q[e0]|)`: a canonical point
/// per parallel line, and the position `x` along it.
fn on_line<const N: usize>(b: &[i128; N], q: &[i128; N], e0: usize) -> ([i128; N], i128) {
    let x = div_floor(b[e0], q[e0]);
    (std::array::from_fn(|e| b[e] - x * q[e]), x)
}

/// Most recent returns used to guess the direction of growth.
pub const RETURNS: usize = 5;

/// Try to prove that a trajectory grows forever, from its latest returns to
/// `pattern` (oldest first; the last one is the current state). `steady` marks
/// coordinates that kept their sign at every event over those returns.
/// `plane` also tries two-parameter families, which is slower.
pub fn prove<const N: usize>(
    w: &Reduced<N>,
    pattern: &[i8; N],
    steady: &[bool; N],
    returns: &[[i64; N]],
    plane: bool,
) -> Option<Lap<N>> {
    explain(w, pattern, steady, returns, plane).ok()
}

/// [`prove`], with the reason when there is no proof.
pub fn explain<const N: usize>(
    w: &Reduced<N>,
    pattern: &[i8; N],
    steady: &[bool; N],
    returns: &[[i64; N]],
    plane: bool,
) -> Result<Lap<N>, &'static str> {
    assert!(N <= 8);
    let r0 = returns.last().ok_or("no returns")?;
    let none = [false; N];
    let mut why = "the returns do not show a direction";
    // Candidate directions: differences between the current and earlier returns.
    let mut tried: [[i128; N]; RETURNS * 2] = [[0; N]; RETURNS * 2];
    let mut n_tried = 0;
    for free in [&none, steady] {
        if n_tried > 0 && *free == none {
            break;
        }
        for r in returns.iter().rev().skip(1) {
            let diff: [i128; N] =
                std::array::from_fn(|e| if free[e] { 0 } else { (r0[e] - r[e]) as i128 });
            let g = diff.iter().fold(0, |g, &x| gcd(g, x));
            if g == 0 {
                continue;
            }
            let q: [i128; N] = std::array::from_fn(|e| diff[e] / g);
            if tried[..n_tried].contains(&q) || n_tried == tried.len() {
                continue;
            }
            tried[n_tried] = q;
            n_tried += 1;
            match prove_with(w, pattern, free, &q, r0) {
                Ok(l) => return Ok(l),
                Err(e) => why = e,
            }
        }
    }
    // Returns that spread over a plane: two growth directions at once.
    if !plane {
        return Err(why);
    }
    match plane::prove(w, pattern, returns) {
        Ok(l) => Ok(l),
        Err(e) if e == plane::NOT_A_PLANE => Err(why),
        Err(e) => Err(e),
    }
}

/// Most (line, residue class) states the induction may visit.
const MAX_STATES: usize = 64;

fn prove_with<const N: usize>(
    w: &Reduced<N>,
    pattern: &[i8; N],
    steady: &[bool; N],
    q: &[i128; N],
    r0: &[i64; N],
) -> Result<Lap<N>, &'static str> {
    let e0 = (0..N).find(|&e| q[e] != 0).ok_or("no direction")?;
    let start: [i128; N] = std::array::from_fn(|e| if steady[e] { 0 } else { r0[e] as i128 });
    let (p0, x0) = on_line(&start, q, e0);

    let mut m: i128 = 1;
    'refine: loop {
        if m > MAX_MODULUS {
            return Err("needs too many residue classes");
        }
        // States (line, residue class), visited breadth first.
        let mut states = [([0i128; N], 0i128); MAX_STATES];
        let mut n_states = 1;
        states[0] = (p0, x0.rem_euclid(m));
        let mut laps: [Option<LapOk<N>>; MAX_STATES] = [const { None }; MAX_STATES];
        // Threshold X: every state laps correctly and outward for x ≥ X.
        let mut threshold: Option<i128> = None;
        let mut i = 0;
        while i < n_states {
            let (p, c) = states[i];
            let ok = match lap(w, pattern, steady, &p, q, m, c) {
                Ok(ok) => ok,
                Err(LapErr::Refine(f)) => {
                    m *= f;
                    continue 'refine;
                }
                Err(LapErr::Fail(why)) => return Err(why),
            };
            // Outward: x' − x = (alpha − m)·y + (beta − c) > 0 for all y ≥ y_req.
            let y_req = if ok.alpha < m {
                return Err("the lap shrinks");
            } else if ok.alpha == m {
                if ok.beta <= c {
                    return Err("the lap does not move outward");
                }
                if !ok.grows_in_time {
                    return Err("the lap is a translation (a helix)");
                }
                ok.y_min
            } else {
                let y_grow = div_floor(c - ok.beta, ok.alpha - m) + 1;
                Some(ok.y_min.map_or(y_grow, |y| y.max(y_grow)))
            };
            if let Some(y) = y_req {
                let x = c + m * y;
                threshold = Some(threshold.map_or(x, |t| t.max(x)));
            }
            // States reached: the end line, with every class of x' = alpha·y + beta.
            let step = gcd(ok.alpha, m);
            let mut r = ok.beta.rem_euclid(step);
            while r < m {
                let next = (ok.line, r);
                if !states[..n_states].contains(&next) {
                    if n_states == MAX_STATES {
                        return Err("needs too many residue classes");
                    }
                    states[n_states] = next;
                    n_states += 1;
                }
                r += step;
            }
            laps[i] = Some(ok);
            i += 1;
        }
        if threshold.is_some_and(|t| x0 < t) {
            return Err("the current return is below the threshold");
        }
        // Steady coordinates pushed towards zero: carry |z| ≥ B(x) through the
        // induction with one bound B for all states (only for m = 1, so x = y).
        for e in (0..N).filter(|&e| steady[e]) {
            let laps = laps[..n_states].iter().flatten();
            if !laps.clone().any(|l| l.pushed_in[e]) {
                continue;
            }
            if m != 1 {
                return Err("a steady coordinate with several residue classes");
            }
            // B: the bound that dominates every state's bound for large x.
            let b = laps
                .clone()
                .map(|l| l.bound[e])
                .max_by_key(|b| (b.a, b.b))
                .unwrap_or(Aff { a: 0, b: 1 });
            let mut x_from = x0;
            for l in laps.clone() {
                if l.bound[e].a < b.a {
                    x_from = x_from.max(div_ceil(l.bound[e].b - b.b, b.a - l.bound[e].a));
                }
            }
            if x_from > x0 {
                return Err("a steady coordinate is too close to zero");
            }
            let z0 = (pattern[e] as i64 * r0[e]) as i128;
            if z0 < b.at(x0) {
                return Err("a steady coordinate is too close to zero");
            }
            // |z'| = |z| + gain(x) ≥ B(x) + gain(x) must be ≥ B(alpha·x + beta) for x ≥ x0.
            for l in laps {
                let slack = Aff {
                    a: b.a + l.gain[e].a - b.a * l.alpha,
                    b: l.gain[e].b - b.a * l.beta,
                };
                if slack.a < 0 || slack.at(x0) < 0 {
                    return Err("a steady coordinate does not keep up");
                }
            }
        }
        let ok = laps[0].as_ref().ok_or("no lap from the current state")?;
        let c0 = x0.rem_euclid(m);
        // x' = alpha·(x − c)/m + beta, written over a common denominator m.
        let (num, off, den) = (ok.alpha, ok.beta * m - ok.alpha * c0, m);
        let gg = gcd(gcd(num, off), den);
        return Ok(Lap {
            rank: 1,
            base: std::array::from_fn(|e| p0[e] as i64),
            dir: [std::array::from_fn(|e| q[e] as i64), [0; N]],
            steady: (0..N).filter(|&e| steady[e]).map(|e| 1 << e).sum(),
            slope: [[(num / gg) as i64, 0], [0, 0]],
            off: [(off / gg) as i64, 0],
            den: (den / gg) as i64,
            modulus: [m as u32, 1],
            classes: n_states as u32,
            events: ok.events,
        });
    }
}

/// Spirals whose returns spread over a plane: two unknowns instead of one.
///
/// Same idea as the line case, with `x ∈ ℤ²`. The symbolic lap follows the
/// branch the current return actually takes, and records every condition that
/// branch needs (signs, which coordinate hits zero first, positive event
/// times) as an affine inequality in `x`. A proof is a cone
/// `R = x₀ + {s·E₁ + t·E₂ : s, t ≥ 0}` with its corner at the current return
/// on which every condition holds, which every lap maps into itself, on which
/// some linear functional grows by at least one per lap, and on which the lap
/// time grows along both edges (so the orbit is not a helix). An affine
/// function is nonnegative on `R` iff it is at the corner and its linear part
/// is along both edges. Candidate edges come from the conditions' boundary
/// lines, the axes, and the lap map's rational eigenvectors.
mod plane {
    use super::{div_ceil, div_floor, gcd, Lap, MAX_LAP_EVENTS, MAX_MODULUS, RETURNS};
    use crate::exact::Reduced;

    pub const NOT_A_PLANE: &str = "the returns do not span a plane";

    /// `a·y + b` with `y ∈ ℤ²`.
    #[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
    struct Aff {
        a: [i128; 2],
        b: i128,
    }

    impl Aff {
        fn at(&self, y: [i128; 2]) -> i128 {
            self.a[0] * y[0] + self.a[1] * y[1] + self.b
        }
    }

    /// Distinct directions of the conditions' linear parts (in `x`).
    const MAX_NORMALS: usize = 48;

    struct Normals {
        dirs: [[i128; 2]; MAX_NORMALS],
        n: usize,
    }

    impl Normals {
        fn add(&mut self, a: [i128; 2]) -> Result<(), Err> {
            let g = gcd(a[0], a[1]);
            if g == 0 {
                return Ok(());
            }
            let d = [a[0] / g, a[1] / g];
            if self.dirs[..self.n].contains(&d) {
                return Ok(());
            }
            if self.n == MAX_NORMALS {
                return Err(Err::Fail("too many different conditions"));
            }
            self.dirs[self.n] = d;
            self.n += 1;
            Ok(())
        }
    }

    /// A symbolic lap from one state, in `x` (forms scaled by `m₀·m₁`): its end
    /// `x' = (map·x + off)/mm` on the plane through `line`, and its length.
    #[derive(Clone, Copy)]
    struct LapOk<const N: usize> {
        map: [[i128; 2]; 2],
        off: [i128; 2],
        time: Aff,
        line: [i128; N],
        /// Residues of `x'`: `beta + map_y·ℤ²` modulo `m`.
        beta: [i128; 2],
        map_y: [[i128; 2]; 2],
        events: u32,
    }

    enum Err {
        Refine(usize, i128),
        Fail(&'static str),
    }

    /// The plane's basis `q` restricted to two coordinates with a nonzero minor.
    struct Frame<const N: usize> {
        q: [[i128; N]; 2],
        rows: (usize, usize),
        det: i128,
    }

    impl<const N: usize> Frame<N> {
        fn new(q: [[i128; N]; 2]) -> Option<Self> {
            for i in 0..N {
                for j in i + 1..N {
                    let det = q[0][i] * q[1][j] - q[1][i] * q[0][j];
                    if det != 0 {
                        return Some(Frame {
                            q,
                            rows: (i, j),
                            det,
                        });
                    }
                }
            }
            None
        }

        /// Exact coordinates of `v` in the basis, if `v` lies in its integer span.
        fn coords(&self, v: &[i128; N]) -> Option<[i128; 2]> {
            let (i, j) = self.rows;
            let n0 = self.q[1][j] * v[i] - self.q[1][i] * v[j];
            let n1 = self.q[0][i] * v[j] - self.q[0][j] * v[i];
            if n0 % self.det != 0 || n1 % self.det != 0 {
                return None;
            }
            let x = [n0 / self.det, n1 / self.det];
            (0..N)
                .all(|e| self.q[0][e] * x[0] + self.q[1][e] * x[1] == v[e])
                .then_some(x)
        }

        /// `v = line + x₁·q₁ + x₂·q₂` with a canonical `line` per parallel plane.
        fn split(&self, v: &[i128; N]) -> ([i128; N], [i128; 2]) {
            let (i, j) = self.rows;
            let n0 = self.q[1][j] * v[i] - self.q[1][i] * v[j];
            let n1 = self.q[0][i] * v[j] - self.q[0][j] * v[i];
            let x = [div_floor(n0, self.det), div_floor(n1, self.det)];
            (
                std::array::from_fn(|e| v[e] - self.q[0][e] * x[0] - self.q[1][e] * x[1]),
                x,
            )
        }
    }

    /// An affine form in `x` (integer coefficients), from one in `y` with `x = c + m∘y`,
    /// scaled by `m₀·m₁ > 0`.
    fn in_x(f: &Aff, m: [i128; 2], c: [i128; 2]) -> Aff {
        Aff {
            a: [f.a[0] * m[1], f.a[1] * m[0]],
            b: f.b * m[0] * m[1] - f.a[0] * m[1] * c[0] - f.a[1] * m[0] * c[1],
        }
    }

    /// Run one lap symbolically from `line + Q·(c + m∘y)`, following the branch
    /// taken at `y_rep`. Every condition must also hold at the corner `x0`.
    #[allow(clippy::too_many_arguments)]
    fn lap<const N: usize>(
        w: &Reduced<N>,
        pattern: &[i8; N],
        f: &Frame<N>,
        line: &[i128; N],
        m: [i128; 2],
        c: [i128; 2],
        y_rep: [i128; 2],
        x0: [i128; 2],
        normals: &mut Normals,
    ) -> Result<LapOk<N>, Err> {
        let q = &f.q;
        let mut v: [Aff; N] = std::array::from_fn(|e| Aff {
            a: [m[0] * q[0][e], m[1] * q[1][e]],
            b: line[e] + c[0] * q[0][e] + c[1] * q[1][e],
        });
        let mut need = |g: Aff| -> Result<(), Err> {
            let g = in_x(&g, m, c);
            if g.a[0] * x0[0] + g.a[1] * x0[1] + g.b < 0 {
                return Err(Err::Fail("a lap condition fails at the current return"));
            }
            normals.add(g.a)
        };
        let mut time = Aff::default();
        for events in 0..=MAX_LAP_EVENTS {
            let mut s = [0i8; N];
            for e in 0..N {
                let sg = v[e].at(y_rep).signum();
                if sg == 0 && v[e].a != [0, 0] {
                    return Err(Err::Fail("a coordinate is zero only sometimes"));
                }
                s[e] = sg as i8;
                if v[e].a != [0, 0] {
                    need(Aff {
                        a: [sg * v[e].a[0], sg * v[e].a[1]],
                        b: sg * v[e].b - 1,
                    })?;
                }
            }
            if events == 0 {
                if s != *pattern {
                    return Err(Err::Fail("the family does not start in the pattern"));
                }
            } else if s == *pattern {
                let cols: [[i128; N]; 2] = [
                    std::array::from_fn(|e| v[e].a[0]),
                    std::array::from_fn(|e| v[e].a[1]),
                ];
                let (Some(k0), Some(k1)) = (f.coords(&cols[0]), f.coords(&cols[1])) else {
                    return Err(Err::Fail("the lap ends off the plane"));
                };
                let b: [i128; N] = std::array::from_fn(|e| v[e].b);
                let (line, beta) = f.split(&b);
                let map_y = [[k0[0], k1[0]], [k0[1], k1[1]]];
                // In x: x' = map_y·(x − c)/m + beta, scaled by m₀·m₁.
                let rows: [Aff; 2] = std::array::from_fn(|i| {
                    in_x(
                        &Aff {
                            a: map_y[i],
                            b: beta[i],
                        },
                        m,
                        c,
                    )
                });
                return Ok(LapOk {
                    map: [rows[0].a, rows[1].a],
                    off: [rows[0].b, rows[1].b],
                    time: in_x(&time, m, c),
                    line,
                    beta: [beta[0].rem_euclid(m[0]), beta[1].rem_euclid(m[1])],
                    map_y,
                    events,
                });
            }
            let d = w.drift(&s);
            if d.iter().all(|&x| x == 0) {
                return Err(Err::Fail("the lap reaches a fixed pattern"));
            }
            if (0..N).all(|e| d[e] == 0 || (s[e] != 0 && d[e].signum() as i8 == s[e])) {
                return Err(Err::Fail("the lap reaches an absorbing pattern"));
            }
            let mut times = [Aff::default(); N];
            let mut n = 0;
            for e in 0..N {
                let de = d[e] as i128;
                if s[e] == 0 {
                    if de != 0 {
                        times[n] = Aff { a: [0, 0], b: 1 };
                        n += 1;
                    }
                } else if de.signum() as i8 == -s[e] {
                    let sg = s[e] as i128;
                    let (a, b) = ([sg * v[e].a[0], sg * v[e].a[1]], sg * v[e].b);
                    let dd = de.abs();
                    for j in 0..2 {
                        if a[j] % dd != 0 {
                            return Err(Err::Refine(j, dd / gcd(a[j], dd)));
                        }
                    }
                    times[n] = Aff {
                        a: [a[0] / dd, a[1] / dd],
                        b: div_ceil(b, dd),
                    };
                    n += 1;
                }
            }
            let times = &times[..n];
            // The branch the representative takes: its first sign change.
            let k = *times
                .iter()
                .min_by_key(|t| (t.at(y_rep), t.a[0], t.a[1]))
                .ok_or(Err::Fail("nothing moves"))?;
            for t in times {
                if *t != k {
                    need(Aff {
                        a: [t.a[0] - k.a[0], t.a[1] - k.a[1]],
                        b: t.b - k.b,
                    })?;
                }
            }
            need(Aff { a: k.a, b: k.b - 1 })?;
            time.a[0] += k.a[0];
            time.a[1] += k.a[1];
            time.b += k.b;
            for e in 0..N {
                let de = d[e] as i128;
                v[e].a[0] += k.a[0] * de;
                v[e].a[1] += k.a[1] * de;
                v[e].b += k.b * de;
                if v[e].a.iter().chain([&v[e].b]).any(|x| x.abs() > 1 << 100) {
                    return Err(Err::Fail("numbers too large"));
                }
            }
        }
        Err(Err::Fail("the lap does not return to the pattern"))
    }

    const MAX_STATES: usize = 64;

    pub fn prove<const N: usize>(
        w: &Reduced<N>,
        pattern: &[i8; N],
        returns: &[[i64; N]],
    ) -> Result<Lap<N>, &'static str> {
        let r0 = returns.last().ok_or(NOT_A_PLANE)?;
        let prim = |d: [i128; N]| {
            let g = d.iter().fold(0, |g, &x| gcd(g, x));
            (g != 0).then(|| std::array::from_fn::<i128, N, _>(|e| d[e] / g))
        };
        // Growth directions: differences between the current and earlier returns.
        let mut diffs = [[0i128; N]; RETURNS];
        let mut n = 0;
        for r in returns.iter().rev().skip(1) {
            if let Some(d) = prim(std::array::from_fn(|e| (r0[e] - r[e]) as i128)) {
                diffs[n] = d;
                n += 1;
            }
        }
        let diffs = &diffs[..n];
        let q1 = *diffs.first().ok_or(NOT_A_PLANE)?;
        let q2 = *diffs
            .iter()
            .find(|d| Frame::new([q1, **d]).is_some())
            .ok_or(NOT_A_PLANE)?;
        // Rank exactly two: every difference lies in the plane (det·d in the integer span).
        let probe = Frame::new([q1, q2]).ok_or(NOT_A_PLANE)?;
        if !diffs.iter().all(|d| {
            let scaled: [i128; N] = std::array::from_fn(|e| d[e] * probe.det);
            probe.coords(&scaled).is_some()
        }) {
            return Err(NOT_A_PLANE);
        }
        let (s1, s2) = saturate(q1, q2).ok_or(NOT_A_PLANE)?;
        let mut why = "no cone is mapped into itself";
        for basis in [reduce(s1, s2), (s1, s2)] {
            match prove_with(w, pattern, [basis.0, basis.1], r0) {
                Ok(l) => return Ok(l),
                Err(e) => why = e,
            }
        }
        Err(why)
    }

    /// A basis of every integer point in the plane spanned by `q1`, `q2` (the
    /// saturation): while all 2×2 minors share a prime `p`, the two vectors are
    /// dependent mod `p`, and dividing that combination by `p` enlarges the lattice.
    fn saturate<const N: usize>(
        mut q1: [i128; N],
        mut q2: [i128; N],
    ) -> Option<([i128; N], [i128; N])> {
        let minors = |q1: &[i128; N], q2: &[i128; N]| {
            let mut g = 0;
            for i in 0..N {
                for j in i + 1..N {
                    g = gcd(g, q1[i] * q2[j] - q1[j] * q2[i]);
                }
            }
            g
        };
        let pow = |mut b: i128, mut e: i128, p: i128| {
            let mut r = 1;
            b = b.rem_euclid(p);
            while e > 0 {
                if e & 1 == 1 {
                    r = r * b % p;
                }
                b = b * b % p;
                e >>= 1;
            }
            r
        };
        for _ in 0..128 {
            let g = minors(&q1, &q2);
            if g == 0 {
                return None;
            }
            if g == 1 {
                return Some((q1, q2));
            }
            // Factoring by trial division stays cheap below 2^32; a lattice that
            // coarse is no family worth trying anyway.
            if g > 1 << 32 {
                return None;
            }
            let p = (2..)
                .take_while(|d| d * d <= g)
                .find(|d| g % d == 0)
                .unwrap_or(g);
            if q1.iter().all(|x| x % p == 0) {
                q1 = std::array::from_fn(|e| q1[e] / p);
                continue;
            }
            // q2 ≡ t·q1 (mod p), with t from a coordinate where q1 is a unit.
            let e = (0..N).find(|&e| q1[e] % p != 0)?;
            let t = q2[e] * pow(q1[e], p - 2, p) % p;
            let next: [i128; N] = std::array::from_fn(|i| q2[i] - t * q1[i]);
            if next.iter().any(|x| x % p != 0) {
                return None;
            }
            q2 = std::array::from_fn(|i| next[i] / p);
        }
        None
    }

    /// Gauss-reduce a basis of a 2-D lattice so plane coordinates stay small.
    fn reduce<const N: usize>(q1: [i128; N], q2: [i128; N]) -> ([i128; N], [i128; N]) {
        let (mut q1, mut q2) = (q1, q2);
        let dot = |a: &[i128; N], b: &[i128; N]| (0..N).map(|e| a[e] * b[e]).sum::<i128>();
        for _ in 0..64 {
            if dot(&q1, &q1) > dot(&q2, &q2) {
                std::mem::swap(&mut q1, &mut q2);
            }
            let t = dot(&q1, &q2);
            let n = dot(&q1, &q1);
            let r = div_floor(2 * t + n, 2 * n);
            if r == 0 {
                break;
            }
            q2 = std::array::from_fn(|e| q2[e] - r * q1[e]);
        }
        (q1, q2)
    }

    fn cross(a: [i128; 2], b: [i128; 2]) -> i128 {
        a[0] * b[1] - a[1] * b[0]
    }

    fn dot(a: [i128; 2], b: [i128; 2]) -> i128 {
        a[0] * b[0] + a[1] * b[1]
    }

    /// `d` in the cone spanned by `e1`, `e2` (with `cross(e1, e2) > 0`).
    fn in_cone(d: [i128; 2], e1: [i128; 2], e2: [i128; 2]) -> bool {
        cross(e1, d) >= 0 && cross(d, e2) >= 0
    }

    /// Rational eigenvector directions of a 2×2 integer matrix.
    fn eigen(a: [[i128; 2]; 2]) -> [[i128; 2]; 4] {
        let (tr, det) = (a[0][0] + a[1][1], a[0][0] * a[1][1] - a[0][1] * a[1][0]);
        let disc = tr * tr - 4 * det;
        let mut out = [[0; 2]; 4];
        if disc < 0 {
            return out;
        }
        let r = (disc as f64).sqrt().round() as i128;
        if r * r != disc {
            return out;
        }
        for (i, sgn) in [1, -1].into_iter().enumerate() {
            // (A − λ)v = 0 with 2λ = tr ± r.
            let two_l = tr + sgn * r;
            let v = if a[0][1] != 0 || 2 * a[0][0] != two_l {
                [2 * a[0][1], two_l - 2 * a[0][0]]
            } else {
                [two_l - 2 * a[1][1], 2 * a[1][0]]
            };
            out[2 * i] = v;
            out[2 * i + 1] = [-v[0], -v[1]];
        }
        out
    }

    fn prove_with<const N: usize>(
        w: &Reduced<N>,
        pattern: &[i8; N],
        q: [[i128; N]; 2],
        r0: &[i64; N],
    ) -> Result<Lap<N>, &'static str> {
        let f = Frame::new(q).ok_or(NOT_A_PLANE)?;
        let start: [i128; N] = std::array::from_fn(|e| r0[e] as i128);
        let (p0, x0) = f.split(&start);
        let mut m = [1i128, 1];
        'refine: loop {
            if m[0] * m[1] > MAX_MODULUS {
                return Err("needs too many residue classes");
            }
            let c0 = [x0[0].rem_euclid(m[0]), x0[1].rem_euclid(m[1])];
            // States (plane, residue class), visited breadth first.
            let mut states = [([0i128; N], [0i128; 2]); MAX_STATES];
            let mut laps = [None::<LapOk<N>>; MAX_STATES];
            let mut normals = Normals {
                dirs: [[0; 2]; MAX_NORMALS],
                n: 0,
            };
            states[0] = (p0, c0);
            let mut n_states = 1;
            let mut i = 0;
            while i < n_states {
                let (line, c) = states[i];
                // A representative near the current return, in this residue class.
                let y_rep = [div_floor(x0[0] - c[0], m[0]), div_floor(x0[1] - c[1], m[1])];
                let ok = match lap(w, pattern, &f, &line, m, c, y_rep, x0, &mut normals) {
                    Ok(ok) => ok,
                    Err(Err::Refine(j, by)) => {
                        m[j] *= by;
                        continue 'refine;
                    }
                    Err(Err::Fail(why)) => return Err(why),
                };
                // Residues of x' = map_y·y + beta over all y: beta plus the column span.
                let mut res = [[0i128; 2]; MAX_STATES];
                res[0] = ok.beta;
                let (mut n_res, mut j) = (1, 0);
                while j < n_res {
                    for col in 0..2 {
                        let r = [
                            (res[j][0] + ok.map_y[0][col]).rem_euclid(m[0]),
                            (res[j][1] + ok.map_y[1][col]).rem_euclid(m[1]),
                        ];
                        if !res[..n_res].contains(&r) {
                            if n_res == MAX_STATES {
                                return Err("needs too many residue classes");
                            }
                            res[n_res] = r;
                            n_res += 1;
                        }
                    }
                    j += 1;
                }
                for &r in &res[..n_res] {
                    if !states[..n_states].contains(&(ok.line, r)) {
                        if n_states == MAX_STATES {
                            return Err("needs too many residue classes");
                        }
                        states[n_states] = (ok.line, r);
                        n_states += 1;
                    }
                }
                laps[i] = Some(ok);
                i += 1;
            }
            let laps: Vec<LapOk<N>> = laps[..n_states].iter().flatten().copied().collect();
            let mm = m[0] * m[1];
            let normals = &normals.dirs[..normals.n];
            // Candidate edges: boundary lines of the conditions, axes, eigenvectors.
            let mut cands: Vec<[i128; 2]> = vec![[1, 0], [-1, 0], [0, 1], [0, -1]];
            for n in normals {
                cands.push([-n[1], n[0]]);
                cands.push([n[1], -n[0]]);
            }
            cands.extend(eigen(laps[0].map).into_iter().filter(|v| *v != [0, 0]));
            cands.retain(|d| normals.iter().all(|n| dot(*n, *d) >= 0));
            cands.sort();
            cands.dedup();
            for (a, &e1) in cands.iter().enumerate() {
                for &e2 in &cands[a + 1..] {
                    let (e1, e2) = match cross(e1, e2) {
                        0 => continue,
                        x if x > 0 => (e1, e2),
                        _ => (e2, e1),
                    };
                    let image = |l: &LapOk<N>, d: [i128; 2]| [dot(l.map[0], d), dot(l.map[1], d)];
                    let ok = laps.iter().all(|l| {
                        // Corner and edges map into the cone.
                        let corner = [
                            dot(l.map[0], x0) + l.off[0] - mm * x0[0],
                            dot(l.map[1], x0) + l.off[1] - mm * x0[1],
                        ];
                        in_cone(corner, e1, e2)
                            && in_cone(image(l, e1), e1, e2)
                            && in_cone(image(l, e2), e1, e2)
                            && dot(l.time.a, e1) > 0
                            && dot(l.time.a, e2) > 0
                    });
                    if !ok {
                        continue;
                    }
                    // A functional that grows by ≥ 1 per lap on the whole cone.
                    let grows = [
                        [1, 0],
                        [-1, 0],
                        [0, 1],
                        [0, -1],
                        [-e2[1], e2[0]],
                        [e1[1], -e1[0]],
                    ]
                    .into_iter()
                    .any(|g: [i128; 2]| {
                        laps.iter().all(|l| {
                            // ℓ(x') − ℓ(x) − 1 ≥ 0, scaled by mm.
                            let lin = [
                                g[0] * l.map[0][0] + g[1] * l.map[1][0] - mm * g[0],
                                g[0] * l.map[0][1] + g[1] * l.map[1][1] - mm * g[1],
                            ];
                            let at0 = dot(lin, x0) + g[0] * l.off[0] + g[1] * l.off[1] - mm;
                            at0 >= 0 && dot(lin, e1) >= 0 && dot(lin, e2) >= 0
                        })
                    });
                    if !grows {
                        continue;
                    }
                    let l = &laps[0];
                    let g = l
                        .map
                        .iter()
                        .flatten()
                        .chain(l.off.iter())
                        .fold(mm, |g, &x| gcd(g, x));
                    return Ok(Lap {
                        rank: 2,
                        base: std::array::from_fn(|e| p0[e] as i64),
                        dir: [
                            std::array::from_fn(|e| q[0][e] as i64),
                            std::array::from_fn(|e| q[1][e] as i64),
                        ],
                        steady: 0,
                        slope: [
                            [(l.map[0][0] / g) as i64, (l.map[0][1] / g) as i64],
                            [(l.map[1][0] / g) as i64, (l.map[1][1] / g) as i64],
                        ],
                        off: [(l.off[0] / g) as i64, (l.off[1] / g) as i64],
                        den: (mm / g) as i64,
                        modulus: [m[0] as u32, m[1] as u32],
                        classes: n_states as u32,
                        events: l.events,
                    });
                }
            }
            return Err("no cone is mapped into itself");
        }
    }
}
