//! Is a universe bounded from *every* start, not just the all-fire one?
//!
//! **Unstable** needs one witness: a start whose fate is proven unbounded (a
//! ray, a helix or a proven spiral). The search tries the absorbing sign
//! patterns, every start in a small box, and starts far out along every face.
//!
//! **Bounded** needs a proof about infinitely many starts. For `N = 2` it comes
//! from the far field. A state far from the origin is either inside a quadrant,
//! where it moves in a straight line until a coordinate crosses zero, or *near
//! an axis*: one coordinate `u` small (`|u| ≤ K`, the largest drift) and the
//! other, `X`, large. Near an axis the small coordinate's motion does not depend
//! on `X` at all (only its sign does), so it can be run once with `X = ∞`. That
//! run ends in one of five ways, valid for every `X` above a threshold:
//!
//! * the drift becomes zero, or `u` cycles with `X` coming back exactly: bounded;
//! * `u` cycles while `X` shrinks: the orbit drains into the box;
//! * `u` cycles while `X` grows, or the pattern is absorbing: an escape;
//! * `u` runs off while `X` heads to zero: when `X` crosses, `u` has become
//!   `X' = λ·(X + δ) + c` with `λ = |d_u| / |d_X|`, a state near the *other*
//!   axis. The rounding of the crossing time only depends on `X mod |d_X|`.
//!
//! The last case gives a finite graph of axis states with affine maps on `X`.
//! If every simple cycle of the graph contracts (`Λ < 1`) or translates
//! non-positively (`Λ = 1`, shift `≤ 0`), then along any walk `X` stays bounded
//! (a closed walk is a product of simple cycles, and a potential `φ·X` exists).
//! So every orbit is periodic, stays bounded in the graph, or enters the box
//! `[-R, R]²`; if every start in the box is a fixed point or a cycle, the
//! universe is bounded from every start. For `N ≥ 3` the far field has states
//! with two independent large coordinates, so only instability is proven.

// Index loops mirror the formulas.
#![allow(clippy::needless_range_loop)]

use crate::exact::{analyze_ws, Budget, Fate, Reduced, Workspace};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Stability<const N: usize> {
    /// The trajectory from `start` is proven unbounded.
    Unstable { start: [i64; N], fate: Fate },
    /// Proven bounded from every start (`N = 2`): the far-field graph has
    /// `nodes` axis states and `cycles` simple cycles, none of them expanding,
    /// and every start in `[-radius, radius]²` is fixed or a cycle.
    Bounded {
        radius: i64,
        nodes: u32,
        cycles: u32,
    },
    /// No start tried escapes, but there is no proof.
    Open { why: &'static str },
}

impl<const N: usize> Stability<N> {
    pub fn kind(&self) -> &'static str {
        match self {
            Stability::Unstable { .. } => "unstable",
            Stability::Bounded { .. } => "bounded",
            Stability::Open { .. } => "open",
        }
    }
}

/// Largest drift magnitude of any coordinate in any sign pattern.
fn max_drift<const N: usize>(w: &Reduced<N>) -> i64 {
    let mut k = 0;
    for_each_pattern::<N>(|s| {
        k = w.drift(s).iter().fold(k, |k, x| k.max(x.abs()));
    });
    k
}

fn for_each_pattern<const N: usize>(mut f: impl FnMut(&[i8; N])) {
    let mut s = [-1i8; N];
    loop {
        f(&s);
        let mut e = 0;
        while e < N && s[e] == 1 {
            s[e] = -1;
            e += 1;
        }
        if e == N {
            return;
        }
        s[e] += 1;
    }
}

fn absorbing<const N: usize>(s: &[i8; N], d: &[i64; N]) -> bool {
    d.iter().any(|&x| x != 0)
        && (0..N).all(|e| d[e] == 0 || (s[e] != 0 && d[e].signum() as i8 == s[e]))
}

/// A sign pattern whose drift never changes a sign and is not zero: its own
/// sign vector is a start that escapes along a ray. Cheap, and it settles most
/// universes before any trajectory is run.
pub fn absorbing_pattern<const N: usize>(w: &Reduced<N>) -> Option<[i8; N]> {
    let mut found = None;
    for_each_pattern::<N>(|s| {
        if found.is_none() && absorbing(s, &w.drift(s)) {
            found = Some(*s);
        }
    });
    found
}

/// How far out the witness search looks.
const FAR: i64 = 1000;

/// Look for a start with a proven unbounded fate: the all-fire start, the box `[-radius, radius]^N`,
/// then starts with some coordinates far out (`±FAR`, `±(FAR + 1)`, and an
/// uneven pair) and the rest small.
pub fn witness<const N: usize>(
    w: &Reduced<N>,
    radius: i64,
    small: i64,
    budget: Budget,
    ws: &mut Workspace<N>,
) -> Option<([i64; N], Fate)> {
    let mut try_start = |v: [i64; N]| {
        let fate = analyze_ws(w, v, 0, budget, ws, &mut |_, _, _| {});
        (fate.bounded() == Some(false)).then_some((v, fate))
    };
    if let Some(s) = absorbing_pattern(w) {
        return try_start(s.map(|x| x as i64));
    }
    if let Some(found) = try_start(w.initial_state()) {
        return Some(found);
    }
    let side = 2 * radius + 1;
    for k in 0..side.pow(N as u32) {
        let mut v = [0i64; N];
        let mut rem = k;
        for x in v.iter_mut() {
            *x = rem % side - radius;
            rem /= side;
        }
        if let Some(found) = try_start(v) {
            return Some(found);
        }
    }
    // Far starts: a nonempty set L of large coordinates with signs, the rest in
    // [-small, small]. Magnitudes: all FAR, then FAR + 1 on the first (parity),
    // then FAR and 7·FAR/4 alternating (an uneven ratio).
    let sside = 2 * small + 1;
    for large in 1u32..(1 << N) {
        let l: Vec<usize> = (0..N).filter(|&e| large >> e & 1 == 1).collect();
        let s_count = (N - l.len()) as u32;
        for signs in 0u32..(1 << l.len()) {
            for shape in 0..3 {
                if shape == 2 && l.len() < 2 {
                    continue;
                }
                for k in 0..sside.pow(s_count) {
                    let mut v = [0i64; N];
                    let mut rem = k;
                    for e in (0..N).filter(|&e| large >> e & 1 == 0) {
                        v[e] = rem % sside - small;
                        rem /= sside;
                    }
                    for (j, &e) in l.iter().enumerate() {
                        let mag = match shape {
                            0 => FAR,
                            1 => FAR + (j == 0) as i64,
                            _ => {
                                if j % 2 == 0 {
                                    FAR
                                } else {
                                    7 * FAR / 4
                                }
                            }
                        };
                        v[e] = if signs >> j & 1 == 1 { mag } else { -mag };
                    }
                    if let Some(found) = try_start(v) {
                        return Some(found);
                    }
                }
            }
        }
    }
    None
}

/// What the near-axis run from one axis state does, for every large `X`.
#[derive(Clone, Copy, Debug)]
enum Run {
    /// Fixed, a cycle, or drains into the box: bounded once the box is.
    Settles,
    /// A ray or a helix along the axis.
    Escapes,
    /// `u` runs off while `X` (offset by `off`) heads to zero at speed `de`;
    /// `u` is the small coordinate's value then and `du` its speed.
    Exit { off: i64, u: i64, de: i64, du: i64 },
}

/// An axis state's run and how far out it is valid.
#[derive(Clone, Copy, Debug)]
struct NodeRun {
    run: Run,
    /// The run is exact for `X ≥ threshold`.
    #[allow(dead_code)]
    threshold: i64,
    /// Largest `|offset of X|` and `|u|` along the run (for the box radius).
    reach: i64,
}

/// Run the small coordinate with the large one (`e`, sign `sigma`) at `X = ∞`.
fn near_axis(w: &Reduced<2>, e: usize, sigma: i8, u0: i64) -> Result<NodeRun, &'static str> {
    let o = 1 - e;
    let (mut off, mut u) = (0i64, u0);
    let (mut min_off, mut reach) = (0i64, u0.abs());
    let mut seen: [(i64, i64); 64] = [(0, 0); 64];
    for n in 0..seen.len() {
        let mut s = [0i8; 2];
        s[e] = sigma;
        s[o] = u.signum() as i8;
        let d = w.drift(&s);
        let threshold = 1 - min_off;
        let done = |run| {
            Ok(NodeRun {
                run,
                threshold,
                reach: reach.max(threshold + off.abs()),
            })
        };
        if d == [0, 0] {
            return done(Run::Settles);
        }
        if absorbing(&s, &d) {
            return done(Run::Escapes);
        }
        // u decides everything that follows, so a repeat of u is a cycle;
        // X moved by the change in offset (positive: away from zero).
        if let Some(&(_, off0)) = seen[..n].iter().find(|(u1, _)| *u1 == u) {
            return done(if off > off0 {
                Run::Escapes
            } else {
                Run::Settles
            });
        }
        seen[n] = (u, off);
        // Next sign change of u, if u is heading for one.
        let k = if s[o] == 0 {
            (d[o] != 0).then_some(1)
        } else if d[o].signum() as i8 == -s[o] {
            Some(u.unsigned_abs().div_ceil(d[o].unsigned_abs()) as i64)
        } else {
            None
        };
        match k {
            Some(k) => {
                u += k * d[o];
                off += k * sigma as i64 * d[e];
                min_off = min_off.min(off);
                reach = reach.max(u.abs()).max(off.abs());
            }
            // Only X can change sign (it heads to zero: the pattern is not absorbing).
            None => {
                return done(Run::Exit {
                    off,
                    u,
                    de: d[e].abs(),
                    du: d[o],
                })
            }
        }
    }
    Err("an axis run is too long")
}

/// `X ↦ (p·X + c) / q` from one axis state to another.
#[derive(Clone, Copy, Debug)]
struct Edge {
    to: usize,
    p: i128,
    c: i128,
    q: i128,
}

fn gcd(mut a: i128, mut b: i128) -> i128 {
    a = a.abs();
    b = b.abs();
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

const MAX_CYCLES: u32 = 200_000;

/// Prove a two-element universe bounded from every start, or say why not.
/// Starts that turn out to escape on the way are returned as `Unstable`.
pub fn bounded_far_field(
    w: &Reduced<2>,
    budget: Budget,
    ws: &mut Workspace<2>,
) -> Result<Stability<2>, &'static str> {
    let k = max_drift(w).max(1);
    let width = (2 * k + 1) as usize;
    // Axis state (large coordinate e, its sign, small coordinate u) ↦ index.
    let id =
        |e: usize, sigma: i8, u: i64| (e * 2 + (sigma > 0) as usize) * width + (u + k) as usize;
    let n = 4 * width;
    let mut edges: Vec<Vec<Edge>> = vec![Vec::new(); n];
    let mut radius = k;
    for e in 0..2 {
        for sigma in [-1i8, 1] {
            for u0 in -k..=k {
                let r = near_axis(w, e, sigma, u0)?;
                radius = radius.max(r.reach + k);
                let Run::Exit { off, u, de, du } = r.run else {
                    if let Run::Escapes = r.run {
                        return Err("an axis run escapes");
                    }
                    continue;
                };
                if du == 0 {
                    continue; // u stops moving once X has crossed: a small state, in the box
                }
                // With (X + off) ≡ res (mod de), X crosses after (X + off + ρ) / de
                // steps, ρ = (de − res) mod de, and lands on −σ·ρ.
                let (p, q) = (du.abs() as i128, de as i128);
                for res in 0..de {
                    let rho = (de - res) % de;
                    edges[id(e, sigma, u0)].push(Edge {
                        to: id(1 - e, du.signum() as i8, -(sigma as i64) * rho),
                        p,
                        c: p * (off + rho) as i128 + q * u.abs() as i128,
                        q,
                    });
                }
            }
        }
    }
    // Every simple cycle must contract, or translate by a non-positive amount.
    let mut cycles = 0u32;
    let mut on_path = vec![false; n];
    for start in 0..n {
        // Depth-first over simple cycles whose smallest state is `start`; maps[i]
        // is the composed map X ↦ (P·X + C) / Q from start to path[i].
        let mut path: Vec<(usize, usize)> = vec![(start, 0)];
        let mut maps: Vec<(i128, i128, i128)> = vec![(1, 0, 1)];
        on_path[start] = true;
        while let Some(&(node, i)) = path.last() {
            if i == edges[node].len() {
                on_path[node] = false;
                path.pop();
                maps.pop();
                continue;
            }
            path.last_mut().unwrap().1 += 1;
            let ed = edges[node][i];
            if ed.to < start {
                continue;
            }
            let (pp, cc, qq) = *maps.last().unwrap();
            let (np, nc, nq) = (ed.p * pp, ed.p * cc + ed.c * qq, ed.q * qq);
            let g = gcd(gcd(np, nc), nq);
            let (np, nc, nq) = (np / g, nc / g, nq / g);
            if ed.to == start {
                cycles += 1;
                if np > nq || (np == nq && nc > 0) {
                    return Err("a far-field cycle expands");
                }
                if cycles > MAX_CYCLES {
                    return Err("too many far-field cycles");
                }
            } else if !on_path[ed.to] {
                on_path[ed.to] = true;
                path.push((ed.to, 0));
                maps.push((np, nc, nq));
            }
        }
    }
    // Every start in the box must settle.
    for x in -radius..=radius {
        for y in -radius..=radius {
            let fate = analyze_ws(w, [x, y], 0, budget, ws, &mut |_, _, _| {});
            match fate.bounded() {
                Some(true) => {}
                Some(false) => {
                    return Ok(Stability::Unstable {
                        start: [x, y],
                        fate,
                    })
                }
                None => return Err("a start in the box is undecided"),
            }
        }
    }
    Ok(Stability::Bounded {
        radius,
        nodes: n as u32,
        cycles,
    })
}

/// Classify one universe: a witness of escape, a proof of boundedness (`N = 2`),
/// or open.
pub fn stability<const N: usize>(
    w: &Reduced<N>,
    budget: Budget,
    ws: &mut Workspace<N>,
) -> Stability<N> {
    let (radius, small) = if N <= 2 { (3, max_drift(w)) } else { (2, 1) };
    if let Some((start, fate)) = witness(w, radius, small, budget, ws) {
        return Stability::Unstable { start, fate };
    }
    if N == 1 {
        // No absorbing pattern: far from zero the state heads straight for it,
        // so every orbit reaches [-k, k], and the box search above covered that.
        return Stability::Bounded {
            radius: radius.max(max_drift(w)),
            nodes: 0,
            cycles: 0,
        };
    }
    if N != 2 {
        return Stability::Open {
            why: "no start escapes; a proof needs two large coordinates at once",
        };
    }
    let w2 = Reduced::<2>::from_index(w.index());
    match bounded_far_field(&w2, budget, &mut Workspace::new()) {
        Ok(Stability::Unstable { start, fate }) => Stability::Unstable {
            start: std::array::from_fn(|e| start[e]),
            fate,
        },
        Ok(Stability::Bounded {
            radius,
            nodes,
            cycles,
        }) => Stability::Bounded {
            radius,
            nodes,
            cycles,
        },
        Ok(Stability::Open { why }) | Err(why) => Stability::Open { why },
    }
}
