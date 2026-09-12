//! Exact analysis of a universe.
//!
//! After evaluation each element carries a single integer `v_e = count + antiCount`
//! (exactly one of the two is nonzero). Which channel fires depends only on
//! `sign(v_e)`, and a firing row adds `M[r][2e] - M[r][2e+1] ∈ {-1, 0, 1}` to
//! `v_e`. So the whole model collapses to
//!
//! ```text
//!     v(t+1) = v(t) + d(sign v(t)),
//!     d(s)   = Σ_{f: s_f = +1} pos[f] + Σ_{f: s_f = -1} neg[f]
//! ```
//!
//! with ternary drift rows `pos[f], neg[f] ∈ {-1,0,1}^N`. Inside one sign
//! pattern (a generalised orthant) the trajectory is a straight line
//! `v(t0) + (t - t0)·d`, so we only need to simulate the *events* where a
//! coordinate crosses zero. That makes three fates provable:
//!
//! * **Fixed**: `d = 0` in the current orthant — nothing will ever change.
//! * **Ray**: the drift points away from every face of the orthant, so the
//!   trajectory stays there forever and grows linearly — a certificate of
//!   divergence, no step limit needed.
//! * **Cycle**: an event state repeats, so the trajectory is eventually periodic.
//!
//! What is left (`Undecided`) is a trajectory that keeps switching orthants
//! with growing norm — the discrete analogue of an expanding spiral.

// Index loops mirror the matrix formulas and read better than zipped iterators here.
#![allow(clippy::needless_range_loop)]

use std::collections::HashMap;

/// Event callback: `(step, state, drift)`.
pub type Trace<'a, const N: usize> = &'a mut dyn FnMut(u64, &[i64; N], &[i64; N]);

/// Drift rows of a universe with `N` elements (its `M` matrix reduced modulo
/// everything that does not affect the dynamics).
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, PartialOrd, Ord)]
pub struct Reduced<const N: usize> {
    /// `pos[f][e]`: change of `v_e` per step while `v_f > 0`.
    pub pos: [[i8; N]; N],
    /// `neg[f][e]`: change of `v_e` per step while `v_f < 0`.
    pub neg: [[i8; N]; N],
}

impl<const N: usize> Reduced<N> {
    /// Number of distinct reduced universes: `3^(2N²)`.
    pub fn count() -> u128 {
        3u128.pow((2 * N * N) as u32)
    }

    /// Reduce a bit-row matrix (`2N` rows, bit `c` of row `r` = `M[r][c]`).
    pub fn from_rows(rows: &[u16]) -> Self {
        assert_eq!(rows.len(), 2 * N);
        let entry = |r: usize, e: usize| {
            ((rows[r] >> (2 * e)) & 1) as i8 - ((rows[r] >> (2 * e + 1)) & 1) as i8
        };
        let mut w = Reduced {
            pos: [[0; N]; N],
            neg: [[0; N]; N],
        };
        for f in 0..N {
            for e in 0..N {
                w.pos[f][e] = entry(2 * f, e);
                w.neg[f][e] = entry(2 * f + 1, e);
            }
        }
        w
    }

    /// A representative bit-row matrix: `+1 → (1,0)`, `-1 → (0,1)`, `0 → (0,0)`.
    pub fn to_rows(&self) -> Vec<u16> {
        let mut rows = vec![0u16; 2 * N];
        for f in 0..N {
            for e in 0..N {
                for (r, w) in [(2 * f, self.pos[f][e]), (2 * f + 1, self.neg[f][e])] {
                    match w {
                        1 => rows[r] |= 1 << (2 * e),
                        -1 => rows[r] |= 1 << (2 * e + 1),
                        _ => {}
                    }
                }
            }
        }
        rows
    }

    /// Decode a base-3 index (`2N²` digits, `pos` rows first, row-major, most significant first).
    pub fn from_index(mut index: u64) -> Self {
        let mut w = Reduced {
            pos: [[0; N]; N],
            neg: [[0; N]; N],
        };
        for f in (0..N).rev() {
            for e in (0..N).rev() {
                w.neg[f][e] = (index % 3) as i8 - 1;
                index /= 3;
            }
        }
        for f in (0..N).rev() {
            for e in (0..N).rev() {
                w.pos[f][e] = (index % 3) as i8 - 1;
                index /= 3;
            }
        }
        w
    }

    /// Base-3 index of this universe (inverse of [`Reduced::from_index`]).
    pub fn index(&self) -> u64 {
        self.pos
            .iter()
            .chain(self.neg.iter())
            .flatten()
            .fold(0u64, |acc, &x| acc * 3 + (x + 1) as u64)
    }

    /// Number of bit matrices `M` that reduce to this: `2^(number of zero entries)`.
    pub fn multiplicity(&self) -> u64 {
        let zeros = self
            .pos
            .iter()
            .chain(self.neg.iter())
            .flatten()
            .filter(|&&w| w == 0)
            .count();
        1u64 << zeros
    }

    /// State after the first step, when every channel fires: `Σ_f pos[f] + neg[f]`.
    pub fn initial_state(&self) -> [i64; N] {
        let mut v = [0i64; N];
        for f in 0..N {
            for e in 0..N {
                v[e] += (self.pos[f][e] + self.neg[f][e]) as i64;
            }
        }
        v
    }

    /// Drift for a sign pattern.
    #[inline]
    pub fn drift(&self, sign: &[i8; N]) -> [i64; N] {
        let mut d = [0i64; N];
        for f in 0..N {
            let row = match sign[f] {
                1 => &self.pos[f],
                -1 => &self.neg[f],
                _ => continue,
            };
            for e in 0..N {
                d[e] += row[e] as i64;
            }
        }
        d
    }

    /// Apply a signed permutation: element `e` goes to position `perm[e]` with
    /// its sign multiplied by `flip[e]` (`+1`/`-1`). Conjugate universes have
    /// conjugate trajectories, so fates are invariant under this action.
    pub fn transform(&self, perm: &[usize; N], flip: &[i8; N]) -> Self {
        let mut w = Reduced {
            pos: [[0; N]; N],
            neg: [[0; N]; N],
        };
        for f in 0..N {
            for e in 0..N {
                let (p, n) = if flip[f] == 1 {
                    (self.pos[f][e], self.neg[f][e])
                } else {
                    (self.neg[f][e], self.pos[f][e])
                };
                w.pos[perm[f]][perm[e]] = flip[e] * p;
                w.neg[perm[f]][perm[e]] = flip[e] * n;
            }
        }
        w
    }

    /// Lexicographically smallest member of this universe's symmetry orbit.
    pub fn canonical(&self) -> Self {
        let mut best = *self;
        for_each_signed_permutation::<N>(|perm, flip| {
            let t = self.transform(perm, flip);
            if t < best {
                best = t;
            }
        });
        best
    }
}

/// Enumerate the hyperoctahedral group `B_N` (all `2^N · N!` signed permutations).
pub fn for_each_signed_permutation<const N: usize>(mut f: impl FnMut(&[usize; N], &[i8; N])) {
    let mut perm: [usize; N] = std::array::from_fn(|i| i);
    loop {
        for mask in 0..(1u32 << N) {
            let flip: [i8; N] = std::array::from_fn(|i| if (mask >> i) & 1 == 1 { -1 } else { 1 });
            f(&perm, &flip);
        }
        if !next_permutation(&mut perm) {
            break;
        }
    }
}

fn next_permutation(p: &mut [usize]) -> bool {
    let n = p.len();
    if n < 2 {
        return false;
    }
    let mut i = n - 1;
    while i > 0 && p[i - 1] >= p[i] {
        i -= 1;
    }
    if i == 0 {
        return false;
    }
    let mut j = n - 1;
    while p[j] <= p[i - 1] {
        j -= 1;
    }
    p.swap(i - 1, j);
    p[i..].reverse();
    true
}

/// Number of symmetry orbits of reduced universes (Burnside's lemma).
pub fn reduced_orbits<const N: usize>() -> u128 {
    // Positions: (f, side, e) with side 0 = pos, 1 = neg. A group element permutes
    // positions and multiplies entries by flip[e]; an entry is fixed by a cycle
    // iff the product of signs along the cycle is +1 (3 choices) else it must be 0.
    let mut total: u128 = 0;
    let mut order: u128 = 0;
    for_each_signed_permutation::<N>(|perm, flip| {
        order += 1;
        let idx = |f: usize, side: usize, e: usize| (f * 2 + side) * N + e;
        let len = 2 * N * N;
        let mut target = vec![0usize; len];
        let mut sign = vec![1i8; len];
        for f in 0..N {
            for side in 0..2 {
                for e in 0..N {
                    let new_side = if flip[f] == 1 { side } else { 1 - side };
                    target[idx(f, side, e)] = idx(perm[f], new_side, perm[e]);
                    sign[idx(f, side, e)] = flip[e];
                }
            }
        }
        let mut seen = vec![false; len];
        let mut fixed: u128 = 1;
        for start in 0..len {
            if seen[start] {
                continue;
            }
            let (mut i, mut s) = (start, 1i8);
            loop {
                seen[i] = true;
                s *= sign[i];
                i = target[i];
                if i == start {
                    break;
                }
            }
            fixed *= if s == 1 { 3 } else { 1 };
        }
        total += fixed;
    });
    total / order
}

/// Number of symmetry orbits of bit matrices `M` (Burnside's lemma; pure permutation of `4N²` bits).
pub fn matrix_orbits<const N: usize>() -> u128 {
    let c = 2 * N;
    let mut total: u128 = 0;
    let mut order: u128 = 0;
    for_each_signed_permutation::<N>(|perm, flip| {
        order += 1;
        // channel map: 2e+k -> 2perm[e] + (k xor (flip[e] == -1))
        let ch = |x: usize| 2 * perm[x / 2] + ((x % 2) ^ (flip[x / 2] == -1) as usize);
        let len = c * c;
        let mut seen = vec![false; len];
        let mut cycles = 0u32;
        for start in 0..len {
            if seen[start] {
                continue;
            }
            cycles += 1;
            let mut i = start;
            loop {
                seen[i] = true;
                i = ch(i / c) * c + ch(i % c);
                if i == start {
                    break;
                }
            }
        }
        total += 1u128 << cycles;
    });
    total / order
}

/// Outcome of an exact analysis.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum Fate {
    /// Drift is zero in the current orthant from `step` on (includes the origin with nothing firing).
    Fixed { step: u64 },
    /// An event state repeated: eventually periodic with the given transient and period.
    Cycle { transient: u64, period: u64 },
    /// From `step` on the trajectory never leaves its orthant and grows linearly.
    Ray { step: u64 },
    /// From `step` on, `v(t + period) = v(t) + shift` forever (`shift != 0`): a
    /// drifting cycle. Proven by replay: an earlier event state with the same
    /// sign pattern differs by `shift`, every coordinate with a nonzero shift kept
    /// a constant sign along the segment and shifts away from zero, so the segment
    /// repeats verbatim. The analogue of Karp–Miller's `M →σ M + L` test.
    Helix {
        step: u64,
        period: u64,
        shift_norm: i64,
    },
    /// Empirical: the trajectory keeps returning to the same orthant with strictly
    /// growing norm (checked over `returns` consecutive returns). Not a proof.
    Spiral { returns: u64, max_norm: i64 },
    /// No certificate within the budget; the norm reached `max_norm` after `events` orthant changes.
    Undecided { events: u64, max_norm: i64 },
}

impl Fate {
    /// Coarse class name.
    pub fn kind(&self) -> &'static str {
        match self {
            Fate::Fixed { .. } => "fixed",
            Fate::Cycle { .. } => "cycle",
            Fate::Ray { .. } => "ray",
            Fate::Helix { .. } => "helix",
            Fate::Spiral { .. } => "spiral",
            Fate::Undecided { .. } => "undecided",
        }
    }

    /// `true` when the trajectory is proven bounded.
    pub fn bounded(&self) -> Option<bool> {
        match self {
            Fate::Fixed { .. } | Fate::Cycle { .. } => Some(true),
            Fate::Ray { .. } | Fate::Helix { .. } => Some(false),
            Fate::Spiral { .. } | Fate::Undecided { .. } => None,
        }
    }
}

/// Limits for the event loop.
#[derive(Clone, Copy, Debug)]
pub struct Budget {
    pub max_events: u64,
    pub max_norm: i64,
    /// Consecutive returns to the section with strictly increasing norm needed to call a spiral (0 disables).
    pub spiral_returns: u64,
}

impl Default for Budget {
    fn default() -> Self {
        Budget {
            max_events: 10_000,
            max_norm: 1 << 56,
            spiral_returns: 8,
        }
    }
}

/// How many past events the helix search looks back over.
const HELIX_WINDOW: usize = 256;

/// Search the event log for an earlier state `u` with the same pattern such that
/// `v = u + D` replays forever. `pattern` is the per-coordinate regime (sign, or
/// support for the clamped model); `allowed(e, D_e, regime_e)` says whether a nonzero
/// shift on coordinate `e` moves it away from every regime boundary.
fn helix_shift<const N: usize>(
    log: &[([i64; N], u64, [i8; N])],
    v: &[i64; N],
    pattern: &[i8; N],
    allowed: impl Fn(usize, i64, i8) -> bool,
) -> Option<(u64, i64)> {
    let start = log.len().saturating_sub(HELIX_WINDOW);
    for j in (start..log.len()).rev() {
        let (u, ustep, upat) = &log[j];
        if upat != pattern {
            continue;
        }
        let shift: [i64; N] = std::array::from_fn(|e| v[e] - u[e]);
        if shift.iter().all(|&x| x == 0) {
            continue;
        }
        let ok = (0..N).all(|e| {
            shift[e] == 0
                || (allowed(e, shift[e], pattern[e])
                    && log[j..].iter().all(|(_, _, p)| p[e] == pattern[e]))
        });
        if ok {
            let norm = shift.iter().map(|x| x.abs()).max().unwrap_or(0);
            return Some((*ustep, norm));
        }
    }
    None
}

/// Track returns to each pattern; report a spiral once some pattern has been
/// re-entered `needed` consecutive times with strictly growing max-norm.
fn spiral_return<const N: usize>(
    returns: &mut HashMap<[i8; N], (i64, u64)>,
    pattern: &[i8; N],
    v: &[i64; N],
    needed: u64,
    max_norm: i64,
) -> Option<Fate> {
    let norm = v.iter().map(|x| x.abs()).max().unwrap_or(0);
    match returns.get_mut(pattern) {
        None => {
            returns.insert(*pattern, (norm, 0));
        }
        Some(entry) => {
            entry.1 = if norm > entry.0 { entry.1 + 1 } else { 0 };
            entry.0 = norm;
            if entry.1 >= needed {
                return Some(Fate::Spiral {
                    returns: entry.1,
                    max_norm: max_norm.max(norm),
                });
            }
        }
    }
    None
}

#[inline]
fn signum<const N: usize>(v: &[i64; N]) -> [i8; N] {
    std::array::from_fn(|e| v[e].signum() as i8)
}

/// Analyse the trajectory starting at `start` (step 1 of the original model).
pub fn analyze<const N: usize>(w: &Reduced<N>, start: [i64; N], budget: Budget) -> Fate {
    analyze_from(w, start, 1, budget)
}

/// Analyse the trajectory of the original model (all channels fire first).
pub fn analyze_universe<const N: usize>(w: &Reduced<N>, budget: Budget) -> Fate {
    analyze(w, w.initial_state(), budget)
}

pub fn analyze_from<const N: usize>(
    w: &Reduced<N>,
    start: [i64; N],
    start_step: u64,
    budget: Budget,
) -> Fate {
    analyze_traced(w, start, start_step, budget, &mut |_, _, _| {})
}

/// Like [`analyze_from`], calling `trace(step, state, drift)` at every event.
pub fn analyze_traced<const N: usize>(
    w: &Reduced<N>,
    start: [i64; N],
    start_step: u64,
    budget: Budget,
    trace: Trace<N>,
) -> Fate {
    let mut v = start;
    let mut step = start_step;
    let mut events = 0u64;
    let mut max_norm = 0i64;
    // Event log: (state, step, sign pattern). Short trajectories dominate, so
    // cycle detection scans the log until it is long enough to justify a map.
    let mut log: Vec<([i64; N], u64, [i8; N])> = Vec::new();
    let mut seen: HashMap<[i64; N], u64> = HashMap::new();
    // Per-pattern Poincaré sections: (last norm, consecutive growing returns).
    let mut returns: HashMap<[i8; N], (i64, u64)> = HashMap::new();
    loop {
        let s = signum(&v);
        let d = w.drift(&s);
        trace(step, &v, &d);
        if d.iter().all(|&x| x == 0) {
            return Fate::Fixed { step };
        }
        // Absorbing orthant: every coordinate keeps (or gains nothing against) its sign.
        let absorbing = (0..N).all(|e| {
            if s[e] == 0 {
                d[e] == 0
            } else {
                d[e] == 0 || d[e].signum() as i8 == s[e]
            }
        });
        if absorbing {
            return Fate::Ray { step };
        }
        let first = if log.len() < 64 {
            log.iter().find(|(u, _, _)| *u == v).map(|&(_, t, _)| t)
        } else {
            seen.get(&v).copied()
        };
        if let Some(first) = first {
            return Fate::Cycle {
                transient: first,
                period: step - first,
            };
        }
        if let Some((ustep, norm)) = helix_shift(&log, &v, &s, |_, shift, sign| {
            sign != 0 && shift.signum() as i8 == sign
        }) {
            return Fate::Helix {
                step: ustep,
                period: step - ustep,
                shift_norm: norm,
            };
        }
        if log.len() == 64 {
            seen.extend(log.iter().map(|&(u, t, _)| (u, t)));
        }
        if log.len() >= 64 {
            seen.insert(v, step);
        }
        log.push((v, step, s));
        if budget.spiral_returns > 0 {
            if let Some(f) = spiral_return(&mut returns, &s, &v, budget.spiral_returns, max_norm) {
                return f;
            }
        }
        // Steps until some coordinate changes sign (or leaves zero).
        let mut k = u64::MAX;
        for e in 0..N {
            if s[e] == 0 {
                if d[e] != 0 {
                    k = 1;
                }
            } else if d[e].signum() as i8 == -s[e] {
                k = k.min((v[e].unsigned_abs()).div_ceil(d[e].unsigned_abs()));
            }
        }
        debug_assert!(k != u64::MAX && k >= 1);
        for e in 0..N {
            v[e] += k as i64 * d[e];
            max_norm = max_norm.max(v[e].abs());
        }
        step += k;
        events += 1;
        if events >= budget.max_events || max_norm >= budget.max_norm {
            return Fate::Undecided { events, max_norm };
        }
    }
}

/// The "create / destroy with positivity" variant: one channel per element,
/// active while `v_e > 0`; a firing row adds `rows[f][e] ∈ {-1,0,1}` and the
/// result is clamped at zero (`v_e = max(0, v_e + drift)`). Firing all rows
/// once gives the start state, as in the original model.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, PartialOrd, Ord)]
pub struct Clamped<const N: usize> {
    pub rows: [[i8; N]; N],
}

impl<const N: usize> Clamped<N> {
    /// Number of distinct clamped universes: `3^(N²)`.
    pub fn count() -> u128 {
        3u128.pow((N * N) as u32)
    }

    pub fn from_index(mut index: u64) -> Self {
        let mut c = Clamped { rows: [[0; N]; N] };
        for f in (0..N).rev() {
            for e in (0..N).rev() {
                c.rows[f][e] = (index % 3) as i8 - 1;
                index /= 3;
            }
        }
        c
    }

    /// Bit matrices (create bit, destroy bit per entry) reducing to this: `2^(zeros)`.
    pub fn multiplicity(&self) -> u64 {
        1u64 << self.rows.iter().flatten().filter(|&&w| w == 0).count()
    }

    pub fn initial_state(&self) -> [i64; N] {
        let mut v = [0i64; N];
        for f in 0..N {
            for e in 0..N {
                v[e] += self.rows[f][e] as i64;
            }
        }
        v.map(|x| x.max(0))
    }

    #[inline]
    pub fn drift(&self, v: &[i64; N]) -> [i64; N] {
        let mut d = [0i64; N];
        for f in 0..N {
            if v[f] > 0 {
                for e in 0..N {
                    d[e] += self.rows[f][e] as i64;
                }
            }
        }
        d
    }
}

pub fn analyze_clamped<const N: usize>(c: &Clamped<N>, start: [i64; N], budget: Budget) -> Fate {
    analyze_clamped_traced(c, start, budget, &mut |_, _, _| {})
}

pub fn analyze_clamped_traced<const N: usize>(
    c: &Clamped<N>,
    start: [i64; N],
    budget: Budget,
    trace: Trace<N>,
) -> Fate {
    let mut v = start;
    let mut step = 1u64;
    let mut events = 0u64;
    let mut max_norm = 0i64;
    let mut log: Vec<([i64; N], u64, [i8; N])> = Vec::new();
    let mut seen: HashMap<[i64; N], u64> = HashMap::new();
    let mut returns: HashMap<[i8; N], (i64, u64)> = HashMap::new();
    loop {
        let d = c.drift(&v);
        trace(step, &v, &d);
        let support: [i8; N] = std::array::from_fn(|e| (v[e] > 0) as i8);
        // Effective change: clamped coordinates cannot go below zero.
        let effective: [i64; N] =
            std::array::from_fn(|e| if support[e] == 1 { d[e] } else { d[e].max(0) });
        if effective.iter().all(|&x| x == 0) {
            return Fate::Fixed { step };
        }
        if (0..N).all(|e| {
            if support[e] == 1 {
                d[e] >= 0
            } else {
                d[e] <= 0
            }
        }) {
            return Fate::Ray { step };
        }
        let first = if log.len() < 64 {
            log.iter().find(|(u, _, _)| *u == v).map(|&(_, t, _)| t)
        } else {
            seen.get(&v).copied()
        };
        if let Some(first) = first {
            return Fate::Cycle {
                transient: first,
                period: step - first,
            };
        }
        if let Some((ustep, norm)) =
            helix_shift(&log, &v, &support, |_, shift, sup| sup == 1 && shift > 0)
        {
            return Fate::Helix {
                step: ustep,
                period: step - ustep,
                shift_norm: norm,
            };
        }
        if log.len() == 64 {
            seen.extend(log.iter().map(|&(u, t, _)| (u, t)));
        }
        if log.len() >= 64 {
            seen.insert(v, step);
        }
        log.push((v, step, support));
        if budget.spiral_returns > 0 {
            if let Some(f) =
                spiral_return(&mut returns, &support, &v, budget.spiral_returns, max_norm)
            {
                return f;
            }
        }
        let mut k = u64::MAX;
        for e in 0..N {
            if support[e] == 0 {
                if d[e] > 0 {
                    k = 1;
                }
            } else if d[e] < 0 {
                k = k.min(v[e].unsigned_abs().div_ceil(d[e].unsigned_abs()));
            }
        }
        debug_assert!(k != u64::MAX && k >= 1);
        for e in 0..N {
            v[e] = (v[e] + k as i64 * effective[e]).max(0);
            max_norm = max_norm.max(v[e]);
        }
        step += k;
        events += 1;
        if events >= budget.max_events || max_norm >= budget.max_norm {
            return Fate::Undecided { events, max_norm };
        }
    }
}
