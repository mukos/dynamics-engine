//! Research experiments over reduced universes.

use std::collections::BTreeMap;

use crate::exact::{
    analyze_clamped_ws, analyze_from, analyze_universe_ws, Budget, Clamped, Fate, Reduced,
    Workspace,
};
use crate::parallel::{par_chunks, par_fold};
use crate::rng::SplitMix64;
use crate::universe::{Rules, Universe};

/// Coarse verdict of the original 10-step heuristic.
fn heuristic_verdict(signature: &str) -> &'static str {
    if signature.split(',').all(|s| s == "|") {
        "terminated"
    } else if signature.split(',').any(|s| s == "1" || s == "-1") {
        "growth"
    } else {
        "bounded"
    }
}

/// Aggregated counts of one experiment.
#[derive(Default, Clone, Debug)]
pub struct Census {
    /// Reduced universes analysed.
    pub reduced: u64,
    /// Bit matrices represented (sum of multiplicities).
    pub matrices: u128,
    /// fate kind -> (reduced count, matrix-weighted count).
    pub fates: BTreeMap<&'static str, (u64, u128)>,
    /// (heuristic verdict, fate kind) -> matrix-weighted count.
    pub confusion: BTreeMap<(&'static str, &'static str), u128>,
    /// Cycle period -> reduced count.
    pub periods: BTreeMap<u64, u64>,
    /// Steps at which fixed points / rays / helices are certified -> reduced count.
    pub certified_at: BTreeMap<u64, u64>,
    /// Largest transient+period seen for a cycle, with its index.
    pub longest_cycle: (u64, u64),
    /// Undecided example indices (first few).
    pub undecided_examples: Vec<u64>,
    /// Spiral example indices (first few).
    pub spiral_examples: Vec<u64>,
    /// Fine classes: (kind code, parameter) -> (count, smallest index). Kind codes: 0 fixed@step,
    /// 1 cycle@period, 2 ray@step, 3 helix@period, 4 spiral, 5 undecided.
    pub classes: BTreeMap<(u8, u64), (u64, u64)>,
}

pub fn class_of(fate: &Fate) -> (u8, u64) {
    match *fate {
        Fate::Fixed { step } => (0, step),
        Fate::Cycle { period, .. } => (1, period),
        Fate::Ray { step } => (2, step),
        Fate::Helix { period, .. } => (3, period),
        Fate::Spiral { .. } => (4, 0),
        Fate::Undecided { .. } => (5, 0),
    }
}

pub const KINDS: [&str; 6] = ["fixed", "cycle", "ray", "helix", "spiral", "undecided"];

pub fn class_name(kind: u8) -> &'static str {
    KINDS[kind as usize]
}

impl Census {
    pub fn heuristic_available(&self) -> bool {
        !self.confusion.is_empty()
    }
}

/// Parameters (steps, periods) below this are tallied in flat arrays; larger ones in a map.
const DENSE: usize = 4096;

/// Per-thread accumulator: flat arrays on the hot path, maps only for rare overflow.
struct Acc {
    reduced: u64,
    matrices: u128,
    fates: [(u64, u128); 6],
    classes: Vec<(u64, u64)>,
    classes_overflow: BTreeMap<(u8, u64), (u64, u64)>,
    certified: Vec<u64>,
    certified_overflow: BTreeMap<u64, u64>,
    longest_cycle: (u64, u64),
    undecided_examples: Vec<u64>,
    spiral_examples: Vec<u64>,
    confusion: BTreeMap<(&'static str, &'static str), u128>,
}

impl Acc {
    fn new() -> Self {
        Acc {
            reduced: 0,
            matrices: 0,
            fates: [(0, 0); 6],
            classes: vec![(0, u64::MAX); 6 * DENSE],
            classes_overflow: BTreeMap::new(),
            certified: vec![0; DENSE],
            certified_overflow: BTreeMap::new(),
            longest_cycle: (0, 0),
            undecided_examples: Vec::new(),
            spiral_examples: Vec::new(),
            confusion: BTreeMap::new(),
        }
    }

    fn merge(&mut self, o: Acc) {
        self.reduced += o.reduced;
        self.matrices += o.matrices;
        for k in 0..6 {
            self.fates[k].0 += o.fates[k].0;
            self.fates[k].1 += o.fates[k].1;
        }
        for (a, b) in self.classes.iter_mut().zip(o.classes) {
            a.0 += b.0;
            a.1 = a.1.min(b.1);
        }
        for (k, (n, idx)) in o.classes_overflow {
            let e = self.classes_overflow.entry(k).or_insert((0, u64::MAX));
            e.0 += n;
            e.1 = e.1.min(idx);
        }
        for (a, b) in self.certified.iter_mut().zip(o.certified) {
            *a += b;
        }
        for (k, n) in o.certified_overflow {
            *self.certified_overflow.entry(k).or_insert(0) += n;
        }
        if o.longest_cycle.0 > self.longest_cycle.0 {
            self.longest_cycle = o.longest_cycle;
        }
        self.undecided_examples.extend(o.undecided_examples);
        self.undecided_examples.truncate(10);
        self.spiral_examples.extend(o.spiral_examples);
        self.spiral_examples.truncate(10);
        for (k, n) in o.confusion {
            *self.confusion.entry(k).or_insert(0) += n;
        }
    }

    #[inline]
    fn record<const N: usize>(&mut self, index: u64, w: &Reduced<N>, fate: Fate, heuristic: bool) {
        self.record_with(
            index,
            w.multiplicity() as u128,
            fate,
            if heuristic { Some(w.to_rows()) } else { None },
        );
    }

    #[inline]
    fn record_with(
        &mut self,
        index: u64,
        mult: u128,
        fate: Fate,
        heuristic_rows: Option<Vec<u16>>,
    ) {
        self.reduced += 1;
        self.matrices += mult;
        let (kind, param) = class_of(&fate);
        let f = &mut self.fates[kind as usize];
        f.0 += 1;
        f.1 += mult;
        if (param as usize) < DENSE {
            let e = &mut self.classes[kind as usize * DENSE + param as usize];
            e.0 += 1;
            e.1 = e.1.min(index);
        } else {
            let e = self
                .classes_overflow
                .entry((kind, param))
                .or_insert((0, u64::MAX));
            e.0 += 1;
            e.1 = e.1.min(index);
        }
        let certified = match fate {
            Fate::Fixed { step } | Fate::Ray { step } => Some(step),
            Fate::Helix { step, period, .. } => Some(step + period),
            Fate::Cycle { transient, period } => {
                if transient + period > self.longest_cycle.0 {
                    self.longest_cycle = (transient + period, index);
                }
                None
            }
            Fate::Undecided { .. } => {
                if self.undecided_examples.len() < 10 {
                    self.undecided_examples.push(index);
                }
                None
            }
            Fate::Spiral { .. } => {
                if self.spiral_examples.len() < 10 {
                    self.spiral_examples.push(index);
                }
                None
            }
        };
        if let Some(step) = certified {
            if (step as usize) < DENSE {
                self.certified[step as usize] += 1;
            } else {
                *self.certified_overflow.entry(step).or_insert(0) += 1;
            }
        }
        if let Some(rows) = heuristic_rows {
            let verdict = with_channels_rows(&rows);
            *self.confusion.entry((verdict, fate.kind())).or_insert(0) += mult;
        }
    }
}

impl From<Acc> for Census {
    fn from(a: Acc) -> Census {
        let mut c = Census {
            reduced: a.reduced,
            matrices: a.matrices,
            confusion: a.confusion,
            longest_cycle: a.longest_cycle,
            undecided_examples: a.undecided_examples,
            spiral_examples: a.spiral_examples,
            ..Census::default()
        };
        for (k, &(r, m)) in a.fates.iter().enumerate() {
            if r > 0 {
                c.fates.insert(KINDS[k], (r, m));
            }
        }
        for (i, &(n, idx)) in a.classes.iter().enumerate() {
            if n > 0 {
                c.classes
                    .insert(((i / DENSE) as u8, (i % DENSE) as u64), (n, idx));
            }
        }
        for (k, v) in a.classes_overflow {
            c.classes.insert(k, v);
        }
        for (&(kind, param), &(n, _)) in &c.classes {
            if kind == 1 {
                c.periods.insert(param, n);
            }
        }
        for (step, &n) in a.certified.iter().enumerate() {
            if n > 0 {
                c.certified_at.insert(step as u64, n);
            }
        }
        for (k, v) in a.certified_overflow {
            c.certified_at.insert(k, v);
        }
        c
    }
}

/// Run the original 10-step heuristic on representative rows and return its coarse verdict.
fn with_channels_rows(rows: &[u16]) -> &'static str {
    macro_rules! go {
        ($c:literal) => {{
            let mut arr = [0u16; $c];
            arr.copy_from_slice(rows);
            let mut u = Universe::<$c>::new(arr, Universe::<$c>::ALL, Rules::default());
            u.run();
            heuristic_verdict(&u.signature())
        }};
    }
    match rows.len() {
        2 => go!(2),
        4 => go!(4),
        6 => go!(6),
        8 => go!(8),
        10 => go!(10),
        12 => go!(12),
        14 => go!(14),
        16 => go!(16),
        _ => unreachable!(),
    }
}

/// Exact fates of every reduced universe with `N` elements (original start: all channels fire).
pub fn census_exhaustive<const N: usize>(
    budget: Budget,
    threads: usize,
    heuristic: bool,
) -> Census {
    let count = Reduced::<N>::count();
    assert!(
        count <= u64::MAX as u128,
        "too many reduced universes to enumerate"
    );
    par_chunks(
        count as u64,
        threads,
        Acc::new,
        |first, len, acc| {
            let mut ws = Workspace::<N>::new();
            let mut w = Reduced::<N>::from_index(first);
            for k in first..first + len {
                acc.record(k, &w, analyze_universe_ws(&w, budget, &mut ws), heuristic);
                w.increment();
            }
        },
        Acc::merge,
    )
    .into()
}

/// Exact fates of `samples` random reduced universes (uniform over ternary entries).
pub fn census_sampled<const N: usize>(
    samples: u64,
    seed: u64,
    budget: Budget,
    threads: usize,
    heuristic: bool,
) -> Census {
    par_chunks(
        samples,
        threads,
        Acc::new,
        |first, len, acc| {
            let mut ws = Workspace::<N>::new();
            let mut rng = SplitMix64::new(seed ^ first.wrapping_mul(0x9E37_79B9_7F4A_7C15));
            for k in first..first + len {
                let w = random_reduced::<N>(&mut rng);
                // Record the universe's own index when it fits, so examples can be replayed with `show`.
                let index = if 2 * N * N <= 40 { w.index() } else { k };
                acc.record(
                    index,
                    &w,
                    analyze_universe_ws(&w, budget, &mut ws),
                    heuristic,
                );
            }
        },
        Acc::merge,
    )
    .into()
}

/// Uniformly random ternary drift rows (two random bits per digit, rejecting `11`).
pub fn random_reduced<const N: usize>(rng: &mut SplitMix64) -> Reduced<N> {
    let mut w = Reduced::<N> {
        pos: [[0; N]; N],
        neg: [[0; N]; N],
    };
    let mut bits = rng.next_u64();
    let mut avail = 32;
    for row in w.pos.iter_mut().chain(w.neg.iter_mut()) {
        for x in row.iter_mut() {
            loop {
                if avail == 0 {
                    bits = rng.next_u64();
                    avail = 32;
                }
                let d = (bits & 3) as i8;
                bits >>= 2;
                avail -= 1;
                if d < 3 {
                    *x = d - 1;
                    break;
                }
            }
        }
    }
    w
}

/// Exact fates of every clamped (create/destroy) universe with `N` elements.
pub fn census_clamped<const N: usize>(budget: Budget, threads: usize) -> Census {
    let count = Clamped::<N>::count();
    assert!(
        count <= u64::MAX as u128,
        "too many clamped universes to enumerate"
    );
    par_chunks(
        count as u64,
        threads,
        Acc::new,
        |first, len, c| {
            let mut ws = Workspace::<N>::new();
            for k in first..first + len {
                let u = Clamped::<N>::from_index(k);
                let fate =
                    analyze_clamped_ws(&u, u.initial_state(), budget, &mut ws, &mut |_, _, _| {});
                c.record_with(k, u.multiplicity() as u128, fate, None);
            }
        },
        Acc::merge,
    )
    .into()
}

/// Does the fate depend on the initial state? For every reduced universe, analyse
/// all starts in `[-radius, radius]^N` and count how many distinct fate kinds occur.
#[derive(Default, Clone, Debug)]
pub struct Basins {
    pub reduced: u64,
    /// number of distinct fate kinds over the start box -> reduced count
    pub kinds_per_universe: BTreeMap<usize, u64>,
    /// set of fate kinds (sorted, joined) -> reduced count
    pub kind_sets: BTreeMap<String, u64>,
    /// fate kind -> number of (universe, start) pairs
    pub by_start: BTreeMap<&'static str, u64>,
    /// fate-kind set -> one representative universe (smallest index) with one start per kind
    pub representatives: BTreeMap<String, Representative>,
}

/// A universe standing for a fate-kind set, with one start state per fate it shows.
#[derive(Clone, Debug)]
pub struct Representative {
    pub index: u64,
    pub pos: Vec<Vec<i8>>,
    pub neg: Vec<Vec<i8>>,
    /// (start state, fate kind, fate debug text)
    pub starts: Vec<(Vec<i64>, &'static str, String)>,
}

/// `samples = None` enumerates every reduced universe; `Some((n, seed))` draws `n` at random.
pub fn basins<const N: usize>(
    radius: i64,
    budget: Budget,
    threads: usize,
    samples: Option<(u64, u64)>,
) -> Basins {
    let count = samples
        .map(|(n, _)| n)
        .unwrap_or(Reduced::<N>::count() as u64);
    let side = (2 * radius + 1) as u64;
    let starts = side.pow(N as u32);
    par_fold(
        count,
        threads,
        Basins::default,
        |k, b| {
            let w = match samples {
                None => Reduced::<N>::from_index(k),
                Some((_, seed)) => random_reduced::<N>(&mut SplitMix64::new(
                    seed ^ k.wrapping_mul(0x9E37_79B9_7F4A_7C15),
                )),
            };
            let mut kinds: Vec<&'static str> = Vec::new();
            let mut examples: Vec<(Vec<i64>, &'static str, String)> = Vec::new();
            for s in 0..starts {
                let mut v = [0i64; N];
                let mut rem = s;
                for x in v.iter_mut() {
                    *x = (rem % side) as i64 - radius;
                    rem /= side;
                }
                if v.iter().all(|&x| x == 0) {
                    continue; // the origin is trivially fixed in every universe
                }
                let fate = analyze_from(&w, v, 1, budget);
                let kind = fate.kind();
                *b.by_start.entry(kind).or_insert(0) += 1;
                if !kinds.contains(&kind) {
                    kinds.push(kind);
                }
                // Keep up to two example starts per fate so single-fate universes still show variety.
                if examples.iter().filter(|(_, k, _)| *k == kind).count() < 2 {
                    examples.push((v.to_vec(), kind, format!("{fate:?}")));
                }
            }
            kinds.sort();
            b.reduced += 1;
            *b.kinds_per_universe.entry(kinds.len()).or_insert(0) += 1;
            let key = kinds.join("+");
            *b.kind_sets.entry(key.clone()).or_insert(0) += 1;
            let index = w.index();
            let better = b.representatives.get(&key).is_none_or(|r| index < r.index);
            if better {
                b.representatives.insert(
                    key,
                    Representative {
                        index,
                        pos: w.pos.iter().map(|r| r.to_vec()).collect(),
                        neg: w.neg.iter().map(|r| r.to_vec()).collect(),
                        starts: examples,
                    },
                );
            }
        },
        |a, o| {
            a.reduced += o.reduced;
            for (k, n) in o.kinds_per_universe {
                *a.kinds_per_universe.entry(k).or_insert(0) += n;
            }
            for (k, n) in o.kind_sets {
                *a.kind_sets.entry(k).or_insert(0) += n;
            }
            for (k, n) in o.by_start {
                *a.by_start.entry(k).or_insert(0) += n;
            }
            for (k, r) in o.representatives {
                let better = a
                    .representatives
                    .get(&k)
                    .is_none_or(|cur| r.index < cur.index);
                if better {
                    a.representatives.insert(k, r);
                }
            }
        },
    )
}

// ---------------------------------------------------------------------------
// Basin maps: which attractor each start state reaches.

use crate::exact::analyze_ws;

/// One cell of a basin map: fate kind code and attractor identity.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Attractor {
    pub kind: u8,
    pub id: Vec<i64>,
}

/// Basin map of one universe over the box `[-radius, radius]^N` (row-major,
/// last coordinate fastest), with the list of distinct attractors and
/// pattern features.
#[derive(Clone, Debug)]
pub struct BasinMap {
    pub index: u64,
    pub radius: i64,
    pub attractors: Vec<Attractor>,
    /// Attractor index per cell (into `attractors`); origin gets the fixed attractor at 0.
    pub cells: Vec<u16>,
    pub features: BasinFeatures,
}

/// Pattern statistics of a basin map.
#[derive(Clone, Debug, Default)]
pub struct BasinFeatures {
    /// Distinct attractors reached from the box.
    pub attractors: usize,
    /// Attractors that are not fixed points (cycles, rays, helices, spirals).
    pub moving: usize,
    /// Fixed points reached from at least one other cell (a real basin, not an inert cell).
    pub sinks: usize,
    /// Cells that are fixed and reached from nowhere else.
    pub inert: usize,
    /// Distinct fate kinds.
    pub kinds: usize,
    /// Fraction of axis-neighbouring cell pairs that reach different attractors.
    pub boundary: f64,
    /// Fraction of directions whose fate is the same at radius r and 2r (1 = conic basins).
    pub conic: f64,
    /// Fraction of cells that differ from a neighbour but agree with the cell two steps on
    /// (alternating stripes: a lattice / residue effect).
    pub stripes: f64,
    /// Fraction of cells per kind code (fixed, cycle, ray, helix, spiral, undecided).
    pub share: [f64; 6],
}

fn cell_coords<const N: usize>(mut k: usize, radius: i64) -> [i64; N] {
    let side = (2 * radius + 1) as usize;
    let mut v = [0i64; N];
    for e in (0..N).rev() {
        v[e] = (k % side) as i64 - radius;
        k /= side;
    }
    v
}

fn cell_index<const N: usize>(v: &[i64; N], radius: i64) -> Option<usize> {
    let side = (2 * radius + 1) as usize;
    let mut k = 0usize;
    for &x in v {
        if x < -radius || x > radius {
            return None;
        }
        k = k * side + (x + radius) as usize;
    }
    Some(k)
}

/// Compute the basin map of `w` over `[-radius, radius]^N`.
pub fn basin_map<const N: usize>(w: &Reduced<N>, radius: i64, budget: Budget) -> BasinMap {
    let side = (2 * radius + 1) as usize;
    let total = side.pow(N as u32);
    let mut ws = Workspace::<N>::new();
    let mut attractors: Vec<Attractor> = Vec::new();
    let mut lookup: BTreeMap<Attractor, u16> = BTreeMap::new();
    let mut cells = vec![0u16; total];
    for (k, cell) in cells.iter_mut().enumerate() {
        let v = cell_coords::<N>(k, radius);
        let fate = analyze_ws(w, v, 1, budget, &mut ws, &mut |_, _, _| {});
        let (kind, _) = class_of(&fate);
        let a = Attractor {
            kind,
            id: ws.attractor.clone(),
        };
        let id = match lookup.get(&a) {
            Some(&i) => i,
            None => {
                let i = attractors.len() as u16;
                attractors.push(a.clone());
                lookup.insert(a, i);
                i
            }
        };
        *cell = id;
    }
    let features = basin_features::<N>(&attractors, &cells, radius);
    BasinMap {
        index: w.index(),
        radius,
        attractors,
        cells,
        features,
    }
}

fn basin_features<const N: usize>(
    attractors: &[Attractor],
    cells: &[u16],
    radius: i64,
) -> BasinFeatures {
    let mut f = BasinFeatures {
        attractors: attractors.len(),
        kinds: {
            let mut ks: Vec<u8> = attractors.iter().map(|a| a.kind).collect();
            ks.sort();
            ks.dedup();
            ks.len()
        },
        ..Default::default()
    };
    let total = cells.len() as f64;
    let mut basin_size = vec![0usize; attractors.len()];
    for &c in cells {
        f.share[attractors[c as usize].kind as usize] += 1.0 / total;
        basin_size[c as usize] += 1;
    }
    // Merge every fixed point that captures only itself into one "inert" class,
    // so a zero-drift universe does not read as hundreds of attractors.
    let inert_class = attractors.len() as u16;
    let class: Vec<u16> = (0..attractors.len())
        .map(|i| {
            if attractors[i].kind == 0 && basin_size[i] == 1 {
                inert_class
            } else {
                i as u16
            }
        })
        .collect();
    f.moving = attractors.iter().filter(|a| a.kind != 0).count();
    f.sinks = (0..attractors.len())
        .filter(|&i| attractors[i].kind == 0 && basin_size[i] > 1)
        .count();
    f.inert = (0..attractors.len())
        .filter(|&i| attractors[i].kind == 0 && basin_size[i] == 1)
        .count();
    let cls = |k: usize| class[cells[k] as usize];
    let (mut pairs, mut diff, mut stripes) = (0u64, 0u64, 0u64);
    let (mut dirs, mut same) = (0u64, 0u64);
    for k in 0..cells.len() {
        let v = cell_coords::<N>(k, radius);
        for e in 0..N {
            let mut u = v;
            u[e] += 1;
            if let Some(j) = cell_index(&u, radius) {
                pairs += 1;
                if cls(j) != cls(k) {
                    diff += 1;
                    let mut u2 = v;
                    u2[e] += 2;
                    if let Some(j2) = cell_index(&u2, radius) {
                        if cls(j2) == cls(k) {
                            stripes += 1;
                        }
                    }
                }
            }
        }
        // Conic test: compare with the doubled point for cells in the middle ring.
        let norm = v.iter().map(|x| x.abs()).max().unwrap_or(0);
        if norm > 0 && norm * 2 <= radius && norm * 2 > radius / 2 {
            let doubled: [i64; N] = std::array::from_fn(|e| 2 * v[e]);
            if let Some(j) = cell_index(&doubled, radius) {
                dirs += 1;
                if cls(j) == cls(k) {
                    same += 1;
                }
            }
        }
    }
    f.boundary = diff as f64 / pairs.max(1) as f64;
    f.stripes = stripes as f64 / pairs.max(1) as f64;
    f.conic = same as f64 / dirs.max(1) as f64;
    f
}

/// Basin maps of every canonical universe (one per symmetry orbit), with features.
/// `samples = None` scans every canonical universe; `Some((n, seed))` draws `n` random
/// universes and maps their canonical representatives (duplicates removed).
pub fn basin_scan<const N: usize>(
    radius: i64,
    budget: Budget,
    threads: usize,
    samples: Option<(u64, u64)>,
) -> Vec<BasinMap> {
    let count = samples
        .map(|(n, _)| n)
        .unwrap_or(Reduced::<N>::count() as u64);
    let mut maps: Vec<BasinMap> = par_chunks(
        count,
        threads,
        Vec::new,
        |first, len, acc: &mut Vec<BasinMap>| match samples {
            None => {
                let mut w = Reduced::<N>::from_index(first);
                for _ in first..first + len {
                    if w.canonical() == w {
                        acc.push(basin_map(&w, radius, budget));
                    }
                    w.increment();
                }
            }
            Some((_, seed)) => {
                let mut rng = SplitMix64::new(seed ^ first.wrapping_mul(0x9E37_79B9_7F4A_7C15));
                for _ in 0..len {
                    let w = random_reduced::<N>(&mut rng).canonical();
                    acc.push(basin_map(&w, radius, budget));
                }
            }
        },
        |a, b| a.extend(b),
    );
    maps.sort_by_key(|m| m.index);
    maps.dedup_by_key(|m| m.index);
    maps
}
