//! Research experiments over reduced universes.

use std::collections::BTreeMap;

use crate::exact::{
    analyze_clamped, analyze_from, analyze_universe, Budget, Clamped, Fate, Reduced,
};
use crate::parallel::par_fold;
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
    /// Steps at which fixed points / rays are certified -> reduced count.
    pub certified_at: BTreeMap<u64, u64>,
    /// Largest transient+period seen for a cycle, with its index.
    pub longest_cycle: (u64, u64),
    /// Undecided example indices (first few).
    pub undecided_examples: Vec<u64>,
    /// Spiral example indices (first few).
    pub spiral_examples: Vec<u64>,
}

impl Census {
    pub fn heuristic_available(&self) -> bool {
        !self.confusion.is_empty()
    }

    fn merge(&mut self, o: Census) {
        self.reduced += o.reduced;
        self.matrices += o.matrices;
        for (k, (a, b)) in o.fates {
            let e = self.fates.entry(k).or_insert((0, 0));
            e.0 += a;
            e.1 += b;
        }
        for (k, n) in o.confusion {
            *self.confusion.entry(k).or_insert(0) += n;
        }
        for (k, n) in o.periods {
            *self.periods.entry(k).or_insert(0) += n;
        }
        for (k, n) in o.certified_at {
            *self.certified_at.entry(k).or_insert(0) += n;
        }
        if o.longest_cycle.0 > self.longest_cycle.0 {
            self.longest_cycle = o.longest_cycle;
        }
        self.undecided_examples.extend(o.undecided_examples);
        self.undecided_examples.truncate(10);
        self.spiral_examples.extend(o.spiral_examples);
        self.spiral_examples.truncate(10);
    }

    fn record<const N: usize>(&mut self, index: u64, w: &Reduced<N>, fate: Fate, heuristic: bool) {
        self.record_with(
            index,
            w.multiplicity() as u128,
            fate,
            if heuristic { Some(w.to_rows()) } else { None },
        );
    }

    fn record_with(
        &mut self,
        index: u64,
        mult: u128,
        fate: Fate,
        heuristic_rows: Option<Vec<u16>>,
    ) {
        self.reduced += 1;
        self.matrices += mult;
        let e = self.fates.entry(fate.kind()).or_insert((0, 0));
        e.0 += 1;
        e.1 += mult;
        match fate {
            Fate::Cycle { transient, period } => {
                *self.periods.entry(period).or_insert(0) += 1;
                if transient + period > self.longest_cycle.0 {
                    self.longest_cycle = (transient + period, index);
                }
            }
            Fate::Fixed { step } | Fate::Ray { step } => {
                *self.certified_at.entry(step).or_insert(0) += 1
            }
            Fate::Helix { step, period, .. } => {
                *self.certified_at.entry(step + period).or_insert(0) += 1
            }
            Fate::Undecided { .. } => {
                if self.undecided_examples.len() < 10 {
                    self.undecided_examples.push(index);
                }
            }
            Fate::Spiral { .. } => {
                if self.spiral_examples.len() < 10 {
                    self.spiral_examples.push(index);
                }
            }
        }
        if let Some(rows) = heuristic_rows {
            let verdict = with_channels_rows(&rows);
            *self.confusion.entry((verdict, fate.kind())).or_insert(0) += mult;
        }
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
    par_fold(
        count as u64,
        threads,
        Census::default,
        |k, c| {
            let w = Reduced::<N>::from_index(k);
            c.record(k, &w, analyze_universe(&w, budget), heuristic);
        },
        Census::merge,
    )
}

/// Exact fates of `samples` random reduced universes (uniform over ternary entries).
pub fn census_sampled<const N: usize>(
    samples: u64,
    seed: u64,
    budget: Budget,
    threads: usize,
    heuristic: bool,
) -> Census {
    par_fold(
        samples,
        threads,
        Census::default,
        |k, c| {
            let w = random_reduced::<N>(&mut SplitMix64::new(
                seed ^ k.wrapping_mul(0x9E37_79B9_7F4A_7C15),
            ));
            // Record the universe's own index when it fits, so examples can be replayed with `show`.
            let index = if 2 * N * N <= 40 { w.index() } else { k };
            c.record(index, &w, analyze_universe(&w, budget), heuristic);
        },
        Census::merge,
    )
}

/// Uniformly random ternary drift rows.
pub fn random_reduced<const N: usize>(rng: &mut SplitMix64) -> Reduced<N> {
    let mut w = Reduced::<N> {
        pos: [[0; N]; N],
        neg: [[0; N]; N],
    };
    for row in w.pos.iter_mut().chain(w.neg.iter_mut()) {
        for x in row.iter_mut() {
            *x = (rng.next_u64() % 3) as i8 - 1;
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
    par_fold(
        count as u64,
        threads,
        Census::default,
        |k, c| {
            let u = Clamped::<N>::from_index(k);
            c.record_with(
                k,
                u.multiplicity() as u128,
                analyze_clamped(&u, u.initial_state(), budget),
                None,
            );
        },
        Census::merge,
    )
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
                let kind = analyze_from(&w, v, 1, budget).kind();
                *b.by_start.entry(kind).or_insert(0) += 1;
                if !kinds.contains(&kind) {
                    kinds.push(kind);
                }
            }
            kinds.sort();
            b.reduced += 1;
            *b.kinds_per_universe.entry(kinds.len()).or_insert(0) += 1;
            *b.kind_sets.entry(kinds.join("+")).or_insert(0) += 1;
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
        },
    )
}
