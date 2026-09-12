//! Run many universes and tally their signatures.

use std::collections::{BTreeMap, HashMap};
use std::sync::atomic::{AtomicU64, Ordering};
use std::thread;
use std::time::{Duration, Instant};

use crate::matrix::{self, enumerable, matrix_count};
use crate::rng::SplitMix64;
use crate::universe::{decode_signature, Rules, Universe};

/// Signature text -> number of matrices that produced it.
pub type Results = BTreeMap<String, u64>;

/// Per-thread histogram over packed signature codes. Dense when the code space
/// is small (≤ 8 channels → 65 536 codes), hashed otherwise.
enum Tally {
    Dense(Vec<u64>),
    Sparse(HashMap<u32, u64>),
}

impl Tally {
    fn new<const C: usize>() -> Self {
        if C <= 8 {
            Tally::Dense(vec![0; 1 << (2 * C)])
        } else {
            Tally::Sparse(HashMap::new())
        }
    }

    #[inline]
    fn add(&mut self, code: u32) {
        match self {
            Tally::Dense(v) => v[code as usize] += 1,
            Tally::Sparse(m) => *m.entry(code).or_insert(0) += 1,
        }
    }

    fn merge_into<const C: usize>(self, results: &mut Results) {
        let entries: Vec<(u32, u64)> = match self {
            Tally::Dense(v) => v
                .into_iter()
                .enumerate()
                .filter(|&(_, n)| n > 0)
                .map(|(c, n)| (c as u32, n))
                .collect(),
            Tally::Sparse(m) => m.into_iter().collect(),
        };
        for (code, n) in entries {
            *results.entry(decode_signature::<C>(code)).or_insert(0) += n;
        }
    }
}

/// How many worker threads to use and how often to report progress.
#[derive(Clone, Copy, Debug)]
pub struct Execution {
    pub threads: usize,
    /// Print progress to stderr at this interval; `None` disables it.
    pub progress: Option<Duration>,
}

impl Default for Execution {
    fn default() -> Self {
        Execution {
            threads: thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1),
            progress: None,
        }
    }
}

/// Exhaustive job: every matrix index `offset + k*stride` below `matrix_count`.
#[derive(Clone, Copy, Debug)]
pub struct Exhaustive {
    pub offset: u128,
    pub stride: u128,
}

impl Default for Exhaustive {
    fn default() -> Self {
        Exhaustive {
            offset: 0,
            stride: 1,
        }
    }
}

impl Exhaustive {
    /// Number of matrices this job covers for `channels`.
    pub fn count(&self, channels: usize) -> u128 {
        let total = matrix_count(channels);
        if self.offset >= total {
            0
        } else {
            (total - self.offset).div_ceil(self.stride)
        }
    }
}

/// Sampled job: `samples` random matrices from a seeded PRNG.
#[derive(Clone, Copy, Debug)]
pub struct Sampled {
    pub samples: u64,
    pub seed: u64,
}

/// Simulate every matrix selected by `job` (all channels initially active).
pub fn exhaustive<const C: usize>(job: Exhaustive, rules: Rules, exec: Execution) -> Results {
    assert!(
        enumerable(C),
        "{} channels cannot be enumerated; use sampling",
        C
    );
    assert!(job.stride >= 1, "stride must be at least 1");
    let count = job.count(C);
    run_parallel::<C>(count, exec, move |k, tally| {
        let mut universe = Universe::<C>::new([0; C], 0, rules);
        let index = job.offset + k * job.stride;
        universe.reset(matrix::decode::<C>(index as u64), Universe::<C>::ALL);
        universe.run();
        tally.add(universe.code());
    })
}

/// Simulate `samples` random matrices. Each worker derives its own stream from `seed`.
pub fn sample<const C: usize>(job: Sampled, rules: Rules, exec: Execution) -> Results {
    let mut per_thread: Vec<Option<(SplitMix64, Universe<C>)>> = Vec::new();
    per_thread.resize_with(exec.threads.max(1), || None);
    // Each worker owns one PRNG; workers are numbered by the chunk they receive.
    run_chunked::<C>(
        job.samples as u128,
        exec,
        move |worker, first, len, tally| {
            let mut rng = SplitMix64::new(
                job.seed ^ (first as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ worker as u64,
            );
            let mut universe = Universe::<C>::new([0; C], 0, rules);
            for _ in 0..len {
                universe.reset(matrix::random::<C>(&mut rng), Universe::<C>::ALL);
                universe.run();
                tally.add(universe.code());
            }
        },
    )
}

/// Split `0..count` into one contiguous chunk per thread and call `f(k, tally)` for each k.
fn run_parallel<const C: usize>(
    count: u128,
    exec: Execution,
    f: impl Fn(u128, &mut Tally) + Sync + Send + Copy,
) -> Results {
    run_chunked::<C>(count, exec, move |_, first, len, tally| {
        for k in first..first + len {
            f(k, tally);
        }
    })
}

fn run_chunked<const C: usize>(
    count: u128,
    exec: Execution,
    f: impl Fn(usize, u128, u128, &mut Tally) + Sync + Send + Copy,
) -> Results {
    let threads = exec.threads.max(1).min(count.max(1) as usize);
    let chunk = count.div_ceil(threads as u128);
    let done = AtomicU64::new(0);
    let remaining = AtomicU64::new(threads as u64);
    let started = Instant::now();

    let mut results = Results::new();
    thread::scope(|scope| {
        let mut handles = Vec::with_capacity(threads);
        for worker in 0..threads {
            let first = worker as u128 * chunk;
            let len = chunk.min(count.saturating_sub(first));
            let (done, remaining) = (&done, &remaining);
            handles.push(scope.spawn(move || {
                let mut tally = Tally::new::<C>();
                // Report in slices so progress stays responsive on long runs.
                const SLICE: u128 = 1 << 16;
                let mut at = first;
                while at < first + len {
                    let n = SLICE.min(first + len - at);
                    f(worker, at, n, &mut tally);
                    done.fetch_add(n as u64, Ordering::Relaxed);
                    at += n;
                }
                remaining.fetch_sub(1, Ordering::Release);
                tally
            }));
        }

        if let Some(every) = exec.progress {
            while remaining.load(Ordering::Acquire) > 0 {
                thread::sleep(every);
                if remaining.load(Ordering::Acquire) == 0 {
                    break;
                }
                report(done.load(Ordering::Relaxed), count, started);
            }
        }

        for handle in handles {
            handle
                .join()
                .expect("worker panicked")
                .merge_into::<C>(&mut results);
        }
    });
    results
}

fn report(done: u64, count: u128, started: Instant) {
    let elapsed = started.elapsed().as_secs_f64();
    let rate = done as f64 / elapsed.max(1e-9);
    let left = (count as f64 - done as f64) / rate.max(1e-9);
    eprintln!(
        "progress {:>6.2}%  {} / {}  {:.0}/s  elapsed {}  eta {}",
        100.0 * done as f64 / count as f64,
        done,
        count,
        rate,
        human(elapsed),
        human(left)
    );
}

fn human(secs: f64) -> String {
    let s = secs.round() as u64;
    if s >= 3600 {
        format!("{}h{:02}m", s / 3600, (s % 3600) / 60)
    } else if s >= 60 {
        format!("{}m{:02}s", s / 60, s % 60)
    } else {
        format!("{s}s")
    }
}

/// Results ordered by count (descending), then signature.
pub fn sorted(results: &Results) -> Vec<(&str, u64)> {
    let mut rows: Vec<(&str, u64)> = results.iter().map(|(s, &n)| (s.as_str(), n)).collect();
    rows.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(b.0)));
    rows
}
