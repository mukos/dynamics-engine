//! Minimal fork-join helper: split `0..count` into contiguous chunks, one per thread.

use std::thread;

/// Run `work(k, &mut acc)` for every `k` in `0..count` across `threads` threads and
/// fold the per-thread accumulators with `merge`.
pub fn par_fold<T: Send>(
    count: u64,
    threads: usize,
    init: impl Fn() -> T + Sync,
    work: impl Fn(u64, &mut T) + Sync,
    mut merge: impl FnMut(&mut T, T),
) -> T {
    let threads = threads.max(1).min(count.max(1) as usize);
    let chunk = count.div_ceil(threads as u64);
    let mut total = init();
    thread::scope(|scope| {
        let handles: Vec<_> = (0..threads)
            .map(|t| {
                let (init, work) = (&init, &work);
                scope.spawn(move || {
                    let first = t as u64 * chunk;
                    let last = (first + chunk).min(count);
                    let mut acc = init();
                    for k in first..last {
                        work(k, &mut acc);
                    }
                    acc
                })
            })
            .collect();
        for h in handles {
            merge(&mut total, h.join().expect("worker panicked"));
        }
    });
    total
}
