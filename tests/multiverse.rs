use std::collections::BTreeMap;
use std::time::Duration;

use system::multiverse::{exhaustive, sample, Execution, Exhaustive, Results, Sampled};
use system::universe::Rules;

/// Output of the original TypeScript implementation (`ts-node real/init.ts`, size 2).
fn original_size_2() -> Results {
    [
        ("|,|,|,|", 8416),
        ("0,0,0,-1", 6884),
        ("0,-1,0,-1", 4543),
        ("1,0,0,-1", 4543),
        ("0,0,0,0", 11412),
        ("0,-1,0,0", 6884),
        ("1,0,0,0", 6884),
        ("0,0,1,0", 6884),
        ("0,-1,1,0", 4543),
        ("1,0,1,0", 4543),
    ]
    .into_iter()
    .map(|(s, n)| (s.to_string(), n))
    .collect()
}

fn quiet(threads: usize) -> Execution {
    Execution {
        threads,
        progress: None,
    }
}

#[test]
fn exhaustive_size_2_reproduces_the_original_results() {
    let r = exhaustive::<4>(Exhaustive::default(), Rules::default(), quiet(1));
    assert_eq!(r, original_size_2());
    let r = exhaustive::<4>(Exhaustive::default(), Rules::default(), quiet(8));
    assert_eq!(r, original_size_2());
}

#[test]
fn shards_partition_the_exhaustive_run() {
    let mut merged: Results = BTreeMap::new();
    for offset in 0..5u128 {
        for (k, v) in exhaustive::<4>(Exhaustive { offset, stride: 5 }, Rules::default(), quiet(3))
        {
            *merged.entry(k).or_insert(0) += v;
        }
    }
    assert_eq!(merged, original_size_2());
}

#[test]
fn exhaustive_size_1_covers_all_16_matrices() {
    let r = exhaustive::<2>(Exhaustive::default(), Rules::default(), quiet(2));
    assert_eq!(r.values().sum::<u64>(), 16);
}

#[test]
fn sampling_is_reproducible_and_counts_every_sample() {
    let job = Sampled {
        samples: 5000,
        seed: 7,
    };
    let exec = Execution {
        threads: 4,
        progress: Some(Duration::from_millis(1)),
    };
    let a = sample::<6>(job, Rules::default(), exec);
    let b = sample::<6>(job, Rules::default(), exec);
    assert_eq!(a, b);
    assert_eq!(a.values().sum::<u64>(), 5000);
    assert_ne!(
        a,
        sample::<6>(
            Sampled {
                samples: 5000,
                seed: 8
            },
            Rules::default(),
            exec
        )
    );
}
