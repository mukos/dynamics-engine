use system::matrix::{decode, matrix_count, random, to_vec};
use system::rng::SplitMix64;
use system::universe::{Rules, Symbol, Universe};

fn rows<const C: usize>(m: [[u8; C]; C]) -> [u16; C] {
    let mut out = [0u16; C];
    for (r, row) in m.iter().enumerate() {
        for (c, &cell) in row.iter().enumerate() {
            out[r] |= (cell as u16) << c;
        }
    }
    out
}

#[test]
fn decode_is_row_major_msb_first_like_the_original() {
    assert_eq!(to_vec(&decode::<2>(0)), vec![vec![0, 0], vec![0, 0]]);
    assert_eq!(to_vec(&decode::<2>(1)), vec![vec![0, 0], vec![0, 1]]);
    assert_eq!(to_vec(&decode::<2>(0b1000)), vec![vec![1, 0], vec![0, 0]]);
    assert_eq!(to_vec(&decode::<2>(0b0110)), vec![vec![0, 1], vec![1, 0]]);
    assert_eq!(
        to_vec(&decode::<2>((matrix_count(2) - 1) as u64)),
        vec![vec![1, 1], vec![1, 1]]
    );
    assert_eq!(
        to_vec(&decode::<4>(1 << 15)),
        vec![vec![1, 0, 0, 0], vec![0; 4], vec![0; 4], vec![0; 4]]
    );
}

#[test]
fn random_rows_are_seed_deterministic_and_masked() {
    let a = random::<6>(&mut SplitMix64::new(42));
    let b = random::<6>(&mut SplitMix64::new(42));
    assert_eq!(a, b);
    assert!(a.iter().all(|&r| r < 1 << 6));
}

#[test]
fn self_reinforcing_positive_channel_grows() {
    let mut u = Universe::<2>::new(rows([[1, 0], [0, 0]]), 0b01, Rules::default());
    assert_eq!(u.run(), 10);
    let counts: Vec<i32> = u.history[..10].iter().map(|s| s[0]).collect();
    assert_eq!(counts, (1..=10).collect::<Vec<_>>());
    assert_eq!(u.symbol(0), Symbol::Growth);
    assert_eq!(u.symbol(1), Symbol::Bounded);
    assert_eq!(u.signature(), "1,0");
}

#[test]
fn dead_end_negative_channel_terminates_early() {
    let mut u = Universe::<2>::new(rows([[0, 1], [0, 0]]), 0b01, Rules::default());
    assert_eq!(u.run(), 1);
    assert_eq!(u.history[0], [0, -1]);
    assert_eq!(u.signature(), "|,|");
}

#[test]
fn empty_matrix_performs_no_step() {
    let mut u = Universe::<2>::new(rows([[0, 0], [0, 0]]), Universe::<2>::ALL, Rules::default());
    assert_eq!(u.run(), 0);
    assert_eq!(u.signature(), "|,|");
}

#[test]
fn classifier_skip_previous_ignores_the_step_before_last() {
    // Channel 0 history 1, 2, 2: the last value only repeats the previous step.
    let mut u = Universe::<2>::new(
        rows([[0, 0], [0, 0]]),
        0,
        Rules {
            max_steps: 3,
            skip_previous: true,
        },
    );
    u.steps = 3;
    u.history[0] = [1, 0];
    u.history[1] = [2, 0];
    u.history[2] = [2, 0];
    assert_eq!(
        u.symbol(0),
        Symbol::Growth,
        "original rule: previous step is not compared"
    );
    u.rules.skip_previous = false;
    assert_eq!(
        u.symbol(0),
        Symbol::Bounded,
        "fixed rule: previous step counts as a repeat"
    );
    assert_eq!(u.symbol(1), Symbol::Bounded);
}

#[test]
fn reset_clears_state_between_runs() {
    let mut u = Universe::<2>::new(rows([[1, 0], [0, 0]]), 0b01, Rules::default());
    u.run();
    u.reset(rows([[0, 0], [0, 0]]), Universe::<2>::ALL);
    assert_eq!(u.run(), 0);
    assert_eq!(u.values, [0, 0]);
}
