use dynamics_engine::exact::{
    analyze, analyze_universe, matrix_orbits, reduced_orbits, Budget, Fate, Reduced,
};
use dynamics_engine::matrix::decode;

fn w<const N: usize>(pos: [[i8; N]; N], neg: [[i8; N]; N]) -> Reduced<N> {
    Reduced { pos, neg }
}

#[test]
fn reduction_matches_bit_rows() {
    // rows: channel 0 feeds count of e0 (bit 0) and anti of e1 (bit 3); channel 1 feeds both channels of e0 (net 0).
    let rows = [0b1001u16, 0b0011, 0, 0];
    let r = Reduced::<2>::from_rows(&rows);
    assert_eq!(r.pos[0], [1, -1]);
    assert_eq!(r.neg[0], [0, 0]);
    assert_eq!(r.multiplicity(), 1 << 6);
    assert_eq!(Reduced::<2>::from_rows(&r.to_rows()), r);
}

#[test]
fn index_roundtrip_and_count() {
    assert_eq!(Reduced::<2>::count(), 6561);
    for i in [0u64, 1, 2, 3, 6560, 1234] {
        let r = Reduced::<2>::from_index(i);
        let mut back = 0u64;
        for row in r.pos.iter().chain(r.neg.iter()) {
            for &x in row {
                back = back * 3 + (x + 1) as u64;
            }
        }
        assert_eq!(back, i);
    }
    assert_eq!(Reduced::<2>::from_index(0).pos, [[-1, -1], [-1, -1]]);
}

#[test]
fn multiplicities_sum_to_all_bit_matrices() {
    let total: u64 = (0..Reduced::<2>::count() as u64)
        .map(|i| Reduced::<2>::from_index(i).multiplicity())
        .sum();
    assert_eq!(total, 1 << 16);
}

#[test]
fn every_bit_matrix_reduces_consistently() {
    // The reduced drift must equal what the bit-level Universe computes, for all 4x4 matrices.
    for index in 0..1u64 << 16 {
        let rows = decode::<4>(index);
        let r = Reduced::<2>::from_rows(&rows);
        let mut u = dynamics_engine::universe::Universe::<4>::new(
            rows,
            dynamics_engine::universe::Universe::<4>::ALL,
            Default::default(),
        );
        u.run();
        if u.steps >= 1 {
            let v = r.initial_state();
            assert_eq!(
                [
                    u.history[0][0] + u.history[0][1],
                    u.history[0][2] + u.history[0][3]
                ],
                [v[0] as i32, v[1] as i32],
                "index {index}"
            );
        }
    }
}

#[test]
fn fates() {
    let b = Budget::default();
    // Self-reinforcing: v grows forever.
    assert_eq!(analyze(&w([[1]], [[0]]), [1], b), Fate::Ray { step: 1 });
    // Self-damping from +3: 3,2,1,0 then stuck at origin.
    assert_eq!(analyze(&w([[-1]], [[0]]), [3], b), Fate::Fixed { step: 4 });
    // Drift of magnitude 1 cannot cross zero: +1 -> 0 and the element goes silent.
    assert_eq!(analyze(&w([[-1]], [[1]]), [1], b), Fate::Fixed { step: 2 });
    // Two elements pushing each other with drift 2 flip sign every step: period 2.
    let flip = w([[-1, -1], [-1, -1]], [[1, 1], [1, 1]]);
    assert_eq!(
        analyze(&flip, [1, 1], b),
        Fate::Cycle {
            transient: 1,
            period: 2
        }
    );
    // Zero drift.
    assert_eq!(analyze(&w([[0]], [[0]]), [5], b), Fate::Fixed { step: 1 });
    // 2D rotation: e0>0 pushes e1 up, e1>0 pushes e0 down, ... The axis states add one
    // extra step per quarter turn, so the integer orbit spirals outward by +1 per quarter.
    let rot = w([[0, 1], [-1, 0]], [[0, -1], [1, 0]]);
    let fate = analyze(&rot, [1, 0], b);
    assert!(matches!(fate, Fate::Spiral { returns: 8, .. }), "{fate:?}");
    // Same system without spiral detection runs out of budget with a growing norm.
    let nb = Budget {
        spiral_returns: 0,
        max_events: 400,
        ..b
    };
    assert!(
        matches!(analyze(&rot, [1, 0], nb), Fate::Undecided { events: 400, max_norm } if max_norm >= 90),
        "{:?}",
        analyze(&rot, [1, 0], nb)
    );
}

#[test]
fn undecided_when_budget_is_exhausted() {
    let b = Budget {
        max_events: 3,
        max_norm: 1 << 40,
        spiral_returns: 0,
    };
    let rot = w([[0, 1], [-1, 0]], [[0, -1], [1, 0]]);
    assert!(matches!(
        analyze(&rot, [1, 0], b),
        Fate::Undecided { events: 3, .. }
    ));
}

#[test]
fn exact_and_heuristic_agree_on_clear_cases() {
    // Every size-2 matrix whose 10-step signature shows growth must not be a proven fixed point.
    for index in 0..1u64 << 16 {
        let rows = decode::<4>(index);
        let mut u = dynamics_engine::universe::Universe::<4>::new(
            rows,
            dynamics_engine::universe::Universe::<4>::ALL,
            Default::default(),
        );
        u.run();
        let fate = analyze_universe(&Reduced::<2>::from_rows(&rows), Budget::default());
        if u.signature().split(',').all(|s| s == "|") {
            assert!(
                matches!(fate, Fate::Fixed { .. }),
                "terminated early but fate {fate:?} at {index}"
            );
        }
    }
}

#[test]
fn symmetry_preserves_fate() {
    let b = Budget::default();
    for index in (0..Reduced::<2>::count() as u64).step_by(7) {
        let r = Reduced::<2>::from_index(index);
        let f = analyze_universe(&r, b).kind();
        dynamics_engine::exact::for_each_signed_permutation::<2>(|perm, flip| {
            assert_eq!(analyze_universe(&r.transform(perm, flip), b).kind(), f);
        });
        assert_eq!(r.canonical(), r.canonical().canonical());
    }
}

#[test]
fn orbit_counts() {
    // N=1: group {id, flip}. Reduced entries (pos, neg) in {-1,0,1}^2: 9 total;
    // flip maps (p,n) -> (-n,-p); fixed iff n = -p: 3. Orbits = (9+3)/2 = 6.
    assert_eq!(reduced_orbits::<1>(), 6);
    // N=1 bit matrices: 16 total, flip swaps rows and columns: fixed = 2^cycles = 2^2 -> (16+4)/2 = 10.
    assert_eq!(matrix_orbits::<1>(), 10);
    assert!(reduced_orbits::<2>() * 8 >= 6561);
}

#[test]
fn clamped_fates() {
    use dynamics_engine::exact::{analyze_clamped, Clamped};
    let b = Budget::default();
    // Self-creating element grows forever.
    assert_eq!(
        analyze_clamped(&Clamped::<1> { rows: [[1]] }, [1], b),
        Fate::Ray { step: 1 }
    );
    // Self-destroying element empties and stops.
    assert_eq!(
        analyze_clamped(&Clamped::<1> { rows: [[-1]] }, [3], b),
        Fate::Fixed { step: 4 }
    );
    // e0 creates e1, e1 destroys e0: e0 drains while e1 grows, then e1 alone has nothing to destroy -> fixed.
    let c = Clamped::<2> {
        rows: [[0, 1], [-1, 0]],
    };
    let fate = analyze_clamped(&c, [2, 0], b);
    assert!(matches!(fate, Fate::Fixed { .. }), "{fate:?}");
    assert_eq!(Clamped::<2>::count(), 81);
    let total: u64 = (0..81u64)
        .map(|i| Clamped::<2>::from_index(i).multiplicity())
        .sum();
    assert_eq!(total, 1 << 8);
}

#[test]
fn helix_is_certified() {
    use dynamics_engine::exact::{analyze_clamped, Clamped};
    let b = Budget::default();
    // Size-3 reduced universe 560: sign pattern alternates [0,-,0] / [+,-,+] while v_1 falls by 4 every 2 steps.
    let w = Reduced::<3>::from_index(560);
    let fate = analyze_universe(&w, b);
    assert!(
        matches!(
            fate,
            Fate::Helix {
                period: 2,
                shift_norm: 4,
                ..
            }
        ),
        "{fate:?}"
    );
    // Clamped: e0 is self-sustaining and creates e1; e1 destroys itself. e1 flips 0/1 while e0 grows.
    let c = Clamped::<2> {
        rows: [[1, 1], [0, -1]],
    };
    let fate = analyze_clamped(&c, [1, 0], b);
    assert!(
        matches!(fate, Fate::Ray { .. } | Fate::Helix { .. }),
        "{fate:?}"
    );
}

#[test]
fn index_is_inverse_of_from_index() {
    for i in [0u64, 1, 6560, 12345, 387420488] {
        assert_eq!(Reduced::<3>::from_index(i).index(), i);
    }
    assert_eq!(Reduced::<2>::from_index(613).index(), 613);
}

#[test]
fn increment_walks_index_order() {
    let mut w = Reduced::<2>::from_index(0);
    for i in 1..Reduced::<2>::count() as u64 {
        assert!(w.increment());
        assert_eq!(w, Reduced::<2>::from_index(i), "at {i}");
    }
    assert!(!w.increment());
    assert_eq!(w, Reduced::<2>::from_index(0));
}

#[test]
fn workspace_reports_attractors() {
    use dynamics_engine::exact::{analyze_ws, Workspace};
    let b = Budget::default();
    let mut ws = Workspace::<2>::new();
    // Self-damping pair from (5,3): stops at the origin.
    let damp = w([[-1, 0], [0, -1]], [[0, 0], [0, 0]]);
    assert!(matches!(
        analyze_ws(&damp, [5, 3], 1, b, &mut ws, &mut |_, _, _| {}),
        Fate::Fixed { .. }
    ));
    assert_eq!(ws.attractor, vec![0, 0]);
    // Ray: attractor is the drift.
    let ray = w([[1, 1], [0, 0]], [[0, 0], [0, 0]]);
    assert!(matches!(
        analyze_ws(&ray, [1, 2], 1, b, &mut ws, &mut |_, _, _| {}),
        Fate::Ray { .. }
    ));
    assert_eq!(ws.attractor, vec![1, 1]);
    // Two starts on the same cycle share the attractor id.
    let cyc = Reduced::<2>::from_index(613);
    analyze_ws(&cyc, cyc.initial_state(), 1, b, &mut ws, &mut |_, _, _| {});
    let a = ws.attractor.clone();
    analyze_ws(&cyc, [1, 0], 1, b, &mut ws, &mut |_, _, _| {});
    assert_eq!(a, ws.attractor);
}

#[test]
fn basin_map_features() {
    use dynamics_engine::experiments::basin_map;
    let b = Budget::default();
    // Universe 560: parity checkerboard between two cycles.
    let m = basin_map(&Reduced::<2>::from_index(560), 6, b);
    assert_eq!(m.features.moving, 2);
    assert!(m.features.stripes > 0.8, "{:?}", m.features);
    assert!(m.features.boundary > 0.95);
    // Zero drift everywhere: every cell is an inert fixed point.
    let m = basin_map(&Reduced::<2>::from_index(3280), 3, b);
    assert_eq!(m.features.inert, 49);
    assert_eq!(m.features.boundary, 0.0);
    assert_eq!(m.features.conic, 1.0);
}
