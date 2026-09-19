# What the system is, and what can be known about it

This note answers the research questions behind the project: can a universe be
classified from its matrix, what does "stable" mean here, which symmetries
exist, whether the element count matters, how to detect infinity, whether the
create/destroy variant behaves differently, and whether there is a formula.
Everything stated as a number was computed with the `exact`, `basins` and
`orbits` commands in this repository; everything stated as a theorem has a
proof sketch here and a test in `tests/exact.rs`.

## 1. The model collapses to a sign-driven integer system

After evaluation each element `e` has `count ≥ 0` and `antiCount ≤ 0` with at
most one of them nonzero, so the element is really one integer
`v_e = count + antiCount`. Which channel of `e` fires next depends only on
`sign(v_e)`. A firing row `r` adds `M[r][2e]` to `count` and subtracts
`M[r][2e+1]` from `antiCount`, i.e. it adds `M[r][2e] − M[r][2e+1] ∈ {−1,0,1}`
to `v_e`. Hence the entire dynamics is

```
v(t+1) = v(t) + d(sign v(t)),      d(s) = Σ_{f: s_f=+1} pos[f] + Σ_{f: s_f=−1} neg[f]
```

with two ternary *drift rows* per element, `pos[f], neg[f] ∈ {−1,0,1}^N`.
The very first step, where every channel fires, only sets the start state
`v(1) = Σ_f (pos[f] + neg[f])`.

Consequences:

* **Reduction.** The `4N²`-bit matrix `M` only matters through its `2N²`
  ternary differences. `2^(4N²)` bit matrices collapse to `3^(2N²)` reduced
  universes (N=3: 68.7 billion → 387 million, a factor 177). A reduced
  universe with `z` zero entries represents `2^z` bit matrices, so bit-matrix
  histograms are recovered by weighting. (`Reduced::from_rows`, `multiplicity`.)
* **Geometry.** The state space `ℤ^N` is cut by the coordinate hyperplanes into
  `3^N` *generalised orthants* (sign patterns, zeros included). Inside one
  orthant the velocity is constant, so a trajectory is a straight line until a
  coordinate hits or crosses zero. This is the discrete-time integer version of
  a *piecewise-constant-derivative (PCD) system*.
* **Event simulation.** Instead of stepping, jump: the number of steps until
  the next sign change is `k = min_e ceil(|v_e| / |d_e|)` over coordinates
  moving towards zero (`k = 1` for coordinates leaving zero). A run of a
  million steps costs as many operations as it has orthant changes.
  (`exact::analyze_traced`.)

## 2. Stability, and how to prove it

"Stable" here means *bounded*. Because the state is an integer vector, bounded
implies eventually periodic (finitely many states in a box), so bounded
trajectories are exactly the fixed points and the cycles. The exact analyser
returns one of six fates:

| fate | meaning | certificate |
|------|---------|-------------|
| `fixed` | the drift in the current orthant is zero | drift computed from the sign pattern |
| `cycle` | an event state recurs | state hashing at events; period and transient are exact |
| `ray` | the drift keeps every coordinate on its side of zero | absorbing-orthant test: `s_e = 0 ⇒ d_e = 0`, `s_e ≠ 0 ⇒ d_e = 0 or sign d_e = s_e` |
| `helix` | `v(t + T) = v(t) + D` forever, `D ≠ 0` | replay lemma below |
| `spiral` | returns to a sign pattern with strictly growing norm, 8 times in a row | **empirical** |
| `undecided` | budget exhausted | none |

**Replay lemma (helix).** Suppose two event states `u` (earlier) and `v = u + D`
have the same sign pattern, `D ≠ 0`, and every coordinate `e` with `D_e ≠ 0`
kept a constant nonzero sign along the segment from `u` to `v` with
`sign D_e` equal to that sign. Then the trajectory from `v` repeats the segment
verbatim and arrives at `v + D`, and so on forever. Proof: drifts depend only on
sign patterns, which agree at `u` and `v`. The jump lengths are decided by the
coordinates that reach zero; those have `D_e = 0` and replay exactly. Coordinates
with `D_e ≠ 0` never reached zero in the segment and are now further from zero
in the same direction, so they still do not. Hence the same events occur with
the same jumps, and the displacement is again `D`. This is the analogue of the
Karp–Miller token-generator criterion for Petri nets (`M →σ M + L` implies
unboundedness) and it is what turned 1.4 % undecided size-3 universes into
proven divergence. With it, the full size-3 census has **no** undecided
universe: 99.93 % are proven, 0.07 % are empirical spirals (9 minutes on 8 cores).

**Is a repeating state the only way to verify stability?** No. Stability
(boundedness) needs the cycle or fixed-point certificate, but *instability* has
three independent certificates (ray, helix, and, empirically, spiral), and all
of them are found without a step limit. What the original 10-step heuristic
called `1`/`-1` mixes rays with cycles of period > 10 and with slow helices.

**Detecting infinity.** Two rigorous mechanisms (ray, helix) plus one empirical
one (spiral). Spirals come in two kinds, both visible with `show`:

* additive: `show --size 3 --index 246774369` — the amplitude grows by 6 per
  loop of 7 events;
* multiplicative: `show --size 3 --index 27411247` — `x ↦ 2x + 13` per loop;
  `show --size 3 --index 171708036` — ratio ≈ 32 per loop, norm 10^12 after
  78 events.

A spiral is a loop through orthants whose *return map* is affine with a
non-identity linear part; proving it diverges needs the eigenvalues of that
map (§6). No spiral has been observed to stop growing, but the label is
honest: it is a strong observation, not a theorem.

## 3. Results

Weights are bit-matrix weighted (what `run` would see); "reduced" counts each
reduced universe once. Start state: all channels fire once.

Sign model (original):

| N | reduced universes | fixed | cycle | ray | helix | spiral | undecided |
|---|---|---|---|---|---|---|---|
| 1 | 9 | 62.5 % | 0 | 37.5 % | 0 | 0 | 0 |
| 2 | 6 561 | 29.82 % | 0.54 % | 69.65 % | 0 | 0 | 0 |
| 3 | 387 420 489 | 10.69 % | 0.53 % | 88.11 % | 0.61 % | 0.07 % | 0 |

(The size-3 line is from `results/exact-size3.txt`; see that file for the
exact counts and the full cycle-period histogram.)

Size-3 facts that no step-limited method can see: cycles with periods up to
1510 and transient+period up to 1824 (`show --size 3 --index 7894996`);
fixed points and rays certified as late as step 668. The 10-step heuristic
misclassifies about 0.2 % of size-3 bit matrices (cycles reported as growth,
helices and rays reported as bounded), see the confusion matrix in the results
file.

Size 4, sampled (2·10^9 uniform reduced universes, 158 s; 95 % half-widths
below 0.001 %; `results/exact-size4-sampled.txt`):

| fixed | cycle | ray | helix | spiral | undecided |
|---|---|---|---|---|---|
| 3.02 % | 0.23 % | 95.25 % | 1.38 % | 0.12 % | 0.0001 % |

The trend of §5 continues (rays 95 %), and the event-driven analyser exposes
time scales no step-by-step simulation could reach: in the sample there is a
cycle of period 989 403 630 (`show --size 4 --index 345302372823031`), a fixed
point first reached at step 805 437 199, and a ray certified at step
1 092 831 065 384 — a trajectory that crosses coordinate axes for a trillion
steps before settling into its final orthant. Certificates cost one event
each, so these cases take microseconds.

Clamped model (create / destroy with positivity, §5):

| N | reduced universes | fixed | cycle | ray | helix | spiral |
|---|---|---|---|---|---|---|
| 1 | 3 | 75.0 % | 0 | 25.0 % | 0 | 0 |
| 2 | 81 | 62.89 % | 0 | 37.11 % | 0 | 0 |
| 3 | 19 683 | 48.30 % | 0.23 % | 51.34 % | 0.12 % | 0 |
| 4 | 43 046 721 | 35.63 % | 0.37 % | 63.54 % | 0.46 % | 0.01 % |

Basins (`basins` command: every start in a box, origin excluded):

| N | starts | universes with a single fate over all starts | globally bounded (only fixed/cycle) | globally divergent (only ray) |
|---|---|---|---|---|
| 2 | 48 (radius 3) | 22.7 % | 13.0 % | 12.4 % |
| 3 | 26 (radius 1, 200 k sampled) | 27.3 % | 2.4 % | 25.7 % |

## 4. Symmetries

Relabelling elements and swapping the roles of positive and negative for an
element (`v_e ↦ −v_e`, `pos[e] ↔ neg[e]`, column `e` negated) conjugates the
dynamics, so fates are invariant under the **hyperoctahedral group** `B_N` of
signed permutations, order `2^N · N!` (`Reduced::transform`, test
`symmetry_preserves_fate`). Burnside counts (`orbits` command):

| N | group | reduced universes | orbits | bit matrices | orbits |
|---|---|---|---|---|---|
| 1 | 2 | 9 | 6 | 16 | 10 |
| 2 | 8 | 6 561 | 873 | 65 536 | 8 548 |
| 3 | 48 | 387 420 489 | 8 093 160 | 2^36 | 1.44 · 10^9 |
| 4 | 384 | 1.85 · 10^15 | 4.83 · 10^12 | 2^64 | 4.8 · 10^16 |

The symmetry is visible in the raw histograms: every signature count in
`results/size3.json` is shared by its whole orbit (e.g. the six single-`1`/`-1`
signatures all occur 3 268 864 628 times).

There is no further linear symmetry: the drift is not equivariant under
general permutations of channels because the "which channel fires" rule ties
channel `2e` to `v_e > 0`.

## 5. Does the number of elements matter?

Yes, monotonically, and in a predictable direction. The ray share rises with
`N` (37 → 70 → 88 % for the sign model, 25 → 37 → 51 → 64 % clamped) and the
fixed share falls, because the drift inside an orthant is a sum of `|s|`
ternary rows: with more elements it is less likely to be exactly zero (fixed)
and, once it is nonzero, a random sign-compatible drift is what a ray needs.
Cycles stay at about half a percent of bit matrices for every `N` computed.
Helices and spirals only exist from `N = 3` (sign) / `N = 3–4` (clamped): a
helix needs a coordinate that drifts while others oscillate; a spiral needs a
loop through at least four orthants (in 2-D, integer effects on the axes make
the loop grow additively — see `tests/exact.rs::fates`).

So there is a repeating pattern in the *kinds* of behaviour (the same six fates
suffice for every `N` tried) but not in their proportions.

## 6. Is there a formula? What theory covers this?

*Inside one orthant* the trajectory is the affine formula `v(t) = v(t0) + (t − t0) d`.
*Along one loop of orthants* the return map is affine, `x ↦ A x + b` on the
set of states that follow that loop, with `A` an integer matrix built from the
jump rules. Its spectral radius decides the loop's behaviour: `A = I` gives a
helix (`b ≠ 0`) or a cycle (`b = 0`); `ρ(A) > 1` gives a multiplicative spiral
(`x ↦ 2x + 13` above); `ρ(A) < 1` is only possible on a bounded set and ends
in a cycle. The global question — *which* loop a given start ends in — is the
hard part, and it is provably hard:

* The continuous-time counterpart, PCD systems, was analysed by Asarin, Maler
  and Pnueli, "Reachability analysis of dynamical systems having
  piecewise-constant derivatives", *Theoretical Computer Science* 138 (1995):
  reachability is decidable in two dimensions and **undecidable in three or
  more** (a 3-D PCD system can simulate a Turing machine). No closed-form
  predictor for `N ≥ 3` can exist in general; the certificates in §2 are the
  right tool, and the residual is exactly the spiral class.
* Glass networks (Glass & Kauffman 1973; survey of periodic orbits in
  [arXiv:2311.11879](https://arxiv.org/abs/2311.11879)) are piecewise-linear
  ODEs whose vector field depends only on which side of thresholds the
  variables are on, i.e. the same "sign pattern decides the drift" structure
  with an added decay term. Their analysis via wall-to-wall return maps is the
  continuous version of §2.
* Vector addition systems / Petri nets: the clamped variant lives here.
  Boundedness is decidable by the Karp–Miller coverability tree (1969); see the
  survey [arXiv:2411.01592](https://arxiv.org/abs/2411.01592). Our helix test is
  that criterion adapted to deterministic maximal-step firing.
* Reachability for piecewise affine maps in one dimension is already open in
  general ([arXiv:1510.04121](https://arxiv.org/abs/1510.04121)), which is
  the honest answer to "algebraic or differential formula": piecewise yes,
  globally no.

## 7. The create/destroy variant

`--model clamped`: one channel per element, active while `v_e > 0`; a row
carries create (+1) / destroy (−1) per target; `v_e ← max(0, v_e + drift)`.
It is not the same generative system: it has far more fixed points (a
destroyed element stays at zero instead of becoming negative and firing its
other channel), no cycles at all for `N ≤ 2`, and helices/spirals appear one
size later. Its state space is `ℕ^N`, which makes it a Petri-net-like object
where boundedness is decidable in principle; the sign model's `ℤ^N` with sign
switching is the PCD-like object where it is not. The two models coincide only
on trajectories that never go negative.

## 8. Efficiency: what mattered and what would not

Done, in order of payoff:

1. **Model reduction** (§1): 177× fewer universes at N=3, more at larger N.
2. **Event jumps + certificates**: cost proportional to orthant changes, no
   step cap, and every fate but spiral is proven.
3. **No allocation in the inner loop.** The static musl build's `malloc` takes
   a global lock, so per-universe `Vec`/`HashMap` allocations made two threads
   *slower* than one. A reusable per-thread `Workspace` removed every
   allocation from the hot path: the exhaustive size-3 census went from
   527 s to 6.6 s on 8 threads (59 million universes per second), and a
   base-3 odometer replaced per-index decoding.
4. **Symmetry**: another 48× available at N=3 (orbit representatives via
   `Reduced::canonical`), 384× at N=4. Not needed at N=3 (seconds), but the
   only route to an exhaustive N=4: 1.85·10^15 reduced universes take about
   a year at 59 M/s, 4.8·10^12 orbit representatives about a day — provided
   they are *generated* canonically (orderly generation) rather than filtered,
   since a canonical check costs 384 transforms. Until then N=4 is sampled:
   2·10^9 samples in a few minutes give every fate to ±0.002 %
   (`results/exact-size4-sampled.txt`).
5. **GPU**: the per-universe work is branchy integer logic with tiny state,
   which GPUs handle well for the 10-step heuristic (2^36 matrices would take
   minutes on the RTX 3070 instead of 40). For the exact analyser the
   unbounded loop and hashing map poorly; and the reductions above already
   beat a GPU brute force by orders of magnitude. If a GPU path is wanted, the
   right target is sampling at `N ≥ 4`: no CUDA toolkit is installed and the
   static musl build cannot `dlopen` the driver, so it needs `build-essential`
   and `nvidia-cuda-toolkit` first.

## 9. Open questions worth attacking next

* Turn the spiral observation into a certificate: extract the affine return
  map of the observed loop, verify it symbolically (residue classes of the
  ceilings), and check `ρ(A) > 1`.
* Orderly generation of `B_N` orbit representatives to make exact N=4 sampling
  unbiased at the orbit level.
* Classify the *phase portrait* rather than one start: which reduced universes
  are globally bounded (all starts fixed/cycle)? At N=2 it is 13 %, at N=3
  about 2.4 %. A structural characterisation (in terms of `pos`/`neg`) of that
  set would be a real theorem about this system.
* Relate the cycle-period histogram (heavy tail to 1510 at N=3) to loop return
  maps with `A` of finite order.
