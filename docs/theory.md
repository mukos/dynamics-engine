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
| `spiral` | returns to a sign pattern grow without bound | lap lemma below (*proven*); otherwise 8 growing returns when the budget runs out (*empirical*) |
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
universe; with the lap lemma below, all but 192 of them are proven.

**Is a repeating state the only way to verify stability?** No. Stability
(boundedness) needs the cycle or fixed-point certificate, but *instability* has
three independent certificates (ray, helix, spiral), and all of them are found
without a step limit. What the original 10-step heuristic
called `1`/`-1` mixes rays with cycles of period > 10 and with slow helices.

**Detecting infinity.** Three rigorous mechanisms: ray, helix, and the lap
lemma for spirals. `show` prints the lap map of every proven spiral:

* additive: `show 3-246774369` — `x ↦ x + 6` per lap of 7 events;
* multiplicative: `show 3-27411247` — `x ↦ 2x + 12` per lap of 10 events;
  `show 3-171708036` — `x ↦ 32x + 20`.

**Growth is not divergence.** The first spiral rule, "8 growing returns to one
sign pattern in a row", was wrong for 2 784 size-3 universes: 2 736 of them
later fall into a cycle and 48 into a ray. Universe 3-7295755 grows for 725
steps and then cycles with period 104. Only the lap lemma separates them.

**Lap lemma (spiral, line case).** Let the last few returns to a pattern `s`
lie on a line `p + x·q` (`q` primitive). Run one lap from the point `p + x·q`
with `x` an unknown, written `x = c + m·y` (`y ≥ 0` free), and carry the state
as affine forms in `y`. Every branch the simulator takes compares an affine form
with zero, so it holds for all large `y` or for none; every jump length
`⌈|v_e| / |d_e|⌉` is affine in `y` once `m` is a multiple of `|d_e|` and of the
coefficient, so the residue class is refined (`m ← m · k`, up to 4096) until all
of them are. The lap ends back in `s` at `p' + x'·q` with `x' = α y + β`. The
family is a finite union of parallel lines, each split into residue classes;
the lap is proven when every member maps into the family, with `x' > x` beyond
a threshold `X` that the current return exceeds. By induction the returns grow
forever, so the orbit is unbounded and never repeats. Two refinements matter:

* coordinates that keep their sign through the whole lap ("steady") are left
  out of the line and carried as an affine lower bound `|z| ≥ B(x)` that must
  be preserved;
* a slope-1 lap whose length does not depend on `x` is a helix in disguise and
  is left to the replay lemma.

**Plane case.** Returns that span two directions are written `p + x·e₁ + y·e₂`
on a saturated, Gauss-reduced basis of the lattice plane. The lap is run once
along the branch actually taken near the current return; every branch
condition becomes an affine inequality in `(x, y)`. The proof is a cone
`R = x₀ + cone(E₁, E₂)` (edges from the inequalities' boundary lines, the axes
and rational eigenvectors of the lap's linear part) such that the conditions hold
on `R`, the corner and both edges map back into `R`, the lap gets longer along
both edges (so it is not a helix), and some linear functional grows by at least
1 per lap. This covers multiplicative growth, two-scale growth and shear maps
whose returns grow polynomially.

On the size-3 census the lap lemma proves 277 152 of 277 344 spirals
(99.93 %). The 192 left (e.g. `show 3-12194199`) have returns that span all
three dimensions; they keep the empirical label. As a soundness check, all
277 152 proven spirals were rerun for 30 000 events with the spiral rule off;
none became fixed, cycle, ray or helix (`tests/exact.rs` repeats this on a
3-million-universe slice).

### 2b. Bounded from every start (N = 2)

A state far from the origin is either inside a quadrant, where it moves in a
straight line until a coordinate crosses zero, or *near an axis*: one
coordinate `u` small (`|u| ≤ K`, the largest drift entry) and the other, `X`,
large. Near an axis the motion of `u` depends only on the sign of `X`, not its
size, so it can be run once with `X = ∞`. Each of the `4(2K + 1)` axis states
then ends in one of four ways, valid for every `X` above a threshold: the drift
vanishes or `u` cycles with `X` returning exactly (bounded); `u` cycles while
`X` shrinks (it drains into a box); `u` cycles while `X` grows, or the pattern is
absorbing (an escape); or `u` runs off while `X` heads to zero. In the last case
the crossing time is `⌈(X + δ)/|d_X|⌉`, so the next axis state is reached with
`X' = λ(X + δ + ρ) + c`, `λ = |d_u| / |d_X|`, where `ρ` depends only on
`X mod |d_X|`.

This gives a finite graph of axis states with affine maps on `X`. If every
simple cycle has `Λ < 1`, or `Λ = 1` with a non-positive shift, then every closed
walk does too (a closed walk is a product of simple cycles, and conjugating a
shift by a path keeps its sign), node weights `φ` with `φ_t λ ≤ φ_s` exist, and
`φ·X` stays bounded along any walk. So every orbit is periodic, stays bounded in
the graph, or enters the box `[−R, R]²` with `R` above every threshold; if every
start in that box is fixed or a cycle, the universe is bounded from every start.
All 6 561 two-element universes are decided this way (box radius at most 13);
2.57 million random starts up to `2^40` in the 856 bounded ones all settle. At
`N = 3` a far state can have two independent large coordinates, so the same
argument needs two unknowns and cones (as in the plane case above); that is
the open part.

**A closed form at N = 2.** The census turns the proof into a formula. A
two-element universe is bounded from every start **if and only if** no sign
pattern `s` has a nonzero drift with `s_e = 0 ⇒ d_e = 0` and `s_e ≠ 0 ⇒ d_e·s_e ≥ 0`
(no absorbing pattern), **unless** it is the sign rotation
`d(v) = ±(sgn v₁, −sgn v₀)` (universes 2800 and 3760). Of the 6 561 universes,
5 703 have an absorbing pattern, 2 are the rotation, and the other 856 are
bounded. The rotation is the one neutral far-field loop: its continuous version
conserves `|v₀| + |v₁|`, and the integer rounding at each axis adds 1, so it
spirals out additively. At N = 3 no such formula holds: 2 392 328 universes
without an absorbing pattern still escape (helices and spirals).

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
1714 and transient+period up to 3771 (`show 3-81866077`);
fixed points and rays proven as late as step 668. The 10-step heuristic
misclassifies about 0.2 % of size-3 bit matrices (cycles reported as growth,
helices and rays reported as bounded); `exact --size 3 --heuristic` prints the
confusion matrix.

Size 4, sampled (2·10^9 uniform reduced universes, 1 371 s; 95 % half-widths
below 0.001 %; `results/exact-size4-sampled.txt`):

| fixed | cycle | ray | helix | spiral | undecided |
|---|---|---|---|---|---|
| 3.02 % | 0.23 % | 95.25 % | 1.38 % | 0.12 % | 0.001 % |

The trend of §5 continues (rays 95 %). 2 194 807 of the 2 594 741 sampled
spirals (85 %) are proven; the 63 961 undecided universes ran out of the
10 000-event budget without a run of growing returns. The event-driven
analyser exposes time scales no step-by-step simulation could reach: universe
`4-86460340171460` falls into a cycle of period 111 147 613 754 153 after 143
steps (about 3 600 events; `show 4-86460340171460`), and the latest proof of a
fixed point, ray or helix in the sample comes at step 8.4·10^16. Before the
spiral proofs, such orbits stopped as "spirals" after 8 growing returns; now
they run until they close. Certificates cost one event each, so these cases
take milliseconds.

Clamped model (create / destroy with positivity, §5):

| N | reduced universes | fixed | cycle | ray | helix | spiral |
|---|---|---|---|---|---|---|
| 1 | 3 | 75.0 % | 0 | 25.0 % | 0 | 0 |
| 2 | 81 | 62.89 % | 0 | 37.11 % | 0 | 0 |
| 3 | 19 683 | 48.30 % | 0.23 % | 51.34 % | 0.12 % | 0 |
| 4 | 43 046 721 | 35.63 % | 0.37 % | 63.54 % | 0.46 % | 0.01 % |

**Stable from every start?** (`stable` command, `results/stable-size*.txt`).
The fates above are for one start. A universe is *globally bounded* when every
start in `ℤ^N` ends fixed or in a cycle. Instability needs one witness: a sign
pattern whose drift keeps every sign (its own sign vector starts a ray), or any
start the analyser proves to be a helix or a spiral. Boundedness for all starts
is proven for `N ≤ 2` by the far-field argument of §2b.

| N | reduced universes | bounded from every start | unstable (proven) | open |
|---|---|---|---|---|
| 1 | 9 | 4 (44.4 %) | 5 | 0 |
| 2 | 6 561 | 856 (13.05 %) | 5 705 | 0 |
| 3 | 387 420 489 | — | 377 886 315 (97.54 %) | 9 534 174 (2.46 %) |

At size 3, 375 493 987 universes are unstable through an absorbing pattern,
2 184 888 through a start that becomes a helix and 207 440 through a proven
spiral. The open 2.46 % (200 038 symmetry classes) have no escaping start in
`[−2, 2]³` or among starts far out along every face; they are an upper bound on
the globally bounded share. The all-fire fate is a poor guide: 37 % of
universes whose all-fire start is a cycle, and 74 % of those whose all-fire
start is a fixed point, escape from some other start. Size 2 has spirals too,
just not from the all-fire start: `show 2-2800 --start -3,-3` gains 4 per lap.

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
(`x ↦ 2x + 12` above); `ρ(A) < 1` is only possible on a bounded set and ends
in a cycle. The lap lemma of §2 is this argument made rigorous: the map is affine
only on residue classes, so the classes are refined until it is. The global question — *which* loop a given start ends in — is the
hard part, and it is provably hard:

* The continuous-time counterpart, PCD systems, was analysed by Asarin, Maler
  and Pnueli, "Reachability analysis of dynamical systems having
  piecewise-constant derivatives", *Theoretical Computer Science* 138 (1995):
  reachability is decidable in two dimensions and **undecidable in three or
  more** (a 3-D PCD system can simulate a Turing machine). No closed-form
  predictor for `N ≥ 3` can exist in general; the certificates in §2 are the
  right tool, and the residual is the 192 spirals no lap on a line or plane
  explains.
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
   step cap, and every fate is proven (spirals by the lap lemma).
3. **No allocation in the inner loop.** The static musl build's `malloc` takes
   a global lock, so per-universe `Vec`/`HashMap` allocations made two threads
   *slower* than one. A reusable per-thread `Workspace` removed every
   allocation from the hot path: the exhaustive size-3 census went from
   527 s to 6.6 s on 8 threads (59 million universes per second; 9.4 s since
   the spiral proofs), and a
   base-3 odometer replaced per-index decoding.
4. **Symmetry**: another 48× available at N=3 (orbit representatives via
   `Reduced::canonical`), 384× at N=4. Not needed at N=3 (seconds), but the
   only route to an exhaustive N=4: 1.85·10^15 reduced universes take about
   a year at 59 M/s, 4.8·10^12 orbit representatives about a day — provided
   they are *generated* canonically (orderly generation) rather than filtered,
   since a canonical check costs 384 transforms. Until then N=4 is sampled:
   2·10^9 samples in about 23 minutes give every fate to ±0.001 %
   (`results/exact-size4-sampled.txt`).
5. **GPU**: the per-universe work is branchy integer logic with tiny state,
   which GPUs handle well for the 10-step heuristic (2^36 matrices would take
   minutes on the RTX 3070 instead of 40). For the exact analyser the
   unbounded loop and hashing map poorly; and the reductions above already
   beat a GPU brute force by orders of magnitude. If a GPU path is wanted, the
   right target is sampling at `N ≥ 4`: no CUDA toolkit is installed and the
   static musl build cannot `dlopen` the driver, so it needs `build-essential`
   and `nvidia-cuda-toolkit` first.

## 9. Basin patterns

`basinmap` colours every start state in a box by the attractor it reaches
(the fixed state, the smallest state of the cycle, the ray's drift, the helix
shift); `basinscan` does it for every canonical universe and scores the map.
Size 2 was scanned completely (873 canonical universes, box `[−12, 12]²`,
0.2 s); size 3 by sampling (2 990 canonical universes, box `[−4, 4]³`).
`docs/basin-atlas.html` shows the maps. Three families account for almost
everything:

1. **Conic sectors** — 55 % of size-2 universes. The fate depends only on
   the direction of the start: basins are unions of angular sectors bounded
   by lines through the origin. This is the continuous PCD picture, where the
   flow is scale-invariant; it survives discretisation whenever no coordinate
   lands exactly on zero in a way that matters.
2. **Stripe families (conservation laws)** — 40 % of size-2 universes have
   ten or more fixed points with real basins, arranged along a line, with
   parallel basins. Cause: a linear functional `L(v) = a·v` with `a·d(s) = 0`
   for every sign pattern the trajectory visits is conserved, so all starts on
   one level line `L = k` fall to the same attractor (universe 398: fixed points
   `(0, k)`, basins are the diagonals `v₀ + v₁ = k`). At size 3 the same
   mechanism indexes *cycles*: universe 1352401 has a family of cycles through
   `(0, 0, −k)`, one per level of the invariant, separated by rays.
3. **Checkerboards (parity)** — 8 of 873 size-2 universes. A coordinate with
   drift of magnitude 2 either lands exactly on zero (goes silent) or jumps
   over it (flips sign) depending on its parity, so `v₀ + v₁ mod 2` decides
   the fate. Universe 560 is a pure checkerboard: two cycles, each owning one
   colour of the board. Mixed forms (stripes with a checkerboard band near an
   axis, universes 38 and 398) come from the same overshoot effect.

Other structure visible in the atlas: whole orthants of **inert cells** (an
orthant whose drift is zero makes every state in it a fixed point; universe
89), and universes with up to 8 distinct ray directions (universe 5520), all
conic.

Consequences for "classify the universe from its matrix": a universe is a
*phase portrait*, and the portrait is determined by (a) the sector structure
of the drifts, (b) the conserved linear functionals of the drift rows, and
(c) the residues of the drift magnitudes. All three are properties of the
`pos`/`neg` rows alone, so they are the right invariants to compute
symbolically next.

## 10. Open questions worth attacking next

* Prove the 192 three-dimensional size-3 spirals: the lap lemma with three
  unknowns, where the invariant region is a polyhedral cone in 3-D.
* Orderly generation of `B_N` orbit representatives to make exact N=4 sampling
  unbiased at the orbit level.
* Prove global boundedness at N=3 (§2b with two large coordinates): which of
  the 2.46 % open universes are bounded from every start? Conserved linear
  functionals (§9) explain many of them.
* Relate the cycle-period histogram (heavy tail to 1714 at N=3) to loop return
  maps with `A` of finite order.
