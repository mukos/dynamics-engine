# system

Exhaustive and sampled simulation of small *activation-matrix* systems, written
in Rust with no dependencies.

```bash
cargo run --release -- --size 2          # every 4x4 matrix (65 536 universes), instant
cargo run --release -- --size 3          # every 6x6 matrix (2^36 universes), ~40 min on 8 cores
cargo run --release -- exact --size 4 --samples 2000000000   # size 4: sampled census, ~3 min
cargo run --release -- --size 4 --samples 100000000 --seed 1   # 2^64 matrices: sample instead
cargo run --release -- exact --size 3    # exact fates of all 3^18 reduced universes, ~7 s
cargo run --release -- basins --size 2   # exact fates over all start states per universe
cargo run --release -- orbits --size 4   # symmetry orbit counts
cargo run --release -- show --size 3 --index 27411247   # one universe's event trajectory
cargo test
```

`cargo run --release -- --help` lists all commands and options. `run` (the
default) prints the original 10-step *signature* histogram (see below) or JSON
with `--json`. `exact`, `basins`, `orbits` and `show` are the research tools
described in [docs/theory.md](docs/theory.md): the model reduces to
`v(t+1) = v(t) + d(sign v(t))`, whose fates (fixed, cycle, ray, helix, spiral)
can be certified without a step limit. Computed results live in `results/`.

## The model

A universe has `N` elements. Element `e` has two **channels**: `2e` (positive,
holds `count ≥ 0`) and `2e+1` (negative, holds `antiCount ≤ 0`). A universe is
defined by a `2N × 2N` 0/1 **activation matrix**: `M[r][c] = 1` means "when
channel `r` fires, channel `c` receives one unit".

Every step:

1. **Activation.** Each active channel `r` pushes its row: a 1 in column `c`
   adds `+1` to channel `c` if `c` is even, `−1` if `c` is odd.
2. **Evaluation.** Each element collapses to its net value `count + antiCount`.
   Positive net → keeps `count`, clears `antiCount`, and its positive channel
   fires next step. Negative net → the mirror. Zero → both cleared, nothing fires.

A run stops after `--steps` steps (default 10), when nothing fires, or when the
firing channels add nothing. All channels fire on the first step.

### Signatures

After a run, each channel gets one symbol, and the universe's signature is the
symbols of all channels joined by commas (`count,anti` per element):

| symbol | meaning |
|--------|---------|
| `\|`   | the run ended before the step limit (dead or converged) |
| `0`    | the final value already appeared earlier (bounded / periodic) |
| `1`    | positive channel never repeated a value (monotone growth) |
| `-1`   | negative channel never repeated a value (monotone decline) |

The multiverse tallies how many matrices produce each signature.

### Scaling

The number of matrices is `2^(4N²)`; that exponent, not the language, is the
cost.

| N | channels | matrices | exhaustive on 8 cores (~36 M universes/s) |
|---|----------|----------|--------------------------------------------|
| 1 | 2 | 16 | instant |
| 2 | 4 | 65 536 | instant |
| 3 | 6 | 2^36 ≈ 6.9·10^10 | ~30 min |
| 4 | 8 | 2^64 ≈ 1.8·10^19 | ~16 000 years — use `--samples` |

For `N ≥ 4` use `--samples <n> --seed <s>` (reproducible random matrices).
Exhaustive runs can be sharded across machines with `--offset i --stride k`
(`i` in `0..k`); sum the JSON outputs to merge.

## What is known (short version)

* Only the differences `M[r][2e] − M[r][2e+1]` matter: `2^(4N²)` matrices
  collapse to `3^(2N²)` reduced universes (177× fewer at N=3).
* Fates can be proven exactly: at N=3, 99.9 % of reduced universes are certified
  fixed / cycle / ray / helix; the rest are empirically diverging spirals.
* Fates are invariant under signed permutations of elements (group of order
  `2^N N!`); the 10-step heuristic misclassifies ~0.2 % of matrices at N=3.
* The fate depends on the start state for ~75 % of universes; only ~13 % (N=2)
  and ~2.4 % (N=3) are bounded from every start.
* Full analysis, tables, literature and open questions: [docs/theory.md](docs/theory.md).

## Known quirks (kept for compatibility)

* **Classifier skips the previous step.** The final snapshot is compared with
  every earlier snapshot *except* the one immediately before it. A channel that
  changed at every step but the last is therefore still `1`/`-1`. This
  reproduces the original TypeScript output (`tests/multiverse.rs` pins the
  size-2 histogram). Pass `--fixed-cycle` to compare against all earlier steps.
* **Signature strings** use the original format, so old and new results can be
  diffed directly.

## Layout

```
src/universe.rs    one universe: activation, evaluation, 10-step classification
src/matrix.rs      index <-> matrix decoding, random matrices
src/multiverse.rs  exhaustive / sampled runs, threading, progress, tallies
src/exact.rs       reduced model, event-driven exact analysis, certificates, symmetry group
src/experiments.rs censuses (exact fates vs heuristic), basins, clamped variant
src/parallel.rs    fork-join helper
src/main.rs        command-line interface (run / exact / basins / orbits / show)
src/rng.rs         SplitMix64 (deterministic sampling)
tests/             unit tests, the size-2 regression histogram, certificate tests
results/           computed histograms and censuses
docs/theory.md     mathematical analysis and findings
docs/entity-sketch.ts  earlier design sketch (Entity/Action/Event/Cause), not built
```

## Toolchain note

`.cargo/config.toml` builds a static musl binary linked with the `rust-lld`
bundled in the Rust toolchain, so no system C compiler is needed. On a machine
that has `cc`, delete that file to use the default target.
