# system

Exhaustive and sampled simulation of small *activation-matrix* systems, written
in Rust with no dependencies.

```bash
cargo run --release -- --size 2          # every 4x4 matrix (65 536 universes), instant
cargo run --release -- --size 3          # every 6x6 matrix (2^36 universes), ~30 min on 8 cores
cargo run --release -- --size 4 --samples 100000000 --seed 1   # 2^64 matrices: sample instead
cargo test
```

`cargo run --release -- --help` lists all options. Results are printed as a
histogram of *signatures* (see below), or as JSON with `--json`.

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
src/universe.rs    one universe: activation, evaluation, classification
src/matrix.rs      index <-> matrix decoding, random matrices
src/multiverse.rs  exhaustive / sampled runs, threading, progress, tallies
src/main.rs        command-line interface
src/rng.rs         SplitMix64 (deterministic sampling)
tests/             unit tests and the size-2 regression histogram
docs/entity-sketch.ts  earlier design sketch (Entity/Action/Event/Cause), not built
```

## Toolchain note

`.cargo/config.toml` builds a static musl binary linked with the `rust-lld`
bundled in the Rust toolchain, so no system C compiler is needed. On a machine
that has `cc`, delete that file to use the default target.
