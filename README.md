# Dynamics Engine

**An exact classifier for small discrete dynamical systems, written in dependency-free Rust.**

Every universe in this project follows one rule,

```
v(t+1) = v(t) + d(sign v(t))
```

an integer state that moves by a velocity chosen only by the signs of its coordinates.
The engine decides where each universe ends up: it stops, loops, escapes in a straight
line, escapes while repeating, or spirals outward. For three elements it settles all
**387,420,489** universes exactly, in about **6.5 seconds** on 8 cores.

**[Live site](https://mukos.github.io/dynamics-engine/)** · **[Theory notes](docs/theory.md)** · **[Results](results/)**

---

## Contents

- [Results](#results)
- [Quick start](#quick-start)
- [The model](#the-model)
- [How it works](#how-it-works)
- [Commands](#commands)
- [The website](#the-website)
- [Project layout](#project-layout)
- [References](#references)

## Results

Each universe starts from the zero state, fires every channel once, and is then followed
until its fate is proven.

| Elements | Universes | Fixed | Cycle | Ray | Helix | Spiral | Undecided |
|---|---|---|---|---|---|---|---|
| N = 2 | all 6,561 | 22.634 % | 1.341 % | 76.025 % | 0 | 0 | 0 |
| N = 3 | all 387,420,489 | 7.176 % | 0.951 % | 90.504 % | 1.297 % | 0.072 % | 0 |
| N = 4 | 2·10⁹ sampled | 1.740 % | 0.316 % | 95.061 % | 2.745 % | 0.137 % | 7,952 |

The N = 4 rows are estimates from a uniform sample of the 1.85·10¹⁵ universes, with 95 %
confidence half-widths below 0.001 percentage points.

Other findings:

- **Every three-element universe receives a fate.** Fixed, cycle, ray and helix come with
  certificates; spirals are the empirical remainder.
- **The longest cycle at N = 3 has period 1,510.** The longest transient plus period is 1,824 steps.
- **Symmetry cuts the space by about 48×.** Fates are invariant under signed permutations of
  the elements, giving 873 classes at N = 2 and 8,093,160 at N = 3.
- **Basins of attraction come in three families**: sectors (55 % of two-element universes),
  stripes from conserved linear functionals (40 %) and parity checkerboards (1 %).

Full tables, proofs and open questions are in [docs/theory.md](docs/theory.md).

## Quick start

The engine is a command-line tool. You need a stable Rust toolchain
([rustup.rs](https://rustup.rs)); there are no crate dependencies.

```bash
cargo install --git https://github.com/mukos/dynamics-engine
```

That puts `dynamics-engine` on your `PATH`. Look at one universe:

```bash
dynamics-engine show 2-613
```

```
reduced universe 613 (multiplicity 8, canonical: true):
  element 0: pos = [  -1   -1 ]  neg = [   0    1 ]
  element 1: pos = [   1    0 ]  neg = [  -1    0 ]
events:
  step      1  v = [  -1    0 ]  drift = [   0    1 ]
  step      2  v = [  -1    1 ]  drift = [   1    1 ]
  ...
  step     16  v = [   1    0 ]  drift = [  -1   -1 ]
fate: cycle of period 10, repeating from step 6
watch it: https://mukos.github.io/dynamics-engine/#u=2-613
```

A universe id is `N-index`, the same id the website uses in its `#u=` links, so any universe
from the site can be replayed in the terminal and the other way round. Then classify every
three-element universe:

```bash
dynamics-engine exact --size 3
```

To work on the code instead:

```bash
git clone https://github.com/mukos/dynamics-engine.git
cd dynamics-engine
cargo run --release -- show 3-27411247
cargo test --release
```

## The model

A universe has `N` elements. Element `e` has a positive channel and a negative channel, and a
`2N × 2N` matrix of zeros and ones says which channels feed which. At every step each active
channel pushes its row, each element collapses to its net value, and the sign of that value
decides which of its channels fires next. On the first step every channel fires.

**Reduction.** Only the difference between the two bits a row sends to an element matters. Each
pair collapses to a value in {−1, 0, +1}, so the `2^(4N²)` bit matrices reduce to `3^(2N²)`
universes, 177× fewer at N = 3. What remains is two `N × N` drift tables `d⁺` and `d⁻`, and the
rule above, where

```
d(s) = Σ_{s_f = +} d⁺_f + Σ_{s_f = −} d⁻_f
```

Inside a quadrant the signs do not change, so the velocity is constant and the orbit is a
straight line. Orbits are therefore chains of straight segments that bend only on the axes.

**The five fates.**

| Fate | What happens | Certificate |
|---|---|---|
| Fixed | the velocity reaches zero | `d(sign v) = 0` |
| Cycle | the state returns to itself | a repeated state |
| Ray | it escapes in a straight line | every coordinate moves away from zero |
| Helix | it escapes while repeating | a segment replays shifted by `D` with no sign change |
| Spiral | it winds outward | none: 8 consecutive growing returns to one sign pattern |

## How it works

- **Event-driven analysis.** Instead of stepping, the analyser jumps straight to the next sign
  change. A four-element universe that takes 805 million steps to stop
  (`dynamics-engine show 4-1337051325212753`) is settled within the 10,000-event budget.
- **Certificates, not step limits.** Cycles, rays and helices are proven. The helix test is the
  deterministic analogue of the Karp–Miller coverability argument for vector addition systems.
- **Symmetry.** Canonical forms under the hyperoctahedral group (order `2^N N!`) and Burnside
  counts of the orbits.
- **Allocation-free hot loop.** Static musl builds serialise threads on `malloc`. One reusable
  workspace per thread and an index odometer took the N = 3 census from 527 s to 6.5 s.
- **Reproducible scale.** Sampling uses a seeded SplitMix64; exhaustive runs shard across
  machines with `--offset` and `--stride` and merge with `scripts/merge-results.py`.
- **Tests that pin the maths.** Certificates are checked against brute-force simulation, the
  odometer against index decoding, and the original TypeScript histogram is kept as a
  regression test.

## Commands

Run `dynamics-engine` with no arguments for the full help.


| Command | What it does |
|---|---|
| `show <N-index>` | one universe's drift rows, event trajectory and fate |
| `exact --size N` | exact fate of every universe of one size; `--samples` for N ≥ 4 |
| `basins --size N` | fates over every start state in a box, per universe |
| `basinmap --size N --index I` | which attractor each start in a box reaches |
| `basinscan --size N` | basin maps and pattern features of every canonical universe |
| `orbits --size N` | symmetry orbit counts (Burnside) |
| `run --size N` | the original 10-step signature histogram (legacy) |

Add `--json` for machine-readable output and `--model clamped` for the create/destroy variant
with positivity. `dynamics-engine --help` lists every option.

<details>
<summary>The original 10-step heuristic and its quirks</summary>

`run` reproduces the project's first TypeScript version. After a fixed number of steps it gives
each channel a symbol and tallies the joined symbols as a *signature*:

| Symbol | Meaning |
|---|---|
| `\|` | the run ended before the step limit |
| `0` | the final value appeared earlier (bounded or periodic) |
| `1` | the positive channel never repeated a value |
| `-1` | the negative channel never repeated a value |

The comparison skips the step immediately before the last one, exactly as the original did, so
old and new histograms can be diffed. Pass `--fixed-cycle` to compare against every earlier step.
The heuristic misclassifies about 0.2 % of matrices at N = 3, which is why `exact` exists.

| N | Matrices | Exhaustive on 8 cores |
|---|---|---|
| 2 | 65,536 | instant |
| 3 | 2³⁶ ≈ 6.9·10¹⁰ | about 30 minutes |
| 4 | 2⁶⁴ ≈ 1.8·10¹⁹ | about 16,000 years, so use `--samples` |

</details>

## The website

[`docs/`](docs/) is the project's landing page, served by GitHub Pages. It animates real
universes in the browser and includes a JavaScript port of the analyser, so visitors can type
any universe index and see its fate.

```bash
python3 scripts/build-dynamics.py   # results/*.json + scripts/dynamics.template.html → docs/index.html
python3 scripts/serve-docs.py       # http://127.0.0.1:8766, caching disabled
```

## Project layout

```
src/exact.rs        reduced model, event-driven analysis, certificates, symmetry group
src/experiments.rs  censuses, basins, basin scans, the clamped variant
src/multiverse.rs   exhaustive and sampled runs of the original model
src/universe.rs     one universe of the original model and its 10-step classification
src/matrix.rs       index ↔ matrix decoding, random matrices
src/parallel.rs     fork-join helpers
src/rng.rs          SplitMix64
src/main.rs         command-line interface
tests/              certificate, decoding and regression tests
results/            computed censuses, basin scans and orbit counts
docs/               landing page and theory notes
scripts/            site builders, local server, shard merging
```

## References

- E. Asarin, O. Maler, A. Pnueli. *Reachability analysis of dynamical systems having
  piecewise-constant derivatives.* Theoretical Computer Science 138 (1995). Reachability is
  decidable in two dimensions and undecidable from three, which is why spirals stay empirical.
- L. Glass, S. Kauffman. *The logical analysis of continuous, non-linear biochemical control
  networks.* Journal of Theoretical Biology 39 (1973).
- R. M. Karp, R. E. Miller. *Parallel program schemata.* Journal of Computer and System
  Sciences 3 (1969).

## License

[MIT](LICENSE)
