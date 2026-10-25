# Dynamics Engine

**An exact classifier for small discrete dynamical systems, written in dependency-free Rust.**

Every universe in this project follows one rule,

```
v(t+1) = v(t) + d(sign v(t))
```

an integer state that moves by a velocity chosen only by the signs of its coordinates.
The engine decides where each universe ends up: it stops, loops, escapes in a straight
line, escapes while repeating, or spirals outward. For three elements it settles all
**387,420,489** universes exactly, in under **10 seconds** on 8 cores, and proves the fate of
all but 192 of them.

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
| N = 3 | all 387,420,489 | 7.176 % | 0.952 % | 90.504 % | 1.297 % | 0.072 % | 0 |
| N = 4 | 2·10⁹ sampled | 1.740 % | 0.318 % | 95.061 % | 2.748 % | 0.130 % | 63,961 |

The N = 4 rows are estimates from a uniform sample of the 1.85·10¹⁵ universes, with 95 %
confidence half-widths below 0.001 percentage points. 85 % of the sampled N = 4 spirals are
proven; the undecided ones ran out of the 10,000-event budget.

Other findings:

- **Every three-element universe receives a fate, and all but 192 are proven.** Fixed, cycle,
  ray and helix have certificates; 277,152 of the 277,344 spirals are proven by a symbolic lap.
  The 192 left are spirals whose returns span all three dimensions.
- **Watching growth is not a proof.** The old rule (8 growing returns in a row) called 2,784
  universes spirals that later settle: 2,736 into cycles, 48 into rays. Universe 3-7295755
  grows for 725 steps, then cycles with period 104.
- **Stable from every start.** Every one- and two-element universe is decided for *all*
  starts in ℤᴺ: 856 of the 6,561 two-element universes are bounded from every start, with a
  far-field proof (`dynamics-engine stable 2-613`); the rest have a start that escapes. At
  N = 3, 97.54 % are proven unstable and 2.46 % have no escaping start found. A cycle from the
  all-fire start often hides instability: 37 % of those universes escape from another start,
  and two-element universes can spiral (`dynamics-engine stable 2-2800`).
- **At N = 2 stability has a closed form.** A two-element universe is bounded from every start
  exactly when no sign pattern is absorbing (drift nonzero and keeping every sign), except
  for the sign rotation `d = ±(sgn v₁, −sgn v₀)`, which spirals out by 1 per axis crossing.
- **Create/destroy is the same system unless values are clamped at zero.** Letting an element
  create and destroy, and letting values go negative, gives a slice of the original universes.
  Clamping at zero is a genuinely different system: fixed points take over, and at N = 3 the
  fate differs from the unclamped rule in 10 % of universes. Its spirals are proven by the
  same lap argument.
- **The longest cycle at N = 3 has period 1,714.** The longest transient plus period is 3,771 steps.
- **Symmetry cuts the space by about 48×.** Fates are invariant under signed permutations of
  the elements, giving 873 classes at N = 2 and 8,093,160 at N = 3.
- **Basins of attraction come in three families**: sectors (55 % of two-element universes),
  stripes from conserved linear functionals (40 %) and parity checkerboards (1 %).

Full tables, proofs and open questions are in [docs/theory.md](docs/theory.md).

## Quick start

The engine is a single command-line program, `dynamics-engine`.

### Download a prebuilt binary

[Releases](https://github.com/mukos/dynamics-engine/releases/latest) has archives for five
platforms, each holding the binary, this README and the license:

| Platform | Archive |
|---|---|
| Linux x86_64 (static) | [`dynamics-engine-x86_64-unknown-linux-musl.tar.gz`](https://github.com/mukos/dynamics-engine/releases/latest/download/dynamics-engine-x86_64-unknown-linux-musl.tar.gz) |
| Linux ARM64 (static) | [`dynamics-engine-aarch64-unknown-linux-musl.tar.gz`](https://github.com/mukos/dynamics-engine/releases/latest/download/dynamics-engine-aarch64-unknown-linux-musl.tar.gz) |
| macOS Apple silicon | [`dynamics-engine-aarch64-apple-darwin.tar.gz`](https://github.com/mukos/dynamics-engine/releases/latest/download/dynamics-engine-aarch64-apple-darwin.tar.gz) |
| macOS Intel | [`dynamics-engine-x86_64-apple-darwin.tar.gz`](https://github.com/mukos/dynamics-engine/releases/latest/download/dynamics-engine-x86_64-apple-darwin.tar.gz) |
| Windows x86_64 | [`dynamics-engine-x86_64-pc-windows-msvc.zip`](https://github.com/mukos/dynamics-engine/releases/latest/download/dynamics-engine-x86_64-pc-windows-msvc.zip) |

On Linux, for example:

```bash
curl -L https://github.com/mukos/dynamics-engine/releases/latest/download/dynamics-engine-x86_64-unknown-linux-musl.tar.gz | tar xz
./dynamics-engine show 2-613
```

The binaries are not code-signed. If macOS refuses to open one downloaded in a browser, run
`xattr -d com.apple.quarantine dynamics-engine`; on Windows, choose *More info → Run anyway*.
`SHA256SUMS` in each release lists the checksums.

### Or build it with Cargo

With a stable Rust toolchain ([rustup.rs](https://rustup.rs)); there are no crate dependencies:

```bash
cargo install --git https://github.com/mukos/dynamics-engine
```

### First commands

Look at one universe:

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

### Working on the code

```bash
git clone https://github.com/mukos/dynamics-engine.git
cd dynamics-engine
cargo run --release -- show 3-27411247
cargo test --release
```

Pushing a version tag (`git tag v2.1.0 && git push origin v2.1.0`) runs
[`.github/workflows/release.yml`](.github/workflows/release.yml), which tests, builds the five
archives and publishes them as a release.

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

| Fate | What happens | Proof |
|---|---|---|
| Fixed | the velocity reaches zero | `d(sign v) = 0` |
| Cycle | the state returns to itself | a repeated state |
| Ray | it escapes in a straight line | every coordinate moves away from zero |
| Helix | it escapes while repeating | a segment replays shifted by `D` with no sign change |
| Spiral | it winds outward | one lap run symbolically: returns lie on a line or plane and the lap map pushes them outward |

## How it works

- **Event-driven analysis.** Instead of stepping, the analyser jumps straight to the next sign
  change. A four-element universe that takes 805 million steps to stop
  (`dynamics-engine show 4-1337051325212753`) is settled within the 10,000-event budget.
- **Proofs, not step limits.** Every cycle, ray and helix verdict comes with a certificate, a
  checkable proof, rather than a guess after N steps. The helix test is the
  deterministic analogue of the Karp–Miller coverability argument for vector addition systems.
- **Spirals by symbolic execution.** The returns to one sign pattern lie on a line
  `p + x·q` (or a plane). The engine runs one lap with `x` unknown, splitting `x` into residue
  classes so every jump length stays affine, and gets a map such as `x ↦ 2x + 12`. If every
  branch pushes `x` outward beyond the current return, the spiral grows forever
  (`src/spiral.rs`, `dynamics-engine show 3-27411247`).
- **Symmetry.** Canonical forms under the hyperoctahedral group (order `2^N N!`) and Burnside
  counts of the orbits.
- **Allocation-free hot loop.** Static musl builds serialise threads on `malloc`. One reusable
  workspace per thread and an index odometer took the N = 3 census from 527 s to under 10 s.
- **Reproducible scale.** Sampling uses a seeded SplitMix64; exhaustive runs shard across
  machines with `--offset` and `--stride` and merge with `scripts/merge-results.py`.
- **Tests that pin the maths.** The proofs are checked against brute-force simulation, the
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
| `stable <N-index>` | is every start bounded? an escaping start, or a proof (N ≤ 2) |
| `stable --size N` | the same for every universe of one size (N ≤ 3) |
| `variants --size N` | create/destroy rules with and without clamping at zero, compared |
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
src/exact.rs        reduced model, event-driven analysis, proofs, symmetry group
src/spiral.rs       spiral proofs: one lap run symbolically on a line or plane of returns
src/global.rs       stability from every start: escape witnesses, far-field proof for N ≤ 2
src/experiments.rs  censuses, basins, basin scans, the clamped variant
src/multiverse.rs   exhaustive and sampled runs of the original model
src/universe.rs     one universe of the original model and its 10-step classification
src/matrix.rs       index ↔ matrix decoding, random matrices
src/parallel.rs     fork-join helpers
src/rng.rs          SplitMix64
src/main.rs         command-line interface
tests/              proof, decoding and regression tests
results/            computed censuses, basin scans and orbit counts
docs/               landing page and theory notes
scripts/            site builders, local server, shard merging
```

## References

- E. Asarin, O. Maler, A. Pnueli. *Reachability analysis of dynamical systems having
  piecewise-constant derivatives.* Theoretical Computer Science 138 (1995). Reachability is
  decidable in two dimensions and undecidable from three, which is why a few spirals still resist proof.
- L. Glass, S. Kauffman. *The logical analysis of continuous, non-linear biochemical control
  networks.* Journal of Theoretical Biology 39 (1973).
- R. M. Karp, R. E. Miller. *Parallel program schemata.* Journal of Computer and System
  Sciences 3 (1969).

## License

[MIT](LICENSE)
