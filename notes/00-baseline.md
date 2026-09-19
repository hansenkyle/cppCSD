# Stage 0 — Baseline

Branch: `ex/baseline`. Date: 2026-09-19.

Everything here describes the code **as found**, before any change in this
session. It exists so later stages can prove they didn't move the physics.

## What was captured

`scripts/run_all_methods.sh` runs every deck under `test/method/`; the six
that exist today all succeed. Their scalar flux was copied to
`notes/reference/*.results` (264 KB total, checked in — `runs/` itself is
gitignored, so this is the only in-repo copy of a known-good answer).

`scripts/compare_reference.py` compares a fresh run against those
references **numerically, not textually**: it pulls every float out of the
`Scalar Flux` block in document order and compares the sequences. That is
deliberate — Stage 2 changes the printed precision from 6 digits to 4 and
reshapes the blocks, so a text diff would go red for reasons that have
nothing to do with the answer. Default tolerance is `rtol = 2e-4`, which is
the floor set by 4-digit output.

```bash
scripts/run_all_methods.sh          # produce runs/latest for every deck
python3 scripts/compare_reference.py   # check them against Stage 0
```

Value counts pinned: 36 / 36 / 100 / 100 / 1920 / 14400 (= 4 corner values
× n_x × G in each case).

## Baseline numbers

| deck | cells | groups | ordinates | wall time | iterations (per group) |
|---|---|---|---|---|---|
| `3g_3x` | 3 | 3 | 6 | 0.0 s | 3 3 3 |
| `3g_3x/scaled` | 3 | 3 | 6 | 0.1 s | 20 18 16 |
| `5g_5x` | 5 | 5 | 6 | 0.0 s | 3 3 3 3 2 |
| `5g_5x/no_scatter` | 5 | 5 | 6 | 0.0 s | 2 2 2 2 2 |
| `12g_40x` | 40 | 12 | 6 | 2.8 s | 36 10 7 6 5 4 3 3 3 2 2 2 |
| `36g_100x` | 100 | 36 | 6 | 24.3 s | 12 10 8 7 7 6 6 5 5 5 5 4 4 4 4 4 3 3 3 3 3 3 3 2 3 3 2 2 2 2 2 2 2 2 2 2 |

Whole suite: **27 s**.

## Output files as found

| deck | `info` | `log` | `results` | `out` |
|---|---|---|---|---|
| `3g_3x` | 9.3 K | 15 K | 681 B | **0** |
| `5g_5x` | 22 K | 23 K | 1.7 K | **0** |
| `12g_40x` | 688 K | 128 K | 29 K | **0** |
| `36g_100x` | **6.1 M** | 226 K | 212 K | **0** |

### Problems this makes obvious

1. **`out` is always empty.** `Logger::print` is the documented "only stdout
   interface", and nothing calls it. All terminal output today comes from
   `log(..., echo=true)`, which bypasses `out` entirely. So the run
   directory has no record of what the user saw. *(Stage 2c)*

2. **`info` is dominated by redundancy.** 6.1 MB for a 36-group, 100-cell
   deck of a *single material*. `InputDeck::echo()` prints one 36×36
   scattering matrix per spatial cell — 100 identical copies — and one
   external-source table per (group, ordinate) pair — 216 tables. The data
   is per-cell because `InputDeck::Xs` is per-cell. *(Stage 2d for the echo,
   Stage 6 for the underlying storage)*

3. **`log` is residual spam.** 3122 lines for 36g/100x, and the great
   majority are the per-iteration `debug_max` term-by-term dump, which
   `sourceIterate` requests unconditionally (`calculateResiduals(..., true)`).
   Useful when hunting a bug, noise otherwise, and there is no switch.
   *(Stage 2b)*

4. **`results` holds only the scalar flux.** `main` calls
   `writeResults(path, phi, {})` — the angular flux argument is always
   empty, and `writeResiduals` is never called at all. The formatter
   supports both; the driver doesn't use them. *(Stage 2d)*

5. **No convergence record survives the run.** Iteration counts and norms
   exist only as log lines. *(Stage 1)*

## Solver behaviour worth writing down before changing it

These are things I had to work out by reading, and that a future change
could silently break.

- **Groups are solved in one downward pass, never revisited.** The group
  loop in `sourceIterate` runs `g = 0 .. G-1` and converges each group
  before moving on. That's correct for pure CSD + downscatter (no upscatter
  in the decks: every scattering matrix is upper-triangular in `from -> to`
  with `from <= to`), but it means the whole method silently assumes no
  upscatter. Nothing validates that assumption.

- **`phi_g` is *not* reset between groups.** Group `g` starts its iteration
  from group `g-1`'s converged scalar flux, and the first thing the loop
  does is `phi.col(g) = phi_g`. So the in-group scattering source on
  iteration 1 of group `g` is built from the previous group's answer. It's a
  warm start, and it works, but it also means the convergence test on the
  first iteration of each group compares against the *previous* group's
  norm.

- **The convergence test is absolute-difference vs. scaled norm, with no
  iteration cap.** `while (abs_diff_norm > phi_g.norm() * epsilon)` with
  `epsilon = 1e-8` hard-coded in `main`. A non-converging problem spins
  forever. `abs_diff_norm` is seeded to `1e10` so the body always runs once.

- **`transportSweep`'s two branches are near-identical.** The `mu > 0` and
  `mu < 0` cases duplicate ~25 lines each, differing only in loop bounds and
  which neighbour supplies the upwind value. This is the single biggest
  readability target in the file. *(Stage 7)*

- **`Kernel::solveDirect` rotates negative-`mu` problems** so the matrix is
  always assembled for `mu > 0`, then un-rotates the answer on return
  (`x({1,0,3,2})`). Easy to miss; anything touching the kernel has to
  respect it.

- **`sigma_sdEprime` is rebuilt per cell, per angle, per iteration.**
  `input_deck.xs.scatter[i].col(g).cwiseProduct(input_deck.energy.dE)` sits
  inside the innermost loop in `transportSweep`, and it depends only on
  `(i, g)`. It is recomputed `M ×` more often than necessary. *(Stage 7
  names it; a later stage can hoist it.)*

## Test/format state at baseline

- `./build/unit_test` — 91 cases, 233 assertions, all passing.
- `scripts/check-format.sh` — clean.
- Build: CMake 4.0, Release `-O3 -march=native`, deps fetched (Eigen 5.0.1,
  Boost 1.87 multiprecision, yaml-cpp 0.8.0, CLI11 2.7.2, doctest 2.4.11).
