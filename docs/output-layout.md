# Output directory layout

Conventions for organizing input decks and the output of solver runs on
disk. This describes the intended structure only — nothing here is
enforced by the code yet.

## Structure

Each input deck lives in its own directory, named after the deck. All
runs of that deck are collected under a `runs/` subdirectory, one folder
per run, named with a monotonic counter.

```
wafer/
    wafer.yaml
    runs/
        run_0001/
            info
            log
            out
        run_0002/
            info
            log
            out
        latest -> run_0002/
shield/
    shield.yaml
    runs/
        run_0001/
            info
            log
            out
        latest -> run_0001/
```

- **`<deck>/<deck>.yaml`** — the input deck, in the format read by
  `Parser::read` (see [../src/parser.h](../src/parser.h)).
- **`<deck>/runs/run_NNNN/`** — one directory per run. `NNNN` is a
  monotonically increasing, zero-padded counter local to that deck (not a
  timestamp), so directory listings always sort in run order regardless
  of when runs happened.
- **`<deck>/runs/latest`** — a symlink to the most recent run directory,
  updated on each run, so the newest results are always reachable at a
  fixed path without scanning `runs/`.

Each run directory holds exactly three files today, all describing the
*run itself* rather than the numerical method or its results — that
separation is deliberate: solver output formats can change independently
of how a run is logged and identified. Additional files (formatted
numerical results, plots, etc.) will be added under the same run
directory later.

## File purposes

### `info`

User-readable run metadata, written once at the start of the run.
Intended to make a run directory self-describing and reproducible on
its own, independent of the current state of the deck file (which may
have been edited since). Includes:

- A hash of the input deck's contents (so a run can always be tied back
  to the exact deck that produced it, even if `<deck>.yaml` changes
  later).
- Basic execution metadata: machine architecture and the date/time the
  run started.
- Any other run-identifying parameters worth surfacing to a human at a
  glance (e.g. solver options in effect for that run), as they're added.

### `log`

Runtime log messages emitted while the solver runs, at adjustable and
sortable verbosity. Intended for a future parsing/filtering script
(e.g. to pull out only warnings, or only messages above a given
verbosity level) rather than for reading top-to-bottom.

### `out`

A plain echo of the terminal output produced during execution (stdout/
stderr as seen by the user running the solver). This is a raw capture,
not a structured result file.

### `<deck>.h5`

The run's complete numerical record, named after the input deck (e.g.
`wafer/runs/run_0003/wafer.h5`), written once after the solve by
`writeH5`. It holds everything in `results.txt` and `residuals.txt`
unaveraged: every corner of every cell, every angle. Cell averages,
spectrum and multigroup slices are left for the reader to compute.

Corner-valued rows follow the solver's layout: row `4i+k` is cell `i`,
`k = (up L, up R, down L, down R)` (`reshape(nx, 2, 2)` in numpy), noted
on each dataset as a `corner_order` attribute. Shapes use `nx` cells, `G`
groups and `M` angles.

```
attrs: execution_datetime, method, ldcsd_version, n_groups, n_cells, n_angles
/input
    mesh/       x_boundary [nx+1], dx [nx], x_center [nx], material [nx] (names)
    energy/     E_boundary [G+1], dE [G]
    angle/      mu [M], w [M]
    bc                      [2G x M]   row 2g+k is group g, k = (up, down)
    source/     values [G x 4nx x M], q0 [G x 4nx], q1 [G x 4nx]
    materials/<name>/       total [G], S [G], S_b [G+1], scatter [G x G] (from, to)
/solution
    scalar_flux             [4nx x G]
    angular_flux            [G x 4nx x M]
    current, reconstructed_scalar       [4nx x G]   (SMM only)
    closures/{F, F+, F-, K+, K-, T+, T-} [G x 4nx]  (SMM only)
/convergence
    iterations [G], group_times [G]
    gNNN/       delta_l2, delta_linf, phi_l2, phi_linf  [iterations in group]
/residuals
    transport               [G x 4nx x M]
    second_moment           [G x 8nx]  (SMM only) row 8i+k is cell i, k = the 8 SM equations
```

## Non-goals (for now)

- No numerical results (fluxes, convergence data, etc.) live in `info`,
  `log`, or `out` — those three are purely about the run's
  parameters/provenance and its console/log output. Formatted result
  files will be added alongside them in the run directory once that
  format is designed.
- Directory/file naming and content here are a convention, not something
  the current codebase creates or reads.
