# Python Bindings

`import ldcsd` gives an input-deck builder and the solvers (the compiled
`ldcsd._core` extension, [../python/src/bindings.cpp](../python/src/bindings.cpp))
alongside the `.h5` reader (`ldcsd.read`, see [output-layout.md](output-layout.md)).
A full example is [../python/examples/slab.py](../python/examples/slab.py).

## Install

Needs what the C++ build needs, plus the Python headers (`python3-dev`):

```bash
pip install scikit-build-core nanobind ninja
pip install --no-build-isolation -e "python/[test]"   # rebuilds on import after C++ edits
```

Without the extension (e.g. `PYTHONPATH=python`), `import ldcsd` still gives
the reader, just not the classes below. The extension is compiled with
`-march=native`, so build it on the machine that runs it.

## Building a deck

A deck is changed only through its setters. Each setter checks its own piece
and raises without changing the deck if that piece is invalid. Checks that
span pieces (e.g. regions covering every cell, bc/source shapes matching
`G`, `M`, `nx`) run in `deck.validate()`, which every solver constructor
calls, so set things in any order.

| Setter | Argument(s) | Checks |
|---|---|---|
| `set_mesh(x_boundary)` | `[nx+1]` | strictly ascending |
| `set_energy(E_boundary)` | `[G+1]` | strictly descending |
| `set_angle(mu, w)` | `[M]`, `[M]` | `mu` strictly ascending; `w` sums to 2 (within 1e-4, then rescaled) |
| `set_materials(materials, regions)` | `{name: Material}`, one name per cell | material sizes agree, values non-negative, every region names a material |
| `set_bc(bc)` | `[G, 2, M]` | |
| `set_source(source)` | `[G, M, nx, 2, 2]` | non-negative |

`ldcsd.Material(total, S, S_b, scatter)` takes `total [G]`, group-average
stopping power `S [G]`, stopping power at the group boundaries `S_b [G+1]`,
and `scatter [G, G]` indexed `(from, to)`. The dict key becomes its name.

Array arguments accept anything `numpy.asarray` does (lists, int arrays, ...).

`nx`, `G` and `M` aren't set directly; they come from the array sizes and
are readable as `deck.mesh.n_x`, `deck.energy.G`, `deck.angle.M`.

## Array layout

Arrays use the same layout as `ldcsd.read()`:

- Group 0 is the highest-energy group.
- A corner field is `[..., nx, 2, 2]`: axis `-2` is the group's energy edge
  `(up, down)`, axis `-1` the cell's spatial edge `(L, R)`.
- `source` is `[G, M, nx, 2, 2]`, a corner field per group and ordinate.
- `bc[g, (up, down), m]` is the incoming flux for ordinate `m` on the face it
  enters through: `mu > 0` at the left boundary, `mu < 0` at the right. The
  other face's value for that ordinate is unused.

## Reading a deck back

`deck.mesh.x_boundary`, `deck.energy.E_boundary`, `deck.angle.mu` / `.w`,
`deck.materials` (dict, in the order given), `deck.regions`, `deck.bc` and
`deck.source` are read-only copies: writing to one raises instead of
silently doing nothing, so the setters are the only way in.

## YAML

```python
deck = ldcsd.InputDeck.read("deck.yaml")   # same format as the ldcsd executable
deck.to_yaml("copy.yaml")                  # dense scattering; reads back exactly
```

## Solving

```python
solver = ldcsd.SecondMoment(deck)          # or ldcsd.SourceIteration; copies + validates the deck
solver.solve(1e-10, max_iterations=1000)   # relative convergence tolerance
solver.scalar_flux                         # [G, nx, 2, 2]
solver.angular_flux                        # [G, M, nx, 2, 2]
solver.current, solver.reconstructed_scalar   # SecondMoment only, [G, nx, 2, 2]
solver.write_h5("run.h5")                  # the full run, for ldcsd.read()
solver.write_results("results.txt")        # the ldcsd executable's results.txt
solver.write_residuals("residuals.txt")    # ... and residuals.txt (all writers replace the file)
```

The corner-field helpers from the reader work on these too, e.g.
`ldcsd.cell_average(solver.scalar_flux)` gives `[G, nx]`.

## Errors and logging

- Invalid data raises `RuntimeError` with the C++ message (e.g.
  `mesh.x_boundary must be strictly ascending`). So does `InputDeck.read`
  for a missing or malformed file.
- An array of the wrong rank or shape raises `ValueError` naming the
  argument, e.g. `expected bc [G, 2, M], got shape (3, 4)`.
- The C++ log is off unless you call
  `ldcsd.configure_logging("run.log", out_path=None)`. `out_path` sets the
  C++ terminal-output file (the `ldcsd` executable's `out.txt`), which the
  solvers don't currently write to.
