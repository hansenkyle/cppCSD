# Python bindings plan

Scope: Python bindings as an input-deck builder, plus the solvers. `InputDeck`
(and its struct members) and `SourceIteration`/`SecondMoment` are bound.

## Decisions

| Area | Choice |
|---|---|
| Library | nanobind |
| Build | scikit-build-core: `pip install -e python/` builds via the root CMakeLists |
| Namespace | compiled `ldcsd._core`, re-exported from `ldcsd` (optional import) |
| C++ glue | `python/src/bindings.cpp` |
| Validation | validating `InputDeck` setters in C++; fields read-only in Python; full `validate()` on solver construction |
| `n_x`/`G`/`M` | derived by setters, read-only |
| Array layout | reader layout: `source [G, M, nx, 2, 2]`, `bc [G, 2, M]` |
| `read()` | raises with the real error message |
| Materials | by name: dict + regions list |
| Solvers | construct, `solve`, `write_h5`, in-memory results as numpy |
| Logging | optional `configure_logging(log, out=None)` |
| Extras | `to_yaml`; no derived fields exposed; GIL held during solve; `-march=native`, local builds only |
| Name clash | builder gets `ldcsd.Material`; reader dataclass becomes `MaterialData` |

## Target API

```python
deck = ldcsd.InputDeck()
deck.set_mesh(x_boundary)                    # n_x derived
deck.set_energy(E_boundary)                  # G derived
deck.set_angle(mu, w)                        # M derived
deck.set_materials({"water": ldcsd.Material(total=..., S=..., S_b=..., scatter=...)},
                   regions=["water", "water", "lead"])
deck.set_bc(bc)                              # [G, 2, M]
deck.set_source(q)                           # [G, M, nx, 2, 2]
deck.to_yaml("deck.yaml")
deck = ldcsd.InputDeck.read("deck.yaml")     # raises RuntimeError with the real message

s = ldcsd.SecondMoment(deck)                 # full validate() runs here
s.solve(1e-12, max_iterations=1000)
s.scalar_flux, s.angular_flux                # + s.current, s.reconstructed_scalar (SMM)
s.write_h5("run.h5")
```

Getters return read-only copies (`deck.mesh.x_boundary`, `deck.energy.G`,
`deck.materials`, `deck.regions`, `deck.bc`, `deck.source`), so nothing hands
out a live numpy view that bypasses validation or silently loses writes.

## Phase 1 -- C++ core changes (`src/`)

1. `InputDeck` setters: `set_mesh`, `set_energy`, `set_angle`,
   `set_materials(vector<Material>, vector<string> regions)`, `set_bc`,
   `set_source`. Each fills its struct, derives the count, runs that struct's
   own `validate()`. `read()` is rewritten on top of them (name->index mapping
   moves from `read()` into `set_materials`). Fields stay public in C++.
2. Split `read()` into throwing `load(path)` + the existing `int read(path)`
   wrapper for `main`. Python binds `load` as `read`.
3. `InputDeck::write(path)` via `YAML::Emitter` (dense scattering, named-corner
   source).
4. `Method` constructor calls `input_deck.validate()`. Fix: `SourceIteration`
   builds `transport_operator` from the unvalidated ctor argument; use the
   validated member instead (same check for `SecondMoment`).
5. Move `make_timestamp()` from `main.cpp` into `ldcsd_core` so `write_h5` can
   default it.
6. Existing doctests pass; add a read -> write -> read round-trip doctest.

## Phase 2 -- build plumbing

- `python/pyproject.toml`: `scikit_build_core.build` backend, requires
  `scikit-build-core` + `nanobind`, `cmake.source-dir = ".."`,
  `cmake.define.LDCSD_BUILD_TESTS = "OFF"`, `wheel.packages = ["ldcsd"]`.
- Root `CMakeLists.txt`, under `if(SKBUILD)`: `CMAKE_POSITION_INDEPENDENT_CODE ON`
  before FetchContent (yaml-cpp, Boost's compiled parts and `ldcsd_core` are
  static and get linked into a shared `.so`); find Python + nanobind;
  `nanobind_add_module(_core python/src/bindings.cpp)` linking `ldcsd_core`
  with the same `-O3 -march=native`; `install(TARGETS _core DESTINATION ldcsd)`.
  Skip the `ldcsd` executable.
- Editable rebuilds: `pip install -e python/ --no-build-isolation` with
  `editable.rebuild = true`.

## Phase 3 -- `python/src/bindings.cpp`

- `Material`: kwargs constructor (`total`, `S`, `S_b`, `scatter`), read-only fields.
- `Mesh`, `Energy`, `Angle`: read-only (inputs + counts only).
- `InputDeck`: ctor, setters, `read` (static, raises), `to_yaml`, `bc`/`source`
  getters (native layout), `materials`/`regions` as dict/list.
- `Method` base -> `SourceIteration`, `SecondMoment`: ctor,
  `solve(epsilon, max_iterations=1000)`, `write_h5(path, timestamp=...)`, raw
  `solution` fields.
- `configure_logging(log_path, out_path=None)` -> `Logger::configure` /
  `Terminal::configure`.

## Phase 4 -- Python layer (`python/ldcsd/`)

- `deck.py`: thin subclasses of the `_core` classes doing reader <-> native
  layout conversion, reusing `results._corners` + `swapaxes`/`reshape`.
- `__init__.py`: `try: from .deck import ...` / `except ImportError: pass`.
- Rename `results.Material` -> `MaterialData`.

## Phase 5 -- tests, CI, docs

- `python/tests/test_deck.py`: Python-built deck == `InputDeck.read(sample_input.yaml)`
  field-by-field; each setter raises on bad input; `to_yaml` round trip.
- `python/tests/test_solve.py`: tiny deck, both methods, in-memory results ==
  `ldcsd.read(write_h5(...))`. Both use `pytest.importorskip("ldcsd._core")`.
- CI: existing `python-tests` switches to `PYTHONPATH=python pytest python` (no
  C++); new `python-bindings` job (CMake 4, libhdf5-dev,
  `pip install python/[test]`, `pytest python`).
- `format.sh`, `check-format.sh`, CI format step also cover `python/src`.
- CLAUDE.md: fix the InputDeck section (setters are the mutation path;
  cross-struct checks at `validate()`/solver construction); add a Python
  bindings build section. `docs/input-deck.md`: mention `to_yaml`.

## Deferred / defaults

- Scattering in Python is 0-indexed dense `[from, to]`; sparse 1-indexed form
  stays YAML-only.
- `Terminal::Print` writes C++ `std::cout`; in Jupyter it lands in the kernel
  terminal, not the cell.
- Three copies of `timestamp()` (main, Logger, Terminal); only main's moves.
- No convergence history / closures / residuals in memory (they're in the h5).
