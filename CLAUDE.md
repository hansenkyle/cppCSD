# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

cppCSD (`ldcsd`) implements method development for solving the Continuous
Slowing-Down (CSD) equation using Lewis and Miller's Second Moment Method
(SMM) with linear-discontinuous (LD) finite elements.

## Build

Requires CMake 4.0+ and the HDF5 C library (`sudo apt install libhdf5-dev`);
the Python extension additionally needs the Python headers
(`sudo apt install python3-dev`).
Everything else (Eigen, Boost.Multiprecision, yaml-cpp, CLI11, HighFive, doctest)
is fetched automatically via `FetchContent`.

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DLDCSD_BUILD_TESTS=ON
cmake --build build -j
```

This produces the `ldcsd` executable and, when `LDCSD_BUILD_TESTS=ON`
(default), the `unit_test` binary.

## Testing

```bash
cmake --build build -j --target unit_test
./build/unit_test
```

Tests use doctest. To run a single test case or suite, use doctest's CLI
filters:

```bash
./build/unit_test --test-case="example"
./build/unit_test --test-suite="example"
```

Unit test sources live in `test/unit/` and link against `ldcsd_core` (the
static library built from `src/module.cpp`), not the `ldcsd` executable
target.

### Coverage

```bash
scripts/coverage.sh
```

Configures a separate `build-coverage/` tree with `LDCSD_ENABLE_COVERAGE=ON`,
runs the tests, and generates an HTML report at
`build-coverage/coverage-report/index.html` via gcovr (requires
`pip install gcovr`). CI extracts a line-coverage percentage from this to
regenerate `badges/coverage.svg` on pushes to `main`.

## Python package

`python/` is the `ldcsd` package. It has two halves:

- a pure-Python reader for a run's `.h5` file (`python/ldcsd/results.py`:
  `ldcsd.read(path)`, plus `cell_average`, `spectrum`, `multigroup`
  helpers; per-material data is `MaterialData`);
- `ldcsd._core`, a nanobind extension (`python/src/bindings.cpp`) binding
  `InputDeck` (through its setters), `Material`, `SourceIteration`,
  `SecondMoment` and `configure_logging`, re-exported from `ldcsd`. The
  import is optional, so the reader works from a plain checkout.

Arrays cross the binding in the reader's layout -- corner fields
`[G, nx, 2 (up, down), 2 (L, R)]`, source `[G, M, nx, 2, 2]`, bc
`[G, 2 (up, down), M]` -- and everything returned is a read-only copy.

`pip install python/` builds the extension with scikit-build-core, which
runs the root `CMakeLists.txt` with `SKBUILD` set (PIC static deps, no
`ldcsd` executable, `_core` built with the same `-O3 -march=native`). For
development, install editable without build isolation; the extension then
rebuilds on import after C++ edits:

```bash
pip install scikit-build-core nanobind ninja
pip install --no-build-isolation -e "python/[test]"
pytest python
```

Reader-only (no C++ build), as the `python-tests` CI job does; tests that
need `ldcsd._core` skip themselves:

```bash
PYTHONPATH=python pytest python
```

## Formatting

Formatting is enforced by clang-format (config in `.clang-format`) and
checked in CI.

```bash
scripts/format.sh          # reformat src/, test/ and python/src/ in place
scripts/check-format.sh    # check only, nonzero exit if reformatting needed
```

Both the scripts and CI's format check cover `*.cpp` and `*.h` files
(CI also `*.hpp`).

## Architecture

- `ldcsd_core` — static library built from `src/module.cpp`; all
  reusable logic should live here so it's linked into both the `ldcsd`
  executable and `unit_test`.
- `ldcsd` executable (`src/main.cpp`) — thin entry point linking against
  `ldcsd_core`.
- `unit_test` — doctest-based test binary linking `ldcsd_core` +
  `doctest::doctest`; built only when `LDCSD_BUILD_TESTS=ON`.
- `_core` — the Python extension (`python/src/bindings.cpp`) linking
  `ldcsd_core`; built only by `pip install python/` (see Python package).
- Eigen is vendored via `FetchContent` and exposed as the `Eigen3::Eigen`
  interface target; link against it rather than assuming a system install.
- Release builds compile with `-O3 -march=native` (see `CMakeLists.txt`) —
  binaries are not portable across differing CPU microarchitectures.
  
## InputDeck: the trust boundary for input validation

`InputDeck` -- its YAML parser (`InputDeck::load()`/`read()`) and its
setters, which are also the only way the Python bindings change a deck --
is the single point where input legality is enforced. Once data lives inside an `InputDeck`, the rest of
the code trusts it completely -- no downstream code should re-validate
mesh sizes, cross-section shapes, quadrature ordering, etc.

A `Solver` should own an `InputDeck` (it's a small struct, cheap to hold
by value). Any method on `Solver` that lets a caller reconfigure part of
the problem (e.g. swap in a different angular quadrature) should just
call `InputDeck`'s own setters (`set_mesh`, `set_energy`, `set_angle`,
`set_materials`, `set_bc`, `set_source`) to change that data, rather than
mutating fields directly or re-deriving validation logic itself. Each
setter validates its own struct before storing it; checks that span
structs (e.g. xs covering `n_x` cells) run in `InputDeck::validate()`,
which `load()`/`read()` and the `Method` constructor both call. So no
other method needs to worry about whether the data it's handed is legal
-- that question is answered once, at the `InputDeck` boundary.

## Repository Interactions

- Claude should *never* create a pull request unless specifically asked to do so.
- Claude should provide a warning to the user if recent changes are commits differ 
from the idea behind the current branch, and suggest switching to a different 
branch or opening a new branch if appropriate.
- Claude should always ask before creating a new branch.
- Claude should not write Doxygen-formatted comments, but should still write descriptive comments before classes and methods (invisible to doxygen)
