# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

cppCSD (`ldcsd`) implements method development for solving the Continuous
Slowing-Down (CSD) equation using Lewis and Miller's Second Moment Method
(SMM) with linear-discontinuous (LD) finite elements.

## Build

Requires CMake 4.0+. Dependencies (Eigen, doctest) are fetched automatically
via `FetchContent` — no manual dependency install needed.

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

## Formatting

Formatting is enforced by clang-format (config in `.clang-format`) and
checked in CI.

```bash
scripts/format.sh          # reformat src/ and test/ in place
scripts/check-format.sh    # check only, nonzero exit if reformatting needed
```

Note: these scripts currently only glob `*.cpp`/`*.hpp` files — header files
using the `.h` extension (e.g. `src/module.h`, `src/input_deck.h`) are not
covered.

## Numerical data conventions

Eigen handles all actual linear algebra, but bare Eigen objects should not be the
public-facing representation of domain quantities — wrap them in named types so
the code says what a value *means*, not just its numeric shape.

- Use Eigen types directly wherever code is actually doing linear algebra: local
  element matrix/vector assembly, the global system matrix/RHS, and the linear
  solve itself.
- Domain quantities that are numerically vector-shaped but not being
  algebraically manipulated where they're used (e.g. a cell's corner values of
  scalar flux, current, etc.) should be their own small named type, not a bare
  `Eigen::VectorXd`/`Vector4d` passed around by convention.
- Compose, don't inherit — hold an Eigen fixed-size type (e.g. `Vector4d`) as a
  private member. This is free: fixed-size Eigen vectors are stack-allocated
  with no overhead versus a `std::array`, so wrapping costs nothing.
- Convert into Eigen at the point code actually becomes linear algebra (e.g.
  inside assembly), not earlier, and not as the domain type's only interface.

## Architecture

- `ldcsd_core` — static library built from `src/module.cpp`; all
  reusable logic should live here so it's linked into both the `ldcsd`
  executable and `unit_test`.
- `ldcsd` executable (`src/main.cpp`) — thin entry point linking against
  `ldcsd_core`.
- `unit_test` — doctest-based test binary linking `ldcsd_core` +
  `doctest::doctest`; built only when `LDCSD_BUILD_TESTS=ON`.
- Eigen is vendored via `FetchContent` and exposed as the `Eigen3::Eigen`
  interface target; link against it rather than assuming a system install.
- Release builds compile with `-O3 -march=native` (see `CMakeLists.txt`) —
  binaries are not portable across differing CPU microarchitectures.
  
## Repository Interactions

- Claude should *never* create a pull request unless specifically asked to do so.
- Claude should provide a warning to the user if recent changes are commits differ 
from the idea behind the current branch, and suggest switching to a different 
branch or opening a new branch if appropriate.
- Claude should always ask before creating a new branch.
- Claude should not write Doxygen-formatted comments, but should still write descriptive comments before classes and methods (invisible to doxygen)
