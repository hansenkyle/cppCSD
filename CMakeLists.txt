cmake_minimum_required(VERSION 4.0)
project(ldcsd LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(LDCSD_BUILD_TESTS "Build ldcsd unit tests" ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

include(FetchContent)

# --- Eigen (header-only, fetched at configure time) ---
# Eigen's own CMakeLists.txt runs its full build (test scripts, an old
# cmake_minimum_required, etc.) if add_subdirectory() is ever called on it.
# Pointing SOURCE_SUBDIR at a path with no CMakeLists.txt tells
# FetchContent_MakeAvailable to populate the source only, and skip
# add_subdirectory entirely -- no deprecated FetchContent_Populate call,
# and Eigen's CMakeLists.txt (and its stale minimum version) is never parsed.
FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        5.0.1
    SOURCE_SUBDIR  cmake-ignore   # deliberately nonexistent -- see comment above
)
FetchContent_MakeAvailable(eigen)

add_library(Eigen3::Eigen INTERFACE IMPORTED)
target_include_directories(Eigen3::Eigen INTERFACE ${eigen_SOURCE_DIR})

# --- Core library (everything except main.cpp) ---
# Both the ldcsd executable and unit_test link against this, so:
#   1. src/ headers can be included by bare name ("solver.h") from
#      anywhere that links ldcsd_core, tests included -- no relative
#      paths like "../../src/solver.h" needed.
#   2. Test files actually get solver.cpp's implementation linked in,
#      instead of only compiling the test files themselves.
add_library(ldcsd_core STATIC
    src/solver.cpp
)
target_include_directories(ldcsd_core PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(ldcsd_core PUBLIC Eigen3::Eigen)

# --- Executable ---
add_executable(ldcsd
    src/main.cpp
)
target_link_libraries(ldcsd PRIVATE ldcsd_core)

target_compile_options(ldcsd PRIVATE
    $<$<CONFIG:Release>:-O3 -march=native>
)

# --- Tests (doctest) ---
# doctest is also single-header and carries the same stale
# cmake_minimum_required as Eigen -- same fix: skip add_subdirectory
# and declare the target ourselves.
if(LDCSD_BUILD_TESTS)
    FetchContent_Declare(
        doctest
        GIT_REPOSITORY https://github.com/doctest/doctest.git
        GIT_TAG        v2.4.11
        SOURCE_SUBDIR  cmake-ignore   # deliberately nonexistent -- see Eigen comment above
    )
    FetchContent_MakeAvailable(doctest)

    add_library(doctest::doctest INTERFACE IMPORTED)
    target_include_directories(doctest::doctest INTERFACE ${doctest_SOURCE_DIR}/doctest)

    add_executable(unit_test
        test/unit/main.test.cpp
        test/unit/solver.test.cpp)
    target_link_libraries(unit_test PRIVATE ldcsd_core doctest::doctest)

    # --- Coverage instrumentation (opt-in, GCC/Clang gcov-compatible) ---
    # Off by default: instrumentation disables optimization and adds
    # profiling overhead, so this should never be on for a normal build.
    option(LDCSD_ENABLE_COVERAGE "Build unit_test with coverage instrumentation" OFF)
    if(LDCSD_ENABLE_COVERAGE)
        target_compile_options(unit_test PRIVATE --coverage -O0 -g)
        target_link_options(unit_test PRIVATE --coverage)
    endif()
endif()