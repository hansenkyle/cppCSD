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

# --- Executable ---
add_executable(ldcsd
    src/main.cpp
    src/solver.cpp
)
target_link_libraries(ldcsd PRIVATE Eigen3::Eigen)

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
    target_link_libraries(unit_test PRIVATE Eigen3::Eigen doctest::doctest)
endif()