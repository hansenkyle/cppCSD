# Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
#
# Funded by CARRE (https://carre-psaapiv.org/)
#
# Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
# or without modification are permitted provided that the terms of the license are met.

cmake_minimum_required(VERSION 4.0)
project(ldcsd LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(LDCSD_BUILD_TESTS "Build ldcsd unit tests" ON)
option(LDCSD_ENABLE_DEBUG_LOGGING "Compile in Logger trace/debug messages" OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

include(FetchContent)

# Eigen and other dependencies; skip eigen tests
FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        5.0.1
    SOURCE_SUBDIR  cmake-ignore   # deliberately nonexistent -- see comment above
)
FetchContent_MakeAvailable(eigen)

add_library(Eigen3::Eigen INTERFACE IMPORTED)
target_include_directories(Eigen3::Eigen INTERFACE ${eigen_SOURCE_DIR})

# Boost.Multiprecision (cpp_bin_float backend), used for evaluating
# residuals in extended precision during verification. Fetched from Boost's
# official CMake-ready release archive (all libraries pre-populated, no git
# submodule step) rather than the boostorg/boost git superproject, which
# needs submodules already checked out before its own CMake can run.
# BOOST_INCLUDE_LIBRARIES limits what actually gets configured/built to
# multiprecision and its dependencies (a couple of which -- container,
# random -- are compiled, not header-only, but small and fast to build).
set(BOOST_INCLUDE_LIBRARIES multiprecision)
set(BOOST_ENABLE_CMAKE ON)
FetchContent_Declare(
    Boost
    URL https://github.com/boostorg/boost/releases/download/boost-1.87.0/boost-1.87.0-cmake.tar.xz
    URL_HASH SHA256=7da75f171837577a52bbf217e17f8ea576c7c246e4594d617bfde7fafd408be5
)
FetchContent_MakeAvailable(Boost)

# yaml-cpp, used by InputDeck to load input decks
set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "" FORCE)
set(YAML_BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
    yaml-cpp
    GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
    GIT_TAG        0.8.0
)
# yaml-cpp 0.8.0's CMakeLists predates CMake's policy-version enforcement.
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
FetchContent_MakeAvailable(yaml-cpp)
unset(CMAKE_POLICY_VERSION_MINIMUM)

# CLI11, used by main to parse command-line arguments. Header-only; only
# cli.cpp's implementation includes it (cli.h doesn't), so it's linked
# PRIVATE below rather than alongside ldcsd_core's other, publicly-exposed
# dependencies.
FetchContent_Declare(
    CLI11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG        v2.7.2
)
FetchContent_MakeAvailable(CLI11)

# --- Core library, included by main and tests
add_library(ldcsd_core STATIC
    src/convergence.cpp
    src/input_deck.cpp
    src/output_manager.cpp
    src/logger.cpp
    src/output.cpp
    src/solver.cpp
    src/solver_formatter.cpp
    src/cli.cpp
)
target_include_directories(ldcsd_core PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(ldcsd_core PUBLIC Eigen3::Eigen yaml-cpp::yaml-cpp Boost::multiprecision)
target_link_libraries(ldcsd_core PRIVATE CLI11::CLI11)
target_compile_options(ldcsd_core PRIVATE
    $<$<CONFIG:Release>:-O3 -march=native>
)
if(LDCSD_ENABLE_DEBUG_LOGGING)
    target_compile_definitions(ldcsd_core PUBLIC LDCSD_ENABLE_DEBUG_LOGGING)
endif()

# --- Executable ---
add_executable(ldcsd
    src/main.cpp
)
target_link_libraries(ldcsd PRIVATE ldcsd_core)

target_compile_options(ldcsd PRIVATE
    $<$<CONFIG:Release>:-O3 -march=native>
)

# Unit tests using doctest
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
        test/unit/convergence.test.cpp
        test/unit/input_deck.test.cpp
        test/unit/output_manager.test.cpp
        test/unit/logger.test.cpp
        test/unit/output.test.cpp
        test/unit/solver.test.cpp
        test/unit/solver_formatter.test.cpp
        test/unit/cli.test.cpp)
    target_link_libraries(unit_test PRIVATE ldcsd_core doctest::doctest)
    target_compile_definitions(unit_test PRIVATE
        TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/test/unit/data")
    # Must match ldcsd_core's ISA flags: Eigen fixed-size types (e.g. Kernel's Matrix4d) get
    # different vectorization alignment depending on -march, and calling across an alignment
    # mismatch between this target and ldcsd_core.a segfaults intermittently.
    target_compile_options(unit_test PRIVATE
        $<$<CONFIG:Release>:-O3 -march=native>)

    # Coverage instrumentation
    option(LDCSD_ENABLE_COVERAGE "Build unit_test with coverage instrumentation" OFF)
    if(LDCSD_ENABLE_COVERAGE)
        target_compile_options(unit_test PRIVATE --coverage -O0 -g)
        target_link_options(unit_test PRIVATE --coverage)
    endif()
endif()
