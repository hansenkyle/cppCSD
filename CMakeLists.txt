cmake_minimum_required(VERSION 4.0)
project(ldcsd LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(LDCSD_BUILD_TESTS "Build ldcsd unit tests" ON)

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

# --- Core library, included by main and tests
add_library(ldcsd_core STATIC
    src/module.cpp
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
        test/unit/module.test.cpp)
    target_link_libraries(unit_test PRIVATE ldcsd_core doctest::doctest)

    # Coverage instrumentation
    option(LDCSD_ENABLE_COVERAGE "Build unit_test with coverage instrumentation" OFF)
    if(LDCSD_ENABLE_COVERAGE)
        target_compile_options(unit_test PRIVATE --coverage -O0 -g)
        target_link_options(unit_test PRIVATE --coverage)
    endif()
endif()
