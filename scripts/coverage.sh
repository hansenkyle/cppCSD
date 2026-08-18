#!/usr/bin/env bash
# Builds with coverage instrumentation, runs the test suite, and generates
# an HTML coverage report at build-coverage/coverage-report/index.html
# using gcovr (a thin, pure-Python wrapper around gcov -- no lcov/perl).
set -euo pipefail

BUILD_DIR="build-coverage"

command -v gcovr >/dev/null || { echo "gcovr not found -- install with: pip install gcovr"; exit 1; }

cmake -B "$BUILD_DIR" \
    -DCMAKE_BUILD_TYPE=Debug \
    -DLDCSD_BUILD_TESTS=ON \
    -DLDCSD_ENABLE_COVERAGE=ON

cmake --build "$BUILD_DIR" -j --target unit_test

# Running the tests populates .gcda files alongside the .gcno files
# already produced by the build.
"$BUILD_DIR/unit_test"

mkdir -p "$BUILD_DIR/coverage-report"

gcovr \
    --root . \
    --object-directory "$BUILD_DIR" \
    --exclude '.*/_deps/.*' \
    --exclude '.*/test/.*' \
    --html --html-details -o "$BUILD_DIR/coverage-report/index.html" \
    --print-summary

echo ""
echo "Full HTML report: $BUILD_DIR/coverage-report/index.html"
