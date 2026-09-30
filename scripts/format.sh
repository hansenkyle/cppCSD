#!/usr/bin/env bash
# Formats all project source files in place.
set -euo pipefail

find src test python/src -type f \( -name '*.cpp' -o -name '*.h' \) \
    -print0 | xargs -0 clang-format -i

echo "Formatted all source files."
