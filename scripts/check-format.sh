#!/usr/bin/env bash
# Checks formatting without modifying files. Exits nonzero if anything
# would be reformatted. Used by CI; also safe to run locally.
set -euo pipefail

FILES=$(find src test -type f \( -name '*.cpp' -o -name '*.hpp' \))

FAILED=0
for f in $FILES; do
    if ! diff -q "$f" <(clang-format "$f") > /dev/null; then
        echo "Needs formatting: $f"
        FAILED=1
    fi
done

if [ "$FAILED" -ne 0 ]; then
    echo ""
    echo "Some files are not formatted. Run scripts/format.sh to fix."
    exit 1
fi

echo "All files are formatted correctly."