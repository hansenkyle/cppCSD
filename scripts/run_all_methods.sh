#!/usr/bin/env bash
# Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
#
# Funded by CARRE (https://carre-psaapiv.org/)
#
# Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
# or without modification are permitted provided that the terms of the license are met.

# Run ldcsd on every input deck under test/method
# This is how results are moved between machines

# Usage:
#   scripts/run_all_methods.sh [-n] [-b BUILD_DIR] [-d DECK_ROOT] [PATTERN]
#
#   -n          dry run: list the decks that would run, run nothing
#   -b DIR      build directory holding the ldcsd binary (default: build)
#   -d DIR      root to search for decks (default: test/method)
#   PATTERN     only run decks whose path contains PATTERN

set -uo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="$repo_root/build"
deck_root="$repo_root/test/method"
dry_run=0

while getopts ":nb:d:h" opt; do
    case "$opt" in
        n) dry_run=1 ;;
        b) build_dir="$OPTARG" ;;
        d) deck_root="$OPTARG" ;;
        h) sed -n '9,22p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) echo "unknown option -$OPTARG" >&2; exit 2 ;;
    esac
done
shift $((OPTIND - 1))
pattern="${1:-}"

binary="$build_dir/ldcsd"
if [[ $dry_run -eq 0 && ! -x "$binary" ]]; then
    echo "error: no ldcsd binary at '$binary'" >&2
    echo "       build first: cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j" >&2
    exit 1
fi

# "Input deck" = any *.yaml directly inside a directory; a directory may hold at
# most one. Exclude supporting data from runs/ and xs_data/
mapfile -t decks < <(find "$deck_root" -name '*.yaml' -not -path '*/runs/*' -not -path '*/xs_data/*' | sort)

if [[ -n "$pattern" ]]; then
    filtered=()
    for deck in "${decks[@]}"; do
        [[ "$deck" == *"$pattern"* ]] && filtered+=("$deck")
    done
    decks=("${filtered[@]+"${filtered[@]}"}")
fi

if [[ ${#decks[@]} -eq 0 ]]; then
    echo "no decks found under '$deck_root'${pattern:+ matching '$pattern'}" >&2
    exit 1
fi

echo "found ${#decks[@]} deck(s) under ${deck_root#"$repo_root"/}"

failed=0
for deck in "${decks[@]}"; do
    rel="${deck#"$repo_root"/}"
    if [[ $dry_run -eq 1 ]]; then
        echo "  would run: $rel"
        continue
    fi

    printf '  %-56s ' "$rel"
    start=$(date +%s.%N)
    if output=$("$binary" "$deck" 2>&1); then
        elapsed=$(awk "BEGIN {printf \"%.1f\", $(date +%s.%N) - $start}")
        run_dir="$(dirname "$deck")/runs/latest"
        echo "ok (${elapsed}s) -> ${run_dir#"$repo_root"/}"
    else
        echo "FAILED"
        echo "$output" | sed 's/^/      /'
        failed=$((failed + 1))
    fi
done

if [[ $failed -gt 0 ]]; then
    echo "$failed deck(s) failed" >&2
    exit 1
fi
