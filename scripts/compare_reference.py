#!/usr/bin/env python3
# Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
#
# Funded by CARRE (https://carre-psaapiv.org/)
#
# Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
# or without modification are permitted provided that the terms of the license are met.

"""Compares a run's scalar flux against the Stage 0 reference results in
notes/reference/.

The point is to survive formatting changes: rather than diffing text, this
pulls every floating-point number out of the "Scalar Flux" block in document
order and compares the two sequences numerically. Ordering (group-major,
then the up/down row pair, then spatial cells) is what's being pinned, not
the layout, the column widths, or the printed precision.

Because the output format is moving to 4 significant decimal digits, the
default tolerance is 2e-4 relative -- tight enough to catch a real change in
the physics, loose enough to ignore a precision change.

Usage:
    scripts/compare_reference.py                  # re-run nothing, compare runs/latest of every deck
    scripts/compare_reference.py --update         # overwrite the references with current results
    scripts/compare_reference.py --rtol 1e-6      # stricter
"""

import argparse
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
REFERENCE_DIR = REPO / "notes" / "reference"
DECK_ROOT = REPO / "test" / "method"

FLOAT = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


def decks() -> list[Path]:
    """Every input deck under test/method/, same rule as run_all_methods.sh."""
    return sorted(
        p
        for p in DECK_ROOT.rglob("*.yaml")
        if "runs" not in p.parts and "xs_data" not in p.parts
    )


def reference_name(deck: Path) -> str:
    return str(deck.parent.relative_to(DECK_ROOT)).replace("/", "_")


def scalar_flux_numbers(text: str) -> list[float]:
    """Every number in the scalar-flux section, in document order.

    Starts at the 'Scalar Flux' banner and stops at the next banner that
    isn't part of it, so angular flux / residual blocks added later don't
    silently join the comparison. Label text ('group 0:', 'cells 0..2') is
    stripped first so only data numbers survive.
    """
    lines = text.splitlines()
    try:
        start = next(i for i, line in enumerate(lines) if "Scalar Flux" in line)
    except StopIteration:
        return []

    numbers: list[float] = []
    for line in lines[start + 1 :]:
        if line.startswith("--- ") and "Scalar Flux" not in line:
            break
        # Drop anything that labels rather than reports: the "columns are
        # spatial cells 0..N" caption, and "group N:" / "ordinate N" labels.
        if ":" in line or "columns are" in line or not line.strip():
            continue
        numbers.extend(float(m.group()) for m in FLOAT.finditer(line))
    return numbers


def compare(name: str, reference: list[float], current: list[float], rtol: float, atol: float):
    """Returns (ok, message)."""
    if len(reference) != len(current):
        return False, f"value count changed: reference {len(reference)}, current {len(current)}"

    worst_rel, worst_at = 0.0, -1
    failures = 0
    for i, (r, c) in enumerate(zip(reference, current)):
        if abs(r - c) <= atol + rtol * abs(r):
            continue
        failures += 1
        rel = abs(r - c) / max(abs(r), 1e-300)
        if rel > worst_rel:
            worst_rel, worst_at = rel, i

    if failures == 0:
        return True, f"{len(current)} values match (rtol={rtol:g})"
    return False, (
        f"{failures}/{len(current)} values differ; worst relative difference "
        f"{worst_rel:.3e} at index {worst_at} "
        f"(reference {reference[worst_at]:.8e}, current {current[worst_at]:.8e})"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--rtol", type=float, default=2e-4, help="relative tolerance (default 2e-4)")
    parser.add_argument("--atol", type=float, default=1e-12, help="absolute floor (default 1e-12)")
    parser.add_argument("--update", action="store_true", help="overwrite references with current results")
    args = parser.parse_args()

    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    failures = 0

    for deck in decks():
        name = reference_name(deck)
        results = deck.parent / "runs" / "latest" / "results"
        reference_path = REFERENCE_DIR / f"{name}.results"

        if not results.exists():
            print(f"  {name:<28} SKIP (no runs/latest/results -- run scripts/run_all_methods.sh)")
            continue

        if args.update:
            reference_path.write_text(results.read_text())
            print(f"  {name:<28} updated")
            continue

        if not reference_path.exists():
            print(f"  {name:<28} SKIP (no reference captured)")
            continue

        ok, message = compare(
            name,
            scalar_flux_numbers(reference_path.read_text()),
            scalar_flux_numbers(results.read_text()),
            args.rtol,
            args.atol,
        )
        print(f"  {name:<28} {'ok  ' if ok else 'FAIL'} {message}")
        failures += not ok

    if failures:
        print(f"\n{failures} deck(s) differ from reference", file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
