# Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
#
# Funded by CARRE (https://carre-psaapiv.org/)
#
# Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
# or without modification are permitted provided that the terms of the license are met.

from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("ldcsd._core")
import ldcsd  # noqa: E402

SAMPLE = Path(__file__).resolve().parents[2] / "test" / "unit" / "data" / "sample_input.yaml"


@pytest.mark.parametrize("method", ["SourceIteration", "SecondMoment"])
def test_in_memory_results_match_the_h5_file(method, tmp_path):
    deck = ldcsd.InputDeck.read(SAMPLE)
    solver = getattr(ldcsd, method)(deck)
    solver.solve(1e-10, max_iterations=200)
    solver.write_h5(tmp_path / "run.h5", timestamp="2026-01-01 00:00:00")

    run = ldcsd.read(tmp_path / "run.h5")
    assert run.method == solver.name
    assert run.execution_datetime == "2026-01-01 00:00:00"
    assert solver.scalar_flux.shape == (3, 3, 2, 2)
    assert solver.angular_flux.shape == (3, 4, 3, 2, 2)
    np.testing.assert_array_equal(solver.scalar_flux, run.scalar_flux)
    np.testing.assert_array_equal(solver.angular_flux, run.angular_flux)
    np.testing.assert_array_equal(deck.bc, run.input.bc)
    np.testing.assert_array_equal(deck.source, run.input.source)
    if method == "SecondMoment":
        np.testing.assert_array_equal(solver.current, run.current)
        np.testing.assert_array_equal(solver.reconstructed_scalar, run.reconstructed_scalar)


def test_solver_validates_the_deck():
    deck = ldcsd.InputDeck()
    deck.set_mesh([0, 1])  # nothing else set
    with pytest.raises(RuntimeError):
        ldcsd.SecondMoment(deck)


def test_configure_logging_writes_the_log(tmp_path):
    ldcsd.configure_logging(tmp_path / "run.log")
    ldcsd.InputDeck.read(SAMPLE)
    assert "read input deck" in (tmp_path / "run.log").read_text()
