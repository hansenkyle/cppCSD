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
G, M, NX = 3, 4, 3


def sample_bc():
    # sample_input.yaml: down[g][m] = 10m + g, up = down + 0.5.
    down = np.array([[10 * m + g for m in range(M)] for g in range(G)], dtype=float)
    return np.stack([down + 0.5, down], axis=1)  # [G, 2 (up, down), M]


def sample_source():
    # sample_input.yaml: base 100m + 10g + c, then up_left/up_right/down_left/down_right add
    # 0.1/0.2/0.3/0.4 -- i.e. 0.1 * (1 + 2 * (up, down) + (L, R)).
    q = np.empty((G, M, NX, 2, 2))
    for g, m, c, e, s in np.ndindex(q.shape):
        # round: the YAML's decimal literal, not 0.1 * 3 = 0.30000000000000004
        q[g, m, c, e, s] = round(100 * m + 10 * g + c + 0.1 * (1 + 2 * e + s), 1)
    return q


def build_sample():
    """sample_input.yaml, built through the setters from its own numbers."""
    deck = ldcsd.InputDeck()
    deck.set_mesh([0, 1, 2, 3])
    deck.set_energy([5.0, 1.0, 0.5, 0.0])
    deck.set_angle([-0.9, -0.3, 0.3, 0.9], [0.5, 0.5, 0.5, 0.5])
    water_scatter = np.zeros((G, G))
    for frm, to, value in [(0, 0, 1.1), (0, 1, 0.05), (1, 1, 0.9), (1, 2, 0.03), (2, 2, 0.7)]:
        water_scatter[frm, to] = value
    deck.set_materials(
        {
            "water": ldcsd.Material(
                total=[1.2, 1.0, 0.8],
                S=[2.0, 1.8, 1.5],
                S_b=[2.2, 1.9, 1.6, 1.3],
                scatter=water_scatter,
            ),
            "lead": ldcsd.Material(
                total=[3.2, 3.0, 2.8],
                S=[5.0, 4.8, 4.5],
                S_b=[5.2, 4.9, 4.6, 4.3],
                scatter=[[2.1, 0.1, 0.0], [0.0, 1.9, 0.08], [0.0, 0.0, 1.7]],
            ),
        },
        regions=["water", "water", "lead"],
    )
    deck.set_bc(sample_bc())
    deck.set_source(sample_source())
    return deck


def assert_same_deck(a, b):
    assert (a.mesh.n_x, a.energy.G, a.angle.M) == (b.mesh.n_x, b.energy.G, b.angle.M)
    np.testing.assert_array_equal(a.mesh.x_boundary, b.mesh.x_boundary)
    np.testing.assert_array_equal(a.energy.E_boundary, b.energy.E_boundary)
    np.testing.assert_array_equal(a.angle.mu, b.angle.mu)
    np.testing.assert_array_equal(a.angle.w, b.angle.w)
    assert a.regions == b.regions
    assert list(a.materials) == list(b.materials)
    for name, mat in a.materials.items():
        other = b.materials[name]
        assert mat.name == other.name == name
        for field in ("total", "S", "S_b", "scatter"):
            np.testing.assert_array_equal(getattr(mat, field), getattr(other, field))
    np.testing.assert_array_equal(a.bc, b.bc)
    np.testing.assert_array_equal(a.source, b.source)


def test_setters_build_the_same_deck_as_read():
    deck = build_sample()
    deck.validate()
    assert_same_deck(deck, ldcsd.InputDeck.read(SAMPLE))


def test_getters_use_reader_layout():
    deck = ldcsd.InputDeck.read(SAMPLE)
    assert deck.bc.shape == (G, 2, M)
    assert deck.source.shape == (G, M, NX, 2, 2)
    np.testing.assert_array_equal(deck.bc, sample_bc())
    np.testing.assert_array_equal(deck.source, sample_source())


def test_to_yaml_round_trips(tmp_path):
    deck = build_sample()
    deck.to_yaml(tmp_path / "deck.yaml")
    assert_same_deck(ldcsd.InputDeck.read(tmp_path / "deck.yaml"), deck)


def test_returned_arrays_are_read_only_copies():
    deck = build_sample()
    with pytest.raises(ValueError):
        deck.bc[0, 0, 0] = -1.0
    with pytest.raises(ValueError):
        deck.mesh.x_boundary[0] = -1.0
    with pytest.raises(AttributeError):
        deck.mesh.n_x = 7


def one_group_material(**overrides):
    fields = dict(total=[1.0], S=[1.0], S_b=[1.0, 1.0], scatter=[[0.5]])
    return ldcsd.Material(**(fields | overrides))


@pytest.mark.parametrize(
    "call, message",
    [
        (lambda d: d.set_mesh([0, 2, 1]), "strictly ascending"),
        (lambda d: d.set_energy([0.0, 1.0]), "strictly descending"),
        (lambda d: d.set_angle([0.5, -0.5], [1, 1]), "strictly ascending"),
        (lambda d: d.set_angle([-0.5, 0.5], [1, 3]), "sum to 2"),
        (lambda d: d.set_materials({"a": one_group_material(S_b=[1.0])}, ["a"]), "sizes"),
        (lambda d: d.set_materials({"a": one_group_material(total=[-1.0])}, ["a"]), "non-negative"),
        (lambda d: d.set_materials({"a": one_group_material()}, ["b"]), "undefined material 'b'"),
        (lambda d: d.set_source(-np.ones((1, 2, 1, 2, 2))), "non-negative"),
    ],
)
def test_setter_rejects_bad_input_and_leaves_deck_unchanged(call, message):
    deck = build_sample()
    with pytest.raises(RuntimeError, match=message):
        call(deck)
    assert_same_deck(deck, build_sample())


@pytest.mark.parametrize(
    "call, message",
    [
        (lambda d: d.set_mesh([[0, 1]]), r"x_boundary \[nx\+1\], got shape \(1, 2\)"),
        (lambda d: d.set_bc(np.zeros((G, M))), r"bc \[G, 2, M\]"),
        (lambda d: d.set_bc(np.zeros((G, 3, M))), r"\[G, 2 \(up, down\), M\]"),
        (lambda d: d.set_source(np.zeros((G, M, NX, 2))), r"source \[G, M, nx, 2, 2\]"),
    ],
)
def test_wrong_array_shape_raises_value_error(call, message):
    with pytest.raises(ValueError, match=message):
        call(ldcsd.InputDeck())


def test_read_raises_with_the_error():
    with pytest.raises(RuntimeError, match="bad file"):
        ldcsd.InputDeck.read("does_not_exist.yaml")


def test_validate_catches_cross_struct_mismatch():
    deck = build_sample()
    deck.set_mesh([0, 1, 2, 3, 4])  # regions still cover 3 cells
    with pytest.raises(RuntimeError, match="expected mesh.n_x = 4"):
        deck.validate()
