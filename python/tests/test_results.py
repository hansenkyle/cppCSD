# Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
#
# Funded by CARRE (https://carre-psaapiv.org/)
#
# Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
# or without modification are permitted provided that the terms of the license are met.

import h5py
import numpy as np
import pytest

import ldcsd

NX, G, M = 3, 2, 4


def code(g, i, k, m=0):
    # Unique value per (group, cell, corner row k, angle), so a wrong reshape shows up as a
    # wrong number.
    return 1000 * g + 100 * m + 10 * i + k


def write_run(path, smm):
    corner = np.array([[code(g, r // 4, r % 4) for r in range(4 * NX)] for g in range(G)])
    angular = np.array(
        [[[code(g, r // 4, r % 4, m) for m in range(M)] for r in range(4 * NX)] for g in range(G)]
    )
    with h5py.File(path, "w") as f:
        f.attrs.update(method="test", execution_datetime="now", ldcsd_version="0.1.0")
        f["input/mesh/x_boundary"] = np.linspace(0, 1, NX + 1)
        f["input/mesh/dx"] = np.full(NX, 1 / NX)
        f["input/mesh/x_center"] = np.arange(NX)
        f["input/mesh/material"] = ["a", "a", "b"]
        f["input/energy/E_boundary"] = np.array([3.0, 2.0, 1.0])
        f["input/energy/dE"] = np.ones(G)
        f["input/angle/mu"] = np.linspace(-1, 1, M)
        f["input/angle/w"] = np.full(M, 2 / M)
        f["input/bc"] = np.arange(2 * G * M).reshape(2 * G, M)
        f["input/source/values"] = angular
        f["input/source/q0"] = corner
        f["input/source/q1"] = corner
        for name in ("a", "b"):
            f[f"input/materials/{name}/total"] = np.ones(G)
            f[f"input/materials/{name}/S"] = np.ones(G)
            f[f"input/materials/{name}/S_b"] = np.ones(G + 1)
            f[f"input/materials/{name}/scatter"] = np.eye(G)
        f["convergence/iterations"] = np.array([3, 2])
        f["convergence/group_times"] = np.array([0.1, 0.2])
        for g, n in enumerate((3, 2)):
            for k in ("delta_l2", "delta_linf", "phi_l2", "phi_linf"):
                f[f"convergence/g{g + 1:03d}/{k}"] = np.arange(n, dtype=float)
        f["solution/scalar_flux"] = corner.T
        f["solution/angular_flux"] = angular
        f["residuals/transport"] = angular
        if smm:
            f["solution/current"] = corner.T
            f["solution/reconstructed_scalar"] = corner.T
            for c in ("F", "F+", "F-", "K+", "K-", "T+", "T-"):
                f[f"solution/closures/{c}"] = corner
            f["residuals/second_moment"] = np.array(
                [[code(g, r // 8, r % 8) for r in range(8 * NX)] for g in range(G)]
            )


@pytest.fixture
def run(tmp_path):
    write_run(tmp_path / "run.h5", smm=True)
    return ldcsd.read(tmp_path / "run.h5")


def test_corner_fields_index_as_group_cell_energy_space(run):
    for g in range(G):
        for i in range(NX):
            for e in range(2):
                for s in range(2):
                    k = 2 * e + s
                    assert run.scalar_flux[g, i, e, s] == code(g, i, k)
                    assert run.closures["K+"][g, i, e, s] == code(g, i, k)
                    for m in range(M):
                        assert run.angular_flux[g, m, i, e, s] == code(g, i, k, m)
    assert run.second_moment_residual[1, 2, 1, 0, 1] == code(1, 2, 4 + 0 + 1)
    assert run.input.bc[1, 0, 2] == 2 * M + 2  # group 1, up, angle 2
    assert run.input.material == ["a", "a", "b"]
    assert [len(d) for d in run.convergence.delta_l2] == [3, 2]


def test_helpers_match_the_cpp_text_output(run):
    phi = run.scalar_flux
    flat = phi.reshape(G, -1)  # rows 4i + (up L, up R, down L, down R), as in the file
    np.testing.assert_allclose(ldcsd.cell_average(phi), flat.reshape(G, NX, 4).mean(-1))
    # MethodResult::spectrum: (up L + up R) / 2 and (down L + down R) / 2
    np.testing.assert_allclose(
        ldcsd.spectrum(phi)[..., 0], (flat[:, 0::4] + flat[:, 1::4]) / 2
    )
    np.testing.assert_allclose(
        ldcsd.spectrum(phi)[..., 1], (flat[:, 2::4] + flat[:, 3::4]) / 2
    )
    # MethodResult::multigroup: (up L + down L) / 2 and (up R + down R) / 2
    np.testing.assert_allclose(
        ldcsd.multigroup(phi)[..., 0], (flat[:, 0::4] + flat[:, 2::4]) / 2
    )
    np.testing.assert_allclose(
        ldcsd.multigroup(phi)[..., 1], (flat[:, 1::4] + flat[:, 3::4]) / 2
    )
    assert ldcsd.cell_average(run.angular_flux).shape == (G, M, NX)


def test_source_iteration_run_has_no_smm_fields(tmp_path):
    write_run(tmp_path / "si.h5", smm=False)
    run = ldcsd.read(tmp_path / "si.h5")
    assert run.current is None and run.closures is None and run.second_moment_residual is None
