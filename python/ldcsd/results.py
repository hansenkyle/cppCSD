# Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
#
# Funded by CARRE (https://carre-psaapiv.org/)
#
# Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
# or without modification are permitted provided that the terms of the license are met.

"""Reads an ldcsd run's .h5 results file (layout in docs/output-layout.md).

Corner-valued fields are reshaped from the file's flat 4nx rows to
[G, nx, 2, 2], indexed [group, cell, energy edge (up, down), space edge (L, R)].
Fields with an angle axis come out as [G, M, nx, 2, 2], so psi[g, m] has the
same shape as phi[g].
"""

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np


@dataclass
class MaterialData:
    total: np.ndarray  # [G]
    S: np.ndarray  # [G], group-average stopping power
    S_b: np.ndarray  # [G+1], stopping power at group boundaries
    scatter: np.ndarray  # [G, G], (from, to)


@dataclass
class Input:
    x_boundary: np.ndarray  # [nx+1]
    dx: np.ndarray  # [nx]
    x_center: np.ndarray  # [nx]
    material: list[str]  # [nx], material name in each cell
    E_boundary: np.ndarray  # [G+1], descending
    dE: np.ndarray  # [G]
    mu: np.ndarray  # [M]
    w: np.ndarray  # [M]
    bc: np.ndarray  # [G, 2 (up, down), M]
    source: np.ndarray  # [G, M, nx, 2, 2]
    q0: np.ndarray  # [G, nx, 2, 2]
    q1: np.ndarray  # [G, nx, 2, 2]
    materials: dict[str, MaterialData]


@dataclass
class Convergence:
    iterations: np.ndarray  # [G]
    group_times: np.ndarray  # [G], seconds
    # One array per group, one entry per iteration.
    delta_l2: list[np.ndarray]
    delta_linf: list[np.ndarray]
    phi_l2: list[np.ndarray]
    phi_linf: list[np.ndarray]


@dataclass
class Run:
    method: str
    execution_datetime: str
    ldcsd_version: str
    input: Input
    convergence: Convergence
    scalar_flux: np.ndarray  # [G, nx, 2, 2]
    angular_flux: np.ndarray  # [G, M, nx, 2, 2]
    transport_residual: np.ndarray  # [G, M, nx, 2, 2]
    # Second moment method only; None for other methods.
    current: np.ndarray | None = None  # [G, nx, 2, 2]
    reconstructed_scalar: np.ndarray | None = None  # [G, nx, 2, 2]
    closures: dict[str, np.ndarray] | None = None  # F, F+, F-, K+, K-, T+, T-: [G, nx, 2, 2]
    # [G, nx, 2 (balance, 1st moment), 2 (up, down), 2 (L, R)]
    second_moment_residual: np.ndarray | None = None


def cell_average(corners: np.ndarray) -> np.ndarray:
    """Average over each space-energy cell: [..., nx, 2, 2] -> [..., nx]."""
    return corners.mean(axis=(-2, -1))


def spectrum(corners: np.ndarray) -> np.ndarray:
    """Average over each spatial cell, keeping each group's upper and lower
    energy edge: [..., nx, 2, 2] -> [..., nx, 2 (up, down)]."""
    return corners.mean(axis=-1)


def multigroup(corners: np.ndarray) -> np.ndarray:
    """Average over each energy group, keeping each cell's left and right
    edge: [..., nx, 2, 2] -> [..., nx, 2 (L, R)]."""
    return corners.mean(axis=-2)


def _corners(flat: np.ndarray) -> np.ndarray:
    # Last axis is the file's corner rows 4i + 2*(up, down) + (L, R).
    return flat.reshape(*flat.shape[:-1], -1, 2, 2)


def _by_group(f: h5py.File, path: str) -> np.ndarray:
    # [4nx, G] -> [G, nx, 2, 2]
    return _corners(f[path][()].T)


def _by_group_angle(f: h5py.File, path: str) -> np.ndarray:
    # [G, 4nx, M] -> [G, M, nx, 2, 2]
    return _corners(np.swapaxes(f[path][()], 1, 2))


def read(path: str | Path) -> Run:
    """Reads a <deck>.h5 file written by ldcsd."""
    with h5py.File(path, "r") as f:
        mats = f["input/materials"]
        inp = Input(
            x_boundary=f["input/mesh/x_boundary"][()],
            dx=f["input/mesh/dx"][()],
            x_center=f["input/mesh/x_center"][()],
            material=list(f["input/mesh/material"].asstr()[()]),
            E_boundary=f["input/energy/E_boundary"][()],
            dE=f["input/energy/dE"][()],
            mu=f["input/angle/mu"][()],
            w=f["input/angle/w"][()],
            bc=f["input/bc"][()].reshape(-1, 2, f["input/bc"].shape[1]),
            source=_by_group_angle(f, "input/source/values"),
            q0=_corners(f["input/source/q0"][()]),
            q1=_corners(f["input/source/q1"][()]),
            materials={
                name: MaterialData(**{k: mats[name][k][()] for k in ("total", "S", "S_b", "scatter")})
                for name in mats
            },
        )

        conv = f["convergence"]
        groups = sorted(k for k, v in conv.items() if isinstance(v, h5py.Group))
        convergence = Convergence(
            iterations=conv["iterations"][()],
            group_times=conv["group_times"][()],
            **{
                k: [conv[g][k][()] for g in groups]
                for k in ("delta_l2", "delta_linf", "phi_l2", "phi_linf")
            },
        )

        run = Run(
            method=f.attrs["method"],
            execution_datetime=f.attrs["execution_datetime"],
            ldcsd_version=f.attrs["ldcsd_version"],
            input=inp,
            convergence=convergence,
            scalar_flux=_by_group(f, "solution/scalar_flux"),
            angular_flux=_by_group_angle(f, "solution/angular_flux"),
            transport_residual=_by_group_angle(f, "residuals/transport"),
        )

        if "solution/current" in f:
            run.current = _by_group(f, "solution/current")
            run.reconstructed_scalar = _by_group(f, "solution/reconstructed_scalar")
            run.closures = {k: _corners(v[()]) for k, v in f["solution/closures"].items()}
            sm = f["residuals/second_moment"][()]
            run.second_moment_residual = sm.reshape(sm.shape[0], -1, 2, 2, 2)

    return run
