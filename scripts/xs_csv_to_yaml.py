#!/usr/bin/env python3
# Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
#
# Funded by CARRE (https://carre-psaapiv.org/)
#
# Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
# or without modification are permitted provided that the terms of the license are met.

"""Loads the CSV cross-section format written by scripts/generate_xs.jl.

YAML conversion isn't implemented yet -- this only parses the CSV into
Python objects (see load_csv()).
"""

import csv
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class MaterialXS:
    name: str
    sigma_t: np.ndarray  # size Ng, cm^-1
    S: np.ndarray  # size Ng, MeV/cm (stopping_power_average)
    S_bound: np.ndarray  # size Ng + 1, MeV/cm (stopping_power_boundary)
    scattering: np.ndarray  # Ng x Ng, cm^-1, [from_group, to_group]


def load_csv(path: Path) -> tuple[np.ndarray, dict[str, MaterialXS]]:
    """Parses a generate_xs.jl CSV. Returns (energy_mesh_MeV, {material_name: MaterialXS})."""
    with open(path, newline="") as f:
        rows = [row for row in csv.reader(f) if row and not row[0].startswith("#")]

    row_iter = iter(rows)

    def next_row() -> list[str]:
        return next(row_iter)

    assert next_row() == ["energy_mesh_MeV"]
    energy_mesh = np.array([float(x) for x in next_row()])
    Ng = len(energy_mesh) - 1

    materials: dict[str, MaterialXS] = {}
    for row in row_iter:
        assert row[0] == "material", f"expected a 'material' row, got {row}"
        name = row[1]

        assert next_row() == ["composition"]
        for row in row_iter:
            if row[0] == "stopping_power_boundary_MeV_cm":
                break

        S_bound = np.array([float(x) for x in next_row()])

        assert next_row() == ["group", "sigma_t_cm-1", "stopping_power_average_MeV_cm"]
        sigma_t = np.empty(Ng)
        S = np.empty(Ng)
        for g in range(Ng):
            _, sigma_t_g, S_g = next_row()
            sigma_t[g] = float(sigma_t_g)
            S[g] = float(S_g)

        next_row()  # "scattering_matrix_cm-1 (...)" title row
        next_row()  # "from\to,1,2,...,Ng" header row
        scattering = np.zeros((Ng, Ng))
        for f in range(Ng):
            values = next_row()[1:]
            for t, value in enumerate(values):
                if value != "":
                    scattering[f, t] = float(value)

        materials[name] = MaterialXS(name, sigma_t, S, S_bound, scattering)

    return energy_mesh, materials


if __name__ == "__main__":
    path = Path(sys.argv[1] if len(sys.argv) > 1 else "xs_data/al_27gcc.csv")
    energy_mesh, materials = load_csv(path)
    print(f"energy_mesh_MeV: {energy_mesh}")
    for name, mat in materials.items():
        print(f"\nmaterial '{name}': {len(mat.sigma_t)} groups")
        print(f"  sigma_t: {mat.sigma_t}")
        print(f"  S: {mat.S}")
        print(f"  S_bound: {mat.S_bound}")
        print(f"  scattering:\n{mat.scattering}")
