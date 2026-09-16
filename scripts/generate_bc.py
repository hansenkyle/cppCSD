#!/usr/bin/env python3
# Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
#
# Funded by CARRE (https://carre-psaapiv.org/)
#
# Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
# or without modification are permitted provided that the terms of the license are met.

"""Generates InputDeck's boundary_conditions.up/down from a continuous function
of energy psi(E).

Each energy group is a linear-discontinuous (LD) element: psi is reconstructed
on the group from its average and average slope (the Legendre P0 and P1
moments of psi over the group, computed by Gauss-Legendre quadrature), then
evaluated at the group's low- and high-energy corners to give down/up.
Plotting requires matplotlib (pip install matplotlib), only if --plot is used.
"""

import argparse

import numpy as np
import yaml

# All floats render as e-notation with 8 decimal digits, matching the rest of the
# codebase's cross-section/quadrature data.
yaml.add_representer(float, lambda dumper, v: dumper.represent_scalar("tag:yaml.org,2002:float", f"{v:.8e}"))

def al_slab_spectrum(E_mev: np.ndarray) -> np.ndarray:
    """Van Allen upper-belt spectrum used as the al_slab test's incoming boundary
    flux (test/method/al_slab/description.md, Fennel et al. 2015, figure 2):
    psi = 1.5e5 * exp(-0.01 * E[keV]) 1/(cm^2 s ster keV), E < 1000 keV; 0 above.
    Converted to 1/(cm^2 s ster MeV) (x1000, since d(keV) = 1e-3 d(MeV)) to match
    the energy mesh's MeV units.
    """
    E_kev = E_mev * 1000.0
    return np.where(E_kev < 1000.0, 1.5e5 * np.exp(-0.01 * E_kev) * 1000.0, 0.0)


def group_corners(psi, E_boundary: np.ndarray, n_quad: int) -> tuple[np.ndarray, np.ndarray]:
    """LD reconstruction of psi(E) on each group [E_boundary[g+1], E_boundary[g]].
    Returns (down, up): the low- and high-energy corner values, length G each.
    """
    xi, w = np.polynomial.legendre.leggauss(n_quad)
    down = np.empty(len(E_boundary) - 1)
    up = np.empty(len(E_boundary) - 1)
    for g in range(len(E_boundary) - 1):
        E_hi, E_lo = E_boundary[g], E_boundary[g + 1]
        E = 0.5 * (E_hi - E_lo) * xi + 0.5 * (E_hi + E_lo)
        psi_vals = psi(E)
        average = 0.5 * np.sum(w * psi_vals)
        slope = 1.5 * np.sum(w * psi_vals * xi)
        down[g] = average - slope
        up[g] = average + slope
    return down, up


def plot_reconstruction(psi, E_boundary: np.ndarray, down: np.ndarray, up: np.ndarray) -> None:
    import matplotlib.pyplot as plt

    E_fine = np.geomspace(E_boundary[-1], E_boundary[0], 2000)
    _, ax = plt.subplots()
    ax.plot(E_fine, psi(E_fine), label="psi(E)")
    for g in range(len(E_boundary) - 1):
        ax.plot([E_boundary[g + 1], E_boundary[g]], [down[g], up[g]], "o-", color="tab:orange",
                label="LD reconstruction" if g == 0 else None)
    ax.set_xscale("log")
    ax.set_xlabel("E (MeV)")
    ax.set_ylabel("psi")
    ax.legend()
    plt.show()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", type=int, default=36,
                         help="number of log-spaced groups from 1 MeV to 1e-8 MeV (al_slab's energy mesh)")
    parser.add_argument("--M", type=int, default=16, help="ordinates to broadcast each group's isotropic value across")
    parser.add_argument("--n-quad", type=int, default=8, help="Gauss-Legendre order for the average/slope moments")
    parser.add_argument("--plot", action="store_true", help="plot the LD reconstruction against psi(E)")
    args = parser.parse_args()

    E_boundary = np.geomspace(1.0, 1e-8, args.groups + 1)
    down, up = group_corners(al_slab_spectrum, E_boundary, args.n_quad)

    bc = {
        "boundary_conditions": {
            "up": [[float(v)] * args.M for v in up],
            "down": [[float(v)] * args.M for v in down],
        }
    }
    print(yaml.dump(bc, sort_keys=False, default_flow_style=None))

    if args.plot:
        plot_reconstruction(al_slab_spectrum, E_boundary, down, up)


if __name__ == "__main__":
    main()
