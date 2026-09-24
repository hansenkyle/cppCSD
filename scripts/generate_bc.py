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


# Claudepierre et al. (2021), Space Sci. Rev. 217, 80, Fig. 17 (background-corrected
# panel, L=5.46, Probe-A/B, 2013-10-18 19:02:30 UT): (E [keV], flux [1/(cm^2 s ster keV)]),
# read off the plot -- accurate to about a factor of 1.5 (test/method/al_slab/description.md).
_CLAUDEPIERRE_2021_FIG17_E_KEV = np.array([
    28.0, 30.5, 32.8, 35.2, 37.8, 40.6, 43.5, 46.7, 50.2, 53.9, 57.9, 62.1, 66.7, 71.6,
    76.9, 82.5, 88.6, 95.1, 102.2, 109.7, 117.8, 126.4, 135.7, 145.7, 156.5, 168.0, 180.4,
    193.7, 207.9, 223.2, 239.7, 257.3, 276.3, 296.6, 318.5, 341.9, 367.1, 394.2, 423.2,
    454.4, 487.8, 523.8, 562.3, 603.8, 648.2, 696.0, 747.2, 802.3, 861.4, 924.8, 992.9,
    1066.1, 1144.6, 1228.9, 1319.4, 1416.6, 1520.9, 1632.9, 1753.2, 1882.3, 2021.0,
    2169.8, 2329.6, 2501.2, 2685.4, 2883.2, 3095.6, 3323.6, 3568.4, 3831.2, 4113.4, 4416.3,
])
_CLAUDEPIERRE_2021_FIG17_FLUX = np.array([
    46500.6, 56311.1, 63164.8, 60792.2, 44753.9, 39897.9, 35568.8, 31709.4, 27206.9,
    27206.9, 26185.0, 19276.8, 16539.7, 14191.2, 14745.0, 13145.1, 10855.0, 10855.0,
    11278.6, 10855.0, 9313.7, 9677.2, 9677.2, 9313.7, 8627.1, 8303.1, 7991.2, 7991.2,
    6856.5, 7124.1, 7124.1, 6599.0, 4675.6, 4168.2, 3716.0, 4011.7, 3861.0, 3312.8,
    2632.9, 2632.9, 2632.9, 2438.8, 1938.3, 2013.9, 1938.3, 1373.3, 1426.9, 1373.3,
    1091.5, 1091.5, 1011.0, 936.5, 591.5, 470.1, 373.6, 333.1, 346.1, 320.6, 236.0,
    218.6, 90.6, 64.2, 51.0, 49.1, 39.0, 11.5, 2.9, 1.8, 1.0, 0.6, 0.6, 0.5,
])


def claudepierre_2021_fig17_spectrum(E_mev: np.ndarray) -> np.ndarray:
    """Incoming boundary flux digitized from Claudepierre et al. (2021), Fig. 17
    (test/method/al_slab/description.md). log-log interpolated between the
    digitized points; clamped flat outside their [28 keV, 4.4 MeV] range.
    Converted to 1/(cm^2 s ster MeV) (x1000) to match the energy mesh's MeV units.
    """
    E_kev = E_mev * 1000.0
    log_flux = np.interp(np.log(E_kev), np.log(_CLAUDEPIERRE_2021_FIG17_E_KEV),
                          np.log(_CLAUDEPIERRE_2021_FIG17_FLUX))
    return np.exp(log_flux) * 1000.0


SPECTRA = {
    "van_allen": al_slab_spectrum,
    "claudepierre": claudepierre_2021_fig17_spectrum,
}


def group_corners(psi, E_boundary: np.ndarray, n_quad: int) -> tuple[np.ndarray, np.ndarray]:
    """LD reconstruction of psi(E) on each group [E_boundary[g+1], E_boundary[g]].
    Returns (down, up): the low- and high-energy corner values, length G each.

    For a group much wider than psi's decay scale, the linear (P0 +/- P1) corner
    extrapolation can overshoot below zero even though psi itself never goes
    negative; corners are clipped to 0 to keep the reconstruction physical.
    """
    down = np.empty(len(E_boundary) - 1)
    up = np.empty(len(E_boundary) - 1)
    for g in range(len(E_boundary) - 1):
        E_hi, E_lo = E_boundary[g], E_boundary[g + 1]
        E = np.linspace(E_lo, E_hi, 201)
        psi_vals = psi(E)
        average = np.trapezoid(psi_vals, E)/(E_hi - E_lo)
        grad = np.gradient(psi_vals[1:], (E[1:]+E[:-1])/2)
        slope = np.trapezoid(grad, (E[1:]+E[:-1])/2)/(E_hi - E_lo)
        down[g] = max(0.0, average - 0.5*slope*(E_hi - E_lo))
        up[g] = max(0.0, average + 0.5*slope*(E_hi - E_lo))
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
                         help="number of log-spaced groups from --e-max to --e-min (al_slab's energy mesh)")
    parser.add_argument("--e-max", type=float, default=1.0, help="highest energy mesh boundary, MeV")
    parser.add_argument("--e-min", type=float, default=1e-8, help="lowest energy mesh boundary, MeV")
    parser.add_argument("--spectrum", choices=SPECTRA, default="van_allen", help="incoming boundary spectrum")
    parser.add_argument("--M", type=int, default=16, help="ordinates to broadcast each group's isotropic value across")
    parser.add_argument("--n-quad", type=int, default=8, help="Gauss-Legendre order for the average/slope moments")
    parser.add_argument("--plot", action="store_true", help="plot the LD reconstruction against psi(E)")
    args = parser.parse_args()

    spectrum = SPECTRA[args.spectrum]
    E_boundary = np.geomspace(args.e_max, args.e_min, args.groups + 1)
    down, up = group_corners(spectrum, E_boundary, args.n_quad)

    bc = {
        "boundary_conditions": {
            "up": [[float(v)] * args.M for v in up],
            "down": [[float(v)] * args.M for v in down],
        }
    }
    print(yaml.dump(bc, sort_keys=False, default_flow_style=None))

    if args.plot:
        plot_reconstruction(spectrum, E_boundary, down, up)


if __name__ == "__main__":
    main()

