# Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
#
# Funded by CARRE (https://carre-psaapiv.org/)
#
# Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
# or without modification are permitted provided that the terms of the license are met.

"""Builds a two-material slab deck in Python, solves it with the second moment method, and
writes the run (slab.h5) and the deck (slab.yaml) to the current directory. See
docs/python-bindings.md.

The numbers are illustrative, not physical data.
"""

import numpy as np

import ldcsd

NX, G = 10, 3

deck = ldcsd.InputDeck()
deck.set_mesh(np.linspace(0.0, 1.0, NX + 1))  # 10 cells over 1 cm
deck.set_energy([3.0, 2.0, 1.0, 0.0])  # 3 groups, descending to 0
mu, w = np.polynomial.legendre.leggauss(8)  # S8; weights already sum to 2
deck.set_angle(mu, w)
M = deck.angle.M


def material(sigma_t, S_b, downscatter):
    # Within-group scattering at half of sigma_t, plus a little into the next group down.
    sigma_t = np.asarray(sigma_t)
    scatter = np.diag(0.5 * sigma_t) + np.diag(np.full(G - 1, downscatter), k=1)
    S_b = np.asarray(S_b)
    return ldcsd.Material(total=sigma_t, S=(S_b[:-1] + S_b[1:]) / 2, S_b=S_b, scatter=scatter)


deck.set_materials(
    {
        "light": material([1.0, 1.2, 1.5], [1.0, 1.3, 1.8, 2.5], downscatter=0.05),
        "heavy": material([3.0, 3.5, 4.0], [2.0, 2.4, 3.0, 4.0], downscatter=0.1),
    },
    regions=["light"] * 5 + ["heavy"] * 5,
)

# bc[g, (up, down), m] is the incoming flux for ordinate m on whichever face it enters through:
# mu > 0 at x = 0, mu < 0 at x = 1. Here, a beam into the top group from the left only.
bc = np.zeros((G, 2, M))
bc[0, :, mu > 0] = 1.0
deck.set_bc(bc)

deck.set_source(np.zeros((G, M, NX, 2, 2)))  # no volumetric source

solver = ldcsd.SecondMoment(deck)  # runs deck.validate()
solver.solve(1e-10, max_iterations=500)

phi = ldcsd.cell_average(solver.scalar_flux)  # [G, nx]
np.set_printoptions(precision=4, linewidth=120)
print("cell-average scalar flux [group, cell]:")
print(phi)

solver.write_h5("slab.h5")  # read back with ldcsd.read()
deck.to_yaml("slab.yaml")  # rerun with the ldcsd executable
