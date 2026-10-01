# Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
#
# Funded by CARRE (https://carre-psaapiv.org/)
#
# Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
# or without modification are permitted provided that the terms of the license are met.

from xs_data import *
import numpy as np
import ldcsd

from pathlib import Path

NX, G = 40, 12

deck = ldcsd.InputDeck()
deck.set_mesh(np.linspace(0.0, 0.25, NX + 1))
deck.set_energy([5.000000E+00,	3.859339E+00,	2.978899E+00,	2.299316E+00,	1.774768E+00,	1.369887E+00,
                1.057371E+00,	8.161508E-01,	6.299605E-01,	4.862462E-01,	3.753178E-01,	2.896957E-01, 2.236068E-01])
mu, w = np.polynomial.legendre.leggauss(128)  # S8; weights already sum to 2
deck.set_angle(mu, w)
M = deck.angle.M


def material(sigma_t, S_b, self_scatter):
    # Within-group scattering at half of sigma_t, plus a little into the next group down.
    sigma_t = np.asarray(sigma_t)
    S_b = np.asarray(S_b)
    return ldcsd.Material(total=sigma_t, S=(S_b[:-1] + S_b[1:]) / 2, S_b=S_b, scatter=np.diag(self_scatter))


deck.set_materials(
    {
        "al_6061": material(sigma_t, stopping_power_bound, sigma_s),
    },
    regions=["al_6061"]*NX,
)

bc = np.zeros((G, 2, M))
bc[0, :, mu > 0] = 1000
deck.set_bc(bc)

deck.set_source(np.zeros((G, M, NX, 2, 2)))  # no volumetric source

solver = ldcsd.SecondMoment(deck)  # runs deck.validate()
solver.solve(1e-6, max_iterations=100)

phi = ldcsd.cell_average(solver.scalar_flux)  # [G, nx]
np.set_printoptions(precision=4, linewidth=120)
print("cell-average scalar flux [group, cell]:")
print(phi)

current_file_path = Path(__file__).resolve()
solver.write_h5(f"{current_file_path}slab.h5")  # read back with ldcsd.read()
deck.to_yaml(f"{current_file_path}slab.yaml")  # rerun with the ldcsd executable

solver.write_results(f"{current_file_path}results.txt")
solver.write_residuals(f"{current_file_path}residuals.txt")
