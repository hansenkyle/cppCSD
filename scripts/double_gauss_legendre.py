#!/usr/bin/env python3
# Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
#
# Funded by CARRE (https://carre-psaapiv.org/)
#
# Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
# or without modification are permitted provided that the terms of the license are met.

"""Generates 1D double Gauss-Legendre quadrature sets on mu in [-1, 1]: n-point
Gauss-Legendre on [-1, 0] and matching n-point Gauss-Legendre on [0, 1], each
weight halved so weights over both halves sum to 2.
"""

import argparse

import numpy as np


def double_gauss_legendre(n: int) -> tuple[np.ndarray, np.ndarray]:
    """n-point Gauss-Legendre per half-range, mirrored onto [-1, 0] and [0, 1].
    Returns (mu, weight), each length 2n, ordered ascending in mu.
    """
    x, w = np.polynomial.legendre.leggauss(n)
    mu = np.concatenate([0.5 * x - 0.5, 0.5 * x + 0.5])
    weight = np.concatenate([0.5 * w, 0.5 * w])
    return mu, weight


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("n", type=int, help="quadrature order per half-range (total points = 2n)")
    args = parser.parse_args()

    mu, weight = double_gauss_legendre(args.n)
    for m, w in zip(mu, weight):
        print(f"{m:.8e} {w:.8e}")

    fmt = lambda values: "[" + ", ".join(f"{v:.8e}" for v in values) + "]"
    print(fmt(mu))
    print(fmt(weight))


if __name__ == "__main__":
    main()
