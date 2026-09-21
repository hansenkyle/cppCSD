// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "transport_operator.h"

#include <doctest.h>

namespace {
InputDeck makeResidualDeck() {
  InputDeck deck;
  deck.mesh.n_x = 2;
  deck.mesh.x_boundary = Eigen::Vector3d(0.0, 0.6, 1.4);
  deck.energy.G = 2;
  deck.energy.E_boundary = Eigen::Vector3d(2.0, 1.2, 0.5);
  deck.angle.M = 2;
  deck.angle.mu = Eigen::Vector2d(-0.5, 0.5);
  deck.angle.w = Eigen::Vector2d(1.0, 1.0);

  // One material per cell, differing in every quantity, so a bug that reads
  // the wrong cell's data shows up in the residual.
  //
  // Nonzero scattering is what makes this test bite: the scattering source is
  // the only term that vanishes when phi == 0, so a bug in its assembly is
  // invisible on a first iteration and only shows up once phi is populated.
  Material left;
  left.name = "left";
  left.total = Eigen::Vector2d(1.3, 0.7);
  left.S = Eigen::Vector2d(0.6, 0.35);
  left.S_b = Eigen::Vector3d(0.7, 0.5, 0.3);
  left.scatter = Eigen::MatrixXd::Zero(2, 2);
  left.scatter(0, 0) = 0.4; // within-group, g0
  left.scatter(0, 1) = 0.3; // downscatter g0 -> g1
  left.scatter(1, 1) = 0.5; // within-group, g1

  Material right;
  right.name = "right";
  right.total = Eigen::Vector2d(0.9, 1.1);
  right.S = Eigen::Vector2d(0.45, 0.55);
  right.S_b = Eigen::Vector3d(0.5, 0.4, 0.6);
  right.scatter = Eigen::MatrixXd::Zero(2, 2);
  right.scatter(0, 0) = 0.5;
  right.scatter(0, 1) = 0.3;
  right.scatter(1, 1) = 0.5;

  deck.xs.set_materials({left, right}, {0, 1});

  deck.bc.values = Eigen::MatrixXd(4, 2);
  deck.bc.values << 1.1, 0.9, 1.3, 0.8, 0.7, 1.2, 0.6, 1.0;

  deck.source.values.clear();
  for (int g = 0; g < 2; ++g) {
    Eigen::MatrixXd q(8, 2);
    q << 0.5, 0.4, 0.6, 0.3, 0.45, 0.55, 0.35, 0.65, 0.7, 0.2, 0.25, 0.75, 0.15, 0.85, 0.8, 0.1;
    deck.source.values.push_back(q * (1.0 + 0.25 * g));
  }

  deck.validate();
  return deck;
}

} // namespace

TEST_CASE("sweep output satisfies the discretized equations to machine precision") {
  TransportOperator solver(makeResidualDeck());

  const int n_rows = 4 * 2;
  const Eigen::MatrixXd zero_psi = Eigen::MatrixXd::Zero(n_rows, 2);

  // An arbitrary, definitely-nonzero scalar flux in both groups.
  Eigen::MatrixXd phi(n_rows, 2);
  phi << 1.7, 0.8, 1.4, 0.95, 1.2, 1.1, 1.05, 1.35, 0.9, 1.6, 1.25, 0.75, 1.5, 1.15, 0.85, 1.45;

  for (int g = 0; g < 2; ++g) {
    // g == 1 also exercises the CSD coupling, by feeding it a nonzero
    // previous-group psi rather than the zero matrix g == 0 gets.
    const Eigen::MatrixXd psi_gm1 =
        (g == 0) ? zero_psi : Eigen::MatrixXd(solver.sweep(0, zero_psi, phi));

    const Eigen::MatrixXd psi = solver.sweep(g, psi_gm1, phi);
    const Eigen::MatrixXd residuals = solver.calculateResiduals(g, psi, psi_gm1, phi);

    CAPTURE(g);
    CHECK(residuals.cwiseAbs().maxCoeff() < 1e-12);
  }
}