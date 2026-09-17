// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "solver.h"

#include <doctest.h>

namespace {
// solveDirect's output vector is energy-major, space-minor: [up_L, up_R, down_L, down_R].
// Swap L<->R within each energy block, e.g. for comparing against a mu-mirrored solve.
Eigen::Vector4d swapLR(const Eigen::Vector4d& x) { return x({1, 0, 3, 2}); }
} // namespace

// Solving with -cosine must reproduce the +cosine solution with L/R swapped within each energy
// block, provided the L/R-structured inputs (psi_in_E, q_up, q_down) are also pre-swapped --
// solveDirect always assembles as if flow travels L->R using |cosine|, so the caller must
// present already-mirrored inputs for a negative cosine to represent the same physical problem.
TEST_CASE("solveDirect mirrors under cosine sign flip with matching L/R-swapped inputs") {
  Solver::Kernel kernel;

  const double dx = 0.7, dE = 0.3, xs = 1.2, S = 0.4, S_up = 0.5, S_down = 0.6;
  const double psi_in_x_up = 0.8, psi_in_x_down = 0.5;
  const Eigen::Vector2d psi_in_E(1.3, 0.9);
  const Eigen::Vector2d q_up(0.2, 0.4);
  const Eigen::Vector2d q_down(0.6, 0.1);
  const Eigen::VectorXd sigma_sdEprime = Eigen::VectorXd::Zero(0);
  const Eigen::MatrixXd phi_up = Eigen::MatrixXd::Zero(2, 0);
  const Eigen::MatrixXd phi_down = Eigen::MatrixXd::Zero(2, 0);

  const double mu = 0.6;

  Eigen::Vector4d x_pos =
      kernel.solveDirect(mu, dx, dE, xs, S, S_up, S_down, psi_in_E, psi_in_x_down, psi_in_x_up,
                         q_up, q_down, sigma_sdEprime, phi_up, phi_down);

  Eigen::Vector4d x_neg =
      kernel.solveDirect(-mu, dx, dE, xs, S, S_up, S_down, Eigen::Vector2d(psi_in_E.reverse()),
                         psi_in_x_down, psi_in_x_up, Eigen::Vector2d(q_up.reverse()),
                         Eigen::Vector2d(q_down.reverse()), sigma_sdEprime, phi_up, phi_down);

  Eigen::Vector4d expected = swapLR(x_pos);
  for (int i = 0; i < 4; ++i) {
    CHECK(x_neg[i] == doctest::Approx(expected[i]));
  }
}

// With zero stopping power (S = S_up = S_down = 0, no CSD coupling) and matched incoming
// conditions between the two energy nodes, the up- and down-energy solutions must come out
// identical -- nothing in the assembled system distinguishes "up" from "down" once CSD is off.
TEST_CASE("solveDirect gives equal up/down blocks when stopping power is zero and incoming "
          "x fluxes match") {
  Solver::Kernel kernel;

  const double dx = 0.7, dE = 0.3, xs = 1.2;
  const double psi_in_x = 0.9;
  const Eigen::Vector2d psi_in_E(1.3, 0.9); // irrelevant: S_up = 0 zeroes its contribution
  const Eigen::Vector2d q(0.3, 0.7);        // same source on both energy nodes
  const Eigen::VectorXd sigma_sdEprime = Eigen::VectorXd::Zero(0);
  const Eigen::MatrixXd phi = Eigen::MatrixXd::Zero(2, 0);

  Eigen::Vector4d x = kernel.solveDirect(0.6, dx, dE, xs, /*S=*/0.0, /*S_up=*/0.0,
                                         /*S_down=*/0.0, psi_in_E, psi_in_x, psi_in_x, q, q,
                                         sigma_sdEprime, phi, phi);

  CHECK(x[0] == doctest::Approx(x[2]));
  CHECK(x[1] == doctest::Approx(x[3]));
}

// With zero scattering source, zero external source, and zero CSD (S = S_up = S_down = 0), the
// cell has no absorption either (xs = 0), so it's pure void streaming: the outgoing flux must
// equal the incoming flux exactly, uniform across the cell (L == R), independent of mu, dx, dE.
TEST_CASE("solveDirect reduces to pure streaming with no absorption, scattering, CSD, or source") {
  Solver::Kernel kernel;

  const double dx = 0.4, dE = 0.2;
  const double psi_in_x_up = 0.8, psi_in_x_down = 0.9;
  const Eigen::Vector2d zero2 = Eigen::Vector2d::Zero();
  const Eigen::VectorXd sigma_sdEprime = Eigen::VectorXd::Zero(0);
  const Eigen::MatrixXd phi = Eigen::MatrixXd::Zero(2, 0);

  for (double mu : {0.3, 0.6, 0.95}) {
    Eigen::Vector4d x = kernel.solveDirect(mu, dx, dE, /*xs=*/0.0, /*S=*/0.0, /*S_up=*/0.0,
                                           /*S_down=*/0.0, zero2, psi_in_x_down, psi_in_x_up, zero2,
                                           zero2, sigma_sdEprime, phi, phi);

    CHECK(x[0] == doctest::Approx(psi_in_x_up));
    CHECK(x[1] == doctest::Approx(psi_in_x_up));
    CHECK(x[2] == doctest::Approx(psi_in_x_down));
    CHECK(x[3] == doctest::Approx(psi_in_x_down));
  }
}
