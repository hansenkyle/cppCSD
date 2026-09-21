// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "solver.h"

#include <fstream>
#include <sstream>

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

namespace {

// A small but fully-populated deck: 2 cells, 2 groups, 2 ordinates (one of
// each sign), nonzero scattering, stopping power, and external source.
// Everything is O(1) so an absolute residual tolerance is meaningful.
InputDeck makeResidualDeck() {
  InputDeck deck;
  deck.mesh.n_x = 2;
  deck.mesh.x_boundary = Eigen::Vector3d(0.0, 0.6, 1.4);
  deck.energy.G = 2;
  deck.energy.E_boundary = Eigen::Vector3d(2.0, 1.2, 0.5);
  deck.angle.M = 2;
  deck.angle.mu = Eigen::Vector2d(-0.5, 0.5);
  deck.angle.w = Eigen::Vector2d(1.0, 1.0);

  deck.xs.total = Eigen::MatrixXd(2, 2);
  deck.xs.total << 1.3, 0.9, 0.7, 1.1;

  // Nonzero scattering is what makes this test bite: the scattering source is
  // the only term that vanishes when phi == 0, so a bug in its assembly is
  // invisible on a first iteration and only shows up once phi is populated.
  deck.xs.scatter.clear();
  for (int i = 0; i < 2; ++i) {
    Eigen::SparseMatrix<double> s(2, 2);
    s.insert(0, 0) = 0.4 + 0.1 * i; // within-group, g0
    s.insert(0, 1) = 0.3;           // downscatter g0 -> g1
    s.insert(1, 1) = 0.5;           // within-group, g1
    s.makeCompressed();
    deck.xs.scatter.push_back(s);
  }

  deck.xs.S = Eigen::MatrixXd(2, 2);
  deck.xs.S << 0.6, 0.45, 0.35, 0.55;
  deck.xs.S_bound = Eigen::MatrixXd(3, 2);
  deck.xs.S_bound << 0.7, 0.5, 0.5, 0.4, 0.3, 0.6;

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

// The real invariant tying Kernel::solveDirect to Solver::cellResidual: whatever
// solveDirect produces must satisfy, to machine precision, the same discretized
// equations that cellResidual independently re-derives in extended precision.
// Any term one of them assembles differently from the other shows up here.
//
// Crucially this runs with a NONZERO scalar flux. With phi == 0 the scattering
// source drops out of both sides and the check passes no matter how that term
// is assembled -- which is exactly how a missing dx factor on solveDirect's
// scattering source stayed hidden through a first iteration.
TEST_CASE("transportSweep output satisfies the discretized equations to machine precision") {
  Solver solver(makeResidualDeck());

  const int n_rows = 4 * 2;
  const Eigen::MatrixXd zero_psi = Eigen::MatrixXd::Zero(n_rows, 2);

  // An arbitrary, definitely-nonzero scalar flux in both groups.
  Eigen::MatrixXd phi(n_rows, 2);
  phi << 1.7, 0.8, 1.4, 0.95, 1.2, 1.1, 1.05, 1.35, 0.9, 1.6, 1.25, 0.75, 1.5, 1.15, 0.85, 1.45;

  for (int g = 0; g < 2; ++g) {
    // g == 1 also exercises the CSD coupling, by feeding it a nonzero
    // previous-group psi rather than the zero matrix g == 0 gets.
    const Eigen::MatrixXd psi_gm1 =
        (g == 0) ? zero_psi : Eigen::MatrixXd(solver.transportSweep(0, zero_psi, phi));

    const Eigen::MatrixXd psi = solver.transportSweep(g, psi_gm1, phi);
    const Eigen::MatrixXd residuals = solver.calculateResiduals(g, psi, psi_gm1, phi);

    CAPTURE(g);
    CHECK(residuals.cwiseAbs().maxCoeff() < 1e-12);
  }
}

namespace {

// A minimally-valid, 1-cell/1-group/1-ordinate deck.
InputDeck makeDeck() {
  InputDeck deck;
  deck.mesh.n_x = 1;
  deck.mesh.x_boundary = Eigen::Vector2d(0.0, 1.0);
  deck.energy.G = 1;
  deck.energy.E_boundary = Eigen::Vector2d(1.0, 0.0);
  deck.angle.M = 1;
  deck.angle.mu = Eigen::VectorXd::Constant(1, 0.5);
  deck.angle.w = Eigen::VectorXd::Constant(1, 2.0);
  deck.xs.total = Eigen::MatrixXd::Constant(1, 1, 1.0);
  deck.xs.scatter = {Eigen::SparseMatrix<double>(1, 1)};
  deck.xs.S = Eigen::MatrixXd::Constant(1, 1, 1.0);
  deck.xs.S_bound = Eigen::MatrixXd::Constant(2, 1, 1.0);
  deck.bc.values = Eigen::MatrixXd::Constant(2, 1, 0.0);
  deck.source.values = {Eigen::MatrixXd::Constant(4, 1, 0.0)};
  deck.validate();
  return deck;
}

// Reads a whole file into a string.
std::string readFile(const std::filesystem::path& path) {
  std::ifstream in(path);
  std::ostringstream oss;
  oss << in.rdbuf();
  return oss.str();
}

// A fresh path with nothing there yet, so each test starts from a clean
// (nonexistent) file regardless of what earlier runs left behind.
std::filesystem::path freshPath(const std::string& name) {
  const std::filesystem::path path = std::filesystem::temp_directory_path() / name;
  std::filesystem::remove(path);
  return path;
}

} // namespace
