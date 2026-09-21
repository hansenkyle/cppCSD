// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "solver.h"
#include "logger.h"
#include "solver_formatter.h"

#include <array>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

Solver::Kernel::Kernel() {
  M << 2.0, 1.0, 1.0, 2.0;
  M *= (1.0 / 6);
  L << 0.5, 0.5, -0.5, -0.5;
  Lb << -1, 0, 0, 1;

  A = Eigen::Matrix4d::Zero();
  b = Eigen::Vector4d::Zero();
}

Eigen::MatrixXd Solver::sourceIterate(double epsilon, int max_iterations) {
  // Solve the transport equation in all groups via source iteration.
  //
  // Groups are solved in a single downward pass and never revisited: that's
  // correct for pure CSD plus downscatter, and silently wrong if a deck ever
  // carries upscatter. The scalar flux is *not* reset between groups, so
  // group g starts from group g-1's converged answer (a warm start).

  LDCSD_LOG_INFO("Begin source iteration");
  convergence_ = ConvergenceHistory{};

  // initial guess (maybe provided)
  Eigen::MatrixXd phi = Eigen::MatrixXd::Zero(4 * input_deck.mesh.n_x, input_deck.energy.G);
  Eigen::VectorXd phi_g = Eigen::VectorXd::Zero(4 * input_deck.mesh.n_x);

  Eigen::MatrixXd psi = Eigen::MatrixXd::Zero(4 * input_deck.mesh.n_x, input_deck.angle.M);
  Eigen::MatrixXd psi_up = Eigen::MatrixXd::Zero(4 * input_deck.mesh.n_x, input_deck.angle.M);

  // for each E:
  for (int g = 0; g < input_deck.energy.G; g++) {
    LDCSD_LOG_INFO("Beginning group " + std::to_string(g));
    const auto group_start = std::chrono::steady_clock::now();

    int iteration = 0;
    double abs_diff = 0.0;
    bool converged = false;

    while (!converged && iteration < max_iterations) {
      iteration++;
      phi.col(g) = phi_g;
      // solve transport using known phi
      psi = transport_operator.sweep(g, psi_up, phi);

      const Eigen::MatrixXd residuals =
          transport_operator.calculateResiduals(g, psi, psi_up, phi, log_residual_terms);
      Eigen::Index max_row, max_col;
      const double max_residual = residuals.cwiseAbs().maxCoeff(&max_row, &max_col);

      // compute new phi
      phi_g = transport_operator.integrateAngle(psi);

      abs_diff = (phi_g - phi.col(g)).norm();
      const double phi_norm = phi_g.norm();
      converged = abs_diff <= phi_norm * epsilon;

      convergence_.record(IterationRecord{g, iteration, phi_norm, abs_diff, 0.0, max_residual});

      const int max_cell = static_cast<int>(max_row) / 4;
      const int max_corner = static_cast<int>(max_row) - 4 * max_cell;
      LDCSD_LOG_INFO("group " + std::to_string(g) + " iteration " + std::to_string(iteration) +
                     ": |dphi| = " + std::format("{:.4e}", abs_diff) +
                     ", |phi| = " + std::format("{:.4e}", phi_norm) +
                     ", max|residual| = " + std::format("{:.4e}", max_residual) + " (cell " +
                     std::to_string(max_cell) + ", corner " + std::to_string(max_corner) +
                     ", ordinate " + std::to_string(max_col) + ")");
    }

    phi.col(g) = phi_g;

    const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - group_start;
    convergence_.finishGroup(g, converged, elapsed.count());

    if (converged) {
      LDCSD_LOG_INFO("Converged with abs. norm = " + std::format("{:.4e}", abs_diff) + " in " +
                     std::to_string(iteration) + " iterations");
    } else {
      LDCSD_LOG_WARN("group " + std::to_string(g) + " did NOT converge: hit the " +
                     std::to_string(max_iterations) +
                     "-iteration cap with abs. norm = " + std::format("{:.4e}", abs_diff));
    }
    psi_up = psi;
  }
  return phi;
}

Eigen::Vector4d
Solver::Kernel::solveDirect(double cosine, double dx, double dE, double xs, double S, double S_up,
                            double S_down, Eigen::Vector2d psi_in_E, double psi_in_x_down,
                            double psi_in_x_up, Eigen::Vector2d q_up, Eigen::Vector2d q_down,
                            const Eigen::VectorXd& sigma_sdEprime,
                            const Eigen::MatrixXd& phi_gprime_up,
                            const Eigen::MatrixXd& phi_gprime_down, bool check_condition) {
  A = Eigen::Matrix4d::Zero();
  b = Eigen::Vector4d::Zero();

  // All inputs are rotated to solve problem for mu>0 to reduce code duplication.
  //
  //
  double mu = std::abs(cosine);

  Eigen::MatrixXd phi_gprime_up_local = phi_gprime_up;
  Eigen::MatrixXd phi_gprime_down_local = phi_gprime_down;
  if (cosine < 0) {
    psi_in_E.reverseInPlace();
    q_up.reverseInPlace();
    q_down.reverseInPlace();
    phi_gprime_up_local = phi_gprime_up_local.colwise().reverse().eval();
    phi_gprime_down_local = phi_gprime_down_local.colwise().reverse().eval();
  }

  /*
  matrix/vector are energy-major, space-minor:
  up:    L
         R

  down:  L
         R
  */

  // Parameters:
  //   cosine             : angle consine (mu)
  //   dx, dE             : cell widhts
  //   xs                 : total xs
  //   S                  : group-average stopping power
  //   S_up, S_down       : S(g-1), S(g)
  //   psi_in_E           : flux at next-higher energy group, L/R    : [2x1]
  //   psi_in_x_up        : upwind flux in same energy group, up     : scalar
  //   psi_in_x_down      : "                              ", down   : scalar
  //   q_up               : external source, up (L/R)                : [2x1]
  //   q_down             : "             ", down (L/R)              : [2x1]
  //   sigma_sdEprime     : sigma_s(g' -> g) * dE_g' for all g'      : [Gx1]
  //   phi_gprime_up/down : scalar flux in all groups                : [2xG]

  // "Up" LHS
  // streaming (L)
  A({0, 1}, {0, 1}) += (mu / 6) * 2 * L;
  A({0, 1}, {2, 3}) += (mu / 6) * L;
  // streaming (Lb)
  A(1, 1) += (mu / 6) * 2;
  A(1, 3) += (mu / 6);
  // absorption + CSD loss
  A({0, 1}, {0, 1}) += dx * (xs / 3 + S / (2 * dE)) * M;
  A({0, 1}, {2, 3}) += dx * (xs / 6 + S / (2 * dE)) * M;

  // "Up" RHS
  // streaming source
  b(0) += (mu / 6) * (2 * psi_in_x_up + psi_in_x_down);
  // CSD source
  b({0, 1}) += (dx / dE) * S_up * M * psi_in_E;
  // Scattering source
  b({0, 1}) += (dx / 8) * M * (phi_gprime_down_local + phi_gprime_up_local) * sigma_sdEprime;
  // External source
  b({0, 1}) += (dx / 6) * M * (2 * q_up + q_down);

  // "Down" LHS
  // streaming (L)
  A({2, 3}, {0, 1}) += (mu / 6) * L;
  A({2, 3}, {2, 3}) += (mu / 6) * 2 * L;
  // Streaming (Lb)
  A(3, 1) += (mu / 6);
  A(3, 3) += (mu / 6) * 2;
  // absorption + CSD loss
  A({2, 3}, {0, 1}) += dx * (xs / 6 - S / (2 * dE)) * M;
  A({2, 3}, {2, 3}) += dx * (xs / 3 + (S_down - S / 2) / dE) * M;

  // "Down" RHS
  // streaming source
  b(2) += (mu / 6) * (psi_in_x_up + 2 * psi_in_x_down);
  // Scattering source
  b({2, 3}) += (dx / 8) * M * (phi_gprime_down_local + phi_gprime_up_local) * sigma_sdEprime;
  // External source
  b({2, 3}) += (dx / 6) * M * (q_up + 2 * q_down);

  if (check_condition) {
    Eigen::JacobiSVD<Eigen::Matrix4d> svd(A);
    const Eigen::Vector4d& singular_values = svd.singularValues();
    double condition_number = singular_values(0) / singular_values(singular_values.size() - 1);
    LDCSD_LOG_INFO("solveDirect: condition number = " + std::format("{:.4e}", condition_number));
  }

  Eigen::Vector4d x = A.partialPivLu().solve(b);

  // Return result in expected order
  //
  if (cosine < 0) {
    return x({1, 0, 3, 2});
  }

  return x;
}

namespace {

void appendToFile(const std::filesystem::path& file_path, const std::string& text) {
  std::ofstream out(file_path, std::ios::app);
  if (!out.is_open()) {
    throw std::runtime_error("Solver: failed to open '" + file_path.string() + "' for writing");
  }
  out << text;
}

} // namespace

void Solver::write_metadata(const std::filesystem::path& results_path,
                            const std::filesystem::path& deck_path) const {
  appendToFile(results_path,
               SolverFormatter::formatRunMetadata(input_deck, deck_path, "source iteration"));
}

void Solver::writeInputDeckEcho(const std::filesystem::path& results_path) const {
  appendToFile(results_path, SolverFormatter::formatInputEcho(input_deck));
}

void Solver::writeResults(const std::filesystem::path& results_path,
                          const Eigen::MatrixXd& scalar_flux,
                          const std::vector<Eigen::MatrixXd>& angular_flux) const {
  appendToFile(results_path, SolverFormatter::formatResults(scalar_flux, angular_flux, input_deck));
}
