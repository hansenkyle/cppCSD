// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "source_iteration.h"
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

void SourceIteration::solve(double epsilon, int max_iterations) {
  // Solve the transport equation in all groups via source iteration.
  //
  // No iteration over energy groups, assume downscatter only

  LDCSD_LOG_INFO("Begin source iteration");

  // initial guess (maybe provided)
  Eigen::MatrixXd phi = Eigen::MatrixXd::Zero(4 * input_deck.mesh.n_x, input_deck.energy.G);
  Eigen::VectorXd phi_g = Eigen::VectorXd::Zero(4 * input_deck.mesh.n_x);

  Eigen::MatrixXd psi = Eigen::MatrixXd::Zero(4 * input_deck.mesh.n_x, input_deck.angle.M);
  Eigen::MatrixXd psi_up = Eigen::MatrixXd::Zero(4 * input_deck.mesh.n_x, input_deck.angle.M);

  result.angular_flux = std::vector<Eigen::MatrixXd>(
      input_deck.energy.G, Eigen::MatrixXd::Zero(4 * input_deck.mesh.n_x, input_deck.angle.M));

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
    result.angular_flux[g] = psi;
  }

  result.scalar_flux = phi;
}

void appendToFile(const std::filesystem::path& file_path, const std::string& text) {
  std::ofstream out(file_path, std::ios::app);
  if (!out.is_open()) {
    throw std::runtime_error("SourceIteration: failed to open '" + file_path.string() +
                             "' for writing");
  }
  out << text;
}

void SourceIteration::writeResults(const std::filesystem::path& results_path,
                                   const Eigen::MatrixXd& scalar_flux,
                                   const std::vector<Eigen::MatrixXd>& angular_flux) const {
  appendToFile(results_path, SolverFormatter::formatResults(scalar_flux, angular_flux, input_deck));
}
