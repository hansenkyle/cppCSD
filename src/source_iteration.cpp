// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "source_iteration.h"
#include "logger.h"
#include "output_block.h"

#include <array>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <string>
#include <vector>

void SourceIteration::solve(double epsilon, int max_iterations) {
  // Solve the transport equation in all groups via source iteration.
  //
  // No iteration over energy groups, assume downscatter only

  int I = input_deck.mesh.n_x;
  int G = input_deck.energy.G;
  int M = input_deck.angle.M;

  LDCSD_LOG_INFO("Begin source iteration");

  Eigen::VectorXd phi_g = Eigen::VectorXd::Zero(4 * I);
  Eigen::MatrixXd psi_up = Eigen::MatrixXd::Zero(4 * I, M);

  // write phi directly into this->solution
  solution.scalar_flux = Eigen::MatrixXd::Zero(4 * I, G);
  auto& phi = solution.scalar_flux;

  solution.angular_flux =
      std::vector<Eigen::MatrixXd>(input_deck.energy.G, Eigen::MatrixXd::Zero(4 * I, M));
  residuals.high_order = solution.angular_flux;

  for (int g = 0; g < G; g++) {
    LDCSD_LOG_INFO("Beginning group " + std::to_string(g));
    const auto group_start = std::chrono::steady_clock::now();

    // write result directly into this->solution
    auto& psi = solution.angular_flux[g];

    int iteration = 0;
    double delta_phi_l2norm = 0.0;
    Eigen::VectorXd absolute_delta_phi;

    while (iteration < max_iterations) {
      iteration++;
      phi.col(g) = phi_g;
      // solve transport using known phi
      psi = transport_operator.sweep(g, psi_up, phi);
      // compute new phi
      phi_g = transport_operator.integrateAngle(psi);

      absolute_delta_phi = (phi_g - phi.col(g));
      delta_phi_l2norm = l2norm(phi_g - phi.col(g));
      const double phi_norm = l2norm(phi_g);

      convergence_.log_group(g, IterationRecord({l2norm(phi_g), linfnorm(phi_g)},
                                                {delta_phi_l2norm, linfnorm(phi_g - phi.col(g))}));

      if (delta_phi_l2norm <= phi_norm * epsilon) {
        LDCSD_LOG_INFO("Converged with abs. norm = " + std::format("{:.4e}", delta_phi_l2norm) +
                       " in " + std::to_string(iteration) + " iterations");
        break; // exit while loop
      }
    }

    residuals.high_order[g] =
        transport_operator.calculateResiduals(g, psi, psi_up, phi, log_residual_terms);
    phi.col(g) = phi_g;
    psi_up = psi;

    const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - group_start;
    convergence_.time_group(g, elapsed.count());

    if (iteration == max_iterations) {
      LDCSD_LOG_WARN("group " + std::to_string(g) + " did NOT converge: hit the " +
                     std::to_string(max_iterations) +
                     "-iteration cap with abs. norm = " + std::format("{:.4e}", delta_phi_l2norm));
    }
  }
}

void SourceIteration::writeResults(const std::filesystem::path& results_path) const {
  appendToFile(results_path, solutionBlock(solution).render_txt());
}

void SourceIteration::writeResiduals(const std::filesystem::path& file_path,
                                     std::string timestamp) const {
  using Eigen::seqN;
  using Eigen::placeholders::all;
  writeMetadata(file_path, timestamp);

  UnitGroup transport("transport residuals");
  int I = input_deck.mesh.n_x;
  std::vector<std::string> x_i;
  std::vector<std::string> mu_m;

  for (int i = 0; i < input_deck.mesh.n_x; i++) {
    x_i.push_back(std::to_string(i + 1));
  }
  for (int m = 0; m < input_deck.angle.M; m++) {
    mu_m.push_back(std::to_string(m + 1));
  }
  for (int g = 0; g < input_deck.energy.G; g++) {
    UnitGroup group("g = " + std::to_string(g + 1));
    for (int m = 0; m < input_deck.angle.M; m++) {
      HorizontalTable angle("m = " + std::to_string(m + 1));
      angle.add_row("i", x_i);
      angle.add_row("up,L", residuals.high_order[g](seqN(0, I, 4), m));
      angle.add_row("up,R", residuals.high_order[g](seqN(1, I, 4), m));
      angle.add_row("down,L", residuals.high_order[g](seqN(2, I, 4), m));
      angle.add_row("down,R", residuals.high_order[g](seqN(3, I, 4), m));
      group.add(angle, "{:.4e}");
    }
    transport.add(group);
  }

  appendToFile(file_path, transport.render_txt());
}
