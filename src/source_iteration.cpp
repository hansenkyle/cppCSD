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

  LDCSD_LOG_INFO("Begin source iteration");

  // initial guess (maybe provided)
  Eigen::MatrixXd phi = Eigen::MatrixXd::Zero(4 * input_deck.mesh.n_x, input_deck.energy.G);
  Eigen::VectorXd phi_g = Eigen::VectorXd::Zero(4 * input_deck.mesh.n_x);

  Eigen::MatrixXd psi = Eigen::MatrixXd::Zero(4 * input_deck.mesh.n_x, input_deck.angle.M);
  Eigen::MatrixXd psi_up = Eigen::MatrixXd::Zero(4 * input_deck.mesh.n_x, input_deck.angle.M);

  solution.angular_flux = std::vector<Eigen::MatrixXd>(
      input_deck.energy.G, Eigen::MatrixXd::Zero(4 * input_deck.mesh.n_x, input_deck.angle.M));

  // for each E:
  for (int g = 0; g < input_deck.energy.G; g++) {
    LDCSD_LOG_INFO("Beginning group " + std::to_string(g));
    const auto group_start = std::chrono::steady_clock::now();

    int iteration = 0;
    double abs_diff_l2 = 0.0;
    Eigen::VectorXd abs_diff;
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

      abs_diff = (phi_g - phi.col(g));
      abs_diff_l2 = (phi_g - phi.col(g)).norm();
      const double phi_norm = phi_g.norm();
      converged = abs_diff_l2 <= phi_norm * epsilon;

      convergence_.log_group(g,
                             IterationRecord(abs_diff.norm(), abs_diff.lpNorm<Eigen::Infinity>()));

      const int max_cell = static_cast<int>(max_row) / 4;
      const int max_corner = static_cast<int>(max_row) - 4 * max_cell;
      LDCSD_LOG_INFO("group " + std::to_string(g) + " iteration " + std::to_string(iteration) +
                     ": |dphi| = " + std::format("{:.4e}", abs_diff_l2) +
                     ", |phi| = " + std::format("{:.4e}", phi_norm) +
                     ", max|residual| = " + std::format("{:.4e}", max_residual) + " (cell " +
                     std::to_string(max_cell) + ", corner " + std::to_string(max_corner) +
                     ", ordinate " + std::to_string(max_col) + ")");
    }

    phi.col(g) = phi_g;

    const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - group_start;
    convergence_.time_group(g, elapsed.count());

    if (converged) {
      LDCSD_LOG_INFO("Converged with abs. norm = " + std::format("{:.4e}", abs_diff_l2) + " in " +
                     std::to_string(iteration) + " iterations");
    } else {
      LDCSD_LOG_WARN("group " + std::to_string(g) + " did NOT converge: hit the " +
                     std::to_string(max_iterations) +
                     "-iteration cap with abs. norm = " + std::format("{:.4e}", abs_diff_l2));
    }
    psi_up = psi;
    solution.angular_flux[g] = psi;
  }

  solution.scalar_flux = phi;
}

namespace {
auto int_label_seq = [](int max) {
  auto intview =
      std::views::iota(1, max + 1) | std::views::transform([](int x) { return std::to_string(x); });
  std::vector<std::string> result(intview.begin(), intview.end());
  return result;
};

auto evdoub_to_string = [](const Eigen::VectorXd& dvec, std::string_view fmt = "{:.2e}") {
  std::vector<std::string> result;
  for (Eigen::Index i = 0; i < dvec.size(); i++) {
    result.push_back(std::vformat(fmt, std::make_format_args(dvec(i))));
  }
  return result;
};
} // namespace

void SourceIteration::writeResults(const std::filesystem::path& results_path) const {

  auto vint_to_vstring = [](std::vector<int> ivec) {
    std::vector<std::string> result = {};
    for (auto i : ivec) {
      result.push_back(std::to_string(i));
    }
    return result;
  };

  auto vdoub_to_string = [](std::vector<double> dvec, std::string fmt = "{:.2e}") {
    std::vector<std::string> result = {};
    for (auto i : dvec) {
      result.push_back(std::vformat(fmt, std::make_format_args(i)));
    }
    return result;
  };

  auto average_space = [&](Eigen::MatrixXd data) {
    Eigen::MatrixXd left =
        data(Eigen::seqN(0, 2 * input_deck.mesh.n_x, 2), Eigen::placeholders::all);
    Eigen::MatrixXd right =
        data(Eigen::seqN(1, 2 * input_deck.mesh.n_x, 2), Eigen::placeholders::all);
    return (left + right) / 2;
  };

  auto interleave = [&](std::vector<std::string> original) {
    std::vector<std::string> doubled(2 * original.size());
    for (int i = 0; i < original.size(); i++) {
      doubled[2 * i] = original[i];
      doubled[2 * i + 1] = "";
    }
    return doubled;
  };

  auto iseq = int_label_seq(input_deck.mesh.n_x);
  auto gseq = int_label_seq(input_deck.energy.G);
  auto mseq = int_label_seq(input_deck.angle.M);

  auto x_center = evdoub_to_string(input_deck.mesh.x_center);
  auto mu = evdoub_to_string(input_deck.angle.mu);

  std::vector<std::string> xb, eb;

  for (int i = 0; i < input_deck.mesh.n_x; i++) {
    auto s = evdoub_to_string(input_deck.mesh.x_boundary(Eigen::seqN(i, 2)));
    xb.insert(xb.end(), s.begin(), s.end());
  }

  for (int i = 0; i < input_deck.energy.G; i++) {
    auto s = evdoub_to_string(input_deck.energy.E_boundary(Eigen::seqN(i, 2)));
    eb.insert(eb.end(), s.begin(), s.end());
  }

  UnitGroup sol_block("solution");

  // cell-average scalar flux
  MatrixTable scalar("cell-average scalar flux", "averaged over each space-energy cell");
  scalar.set_data(solution.cell_average_scalar().transpose(), x_center, gseq);
  scalar.add_column_label(iseq);

  std::vector<std::vector<std::string>> labels = {{"x_i"}, {"g \\ i"}};
  scalar.set_corner_grid(labels);

  // cell-average angular flux
  labels = {{"", "x_i"}, {"mu", "m \\i"}};
  UnitGroup angular("cell-average angular flux", "averaged over each space-energy cell");
  for (int g = 0; g < input_deck.energy.G; g++) {
    std::cout << g << "\n";
    int gplusone = g + 1;
    MatrixTable group("g = " + std::to_string(gplusone));
    group.set_data(solution.cell_average_angular()[g].transpose(), x_center, mu);
    group.add_column_label(iseq);
    group.add_row_label(mseq);
    group.set_corner_grid(labels);
    angular.add(group, "{:.4e}");
  }

  MatrixTable multigroup("multigroup scalar flux", "averaged over each energy group, not space");
  multigroup.set_data(solution.multigroup(), xb, gseq);
  multigroup.add_column_label(interleave(iseq));
  multigroup.set_corner_grid({{"x_boundary"}, {"g \\ i"}});

  MatrixTable spectrum("energy spectrum", "averaged over each spatial cell, not energy");
  spectrum.set_data(solution.spectrum(), eb, x_center);
  spectrum.add_row_label(iseq);
  spectrum.add_column_label(interleave(gseq));
  spectrum.set_corner_grid({{"", "E_bound"}, {"x_i", "i \\g"}});

  sol_block.add(scalar, "{:.4e}");
  sol_block.add(angular);
  sol_block.add(multigroup, "{:.4e}");
  sol_block.add(spectrum, "{:.4e}");

  appendToFile(results_path, sol_block.render_txt());
}
void SourceIteration::writeConvergence(const std::filesystem::path& results_path) const {
  auto gseq = int_label_seq(input_deck.energy.G);
  VerticalTable summary("iteration summary");
  summary.add_column("g", gseq);
  summary.add_column("# iterations", convergence_.iterations);

  UnitGroup convergence("per-group convergence history",
                        "delta = (phi_n - phi_n-1). absolute change.");
  for (int g = 0; g < input_deck.energy.G; g++) {
    VerticalTable group("g = " + gseq[g]);
    group.add_column("iteration", int_label_seq(convergence_.records[g].size()));

    std::vector<double> l2, li;
    for (int i = 0; i < convergence_.records[g].size(); i++) {
      l2.push_back(convergence_.records[g][i].norm2);
      li.push_back(convergence_.records[g][i].norminf);
    }

    group.add_column("|delta|_2", l2);
    group.add_column("|delta|_infty", li);
    convergence.add(group, "{:.4e}");
  }

  UnitGroup result("convergence");
  result.add(summary);
  result.add(convergence);
  appendToFile(results_path, result.render_txt());
}