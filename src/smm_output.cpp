// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "logger.h"
#include "output_block.h"
#include "smm.h"

#include <format>
#include <highfive/eigen.hpp>
#include <highfive/highfive.hpp>

namespace {
// Every closure in SMClosures, with the name it's written under.
const std::array<std::pair<const char*, Eigen::VectorXd SMClosures::*>, 7> kClosures = {{
    {"F", &SMClosures::F},
    {"F+", &SMClosures::F_pos},
    {"F-", &SMClosures::F_neg},
    {"K+", &SMClosures::K_pos},
    {"K-", &SMClosures::K_neg},
    {"T+", &SMClosures::T_pos},
    {"T-", &SMClosures::T_neg},
}};
} // namespace

void SecondMoment::writeResults(const std::filesystem::path& results_path) const {
  using Eigen::seqN;
  using Eigen::placeholders::all;
  const int I = input_deck.mesh.n_x;
  const int G = input_deck.energy.G;
  const int M = input_deck.angle.M;

  std::vector<std::string> iseq, x_center, x_bound_LD;
  for (int i = 0; i < I; i++) {
    iseq.push_back(std::to_string(i + 1));
    x_center.push_back(std::format("{:.6e}", input_deck.mesh.x_center(i)));
    x_bound_LD.push_back(std::format("{:.6e}", input_deck.mesh.x_boundary(i)));
    x_bound_LD.push_back(std::format("{:.6e}", input_deck.mesh.x_boundary(i + 1)));
  }

  std::vector<std::string> mseq, mu;
  for (int m = 0; m < M; m++) {
    mseq.push_back(std::to_string(m + 1));
    mu.push_back(std::format("{:.6e}", input_deck.angle.mu(m)));
  }

  std::vector<std::string> gseq;
  std::vector<std::string> E_bound_LD;
  for (int g = 0; g < G; g++) {
    gseq.push_back(std::to_string(g + 1));
    E_bound_LD.push_back(std::format("{:.6e}", input_deck.energy.E_boundary(g)));
    E_bound_LD.push_back(std::format("{:.6e}", input_deck.energy.E_boundary(g + 1)));
  }

  // CELL-AVERAGED
  UnitGroup cell_ave_scalar("cell-averaged scalar flux");

  // low-order scalar flux
  MatrixTable low_order_scalar("low-order scalar flux", "cell-averaged");
  low_order_scalar.set_data(solution.cell_average_scalar_flux().transpose(), x_center, gseq);
  low_order_scalar.add_column_label(iseq);
  low_order_scalar.set_corner_grid({{"x_i"}, {"g \\ i"}});

  // high-order scalar flux
  MatrixTable high_order_scalar("high-order scalar flux", "cell-averaged");
  high_order_scalar.set_data(solution.cell_average_scalar_flux(true).transpose(), x_center, gseq);
  high_order_scalar.add_column_label(iseq);
  high_order_scalar.set_corner_grid({{"x_i"}, {"g \\ i"}});

  // compute relative difference
  Eigen::MatrixXd scalar_diff = Eigen::MatrixXd::Ones(I, G);
  scalar_diff -= solution.cell_average_scalar_flux(true).cwiseProduct(
      solution.cell_average_scalar_flux().cwiseInverse());

  // format relative difference
  MatrixTable scalar_difference_unit(
      "1 - (HO/LO)",
      std::format("cell-averaged first, then compared. maximum: {:.4e}", scalar_diff.maxCoeff()));
  scalar_difference_unit.set_data(scalar_diff.transpose(), x_center, gseq);
  scalar_difference_unit.add_column_label(iseq);
  scalar_difference_unit.set_corner_grid({{"x_i"}, {"g \\ i"}});

  cell_ave_scalar.add(low_order_scalar, "{:.6e}");
  cell_ave_scalar.add(high_order_scalar, "{:.6e}");
  cell_ave_scalar.add(scalar_difference_unit, "{:.4e}");

  UnitGroup cell_ave_current("cell-averaged current");
  // low-order current
  MatrixTable low_order_current("low-order current", "cell-averaged");
  low_order_current.set_data(solution.cell_average_current().transpose(), x_center, gseq);
  low_order_current.add_column_label(iseq);
  low_order_current.set_corner_grid({{"x_i"}, {"g \\ i"}});
  // high-order current

  MatrixTable high_order_current("high-order current", "cell-averaged");
  high_order_current.set_data(solution.cell_average_current(true).transpose(), x_center, gseq);
  high_order_current.add_column_label(iseq);
  high_order_current.set_corner_grid({{"x_i"}, {"g \\ i"}});

  // difference
  Eigen::MatrixXd current_diff = Eigen::MatrixXd::Ones(I, G);
  current_diff -= solution.cell_average_current(true).cwiseProduct(
      solution.cell_average_current().cwiseInverse());

  MatrixTable current_difference_unit(
      "1 - (HO/LO)",
      std::format("cell-averaged first, then compared. maximum: {:.4e}", current_diff.maxCoeff()));
  current_difference_unit.set_data(current_diff.transpose(), x_center, gseq);
  current_difference_unit.add_column_label(iseq);
  current_difference_unit.set_corner_grid({{"x_i"}, {"g \\ i"}});

  cell_ave_current.add(low_order_current, "{:.6e}");
  cell_ave_current.add(high_order_current, "{:.6e}");
  cell_ave_current.add(current_difference_unit, "{:.4e}");

  // CORNER VALUES

  UnitGroup cv_scalar("scalar flux corner values");

  Eigen::MatrixXd scalar_cv_diff = Eigen::MatrixXd::Ones(4 * I, G);
  scalar_cv_diff -=
      solution.high_order_scalar_flux.cwiseProduct(solution.scalar_flux.cwiseInverse());

  for (int g = 0; g < G; g++) {
    UnitGroup group("scalar flux corner values, g=" + std::to_string(g + 1));
    HorizontalTable lo_group("low-order");
    HorizontalTable ho_group("high-order");
    HorizontalTable diff_group("1-HO/LO");

    lo_group.add_row("x_i", x_center);
    lo_group.add_row("i", iseq);
    lo_group.add_row("L, up", solution.scalar_flux(seqN(0, I, 4), g));
    lo_group.add_row("R, up", solution.scalar_flux(seqN(1, I, 4), g));
    lo_group.add_row("L, down", solution.scalar_flux(seqN(2, I, 4), g));
    lo_group.add_row("R, down", solution.scalar_flux(seqN(3, I, 4), g));

    ho_group.add_row("x_i", x_center);
    ho_group.add_row("i", iseq);
    ho_group.add_row("L, up", solution.high_order_scalar_flux(seqN(0, I, 4), g));
    ho_group.add_row("R, up", solution.high_order_scalar_flux(seqN(1, I, 4), g));
    ho_group.add_row("L, down", solution.high_order_scalar_flux(seqN(2, I, 4), g));
    ho_group.add_row("R, down", solution.high_order_scalar_flux(seqN(3, I, 4), g));

    diff_group.add_row("x_i", x_center);
    diff_group.add_row("i", iseq);
    diff_group.add_row("L, up", scalar_cv_diff(seqN(0, I, 4), g));
    diff_group.add_row("R, up", scalar_cv_diff(seqN(1, I, 4), g));
    diff_group.add_row("L, down", scalar_cv_diff(seqN(2, I, 4), g));
    diff_group.add_row("R, down", scalar_cv_diff(seqN(3, I, 4), g));

    group.add(lo_group, "{:.6e}");
    group.add(ho_group, "{:.6e}");
    group.add(diff_group, "{:.4e}");

    cv_scalar.add(group);
  }

  // low-order current
  // high-order current
  // difference

  // ANGULAR FLUX

  // cell-averaged
  // corner values

  // CLOSURES

  UnitGroup close("closures");
  for (int g = 0; g < input_deck.energy.G; g++) {
    UnitGroup group("g = " + std::to_string(g + 1));
    for (const auto& [name, field] : kClosures) {
      const Eigen::VectorXd& c = closures[g].*field;
      HorizontalTable table(name);
      table.add_row("x_i", x_center);
      table.add_row("i", iseq);
      table.add_row("up,left", c(seqN(0, I, 4)));
      table.add_row("up,right", c(seqN(1, I, 4)));
      table.add_row("down,left", c(seqN(2, I, 4)));
      table.add_row("down,right", c(seqN(3, I, 4)));
      group.add(table, "{:.6e}");
    }
    close.add(group);
  }

  //   UnitGroup block = solutionBlock(solution);
  //   block.add(close);
  //   block.add(cellAverageTable("high-order cell-average scalar flux",
  //                              MethodResult::cell_average(solution.high_order_scalar)),
  //             "{:.6e}");
  appendToFile(results_path, cell_ave_scalar.render_txt());
  appendToFile(results_path, cell_ave_current.render_txt());

  appendToFile(results_path, cv_scalar.render_txt());

  appendToFile(results_path, close.render_txt());
}

void SecondMoment::writeResiduals(const std::filesystem::path& file_path, std::string timestamp) {
  using Eigen::seqN;
  using Eigen::placeholders::all;
  writeMetadata(file_path, timestamp);

  static constexpr std::array<const char*, 4> kEdgeLabels = {"up,L", "up,R", "down,L", "down,R"};
  static constexpr std::array<const char*, 8> kSMEqLabels = {
      "(balance, up L)",    "(balance, up R)",    "(balance, down L)",    "(balance, down R)",
      "(1st moment, up L)", "(1st moment, up R)", "(1st moment, down L)", "(1st moment, down R)"};

  // Largest |residual| over every group, its (cell, angle, edge), from a per-group Eigen
  // container whose rows are laid out cell-major in blocks of block_size (4 for transport, one
  // edge per row; 8 for the SM equations, one equation per row).
  auto peakResidual = [](const auto& per_group, int block_size) {
    struct Peak {
      double value = -1.0;
      int g = -1, i = -1, sub = -1, m = -1;
    } peak;
    for (int g = 0; g < static_cast<int>(per_group.size()); g++) {
      Eigen::Index row, col;
      const double gmax = per_group[g].cwiseAbs().maxCoeff(&row, &col);
      if (gmax > peak.value) {
        peak = {gmax, g, static_cast<int>(row) / block_size, static_cast<int>(row) % block_size,
                static_cast<int>(col)};
      }
    }
    return peak;
  };

  const auto transport_peak = peakResidual(residuals.high_order, 4);
  const std::string transport_summary =
      std::format("max |residual| = {:.4e} at g= {}, i= {}, angle {}, {}", transport_peak.value,
                  transport_peak.g + 1, transport_peak.i + 1, transport_peak.m + 1,
                  kEdgeLabels[transport_peak.sub]);
  LDCSD_LOG_INFO("transport residuals: " + transport_summary);

  // Re-evaluate the true peak (found above, across all groups) term-by-term.
  const Eigen::MatrixXd zero_psi_gm1 =
      Eigen::MatrixXd::Zero(4 * input_deck.mesh.n_x, input_deck.angle.M);
  transport_operator.calculateResiduals(
      transport_peak.g, solution.angular_flux[transport_peak.g],
      transport_peak.g == 0 ? zero_psi_gm1 : solution.angular_flux[transport_peak.g - 1],
      solution.scalar_flux, /*debug_max=*/true);

  UnitGroup transport("transport residuals", transport_summary);
  int I = input_deck.mesh.n_x;
  std::vector<std::string> x_i, mu_m, x_center;

  for (int i = 0; i < input_deck.mesh.n_x; i++) {
    x_i.push_back(std::to_string(i + 1));
    x_center.push_back(std::format("{:.5e}", input_deck.mesh.x_center(i)));
  }
  for (int m = 0; m < input_deck.angle.M; m++) {
    mu_m.push_back(std::to_string(m + 1));
  }
  for (int g = 0; g < input_deck.energy.G; g++) {
    UnitGroup group("g = " + std::to_string(g + 1));
    for (int m = 0; m < input_deck.angle.M; m++) {
      HorizontalTable angle("m = " + std::to_string(m + 1));
      angle.add_row("x_i", x_center);
      angle.add_row("i", x_i);
      angle.add_row("up,L", residuals.high_order[g](seqN(0, I, 4), m));
      angle.add_row("up,R", residuals.high_order[g](seqN(1, I, 4), m));
      angle.add_row("down,L", residuals.high_order[g](seqN(2, I, 4), m));
      angle.add_row("down,R", residuals.high_order[g](seqN(3, I, 4), m));
      group.add(angle, "{:.6e}");
    }
    transport.add(group);
  }

  appendToFile(file_path, transport.render_txt());

  const auto sm_peak = peakResidual(residuals.low_order, 8);
  const std::string sm_summary =
      std::format("max |residual| = {:.4e} at g= {}, i= {}, equation {}", sm_peak.value,
                  sm_peak.g + 1, sm_peak.i + 1, kSMEqLabels[sm_peak.sub]);
  LDCSD_LOG_INFO("second moment equation residuals: " + sm_summary);

  // Re-evaluate the true peak (found above, across all groups) term-by-term.
  calculateResiduals(sm_peak.g, solution.scalar_flux, solution.current,
                     solution.angular_flux[sm_peak.g], /*debug_max=*/true);

  UnitGroup low_order("second moment equation residuals", sm_summary);
  for (int g = 0; g < input_deck.energy.G; g++) {
    HorizontalTable group("g = " + std::to_string(g + 1));
    group.add_row("x_i", x_center);
    group.add_row("i", x_i);
    group.add_row("(balance, up L)", residuals.low_order[g](seqN(0, I, 8)));
    group.add_row("(balance, up R)", residuals.low_order[g](seqN(1, I, 8)));
    group.add_row("(balance, down L)", residuals.low_order[g](seqN(2, I, 8)));
    group.add_row("(balance, down R)", residuals.low_order[g](seqN(3, I, 8)));
    group.add_row("(1st moment, up L)", residuals.low_order[g](seqN(4, I, 8)));
    group.add_row("(1st moment, up R)", residuals.low_order[g](seqN(5, I, 8)));
    group.add_row("(1st moment, down L)", residuals.low_order[g](seqN(6, I, 8)));
    group.add_row("(1st moment, down R)", residuals.low_order[g](seqN(7, I, 8)));
    low_order.add(group, "{:.6e}");
  }

  appendToFile(file_path, low_order.render_txt());
}

void SecondMoment::writeH5(const std::filesystem::path& file_path,
                           const std::string& timestamp) const {
  HighFive::File file(file_path.string(), HighFive::File::Truncate);
  writeH5Common(file, timestamp);

  const std::string corner_order(kCornerOrder);
  file.createDataSet("/solution/scalar_flux", solution.scalar_flux)
      .createAttribute("corner_order", corner_order);
  file.createDataSet("/solution/angular_flux", stackGroups(solution.angular_flux))
      .createAttribute("corner_order", corner_order);
  file.createDataSet("/solution/current", solution.current)
      .createAttribute("corner_order", corner_order);
  file.createDataSet("/solution/high_order_scalar", solution.high_order_scalar_flux)
      .createAttribute("corner_order", corner_order);

  for (const auto& [closure_name, field] : kClosures) {
    std::vector<Eigen::VectorXd> per_group;
    for (const SMClosures& c : closures) {
      per_group.push_back(c.*field);
    }
    file.createDataSet(std::string("/solution/closures/") + closure_name, stackGroups(per_group))
        .createAttribute("corner_order", corner_order);
  }

  file.createDataSet("/residuals/transport", stackGroups(residuals.high_order))
      .createAttribute("corner_order", corner_order);
  file.createDataSet("/residuals/second_moment", stackGroups(residuals.low_order))
      .createAttribute("row_order",
                       std::string("row 8i+k is cell i, k = (balance up L, up R, down L, down R, "
                                   "1st moment up L, up R, down L, down R)"));
}