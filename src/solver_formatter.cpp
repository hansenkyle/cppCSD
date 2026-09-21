// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "solver_formatter.h"

#include <array>
#include <ctime>
#include <format>
#include <iomanip>
#include <sstream>

#include "output_block.h"

namespace {

constexpr int kFluxPrecision = 6;

std::string formatSci(double value) {
  std::ostringstream oss;
  oss << std::scientific << std::setprecision(kFluxPrecision) << value;
  return oss.str();
}

std::string timestamp() {
  const std::time_t now = std::time(nullptr);
  const std::tm* tm = std::localtime(&now);
  std::ostringstream oss;
  oss << std::put_time(tm, "%Y-%m-%d %H:%M:%S");
  return oss.str();
}

// Renders a (4 * n_x) x labels.size() matrix as a grid: one column per
// spatial cell, one row per labels entry, with each cell's four values
// (rows 4c..4c+3 = up_left, up_right, down_left, down_right) shown as a 2x2
// cluster:
//   up_left    up_right
//   down_left  down_right
std::string formatCornerGrid(const std::string& title, const Eigen::MatrixXd& data, int n_x,
                             const std::vector<std::string>& labels) {
  const auto count = static_cast<int>(labels.size());
  std::vector<std::vector<std::array<std::string, 4>>> formatted(
      count, std::vector<std::array<std::string, 4>>(n_x));
  std::size_t width = 0;
  for (int i = 0; i < count; ++i) {
    for (int c = 0; c < n_x; ++c) {
      for (int k = 0; k < 4; ++k) {
        std::string text = formatSci(data(4 * c + k, i));
        width = std::max(width, text.size());
        formatted[i][c][k] = std::move(text);
      }
    }
  }
  const int col_width = static_cast<int>(width) + 2;

  std::ostringstream out;
  out << "--- " << title << " ---\n";
  out << "columns are spatial cells 0.." << (n_x - 1)
      << "; each entry is [up_left up_right / down_left down_right]\n";
  for (int i = 0; i < count; ++i) {
    out << labels[i] << ":\n";
    for (int row = 0; row < 2; ++row) {
      for (int c = 0; c < n_x; ++c) {
        const std::array<std::string, 4>& cell = formatted[i][c];
        out << std::setw(col_width) << cell[2 * row] << std::setw(col_width) << cell[2 * row + 1];
      }
      out << "\n";
    }
  }
  return out.str();
}

// "group 0", "group 1", ... -- row labels for a G-column matrix (scalar flux).
std::vector<std::string> groupLabels(int G) {
  std::vector<std::string> labels(G);
  for (int g = 0; g < G; ++g) {
    labels[g] = "group " + std::to_string(g);
  }
  return labels;
}

// "ordinate 0 (mu=.., w=..)", ... -- row labels for an M-column matrix
// (angular flux/residuals within one energy group).
std::vector<std::string> ordinateLabels(const InputDeck& deck) {
  std::vector<std::string> labels(deck.angle.M);
  for (int m = 0; m < deck.angle.M; ++m) {
    labels[m] = "ordinate " + std::to_string(m) + " (mu=" + formatSci(deck.angle.mu[m]) +
                ", w=" + formatSci(deck.angle.w[m]) + ")";
  }
  return labels;
}

} // namespace

namespace SolverFormatter {

std::string formatRunMetadata(const InputDeck& deck, const std::filesystem::path& deck_path,
                              std::string method_name) {
  KeyValueOutput metadata("Run Info");
  metadata.add("execution date/time", timestamp());
  metadata.add("input deck path", deck_path.string());
  metadata.add("n_groups", deck.energy.G);
  metadata.add("n_cells", deck.mesh.n_x);
  metadata.add("n_angles", deck.angle.M);
  metadata.add("method", method_name);
  return metadata.render_txt() + "\n";
}

std::string formatInputEcho(const InputDeck& deck) {
  std::string result = "[Input Echo]\n\n";
  HorizontalTable spatialdata;
  spatialdata.add_row("cell boundaries", deck.mesh.x_boundary);
  spatialdata.add_row("dx", deck.mesh.dx);

  HorizontalTable energydata;
  energydata.add_row("group boundaries", deck.energy.E_boundary);
  energydata.add_row("dE", deck.energy.dE);

  result += spatialdata.render_txt("{:.2e}", 1) + "\n";
  result += energydata.render_txt("{:.2e}", 1) + "\n";

  return result;
}

std::string formatResults(const Eigen::MatrixXd& scalar_flux,
                          const std::vector<Eigen::MatrixXd>& angular_flux, const InputDeck& deck) {
  std::ostringstream out;
  out << formatCornerGrid("Scalar Flux", scalar_flux, deck.mesh.n_x, groupLabels(deck.energy.G));
  for (std::size_t g = 0; g < angular_flux.size(); ++g) {
    out << formatCornerGrid("Angular Flux - Group " + std::to_string(g), angular_flux[g],
                            deck.mesh.n_x, ordinateLabels(deck));
  }
  return out.str();
}

std::string formatConvergence(const ConvergenceHistory& history) {
  if (history.empty()) {
    return "[Convergence Summary]\n\nno iterations recorded\n";
  }

  std::ostringstream out;

  const std::vector<GroupSummary>& groups = history.groups();
  const auto n_groups = static_cast<int>(groups.size());

  std::vector<double> group_index(n_groups), iterations(n_groups), phi_norm(n_groups),
      abs_diff(n_groups), rel_diff(n_groups), max_residual(n_groups), seconds(n_groups),
      converged(n_groups);
  for (int i = 0; i < n_groups; ++i) {
    const GroupSummary& summary = groups[i];
    group_index[i] = summary.group;
    iterations[i] = summary.iterations;
    phi_norm[i] = summary.phi_norm;
    abs_diff[i] = summary.abs_diff;
    rel_diff[i] = summary.rel_diff;
    max_residual[i] = summary.max_residual;
    seconds[i] = summary.seconds;
    converged[i] = summary.converged ? 1.0 : 0.0;
  }

  VerticalTable summary_table("Convergence Summary");
  summary_table.add_column("group", group_index);
  summary_table.add_column("iters", iterations);
  summary_table.add_column("|phi|", phi_norm);
  summary_table.add_column("|dphi|", abs_diff);
  summary_table.add_column("|dphi|/|phi|", rel_diff);
  summary_table.add_column("max|residual|", max_residual);
  summary_table.add_column("seconds", seconds);
  summary_table.add_column("converged", converged);
  out << summary_table.render_txt();

  KeyValueOutput totals("Convergence Totals");
  totals.add("Groups", n_groups);
  totals.add("Total iterations", history.totalIterations());
  totals.add("Total solve time (s)", std::format("{:.3f}", history.totalSeconds()));
  totals.add("Worst residual", formatSci(history.worstResidual()));
  totals.add("All groups converged", history.allConverged() ? "yes" : "NO");
  if (!history.allConverged()) {
    std::string unconverged;
    for (int group : history.unconvergedGroups()) {
      unconverged += (unconverged.empty() ? "" : ", ") + std::to_string(group);
    }
    totals.add("Unconverged groups", unconverged);
  }
  out << totals.render_txt();

  const std::vector<IterationRecord>& records = history.records();
  const auto n_records = static_cast<int>(records.size());
  std::vector<double> record_group(n_records), record_iteration(n_records), record_phi(n_records),
      record_abs(n_records), record_rel(n_records), record_res(n_records);
  for (int i = 0; i < n_records; ++i) {
    record_group[i] = records[i].group;
    record_iteration[i] = records[i].iteration;
    record_phi[i] = records[i].phi_norm;
    record_abs[i] = records[i].abs_diff;
    record_rel[i] = records[i].rel_diff;
    record_res[i] = records[i].max_residual;
  }

  VerticalTable history_table("Iteration History");
  history_table.add_column("group", record_group);
  history_table.add_column("iter", record_iteration);
  history_table.add_column("|phi|", record_phi);
  history_table.add_column("|dphi|", record_abs);
  history_table.add_column("|dphi|/|phi|", record_rel);
  history_table.add_column("max|residual|", record_res);
  out << history_table.render_txt();

  return out.str();
}

std::string formatResiduals(const std::vector<Eigen::MatrixXd>& residuals, const InputDeck& deck) {
  std::ostringstream out;
  for (std::size_t g = 0; g < residuals.size(); ++g) {
    out << formatCornerGrid("Angular Residual - Group " + std::to_string(g), residuals[g],
                            deck.mesh.n_x, ordinateLabels(deck));
  }
  return out.str();
}

} // namespace SolverFormatter
