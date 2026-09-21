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

std::string formatRunMetadata() {
  OutputMetadata metadata("Run Metadata");
  metadata.addEntry("Run time", timestamp());
  return metadata.txt();
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
    return "--- Convergence Summary ---\nno iterations recorded\n";
  }

  constexpr OutputTable::Format kIndex{OutputTable::Notation::Integer, 0};
  constexpr OutputTable::Format kNorm{OutputTable::Notation::Scientific, 4};
  constexpr OutputTable::Format kSeconds{OutputTable::Notation::Fixed, 3};

  std::ostringstream out;

  const std::vector<GroupSummary>& groups = history.groups();
  const auto n_groups = static_cast<int>(groups.size());

  std::vector<double> group_index(n_groups), iterations(n_groups), phi_norm(n_groups),
      abs_diff(n_groups), rel_diff(n_groups), max_residual(n_groups), seconds(n_groups);
  std::vector<std::string> status(n_groups);
  for (int i = 0; i < n_groups; ++i) {
    const GroupSummary& summary = groups[i];
    group_index[i] = summary.group;
    iterations[i] = summary.iterations;
    phi_norm[i] = summary.phi_norm;
    abs_diff[i] = summary.abs_diff;
    rel_diff[i] = summary.rel_diff;
    max_residual[i] = summary.max_residual;
    seconds[i] = summary.seconds;
    status[i] = summary.converged ? "converged" : "MAX_ITER";
  }

  OutputTable summary_table("Convergence Summary", n_groups);
  summary_table.addColumn("group", kIndex, group_index);
  summary_table.addColumn("iters", kIndex, iterations);
  summary_table.addColumn("|phi|", kNorm, phi_norm);
  summary_table.addColumn("|dphi|", kNorm, abs_diff);
  summary_table.addColumn("|dphi|/|phi|", kNorm, rel_diff);
  summary_table.addColumn("max|residual|", kNorm, max_residual);
  summary_table.addColumn("seconds", kSeconds, seconds);
  summary_table.addTextColumn("status", status);
  out << summary_table.txt();

  OutputMetadata totals("Convergence Totals");
  totals.addEntry("Groups", std::to_string(n_groups));
  totals.addEntry("Total iterations", std::to_string(history.totalIterations()));
  totals.addEntry("Total solve time (s)", std::format("{:.3f}", history.totalSeconds()));
  totals.addEntry("Worst residual", formatSci(history.worstResidual()));
  totals.addEntry("All groups converged", history.allConverged() ? "yes" : "NO");
  if (!history.allConverged()) {
    std::string unconverged;
    for (int group : history.unconvergedGroups()) {
      unconverged += (unconverged.empty() ? "" : ", ") + std::to_string(group);
    }
    totals.addEntry("Unconverged groups", unconverged);
  }
  out << totals.txt();

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

  OutputTable history_table("Iteration History", n_records);
  history_table.addColumn("group", kIndex, record_group);
  history_table.addColumn("iter", kIndex, record_iteration);
  history_table.addColumn("|phi|", kNorm, record_phi);
  history_table.addColumn("|dphi|", kNorm, record_abs);
  history_table.addColumn("|dphi|/|phi|", kNorm, record_rel);
  history_table.addColumn("max|residual|", kNorm, record_res);
  out << history_table.txt();

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
