// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "solver_formatter.h"

#include <array>
#include <ctime>
#include <iomanip>
#include <sstream>

#include "output.h"

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

std::string formatResiduals(const std::vector<Eigen::MatrixXd>& residuals, const InputDeck& deck) {
  std::ostringstream out;
  for (std::size_t g = 0; g < residuals.size(); ++g) {
    out << formatCornerGrid("Angular Residual - Group " + std::to_string(g), residuals[g],
                            deck.mesh.n_x, ordinateLabels(deck));
  }
  return out.str();
}

} // namespace SolverFormatter
