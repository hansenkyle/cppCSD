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

// Renders a (4 * n_x) x G matrix as a grid: one column per spatial cell,
// one row per energy group, with each cell's four values (rows
// 4c..4c+3 = up_left, up_right, down_left, down_right) shown as a 2x2
// cluster:
//   up_left    up_right
//   down_left  down_right
std::string formatCornerGrid(const std::string& title, const Eigen::MatrixXd& data, int n_x,
                             int G) {
  std::vector<std::vector<std::array<std::string, 4>>> formatted(
      G, std::vector<std::array<std::string, 4>>(n_x));
  std::size_t width = 0;
  for (int g = 0; g < G; ++g) {
    for (int c = 0; c < n_x; ++c) {
      for (int k = 0; k < 4; ++k) {
        std::string text = formatSci(data(4 * c + k, g));
        width = std::max(width, text.size());
        formatted[g][c][k] = std::move(text);
      }
    }
  }
  const int col_width = static_cast<int>(width) + 2;

  std::ostringstream out;
  out << "--- " << title << " ---\n";
  out << "columns are spatial cells 0.." << (n_x - 1)
      << "; each entry is [up_left up_right / down_left down_right]\n";
  for (int g = 0; g < G; ++g) {
    out << "group " << g << ":\n";
    for (int row = 0; row < 2; ++row) {
      for (int c = 0; c < n_x; ++c) {
        const std::array<std::string, 4>& cell = formatted[g][c];
        out << std::setw(col_width) << cell[2 * row] << std::setw(col_width) << cell[2 * row + 1];
      }
      out << "\n";
    }
  }
  return out.str();
}

// Per-ordinate title shared by the angular flux and residual blocks.
std::string ordinateTitle(const std::string& label, int m, const InputDeck& deck) {
  std::ostringstream title;
  title << label << " - Ordinate " << m << " (mu=" << formatSci(deck.angle.mu[m])
        << ", w=" << formatSci(deck.angle.w[m]) << ")";
  return title.str();
}

} // namespace

namespace SolverFormatter {

std::string formatRunMetadata() {
  OutputMetadata metadata("Run Metadata");
  metadata.addEntry("Run time", timestamp());
  return metadata.txt();
}

std::string formatResults(const Solver::Results& results, const InputDeck& deck) {
  std::ostringstream out;
  out << formatCornerGrid("Scalar Flux", results.scalar_flux, deck.mesh.n_x, deck.energy.G);
  for (std::size_t m = 0; m < results.angular_flux.size(); ++m) {
    out << formatCornerGrid(ordinateTitle("Angular Flux", static_cast<int>(m), deck),
                            results.angular_flux[m], deck.mesh.n_x, deck.energy.G);
  }
  return out.str();
}

std::string formatResiduals(const std::vector<Eigen::MatrixXd>& residuals, const InputDeck& deck) {
  std::ostringstream out;
  for (std::size_t m = 0; m < residuals.size(); ++m) {
    out << formatCornerGrid(ordinateTitle("Angular Residual", static_cast<int>(m), deck),
                            residuals[m], deck.mesh.n_x, deck.energy.G);
  }
  return out.str();
}

} // namespace SolverFormatter
