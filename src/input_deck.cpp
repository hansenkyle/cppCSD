// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "input_deck.h"

#include <cmath>
#include <stdexcept>
#include <string>

#include "logger.h"

namespace {
constexpr double kAngleWeightRelTol = 1e-4;

void requireNonNegative(const Eigen::MatrixXd& values, const std::string& name) {
  if (values.size() > 0 && values.minCoeff() < 0.0) {
    throw std::runtime_error(name + " must be non-negative");
  }
}
} // namespace

void InputDeck::Mesh::validate() {
  if (n_x <= 0) {
    throw std::runtime_error("mesh.n_x must be positive");
  }
  if (x_boundary.size() != n_x + 1) {
    throw std::runtime_error("mesh.x_boundary has size " + std::to_string(x_boundary.size()) +
                             ", expected n_x + 1 = " + std::to_string(n_x + 1));
  }
  for (int i = 1; i <= n_x; ++i) {
    if (x_boundary[i] <= x_boundary[i - 1]) {
      throw std::runtime_error("mesh.x_boundary must be strictly ascending");
    }
  }

  dx = x_boundary.tail(n_x) - x_boundary.head(n_x);
}

void InputDeck::Energy::validate() {
  if (G <= 0) {
    throw std::runtime_error("energy.G must be positive");
  }
  if (E_boundary.size() != G + 1) {
    throw std::runtime_error("energy.E_boundary has size " + std::to_string(E_boundary.size()) +
                             ", expected G + 1 = " + std::to_string(G + 1));
  }
  for (int i = 1; i <= G; ++i) {
    if (E_boundary[i] >= E_boundary[i - 1]) {
      throw std::runtime_error("energy.E_boundary must be strictly descending");
    }
  }

  dE = E_boundary.head(G) - E_boundary.tail(G);
  for (int i = 0; i < G; ++i) {
    if (dE[i] <= 0.0) {
      throw std::runtime_error("energy.dE[" + std::to_string(i) + "] must be positive");
    }
  }
}

void InputDeck::Angle::validate() {
  if (M <= 0) {
    throw std::runtime_error("angle.M must be positive");
  }
  if (mu.size() != M) {
    throw std::runtime_error("angle.mu has size " + std::to_string(mu.size()) +
                             ", expected M = " + std::to_string(M));
  }
  if (w.size() != M) {
    throw std::runtime_error("angle.w has size " + std::to_string(w.size()) +
                             ", expected M = " + std::to_string(M));
  }
  for (int m = 1; m < M; ++m) {
    if (mu[m] <= mu[m - 1]) {
      throw std::runtime_error("angle.mu must be strictly ascending");
    }
  }

  const double sum = w.sum();
  const double rel_diff = std::abs(sum - 2.0) / 2.0;
  if (rel_diff > kAngleWeightRelTol) {
    LDCSD_LOG_ERROR("angle.w sums to " + std::to_string(sum) + ", a relative difference of " +
                    std::to_string(rel_diff) + " from 2 (exceeds tolerance " +
                    std::to_string(kAngleWeightRelTol) + ")");
    throw std::runtime_error("angle.w must sum to 2 within a relative tolerance of " +
                             std::to_string(kAngleWeightRelTol));
  }
  if (rel_diff > 0.0) {
    const double scale = 2.0 / sum;
    LDCSD_LOG_WARN("angle.w summed to " + std::to_string(sum) + " (relative difference " +
                   std::to_string(rel_diff) + "); normalizing by " + std::to_string(scale) +
                   " to sum to 2");
    w *= scale;
  }
}

void InputDeck::Xs::validate() const {
  requireNonNegative(total, "xs.total");
  requireNonNegative(scatter, "xs.scatter");
  requireNonNegative(S, "xs.S");
  requireNonNegative(S_bound, "xs.S_bound");
}

void InputDeck::BoundaryConditions::validate() const { requireNonNegative(values, "bc.values"); }

void InputDeck::validate() {}

int InputDeck::read(const std::filesystem::path& path_to_yaml) {
  (void)path_to_yaml;
  return 0;
}
