// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#ifndef INPUT_DECK_H
#define INPUT_DECK_H

#include <filesystem>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

/// @brief Holds all problem data, cell-by cell. Energy structure, spatial mesh, cross sections,
/// angular quadrature, boundary conditions, external source
/// @details Includes validation functions on every member struct. Callers are expected to validate
/// all data before using it (i.e. Solver constructor). Validation checks for non-negativity,
/// ascending mesh, etc.
class InputDeck {
public:
  // Spatial grid: cell boundaries only. Energy structure lives in Energy.
  struct Mesh {
    int n_x = 0;                // number of spatial cells
    Eigen::VectorXd x_boundary; // cell boundary locations, strictly ascending, size n_x + 1
    Eigen::VectorXd dx;         // cell widths, derived from x_boundary by validate(), size n_x

    // Checks n_x > 0, x_boundary.size() == n_x + 1, and x_boundary strictly
    // ascending. Populates dx as the consecutive differences of x_boundary.
    void validate();
  };

  // Energy group structure.
  struct Energy {
    int G = 0;                  // number of energy groups
    Eigen::VectorXd E_boundary; // group boundary locations, strictly descending, size G + 1
    Eigen::VectorXd dE;         // group widths, derived from E_boundary by validate(), size G

    // Checks G > 0, E_boundary.size() == G + 1, and E_boundary strictly
    // descending. Populates dE as the consecutive (positive) differences of
    // E_boundary.
    void validate();
  };

  // Discrete-ordinates angular quadrature: direction cosines and their
  // integration weights, one entry per ordinate.
  struct Angle {
    int M = 0;          // number of ordinates
    Eigen::VectorXd mu; // direction cosines, strictly ascending, size M
    Eigen::VectorXd w;  // quadrature weights, normalized to sum to 2, size M

    // Checks M > 0, mu/w sizes == M, and mu strictly ascending. w's sum is
    // checked against 2: if the relative difference is within 1e-4, w is
    // rescaled in place to sum to exactly 2 (logging a warning); beyond
    // that tolerance, logs an error and throws.
    void validate();
  };

  /// @brief Cross sections and stopping power, per cell, per group.
  struct Xs {
    Eigen::MatrixXd total;                            // group total xs, rows=G, cols=n_x
    std::vector<Eigen::SparseMatrix<double>> scatter; // [cell], each G x G
    Eigen::MatrixXd S;       // group-average stopping power, rows=G, cols=n_x
    Eigen::MatrixXd S_bound; // stopping power at group boundaries, rows=G+1, cols=n_x

    // Checks total/S/S_bound values are non-negative, and every scatter
    // matrix entry is non-negative. Shape against mesh.n_x/energy.G is a
    // cross-struct concern, checked by InputDeck::validate() instead.
    void validate() const;
  };

  /// @brief Incoming angular flux at boundaries. Indexed by mu without knowledge of x
  /// Indexed by group using [int g], which returns [2 by M] matrix of values (U/D moments for all
  /// angles)
  struct BoundaryConditions {
    Eigen::MatrixXd values; // rows = 2 * G, cols = M

    // Views group g's two rows (up, then down) across all M ordinates.
    auto operator[](int g) {
      return values(Eigen::seq(2 * g, 2 * g + 1), Eigen::placeholders::all);
    }
    auto operator[](int g) const {
      return values(Eigen::seq(2 * g, 2 * g + 1), Eigen::placeholders::all);
    }
  };

  /// @brief External source, stored in the same format as angular flux: vector<Eigen Matrix>. Index
  /// using [group](cell/moment, angle)
  struct Source {
    std::vector<Eigen::MatrixXd> values; // size G, each 4 * mesh.n_x rows x angle.M cols

    // Checks every entry is non-negative. Shape against angle.M/mesh.n_x/energy.G is a
    // cross-struct concern, checked by InputDeck::validate() instead.
    void validate() const;
  };

  // Reads and validates path_to_yaml, populating this deck's members.
  // Returns 0 if the file is valid; returns 1 early on the first error
  // found in the file (missing/malformed keys, undefined material
  // references, mismatched sizes, etc.). Material names are used only to
  // expand per-material YAML into the per-cell tables above; they are not
  // retained on the deck.
  int read(const std::filesystem::path& path_to_yaml);

  /// @brief Calls each data member's own validate() and performs cross-struct checks (i.e. ensure
  /// xs has correct G). Call before using/copying data
  void validate();

  // std::string echo() const;

  Mesh mesh;
  Energy energy;
  Angle angle;
  Xs xs;
  BoundaryConditions bc;
  Source source;
};

#endif
