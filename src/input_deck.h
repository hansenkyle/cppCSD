// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#ifndef INPUT_DECK_H
#define INPUT_DECK_H

#include <filesystem>

#include <Eigen/Dense>

// Holds all physical problem data (geometry, energy structure, cross
// sections, angular quadrature, boundary conditions) and is the trust
// boundary for it. Data is grouped into small structs by purpose
// (mesh/energy/angle/xs/bc), each with its own validate(). InputDeck::validate() runs every
// struct's own check plus the cross-struct shape checks (e.g. xs sized against mesh.n_x and
// energy.G) and is the single place that decides whether a deck's data is legal. It is not cached:
// fields are public and meant to be edited directly (e.g. deck.energy.G = ...), so validity is only
// meaningful at the moment validate() is called, which happens at read() and again wherever a deck
// is handed to a Solver.
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

  // Cross sections and stopping power, per energy group and spatial cell.
  // Row g of a table is that group's space-dependent vector; (g, c) is the
  // value for group g in cell c.
  struct Xs {
    Eigen::MatrixXd total;   // group total xs, rows=G, cols=n_x
    Eigen::MatrixXd scatter; // group isotropic scattering xs, rows=G, cols=n_x
    Eigen::MatrixXd S;       // group-average stopping power, rows=G, cols=n_x
    Eigen::MatrixXd S_bound; // stopping power at group boundaries, rows=G+1, cols=n_x

    // Checks every table's values are non-negative. Shape against
    // mesh.n_x/energy.G is a cross-struct concern, checked by
    // InputDeck::validate() instead.
    void validate() const;
  };

  // Incoming boundary angular flux, a (2*G) x M matrix: rows are indexed by
  // group first, then up/down within a group (up before down), columns by
  // ordinate. Row 2*g is group g's up value across ordinates, row 2*g + 1
  // is group g's down value. Which physical boundary (left/right) a value
  // applies to is determined elsewhere by the sign of angle.mu[m], not by
  // anything stored or validated here.
  struct BoundaryConditions {
    Eigen::MatrixXd values; // rows = 2 * G, cols = M

    // Views group g's two rows (up, then down) across all M ordinates.
    auto operator[](int g) {
      return values(Eigen::seq(2 * g, 2 * g + 1), Eigen::placeholders::all);
    }
    auto operator[](int g) const {
      return values(Eigen::seq(2 * g, 2 * g + 1), Eigen::placeholders::all);
    }

    // Checks every entry is non-negative. Shape against angle.M/energy.G is
    // a cross-struct concern, checked by InputDeck::validate() instead.
    void validate() const;
  };

  // Reads and validates path_to_yaml, populating this deck's members.
  // Returns 0 if the file is valid; returns 1 early on the first error
  // found in the file (missing/malformed keys, undefined material
  // references, mismatched sizes, etc.). Material names are used only to
  // expand per-material YAML into the per-cell tables above; they are not
  // retained on the deck.
  int read(const std::filesystem::path& path_to_yaml);

  // Runs every member struct's own validate(), plus cross-struct shape
  // checks (xs against mesh.n_x/energy.G, bc against angle.M/energy.G).
  // Throws std::runtime_error on the first violation found. May also apply
  // small in-place corrections as it validates (e.g. renormalizing
  // angle.w). Not cached -- call this wherever a deck's data is about to be
  // trusted (read() calls it already; Solver's constructor should call it
  // too).
  void validate();

  Mesh mesh;
  Energy energy;
  Angle angle;
  Xs xs;
  BoundaryConditions bc;
};

#endif
