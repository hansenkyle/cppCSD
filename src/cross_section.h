// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#ifndef CROSS_SECTION_H
#define CROSS_SECTION_H

#include <string>
#include <vector>

#include "mesh.h"

/// @class CrossSection
/// @brief Per-cell, per-group cross sections and stopping power
/// @details Per-cell, per-group cross sections and stopping power, parameterized using Mesh.
/// Indexing by [group] alone yields a space-dependent vector of cross sections, indexing
/// [group][cell] yields a single value for that cell-group pair
class CrossSection {
public:
  // Constructs a CrossSection over `mesh` from already-expanded per-cell,
  // per-group data. `total`, `scattering`, and `stop_power` must each have
  // mesh.G rows of mesh.n_x entries; `stop_power_boundary` must have
  // mesh.G + 1 rows of mesh.n_x entries (S at group boundaries); `material`
  // must have mesh.n_x entries. Every numeric value must be non-negative.
  // `mesh` must outlive this CrossSection.
  CrossSection(const Mesh& mesh, std::vector<std::vector<double>> total,
               std::vector<std::vector<double>> scattering,
               std::vector<std::vector<double>> stop_power,
               std::vector<std::vector<double>> stop_power_boundary,
               std::vector<std::string> material);

  const Mesh& mesh;

  std::vector<std::vector<double>> total;      // [group][cell], size G x n_x
  std::vector<std::vector<double>> scattering; // [group][cell], size G x n_x
  std::vector<std::vector<double>> stop_power; // group-average S, [group][cell], size G x n_x
  std::vector<std::vector<double>>
      stop_power_boundary; // S at group boundaries, [group boundary][cell], size (G+1) x n_x
  std::vector<std::string> material; // material name, [cell], size n_x
};

#endif
