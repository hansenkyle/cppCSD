// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#ifndef INPUT_DECK_H
#define INPUT_DECK_H

#include <filesystem>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "cross_section.h"
#include "mesh.h"

/// @struct MaterialData
/// @brief Cross sections and stopping power for a single material, indexed by energy group
struct MaterialData {
  std::vector<double> sigma_t; // group total xs, size == num_groups
  std::vector<ScatterEntry>
      scattering; // sparse group-to-group scattering matrix, (from,to) in [0, num_groups)
  std::vector<double> stopping_power_average;  // group-average S, size == num_groups
  std::vector<double> stopping_power_boundary; // S at group boundaries, size == num_groups + 1
};

// Placeholder convergence criteria; more sophisticated criteria to follow.
struct ConvergenceCriteria {
  int max_iters = 0;
  double epsilon = 0.0;
};

// Discrete ordinates for angular quadrature: direction cosines and their
// integration weights, one entry per ordinate. Validates and normalizes
// itself at construction: mu must be strictly ascending, and w is always
// rescaled so its entries sum to 2 (logging the applied scale factor).
struct AngularQuadrature {
  AngularQuadrature(std::vector<double> mu, std::vector<double> w);

  std::vector<double> mu; // direction cosines, strictly ascending
  std::vector<double> w;  // quadrature weights, normalized to sum to 2
};

// A down/up pair of incoming angular flux values at one boundary, for one
// ordinate and one energy group -- down is the lower-energy edge, up the
// higher-energy edge (matching the same down/up convention used elsewhere,
// e.g. fe_space.h's Corner).
struct DownUp {
  double down;
  double up;
};

// Incoming angular flux at the domain's two spatial boundaries, indexed
// [ordinate][group]. Validates at construction that left and right are each
// exactly num_ordinates x num_groups.
struct BoundaryConditions {
  BoundaryConditions(std::vector<std::vector<DownUp>> left, std::vector<std::vector<DownUp>> right,
                     int num_ordinates, int num_groups);

  std::vector<std::vector<DownUp>> left;
  std::vector<std::vector<DownUp>> right;
};

/// @class InputDeck
/// @brief Holds formatted input data (via helper classes), read from YAML
class InputDeck {
public:
  // Reads and validates path_to_yaml, populating this deck's members.
  // Returns 0 if the file is valid; returns 1 early on the first error
  // found in the file (missing/malformed keys, undefined material
  // references, mismatched sizes, etc.).
  int read(const std::filesystem::path& path_to_yaml);

  // Individually reconfigure this deck after construction/read(), each
  // taking an already-validated domain object (its own constructor is what
  // enforces legality) so this is the only place callers need to look --
  // nothing downstream has to re-check what it's handed.
  //
  // Replacing the mesh invalidates xs (which was expanded against the old
  // one): setMesh clears it and logs that it did, rather than leaving a
  // stale cross section in place silently.
  void setMesh(Mesh new_mesh);
  void setAngularQuadrature(AngularQuadrature new_angular_quadrature);
  void setBoundaryConditions(BoundaryConditions new_boundary_conditions);
  void setConvergence(ConvergenceCriteria new_convergence);

  std::optional<Mesh> mesh;                      // set once read() succeeds
  std::vector<std::string> region_materials;     // material name per cell, size == mesh->n_x
  std::map<std::string, MaterialData> materials; // keyed by material name
  std::optional<CrossSection> xs; // per-cell expansion of materials; set once read() succeeds
  ConvergenceCriteria convergence;
  std::optional<AngularQuadrature> angular_quadrature;   // set once read() succeeds
  std::optional<BoundaryConditions> boundary_conditions; // set once read() succeeds
};

#endif
