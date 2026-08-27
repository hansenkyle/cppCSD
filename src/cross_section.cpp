// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "cross_section.h"

#include <stdexcept>

namespace {

void requireShape(const std::vector<std::vector<double>>& table, std::size_t expected_groups,
                  std::size_t expected_cells, const std::string& name) {
  if (table.size() != expected_groups) {
    throw std::invalid_argument("CrossSection: " + name + " has " + std::to_string(table.size()) +
                                " group row(s), expected " + std::to_string(expected_groups));
  }
  for (const std::vector<double>& row : table) {
    if (row.size() != expected_cells) {
      throw std::invalid_argument("CrossSection: " + name + " row has size " +
                                  std::to_string(row.size()) + ", expected " +
                                  std::to_string(expected_cells));
    }
  }
}

void requireNonNegative(const std::vector<std::vector<double>>& table, const std::string& name) {
  for (const std::vector<double>& row : table) {
    for (double value : row) {
      if (value < 0.0) {
        throw std::invalid_argument("CrossSection: " + name + " must be non-negative");
      }
    }
  }
}

void requireValidScattering(const std::vector<std::vector<ScatterEntry>>& scattering,
                            std::size_t expected_cells, int num_groups) {
  if (scattering.size() != expected_cells) {
    throw std::invalid_argument("CrossSection: scattering has " +
                                std::to_string(scattering.size()) + " cell(s), expected " +
                                std::to_string(expected_cells));
  }
  for (const std::vector<ScatterEntry>& cell_entries : scattering) {
    for (const ScatterEntry& entry : cell_entries) {
      if (entry.from < 0 || entry.from >= num_groups || entry.to < 0 || entry.to >= num_groups) {
        throw std::invalid_argument(
            "CrossSection: scattering entry (from=" + std::to_string(entry.from) + ", to=" +
            std::to_string(entry.to) + ") out of range [0, " + std::to_string(num_groups) + ")");
      }
      if (entry.value < 0.0) {
        throw std::invalid_argument(
            "CrossSection: scattering entry (from=" + std::to_string(entry.from) +
            ", to=" + std::to_string(entry.to) + ") must be non-negative");
      }
    }
  }
}

} // namespace

CrossSection::CrossSection(const Mesh& mesh_in, std::vector<std::vector<double>> total_in,
                           std::vector<std::vector<ScatterEntry>> scattering_in,
                           std::vector<std::vector<double>> stop_power_in,
                           std::vector<std::vector<double>> stop_power_boundary_in,
                           std::vector<std::string> material_in)
    : mesh(mesh_in), total(std::move(total_in)), scattering(std::move(scattering_in)),
      stop_power(std::move(stop_power_in)), stop_power_boundary(std::move(stop_power_boundary_in)),
      material(std::move(material_in)) {
  const auto num_groups = static_cast<std::size_t>(mesh.G);
  const auto num_cells = static_cast<std::size_t>(mesh.n_x);

  requireShape(total, num_groups, num_cells, "total");
  requireShape(stop_power, num_groups, num_cells, "stop_power");
  requireShape(stop_power_boundary, num_groups + 1, num_cells, "stop_power_boundary");
  requireValidScattering(scattering, num_cells, mesh.G);

  requireNonNegative(total, "total");
  requireNonNegative(stop_power, "stop_power");
  requireNonNegative(stop_power_boundary, "stop_power_boundary");

  if (material.size() != num_cells) {
    throw std::invalid_argument("CrossSection: material has size " +
                                std::to_string(material.size()) + ", expected " +
                                std::to_string(num_cells));
  }
}
