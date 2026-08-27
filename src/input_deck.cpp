// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "input_deck.h"

#include <stdexcept>

#include <yaml-cpp/yaml.h>

#include "logger.h"

namespace {

YAML::Node requireNode(const YAML::Node& parent, const std::string& key) {
  const YAML::Node node = parent[key];
  if (!node) {
    throw std::runtime_error("input file is missing required key '" + key + "'");
  }
  return node;
}

void requireSize(const std::vector<double>& values, std::size_t expected, const std::string& name) {
  if (values.size() != expected) {
    throw std::runtime_error(name + " has size " + std::to_string(values.size()) + ", expected " +
                             std::to_string(expected));
  }
}

// Expands a per-material, per-group field (selected via `field`) into a
// per-group, per-cell table by looking up each cell's material.
std::vector<std::vector<double>>
expandByRegion(const std::vector<std::string>& region_materials,
               const std::map<std::string, MaterialData>& materials, std::size_t num_groups,
               std::vector<double> MaterialData::* field) {
  std::vector<std::vector<double>> table(num_groups, std::vector<double>(region_materials.size()));
  for (std::size_t cell = 0; cell < region_materials.size(); ++cell) {
    const std::vector<double>& values = materials.at(region_materials[cell]).*field;
    for (std::size_t g = 0; g < num_groups; ++g) {
      table[g][cell] = values[g];
    }
  }
  return table;
}

// Expands each cell's material into its (shared, unmodified) sparse
// scattering-matrix entry list -- unlike expandByRegion, there's no
// per-group axis to transpose into, so this can't reuse it.
std::vector<std::vector<ScatterEntry>>
expandScatteringByRegion(const std::vector<std::string>& region_materials,
                         const std::map<std::string, MaterialData>& materials) {
  std::vector<std::vector<ScatterEntry>> table(region_materials.size());
  for (std::size_t cell = 0; cell < region_materials.size(); ++cell) {
    table[cell] = materials.at(region_materials[cell]).scattering;
  }
  return table;
}

// Parses a material's `scattering:` block: a list of {from, to, value}
// entries. Only nonzero entries need to appear -- see docs/input-deck.md.
std::vector<ScatterEntry> parseScattering(const YAML::Node& node, int num_groups,
                                          const std::string& name) {
  std::vector<ScatterEntry> entries;
  entries.reserve(node.size());
  for (const YAML::Node& entry_node : node) {
    ScatterEntry entry;
    entry.from = requireNode(entry_node, "from").as<int>();
    entry.to = requireNode(entry_node, "to").as<int>();
    entry.value = requireNode(entry_node, "value").as<double>();

    if (entry.from < 0 || entry.from >= num_groups || entry.to < 0 || entry.to >= num_groups) {
      throw std::runtime_error("material '" + name + "' scattering entry (from=" +
                               std::to_string(entry.from) + ", to=" + std::to_string(entry.to) +
                               ") out of range [0, " + std::to_string(num_groups) + ")");
    }
    entries.push_back(entry);
  }

  for (std::size_t a = 0; a < entries.size(); ++a) {
    for (std::size_t b = a + 1; b < entries.size(); ++b) {
      if (entries[a].from == entries[b].from && entries[a].to == entries[b].to) {
        throw std::runtime_error("material '" + name +
                                 "' has a duplicate scattering entry for "
                                 "(from=" +
                                 std::to_string(entries[a].from) +
                                 ", to=" + std::to_string(entries[a].to) + ")");
      }
    }
  }

  return entries;
}

MaterialData parseMaterial(const YAML::Node& node, const std::string& name,
                           std::size_t num_groups) {
  MaterialData material;
  material.sigma_t = requireNode(node, "sigma_t").as<std::vector<double>>();
  material.scattering =
      parseScattering(requireNode(node, "scattering"), static_cast<int>(num_groups), name);

  const YAML::Node stopping_power = requireNode(node, "stopping_power");
  material.stopping_power_average =
      requireNode(stopping_power, "group_average").as<std::vector<double>>();
  material.stopping_power_boundary =
      requireNode(stopping_power, "group_boundary").as<std::vector<double>>();

  requireSize(material.sigma_t, num_groups, "material '" + name + "' sigma_t");
  requireSize(material.stopping_power_average, num_groups,
              "material '" + name + "' stopping_power.group_average");
  requireSize(material.stopping_power_boundary, num_groups + 1,
              "material '" + name + "' stopping_power.group_boundary");

  return material;
}

// Parses boundary_conditions.<side>.down/up into a [ordinate][group] table
// of DownUp pairs. Only checks that down and up agree with each other in
// shape -- whether that shape matches the deck's actual ordinate/group
// counts is BoundaryConditions' own constructor's job, not this function's.
std::vector<std::vector<DownUp>> parseBoundarySide(const YAML::Node& side_node) {
  const std::vector<std::vector<double>> down =
      requireNode(side_node, "down").as<std::vector<std::vector<double>>>();
  const std::vector<std::vector<double>> up =
      requireNode(side_node, "up").as<std::vector<std::vector<double>>>();

  if (down.size() != up.size()) {
    throw std::runtime_error("boundary_conditions: down and up have different ordinate counts");
  }
  std::vector<std::vector<DownUp>> side(down.size());
  for (std::size_t m = 0; m < down.size(); ++m) {
    if (down[m].size() != up[m].size()) {
      throw std::runtime_error("boundary_conditions: down and up have different group counts");
    }
    side[m].resize(down[m].size());
    for (std::size_t g = 0; g < down[m].size(); ++g) {
      side[m][g] = DownUp{down[m][g], up[m][g]};
    }
  }
  return side;
}

void requireShapeDownUp(const std::vector<std::vector<DownUp>>& table, int expected_rows,
                        int expected_cols, const std::string& name) {
  if (static_cast<int>(table.size()) != expected_rows) {
    throw std::invalid_argument(name + " has " + std::to_string(table.size()) +
                                " row(s), expected " + std::to_string(expected_rows));
  }
  for (const std::vector<DownUp>& row : table) {
    if (static_cast<int>(row.size()) != expected_cols) {
      throw std::invalid_argument(name + " row has size " + std::to_string(row.size()) +
                                  ", expected " + std::to_string(expected_cols));
    }
  }
}

} // namespace

AngularQuadrature::AngularQuadrature(std::vector<double> mu_in, std::vector<double> w_in)
    : mu(std::move(mu_in)), w(std::move(w_in)) {
  if (w.size() != mu.size()) {
    throw std::invalid_argument("AngularQuadrature: w has size " + std::to_string(w.size()) +
                                ", expected " + std::to_string(mu.size()));
  }

  for (std::size_t m = 1; m < mu.size(); ++m) {
    if (mu[m] <= mu[m - 1]) {
      throw std::invalid_argument("AngularQuadrature: mu must be strictly ascending");
    }
  }

  double sum = 0.0;
  for (double weight : w) {
    sum += weight;
  }
  if (sum <= 0.0) {
    throw std::invalid_argument("AngularQuadrature: w must sum to a positive value");
  }

  const double scale = 2.0 / sum;
  for (double& weight : w) {
    weight *= scale;
  }
  LDCSD_LOG_INFO("normalized AngularQuadrature.w: sum was " + std::to_string(sum) + ", scaled by " +
                 std::to_string(scale) + " to sum to 2");
}

BoundaryConditions::BoundaryConditions(std::vector<std::vector<DownUp>> left_in,
                                       std::vector<std::vector<DownUp>> right_in, int num_ordinates,
                                       int num_groups)
    : left(std::move(left_in)), right(std::move(right_in)) {
  requireShapeDownUp(left, num_ordinates, num_groups, "BoundaryConditions: left");
  requireShapeDownUp(right, num_ordinates, num_groups, "BoundaryConditions: right");
}

void InputDeck::setMesh(Mesh new_mesh) {
  mesh.emplace(std::move(new_mesh));
  if (xs.has_value()) {
    xs.reset();
    LDCSD_LOG_INFO("cleared xs: mesh was replaced, so the existing cross-section expansion is no "
                   "longer valid against it");
  }
}

void InputDeck::setAngularQuadrature(AngularQuadrature new_angular_quadrature) {
  angular_quadrature.emplace(std::move(new_angular_quadrature));
}

void InputDeck::setBoundaryConditions(BoundaryConditions new_boundary_conditions) {
  boundary_conditions.emplace(std::move(new_boundary_conditions));
}

void InputDeck::setConvergence(ConvergenceCriteria new_convergence) {
  convergence = new_convergence;
}

int InputDeck::read(const std::filesystem::path& path_to_yaml) {
  try {
    const YAML::Node root = YAML::LoadFile(path_to_yaml.string());
    LDCSD_LOG_TRACE("parsed '" + path_to_yaml.string() + "' as YAML");

    std::vector<double> x_boundary = requireNode(root, "spatial_mesh").as<std::vector<double>>();
    std::vector<double> E_boundary = requireNode(root, "energy_mesh").as<std::vector<double>>();
    setMesh(Mesh(std::move(x_boundary), std::move(E_boundary)));

    const YAML::Node regions = requireNode(root, "regions");
    region_materials = requireNode(regions, "materials").as<std::vector<std::string>>();
    if (static_cast<int>(region_materials.size()) != mesh->n_x) {
      throw std::runtime_error("regions.materials has size " +
                               std::to_string(region_materials.size()) + ", expected " +
                               std::to_string(mesh->n_x) + " (one per spatial cell)");
    }

    const YAML::Node materials_node = requireNode(root, "materials");
    for (const auto& entry : materials_node) {
      const std::string name = entry.first.as<std::string>();
      materials[name] = parseMaterial(entry.second, name, static_cast<std::size_t>(mesh->G));
      LDCSD_LOG_DEBUG("parsed material '" + name + "'");
    }

    for (const std::string& name : region_materials) {
      if (!materials.contains(name)) {
        throw std::runtime_error("region references undefined material '" + name + "'");
      }
    }

    const auto num_groups = static_cast<std::size_t>(mesh->G);
    xs.emplace(*mesh,
               expandByRegion(region_materials, materials, num_groups, &MaterialData::sigma_t),
               expandScatteringByRegion(region_materials, materials),
               expandByRegion(region_materials, materials, num_groups,
                              &MaterialData::stopping_power_average),
               expandByRegion(region_materials, materials, num_groups + 1,
                              &MaterialData::stopping_power_boundary),
               region_materials);

    const YAML::Node angular_quadrature_node = requireNode(root, "angular_quadrature");
    std::vector<double> mu = requireNode(angular_quadrature_node, "mu").as<std::vector<double>>();
    std::vector<double> w = requireNode(angular_quadrature_node, "w").as<std::vector<double>>();
    setAngularQuadrature(AngularQuadrature(std::move(mu), std::move(w)));

    const YAML::Node boundary_conditions_node = requireNode(root, "boundary_conditions");
    std::vector<std::vector<DownUp>> left =
        parseBoundarySide(requireNode(boundary_conditions_node, "left"));
    std::vector<std::vector<DownUp>> right =
        parseBoundarySide(requireNode(boundary_conditions_node, "right"));
    setBoundaryConditions(BoundaryConditions(std::move(left), std::move(right),
                                             static_cast<int>(angular_quadrature->mu.size()),
                                             mesh->G));

    const YAML::Node convergence_node = requireNode(root, "convergence");
    ConvergenceCriteria new_convergence;
    new_convergence.max_iters = requireNode(convergence_node, "max_iters").as<int>();
    new_convergence.epsilon = requireNode(convergence_node, "epsilon").as<double>();
    setConvergence(new_convergence);
  } catch (const std::exception& e) {
    LDCSD_LOG_ERROR(std::string("failed to read input deck '") + path_to_yaml.string() +
                    "': " + e.what());
    return 1;
  }

  LDCSD_LOG_INFO("read input deck '" + path_to_yaml.string() + "': " + std::to_string(mesh->n_x) +
                 " cells, " + std::to_string(mesh->G) + " groups, " +
                 std::to_string(materials.size()) + " materials, " +
                 std::to_string(angular_quadrature->mu.size()) + " ordinates");
  return 0;
}
