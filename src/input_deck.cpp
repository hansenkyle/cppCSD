// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "input_deck.h"

#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>

#include "logger.h"

namespace {
constexpr double kAngleWeightRelTol = 1e-4;

void requireNonNegative(const Eigen::MatrixXd& values, const std::string& name) {
  if (values.size() > 0 && values.minCoeff() < 0.0) {
    throw std::runtime_error(name + " must be non-negative");
  }
}

YAML::Node requireNode(const YAML::Node& parent, const std::string& key) {
  const YAML::Node node = parent[key];
  if (!node) {
    throw std::runtime_error("input file is missing required key '" + key + "'");
  }
  return node;
}

Eigen::VectorXd toVector(const std::vector<double>& values) {
  return Eigen::Map<const Eigen::VectorXd>(values.data(), static_cast<Eigen::Index>(values.size()));
}

void requireSize(Eigen::Index actual, Eigen::Index expected, const std::string& name) {
  if (actual != expected) {
    throw std::runtime_error(name + " has size " + std::to_string(actual) + ", expected " +
                             std::to_string(expected));
  }
}

// Per-material cross sections and stopping power, indexed by energy group.
// Parsing scratch only -- material names and this data are never retained
// on the InputDeck, only the per-cell tables/lists they expand into.
struct Material {
  Eigen::VectorXd sigma_t;                             // size G
  std::vector<InputDeck::Xs::ScatterEntry> scattering; // sparse group-to-group entries
  Eigen::VectorXd stopping_power_avg;                  // size G
  Eigen::VectorXd stopping_power_bnd;                  // size G + 1
};

// Group indices in the scattering: list are 1-indexed in YAML (matching how
// a person would naturally refer to "group 1"); stored 0-indexed internally
// like everything else. Checks from/to fall within [1, G] and that no
// (from, to) pair repeats -- a duplicate is almost certainly a typo, so it's
// rejected rather than summed.
std::vector<InputDeck::Xs::ScatterEntry> parseScattering(const YAML::Node& node,
                                                         const std::string& material_name, int G) {
  std::vector<InputDeck::Xs::ScatterEntry> entries;
  for (const auto& entry : node) {
    const int from_1indexed = requireNode(entry, "from").as<int>();
    const int to_1indexed = requireNode(entry, "to").as<int>();
    const double value = requireNode(entry, "value").as<double>();

    if (from_1indexed < 1 || from_1indexed > G) {
      throw std::runtime_error("material '" + material_name + "' scattering 'from' group " +
                               std::to_string(from_1indexed) + " is out of range [1, " +
                               std::to_string(G) + "]");
    }
    if (to_1indexed < 1 || to_1indexed > G) {
      throw std::runtime_error("material '" + material_name + "' scattering 'to' group " +
                               std::to_string(to_1indexed) + " is out of range [1, " +
                               std::to_string(G) + "]");
    }

    const int from = from_1indexed - 1;
    const int to = to_1indexed - 1;
    for (const InputDeck::Xs::ScatterEntry& existing : entries) {
      if (existing.from == from && existing.to == to) {
        throw std::runtime_error(
            "material '" + material_name + "' has a duplicate scattering entry (from=" +
            std::to_string(from_1indexed) + ", to=" + std::to_string(to_1indexed) + ")");
      }
    }
    entries.push_back(InputDeck::Xs::ScatterEntry{from, to, value});
  }
  return entries;
}

Material parseMaterial(const YAML::Node& node, const std::string& name, int G) {
  Material material;
  material.sigma_t = toVector(requireNode(node, "sigma_t").as<std::vector<double>>());
  material.scattering = parseScattering(requireNode(node, "scattering"), name, G);

  const YAML::Node stopping_power = requireNode(node, "stopping_power");
  material.stopping_power_avg =
      toVector(requireNode(stopping_power, "group_average").as<std::vector<double>>());
  material.stopping_power_bnd =
      toVector(requireNode(stopping_power, "group_boundary").as<std::vector<double>>());

  requireSize(material.sigma_t.size(), G, "material '" + name + "' sigma_t");
  requireSize(material.stopping_power_avg.size(), G,
              "material '" + name + "' stopping_power.group_average");
  requireSize(material.stopping_power_bnd.size(), G + 1,
              "material '" + name + "' stopping_power.group_boundary");
  return material;
}

// Expands a per-material, per-group field (selected via `field`) into a
// per-group, per-cell table by looking up each cell's material.
Eigen::MatrixXd expandByRegion(const std::vector<std::string>& region_materials,
                               const std::map<std::string, Material>& materials, Eigen::Index rows,
                               Eigen::VectorXd Material::* field) {
  Eigen::MatrixXd table(rows, static_cast<Eigen::Index>(region_materials.size()));
  for (std::size_t cell = 0; cell < region_materials.size(); ++cell) {
    table.col(static_cast<Eigen::Index>(cell)) = materials.at(region_materials[cell]).*field;
  }
  return table;
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
  requireNonNegative(S, "xs.S");
  requireNonNegative(S_bound, "xs.S_bound");

  for (const std::vector<ScatterEntry>& cell_entries : scatter) {
    for (const ScatterEntry& entry : cell_entries) {
      if (entry.value < 0.0) {
        throw std::runtime_error("xs.scatter value must be non-negative");
      }
    }
  }
}

void InputDeck::BoundaryConditions::validate() const { requireNonNegative(values, "bc.values"); }

void InputDeck::validate() {
  mesh.validate();
  energy.validate();
  angle.validate();
  xs.validate();
  bc.validate();

  if (xs.total.rows() != energy.G || xs.S.rows() != energy.G) {
    throw std::runtime_error("xs.total/S must have energy.G = " + std::to_string(energy.G) +
                             " rows");
  }
  if (xs.S_bound.rows() != energy.G + 1) {
    throw std::runtime_error("xs.S_bound must have energy.G + 1 = " + std::to_string(energy.G + 1) +
                             " rows");
  }
  if (xs.total.cols() != mesh.n_x || xs.S.cols() != mesh.n_x || xs.S_bound.cols() != mesh.n_x) {
    throw std::runtime_error("xs.total/S/S_bound must have mesh.n_x = " + std::to_string(mesh.n_x) +
                             " columns");
  }
  if (static_cast<int>(xs.scatter.size()) != mesh.n_x) {
    throw std::runtime_error("xs.scatter must have mesh.n_x = " + std::to_string(mesh.n_x) +
                             " cells");
  }

  if (bc.values.rows() != 2 * energy.G) {
    throw std::runtime_error("bc.values must have 2 * energy.G = " + std::to_string(2 * energy.G) +
                             " rows");
  }
  if (bc.values.cols() != angle.M) {
    throw std::runtime_error("bc.values must have angle.M = " + std::to_string(angle.M) +
                             " columns");
  }
}

int InputDeck::read(const std::filesystem::path& path_to_yaml) {
  try {
    const YAML::Node root = YAML::LoadFile(path_to_yaml.string());
    LDCSD_LOG_TRACE("parsed '" + path_to_yaml.string() + "' as YAML");

    const std::vector<double> x_boundary_raw =
        requireNode(root, "spatial_mesh").as<std::vector<double>>();
    mesh.x_boundary = toVector(x_boundary_raw);
    mesh.n_x = static_cast<int>(x_boundary_raw.size()) - 1;

    const std::vector<double> E_boundary_raw =
        requireNode(root, "energy_mesh").as<std::vector<double>>();
    energy.E_boundary = toVector(E_boundary_raw);
    energy.G = static_cast<int>(E_boundary_raw.size()) - 1;

    const YAML::Node regions = requireNode(root, "regions");
    const std::vector<std::string> region_materials =
        requireNode(regions, "materials").as<std::vector<std::string>>();
    if (static_cast<int>(region_materials.size()) != mesh.n_x) {
      throw std::runtime_error("regions.materials has size " +
                               std::to_string(region_materials.size()) + ", expected " +
                               std::to_string(mesh.n_x) + " (one per spatial cell)");
    }

    std::map<std::string, Material> materials;
    const YAML::Node materials_node = requireNode(root, "materials");
    for (const auto& entry : materials_node) {
      const std::string name = entry.first.as<std::string>();
      materials[name] = parseMaterial(entry.second, name, energy.G);
      LDCSD_LOG_DEBUG("parsed material '" + name + "'");
    }
    for (const std::string& name : region_materials) {
      if (!materials.contains(name)) {
        throw std::runtime_error("region references undefined material '" + name + "'");
      }
    }

    xs.total = expandByRegion(region_materials, materials, energy.G, &Material::sigma_t);
    xs.S = expandByRegion(region_materials, materials, energy.G, &Material::stopping_power_avg);
    xs.S_bound =
        expandByRegion(region_materials, materials, energy.G + 1, &Material::stopping_power_bnd);

    xs.scatter.resize(region_materials.size());
    for (std::size_t cell = 0; cell < region_materials.size(); ++cell) {
      xs.scatter[cell] = materials.at(region_materials[cell]).scattering;
    }

    const YAML::Node angular_quadrature_node = requireNode(root, "angular_quadrature");
    const std::vector<double> mu_raw =
        requireNode(angular_quadrature_node, "mu").as<std::vector<double>>();
    angle.mu = toVector(mu_raw);
    angle.w = toVector(requireNode(angular_quadrature_node, "w").as<std::vector<double>>());
    angle.M = static_cast<int>(mu_raw.size());

    const YAML::Node bc_node = requireNode(root, "boundary_conditions");
    const std::vector<std::vector<double>> up =
        requireNode(bc_node, "up").as<std::vector<std::vector<double>>>();
    const std::vector<std::vector<double>> down =
        requireNode(bc_node, "down").as<std::vector<std::vector<double>>>();
    requireSize(static_cast<Eigen::Index>(down.size()), static_cast<Eigen::Index>(up.size()),
                "boundary_conditions.down");

    const auto bc_G = static_cast<Eigen::Index>(up.size());
    const auto bc_M = bc_G > 0 ? static_cast<Eigen::Index>(up[0].size()) : 0;
    bc.values = Eigen::MatrixXd(2 * bc_G, bc_M);
    for (Eigen::Index g = 0; g < bc_G; ++g) {
      requireSize(static_cast<Eigen::Index>(up[g].size()), bc_M,
                  "boundary_conditions.up row " + std::to_string(g));
      requireSize(static_cast<Eigen::Index>(down[g].size()), bc_M,
                  "boundary_conditions.down row " + std::to_string(g));
      for (Eigen::Index m = 0; m < bc_M; ++m) {
        bc.values(2 * g, m) = up[g][m];
        bc.values(2 * g + 1, m) = down[g][m];
      }
    }

    validate();
  } catch (const std::exception& e) {
    LDCSD_LOG_ERROR(std::string("failed to read input deck '") + path_to_yaml.string() +
                    "': " + e.what());
    return 1;
  }

  LDCSD_LOG_INFO("read input deck '" + path_to_yaml.string() + "': " + std::to_string(mesh.n_x) +
                 " cells, " + std::to_string(energy.G) + " groups, " + std::to_string(angle.M) +
                 " ordinates");
  return 0;
}
