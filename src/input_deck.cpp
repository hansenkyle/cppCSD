// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "input_deck.h"

#include <cmath>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <yaml-cpp/yaml.h>

#include "logger.h"

namespace {
constexpr double kAngleWeightRelTol = 1e-4;
void requireNonNegative(const Eigen::Ref<const Eigen::MatrixXd>& values, const std::string& name) {
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

// One entry in a group-to-group scattering matrix: group `from` scatters
// into group `to` with the given macroscopic cross section. Parsing scratch
// only -- collected into a dense matrix (via buildScatterMatrix) before
// being stored anywhere.
struct ScatterEntry {
  int from;
  int to;
  double value;
};

// Builds a G x G matrix from a flat entry list, entry (from, to) -> value;
// every pair the list doesn't mention stays zero.
Eigen::MatrixXd buildScatterMatrix(const std::vector<ScatterEntry>& entries, int G) {
  Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(G, G);
  for (const ScatterEntry& entry : entries) {
    matrix(entry.from, entry.to) = entry.value;
  }
  return matrix;
}

// Dense form: scattering is a list of G rows, each a list of G columns --
// row = from group, column = to group. Zero entries are dropped so the
// result is the same sparse entry list a from/to reader would have produced.
std::vector<ScatterEntry> parseDenseScattering(const YAML::Node& node,
                                               const std::string& material_name, int G) {
  std::vector<ScatterEntry> entries;
  requireSize(static_cast<Eigen::Index>(node.size()), G,
              "material '" + material_name + "' scattering row count");
  for (int from = 0; from < G; ++from) {
    const YAML::Node row = node[from];
    requireSize(static_cast<Eigen::Index>(row.size()), G,
                "material '" + material_name + "' scattering row " + std::to_string(from + 1));
    for (int to = 0; to < G; ++to) {
      const double value = row[to].as<double>();
      if (value != 0.0) {
        entries.push_back(ScatterEntry{from, to, value});
      }
    }
  }
  return entries;
}

// Group indices in the scattering: list are 1-indexed in YAML (matching how
// a person would naturally refer to "group 1"); stored 0-indexed internally
// like everything else. Checks from/to fall within [1, G] and that no
// (from, to) pair repeats -- a duplicate is almost certainly a typo, so it's
// rejected rather than summed.
//
// Accepts two forms: sparse, a list of {from, to, value} maps for a matrix
// that's mostly zero; or dense (see parseDenseScattering), a list of G rows
// of G values each, for a matrix with few zeros. The two are told apart by
// the type of the list's first element (map vs. sequence).
std::vector<ScatterEntry> parseScattering(const YAML::Node& node, const std::string& material_name,
                                          int G) {
  if (node.size() > 0 && node[0].IsSequence()) {
    return parseDenseScattering(node, material_name, G);
  }

  std::vector<ScatterEntry> entries;
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
    for (const ScatterEntry& existing : entries) {
      if (existing.from == from && existing.to == to) {
        throw std::runtime_error(
            "material '" + material_name + "' has a duplicate scattering entry (from=" +
            std::to_string(from_1indexed) + ", to=" + std::to_string(to_1indexed) + ")");
      }
    }
    entries.push_back(ScatterEntry{from, to, value});
  }
  return entries;
}

Material parseMaterial(const YAML::Node& node, const std::string& name, int G) {
  Material material;
  material.name = name;
  material.total = toVector(requireNode(node, "sigma_t").as<std::vector<double>>());
  material.scatter =
      buildScatterMatrix(parseScattering(requireNode(node, "scattering"), name, G), G);

  const YAML::Node stopping_power = requireNode(node, "stopping_power");
  material.S = toVector(requireNode(stopping_power, "group_average").as<std::vector<double>>());
  material.S_b = toVector(requireNode(stopping_power, "group_boundary").as<std::vector<double>>());

  requireSize(material.total.size(), G, "material '" + name + "' sigma_t");
  requireSize(material.S.size(), G, "material '" + name + "' stopping_power.group_average");
  requireSize(material.S_b.size(), G + 1, "material '" + name + "' stopping_power.group_boundary");
  return material;
}

// Parses one group's `source` entry: a list of M ordinate entries, each a
// list of n_x per-cell {up_left, up_right, down_left, down_right} maps.
// Returned as a (4 * n_x) x M matrix, rows 4*c..4*c+3 = that cell's four
// corner values, matching Kernel::solveDirect's per-cell q_up/q_down layout.
Eigen::MatrixXd parseSourceGroup(const YAML::Node& group_node, const std::string& name, int M,
                                 int n_x) {
  requireSize(static_cast<Eigen::Index>(group_node.size()), M, name + " ordinate count");
  Eigen::MatrixXd values(4 * n_x, M);
  for (int m = 0; m < M; ++m) {
    const YAML::Node ordinate_node = group_node[m];
    const std::string ordinate_name = name + " ordinate " + std::to_string(m + 1);
    requireSize(static_cast<Eigen::Index>(ordinate_node.size()), n_x,
                ordinate_name + " cell count");
    for (int c = 0; c < n_x; ++c) {
      const YAML::Node cell_node = ordinate_node[c];
      values(4 * c, m) = requireNode(cell_node, "up_left").as<double>();
      values(4 * c + 1, m) = requireNode(cell_node, "up_right").as<double>();
      values(4 * c + 2, m) = requireNode(cell_node, "down_left").as<double>();
      values(4 * c + 3, m) = requireNode(cell_node, "down_right").as<double>();
    }
  }
  return values;
}
} // namespace

void Material::validate() {
  // check that name exists
  if (name == "") {
    throw std::runtime_error("Material name must be provided");
  }
  // check that all vectors are same size
  int t_size = total.size();
  int s_size = S.size();
  int sb_size = S_b.size();
  int scatter_rows = scatter.rows();
  int scatter_cols = scatter.cols();

  if (!(t_size == s_size and s_size == sb_size - 1 and s_size == scatter_rows and
        scatter_rows == scatter_cols)) {
    throw std::runtime_error("Material data sizes do not match. (tot, S, S_b, scat) = (" +
                             std::to_string(t_size) + ", " + std::to_string(s_size) + ", " +
                             std::to_string(sb_size) + ", " + std::to_string(scatter_rows) + "x" +
                             std::to_string(scatter_cols) + ").");
  }
}

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
  x_center = x_boundary.head(n_x) + dx;
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

void InputDeck::Xs::set_materials(std::vector<Material> materials, std::vector<int> indices) {
  for (Material& material : materials) {
    material.validate();
  }
  for (int index : indices) {
    if (index < 0 || index >= static_cast<int>(materials.size())) {
      throw std::runtime_error("material index " + std::to_string(index) +
                               " is out of range for a material list of size " +
                               std::to_string(materials.size()));
    }
  }

  material_list = std::move(materials);
  material_indices = std::move(indices);
}

double InputDeck::Xs::total(int g, int i) { return material_list[material_indices[i]].total(g); }

double InputDeck::Xs::S(int g, int i) { return material_list[material_indices[i]].S(g); }

double InputDeck::Xs::S_b(int g, int i) { return material_list[material_indices[i]].S_b(g); }

double InputDeck::Xs::S_up(int g, int i) { return material_list[material_indices[i]].S_b(g); }

double InputDeck::Xs::S_down(int g, int i) { return material_list[material_indices[i]].S_b(g + 1); }

double InputDeck::Xs::scatter(int from, int to, int i) {
  return material_list[material_indices[i]].scatter(from, to);
}

Eigen::MatrixXd& InputDeck::Xs::scatter(int i) {
  return material_list[material_indices[i]].scatter;
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

  // create positivity mask
  auto nonneg_indicator = [](const Eigen::VectorXd& v) {
    return (v.array() >= 0.0).select(v, Eigen::VectorXd::Zero(v.size()));
  };
  w_positive = nonneg_indicator(mu);

  auto negative_indicator = [](const Eigen::VectorXd& v) {
    return (v.array() < 0.0).select(v, Eigen::VectorXd::Zero(v.size()));
  };
  w_negative = negative_indicator(mu);
}

void InputDeck::Xs::validate() const {
  for (const Material& material : material_list) {
    const std::string where = " for material '" + material.name + "'";
    requireNonNegative(material.total, "xs.total" + where);
    requireNonNegative(material.S, "xs.S" + where);
    requireNonNegative(material.S_b, "xs.S_b" + where);
    requireNonNegative(material.scatter, "xs.scatter" + where);
  }
}

void InputDeck::Xs::validateShape(int G, int n_x) const {
  if (static_cast<int>(material_indices.size()) != n_x) {
    throw std::runtime_error("xs covers " + std::to_string(material_indices.size()) +
                             " cells, expected mesh.n_x = " + std::to_string(n_x));
  }
  for (const Material& material : material_list) {
    if (material.total.size() != G) {
      throw std::runtime_error("material '" + material.name + "' has " +
                               std::to_string(material.total.size()) +
                               " groups, expected energy.G = " + std::to_string(G));
    }
  }
}

void InputDeck::Source::validate() const {
  for (const Eigen::MatrixXd& ordinate : values) {
    requireNonNegative(ordinate, "source");
  }
}

void InputDeck::validate() {
  mesh.validate();
  energy.validate();
  angle.validate();
  xs.validate();
  source.validate();

  xs.validateShape(energy.G, mesh.n_x);

  if (bc.values.rows() != 2 * energy.G) {
    throw std::runtime_error("bc.values must have 2 * energy.G = " + std::to_string(2 * energy.G) +
                             " rows");
  }
  if (bc.values.cols() != angle.M) {
    throw std::runtime_error("bc.values must have angle.M = " + std::to_string(angle.M) +
                             " columns");
  }

  if (static_cast<int>(source.values.size()) != energy.G) {
    throw std::runtime_error("source must have energy.G = " + std::to_string(energy.G) + " groups");
  }

  // Angular moments of the source, integrated over the quadrature. They pair
  // the source with the angle struct, so they're derived here rather than in
  // Source::validate(), which sees only its own values.
  const Eigen::VectorXd w_mu = angle.w.cwiseProduct(angle.mu);
  source.q0.resize(energy.G);
  source.q1.resize(energy.G);
  for (int g = 0; g < energy.G; ++g) {
    if (source.values[g].rows() != 4 * mesh.n_x || source.values[g].cols() != angle.M) {
      throw std::runtime_error("source group " + std::to_string(g) +
                               " must be shaped 4 * mesh.n_x = " + std::to_string(4 * mesh.n_x) +
                               " rows x angle.M = " + std::to_string(angle.M) + " columns");
    }
    source.q0[g] = source.values[g] * angle.w;
    source.q1[g] = source.values[g] * w_mu;
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

    // Materials keep the order they appear in the file; index_of turns the
    // names in `regions` into positions in that list.
    std::vector<Material> material_list;
    std::map<std::string, int> index_of;
    const YAML::Node materials_node = requireNode(root, "materials");
    for (const auto& entry : materials_node) {
      const std::string name = entry.first.as<std::string>();
      index_of[name] = static_cast<int>(material_list.size());
      material_list.push_back(parseMaterial(entry.second, name, energy.G));
      LDCSD_LOG_DEBUG("parsed material '" + name + "'");
    }

    std::vector<int> material_indices;
    material_indices.reserve(region_materials.size());
    for (const std::string& name : region_materials) {
      const auto found = index_of.find(name);
      if (found == index_of.end()) {
        throw std::runtime_error("region references undefined material '" + name + "'");
      }
      material_indices.push_back(found->second);
    }

    xs.set_materials(std::move(material_list), std::move(material_indices));

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

    const YAML::Node source_node = requireNode(root, "source");
    requireSize(static_cast<Eigen::Index>(source_node.size()), energy.G, "source group count");
    source.values.resize(energy.G);
    for (int g = 0; g < energy.G; ++g) {
      source.values[g] = parseSourceGroup(source_node[g], "source group " + std::to_string(g + 1),
                                          angle.M, mesh.n_x);
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
