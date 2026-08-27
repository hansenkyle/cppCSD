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

MaterialData parseMaterial(const YAML::Node& node, const std::string& name,
                           std::size_t num_groups) {
  MaterialData material;
  material.sigma_t = requireNode(node, "sigma_t").as<std::vector<double>>();
  material.sigma_s = requireNode(node, "sigma_s").as<std::vector<double>>();

  const YAML::Node stopping_power = requireNode(node, "stopping_power");
  material.stopping_power_average =
      requireNode(stopping_power, "group_average").as<std::vector<double>>();
  material.stopping_power_boundary =
      requireNode(stopping_power, "group_boundary").as<std::vector<double>>();

  requireSize(material.sigma_t, num_groups, "material '" + name + "' sigma_t");
  requireSize(material.sigma_s, num_groups, "material '" + name + "' sigma_s");
  requireSize(material.stopping_power_average, num_groups,
              "material '" + name + "' stopping_power.group_average");
  requireSize(material.stopping_power_boundary, num_groups + 1,
              "material '" + name + "' stopping_power.group_boundary");

  return material;
}

AngularQuadrature parseAngularQuadrature(const YAML::Node& node) {
  AngularQuadrature quadrature;
  quadrature.mu = requireNode(node, "mu").as<std::vector<double>>();
  quadrature.w = requireNode(node, "w").as<std::vector<double>>();

  requireSize(quadrature.w, quadrature.mu.size(), "angular_quadrature.w");

  for (std::size_t m = 1; m < quadrature.mu.size(); ++m) {
    if (quadrature.mu[m] <= quadrature.mu[m - 1]) {
      throw std::runtime_error("angular_quadrature.mu must be strictly ascending");
    }
  }

  double sum = 0.0;
  for (double weight : quadrature.w) {
    sum += weight;
  }
  if (sum <= 0.0) {
    throw std::runtime_error("angular_quadrature.w must sum to a positive value");
  }

  const double scale = 2.0 / sum;
  for (double& weight : quadrature.w) {
    weight *= scale;
  }
  LDCSD_LOG_INFO("normalized angular_quadrature.w: sum was " + std::to_string(sum) +
                 ", scaled by " + std::to_string(scale) + " to sum to 2");

  return quadrature;
}

} // namespace

int InputDeck::read(const std::filesystem::path& path_to_yaml) {
  try {
    const YAML::Node root = YAML::LoadFile(path_to_yaml.string());
    LDCSD_LOG_TRACE("parsed '" + path_to_yaml.string() + "' as YAML");

    std::vector<double> x_boundary = requireNode(root, "spatial_mesh").as<std::vector<double>>();
    std::vector<double> E_boundary = requireNode(root, "energy_mesh").as<std::vector<double>>();
    mesh.emplace(std::move(x_boundary), std::move(E_boundary));

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
               expandByRegion(region_materials, materials, num_groups, &MaterialData::sigma_s),
               expandByRegion(region_materials, materials, num_groups,
                              &MaterialData::stopping_power_average),
               expandByRegion(region_materials, materials, num_groups + 1,
                              &MaterialData::stopping_power_boundary),
               region_materials);

    const YAML::Node angular_quadrature_node = requireNode(root, "angular_quadrature");
    angular_quadrature.emplace(parseAngularQuadrature(angular_quadrature_node));

    const YAML::Node convergence_node = requireNode(root, "convergence");
    convergence.max_iters = requireNode(convergence_node, "max_iters").as<int>();
    convergence.epsilon = requireNode(convergence_node, "epsilon").as<double>();
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
