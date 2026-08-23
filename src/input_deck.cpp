#include "input_deck.h"

#include <stdexcept>

#include <yaml-cpp/yaml.h>

#include "logger.h"

namespace {

void requireStrictlyAscending(const std::vector<double>& values, const std::string& name) {
  if (values.size() < 2) {
    throw std::runtime_error(name + " must have at least two entries");
  }
  for (std::size_t i = 1; i < values.size(); ++i) {
    if (values[i] <= values[i - 1]) {
      throw std::runtime_error(name + " must be strictly ascending");
    }
  }
}

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

} // namespace

int InputDeck::read(const std::filesystem::path& path_to_yaml) {
  try {
    const YAML::Node root = YAML::LoadFile(path_to_yaml.string());
    LDCSD_LOG_TRACE("parsed '" + path_to_yaml.string() + "' as YAML");

    spatial_mesh = requireNode(root, "spatial_mesh").as<std::vector<double>>();
    requireStrictlyAscending(spatial_mesh, "spatial_mesh");
    const std::size_t num_cells = spatial_mesh.size() - 1;

    const YAML::Node regions = requireNode(root, "regions");
    region_materials = requireNode(regions, "materials").as<std::vector<std::string>>();
    if (region_materials.size() != num_cells) {
      throw std::runtime_error("regions.materials has size " +
                               std::to_string(region_materials.size()) + ", expected " +
                               std::to_string(num_cells) + " (one per spatial cell)");
    }

    energy_mesh = requireNode(root, "energy_mesh").as<std::vector<double>>();
    requireStrictlyAscending(energy_mesh, "energy_mesh");
    const std::size_t num_groups = energy_mesh.size() - 1;

    const YAML::Node materials_node = requireNode(root, "materials");
    for (const auto& entry : materials_node) {
      const std::string name = entry.first.as<std::string>();
      materials[name] = parseMaterial(entry.second, name, num_groups);
      LDCSD_LOG_DEBUG("parsed material '" + name + "'");
    }

    for (const std::string& name : region_materials) {
      if (!materials.contains(name)) {
        throw std::runtime_error("region references undefined material '" + name + "'");
      }
    }

    const YAML::Node convergence_node = requireNode(root, "convergence");
    convergence.max_iters = requireNode(convergence_node, "max_iters").as<int>();
    convergence.epsilon = requireNode(convergence_node, "epsilon").as<double>();
  } catch (const std::exception& e) {
    LDCSD_LOG_ERROR(std::string("failed to read input deck '") + path_to_yaml.string() +
                    "': " + e.what());
    return 1;
  }

  LDCSD_LOG_INFO("read input deck '" + path_to_yaml.string() +
                 "': " + std::to_string(spatial_mesh.size() - 1) + " cells, " +
                 std::to_string(energy_mesh.size() - 1) + " groups, " +
                 std::to_string(materials.size()) + " materials");
  return 0;
}
