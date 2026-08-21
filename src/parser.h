#ifndef PARSER_H
#define PARSER_H

#include <filesystem>
#include <map>
#include <string>
#include <vector>

// Cross sections and stopping power for a single material, indexed by
// energy group. Group g spans energy_mesh[g] to energy_mesh[g + 1].
struct MaterialData {
  std::vector<double> sigma_t;                 // group total xs, size == num_groups
  std::vector<double> sigma_s;                  // group isotropic scattering xs, size == num_groups
  std::vector<double> stopping_power_average;   // group-average S, size == num_groups
  std::vector<double> stopping_power_boundary;  // S at group boundaries, size == num_groups + 1
};

// Placeholder convergence criteria; more sophisticated criteria to follow.
struct ConvergenceCriteria {
  int max_iters = 0;
  double epsilon = 0.0;
};

struct InputDeck {
  std::vector<double> spatial_mesh;             // strictly ascending cell-boundary locations
  std::vector<std::string> region_materials;    // material name per cell, size == spatial_mesh.size() - 1
  std::vector<double> energy_mesh;              // strictly ascending group boundaries, MeV
  std::map<std::string, MaterialData> materials;  // keyed by material name
  ConvergenceCriteria convergence;
};

class Parser {
public:
  static InputDeck read(const std::filesystem::path &path_to_yaml);
};

#endif
