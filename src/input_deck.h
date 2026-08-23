#ifndef INPUT_DECK_H
#define INPUT_DECK_H

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

// Holds a run's input parameters, read in from a YAML input file.
class InputDeck {
public:
  // Reads and validates path_to_yaml, populating this deck's members.
  // Returns 0 if the file is valid; returns 1 early on the first error
  // found in the file (missing/malformed keys, undefined material
  // references, mismatched sizes, etc.).
  int read(const std::filesystem::path &path_to_yaml);

  std::vector<double> spatial_mesh;             // strictly ascending cell-boundary locations
  std::vector<std::string> region_materials;    // material name per cell, size == spatial_mesh.size() - 1
  std::vector<double> energy_mesh;              // strictly ascending group boundaries, MeV
  std::map<std::string, MaterialData> materials;  // keyed by material name
  ConvergenceCriteria convergence;
};

#endif
