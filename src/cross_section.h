#ifndef CROSS_SECTION_H
#define CROSS_SECTION_H

#include <string>
#include <vector>

#include "mesh.h"

// Per-cell, per-group cross sections and stopping power over a Mesh's
// spatial/energy grid. total/scattering/stop_power/stop_power_boundary are
// indexed [group][cell], so indexing by group alone yields the
// space-dependent vector for that group; material is indexed [cell] only.
class CrossSection {
public:
  // Constructs a CrossSection over `mesh` from already-expanded per-cell,
  // per-group data. `total`, `scattering`, and `stop_power` must each have
  // mesh.G rows of mesh.n_x entries; `stop_power_boundary` must have
  // mesh.G + 1 rows of mesh.n_x entries (S at group boundaries); `material`
  // must have mesh.n_x entries. Every numeric value must be non-negative.
  // `mesh` must outlive this CrossSection.
  CrossSection(const Mesh &mesh, std::vector<std::vector<double>> total,
               std::vector<std::vector<double>> scattering,
               std::vector<std::vector<double>> stop_power,
               std::vector<std::vector<double>> stop_power_boundary,
               std::vector<std::string> material);

  const Mesh &mesh;

  std::vector<std::vector<double>> total;               // [group][cell], size G x n_x
  std::vector<std::vector<double>> scattering;           // [group][cell], size G x n_x
  std::vector<std::vector<double>> stop_power;            // group-average S, [group][cell], size G x n_x
  std::vector<std::vector<double>> stop_power_boundary;   // S at group boundaries, [group boundary][cell], size (G+1) x n_x
  std::vector<std::string> material;                      // material name, [cell], size n_x
};

#endif
