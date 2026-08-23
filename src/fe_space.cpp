#include "fe_space.h"

#include <stdexcept>

FESpace::FESpace(int spatial_degree, int energy_degree)
    : spatial_degree(spatial_degree), energy_degree(energy_degree) {
  if (spatial_degree < 0) {
    throw std::invalid_argument("FESpace: spatial_degree must be non-negative");
  }
  if (energy_degree < 0) {
    throw std::invalid_argument("FESpace: energy_degree must be non-negative");
  }
}
