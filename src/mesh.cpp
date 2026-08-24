#include "mesh.h"

#include <stdexcept>
#include <string>

namespace {

void requireStrictlyAscending(const std::vector<double>& values, const std::string& name) {
  if (values.size() < 2) {
    throw std::invalid_argument("Mesh: " + name + " must have at least two entries");
  }
  for (std::size_t i = 1; i < values.size(); ++i) {
    if (values[i] <= values[i - 1]) {
      throw std::invalid_argument("Mesh: " + name + " must be strictly ascending");
    }
  }
}

void requireStrictlyDescendingEndingAtZero(const std::vector<double>& values,
                                           const std::string& name) {
  if (values.size() < 2) {
    throw std::invalid_argument("Mesh: " + name + " must have at least two entries");
  }
  for (double value : values) {
    if (value < 0.0) {
      throw std::invalid_argument("Mesh: " + name + " must be non-negative");
    }
  }
  for (std::size_t i = 1; i < values.size(); ++i) {
    if (values[i] >= values[i - 1]) {
      throw std::invalid_argument("Mesh: " + name + " must be strictly descending");
    }
  }
  if (values.back() != 0.0) {
    throw std::invalid_argument("Mesh: " + name + " must end at 0");
  }
}

} // namespace

Mesh::Mesh(std::vector<double> x_boundary_in, std::vector<double> E_boundary_in)
    : x_boundary(std::move(x_boundary_in)), E_boundary(std::move(E_boundary_in)) {
  requireStrictlyAscending(x_boundary, "x_boundary");
  requireStrictlyDescendingEndingAtZero(E_boundary, "E_boundary");

  n_x = static_cast<int>(x_boundary.size()) - 1;
  G = static_cast<int>(E_boundary.size()) - 1;

  dx.resize(n_x);
  x_center.resize(n_x);
  for (int i = 0; i < n_x; ++i) {
    dx[i] = x_boundary[i + 1] - x_boundary[i];
    x_center[i] = 0.5 * (x_boundary[i] + x_boundary[i + 1]);
  }

  dE.resize(G);
  E_center.resize(G);
  for (int g = 0; g < G; ++g) {
    dE[g] = E_boundary[g] - E_boundary[g + 1];
    E_center[g] = 0.5 * (E_boundary[g] + E_boundary[g + 1]);
  }
}
