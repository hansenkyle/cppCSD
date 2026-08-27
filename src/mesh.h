#ifndef MESH_H
#define MESH_H

#include <vector>

/// @class Mesh
/// @brief Stores space-energy grid parameters, independent of the finite element discretization
/// scheme
/// @details Construct using x- and E- boundary values, all other fields are derived
class Mesh {
public:
  // Constructs a Mesh from cell and group boundary locations. x_boundary
  // must have at least two strictly ascending entries. E_boundary must have
  // at least two entries, all non-negative, strictly descending, ending at
  // exactly 0.
  Mesh(std::vector<double> x_boundary, std::vector<double> E_boundary);

  int n_x; // number of spatial cells, == x_boundary.size() - 1
  int G;   // number of energy groups, == E_boundary.size() - 1

  std::vector<double> dx; // cell widths, size n_x
  std::vector<double> dE; // group widths, size G

  std::vector<double> x_boundary; // strictly ascending, size n_x + 1
  std::vector<double> E_boundary; // strictly descending, ends at 0, size G + 1

  std::vector<double> x_center; // cell midpoints, size n_x
  std::vector<double> E_center; // group midpoints, size G
};

#endif
