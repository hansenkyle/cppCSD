#ifndef FE_SPACE_H
#define FE_SPACE_H

// Minimal description of the discretization degree used for a solve: the
// polynomial degree in space and, independently, in energy. Everything that
// depends on these degrees -- dofs per cell, basis evaluation, dof maps --
// is expected to be derived from this, not duplicated elsewhere.
class FESpace {
public:
  // Constructs an FESpace for the given spatial and energy polynomial
  // degrees. Both must be non-negative.
  FESpace(int spatial_degree, int energy_degree);

  int spatial_degree;
  int energy_degree;
};

#endif
