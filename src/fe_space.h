#ifndef FE_SPACE_H
#define FE_SPACE_H

#include <Eigen/Core>

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

// Non-owning view of the four corner values of a single linear-discontinuous
// (cell, group) element -- low/high in x crossed with low/high in energy.
// Backed by a Map into a contiguous segment of a Field's flat storage, so
// it's cheap to construct (pointer arithmetic only) and writes through it
// land directly in that storage -- no separate copy to keep in sync.
class CornerValues {
public:
  explicit CornerValues(Eigen::Map<Eigen::Vector4d> values) : values_(values) {}

  double &leftDown() { return values_(0); }
  double &leftUp() { return values_(1); }
  double &rightDown() { return values_(2); }
  double &rightUp() { return values_(3); }

  double leftDown() const { return values_(0); }
  double leftUp() const { return values_(1); }
  double rightDown() const { return values_(2); }
  double rightUp() const { return values_(3); }

private:
  Eigen::Map<Eigen::Vector4d> values_;
};

// A scalar field over an n_x-cell, G-group grid, with 4 linear-discontinuous
// corner values per (group, cell). Storage is a single flat Eigen::VectorXd
// laid out as (group * n_x + cell) * 4 + corner, matching the layout a 4GIx4GI
// stiffness matrix expects -- so the field can be handed directly to a linear
// solve with no packing/unpacking step. `field[group][cell]` returns a
// CornerValues view onto that segment for named, semantic access.
class Field {
public:
  // Constructs a Field over n_x cells and G groups, zero-initialized. Both
  // must be positive.
  Field(int n_x, int G);

  int numCells() const { return n_x_; }
  int numGroups() const { return G_; }

  // Non-owning accessor for one group's row of cells, returned by
  // Field::operator[] so that `field[group][cell]` works.
  class Row {
  public:
    CornerValues operator[](int cell);

  private:
    Row(double *data, int group, int n_x) : data_(data), group_(group), n_x_(n_x) {}
    double *data_;
    int group_;
    int n_x_;
    friend class Field;
  };

  Row operator[](int group);

  // Underlying flat storage, laid out for direct use in a linear solve.
  Eigen::VectorXd &values() { return values_; }
  const Eigen::VectorXd &values() const { return values_; }

private:
  const int n_x_;
  const int G_;
  Eigen::VectorXd values_;
};

#endif
