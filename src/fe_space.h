#ifndef FE_SPACE_H
#define FE_SPACE_H

#include <Eigen/Core>

/// @class FESpace
/// @brief container containing discretization degree in energy and in space
/// @details Everything that depends on discretization degree should use this instead of hard-coded
/// quantities
class FESpace {
public:
  // Constructs an FESpace for the given spatial and energy polynomial
  // degrees. Both must be non-negative.
  FESpace(int spatial_degree, int energy_degree);

  int spatial_degree;
  int energy_degree;
};

/// @class CornerValues
/// @brief References 4 values in a (Field) Array using Eigen::Map
/// @details Non-constant references to four contiguous values in a Field; can be used to directly
/// modify those values with minimal overhead
class CornerValues {
public:
  explicit CornerValues(Eigen::Map<Eigen::Vector4d> values) : values_(values) {}

  double& leftDown() { return values_(0); }
  double& leftUp() { return values_(1); }
  double& rightDown() { return values_(2); }
  double& rightUp() { return values_(3); }

  double leftDown() const { return values_(0); }
  double leftUp() const { return values_(1); }
  double rightDown() const { return values_(2); }
  double rightUp() const { return values_(3); }

private:
  Eigen::Map<Eigen::Vector4d> values_;
};

/// @class Field
/// @brief Scalar field: holds 4 corner values in a n_x by G rectangular grid
/// @details Indexing built-in as (group * n_x + cell) * 4 + corner -- block  all values of the same
/// energy group together; spatially adjacent cells are near each other in the object.
/// field[group][cell] returns a CornerValues view into the EigenXd array.
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
    Row(double* data, int group, int n_x) : data_(data), group_(group), n_x_(n_x) {}
    double* data_;
    int group_;
    int n_x_;
    friend class Field;
  };

  Row operator[](int group);

  // Underlying flat storage, laid out for direct use in a linear solve.
  Eigen::VectorXd& values() { return values_; }
  const Eigen::VectorXd& values() const { return values_; }

private:
  const int n_x_;
  const int G_;
  Eigen::VectorXd values_;
};

#endif
