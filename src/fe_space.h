#ifndef FE_SPACE_H
#define FE_SPACE_H

#include <Eigen/Core>

// Which local mass matrix a linear-discontinuous FESpace uses: the standard
// (consistent) Galerkin mass matrix, or its row-summed (lumped) variant.
enum class MassMatrixKind { Consistent, Lumped };

/// @class FESpace
/// @brief container containing discretization degree in energy and in space
/// @details Everything that depends on discretization degree should use this instead of hard-coded
/// quantities
class FESpace {
public:
  // Constructs an FESpace for the given spatial and energy polynomial
  // degrees. Both must be non-negative. mass_matrix_kind selects the local
  // 2x2 mass matrix M; L and Lb (the local stiffness and boundary-flux
  // matrices) are fixed regardless of it.
  FESpace(int spatial_degree, int energy_degree,
          MassMatrixKind mass_matrix_kind = MassMatrixKind::Consistent);

  int spatial_degree;
  int energy_degree;
  MassMatrixKind mass_matrix_kind;

  // Local 1D linear-basis matrices on the reference element [-1, 1], shared
  // by both the x and E directions (both use the same linear-discontinuous
  // basis). M is the mass matrix (per mass_matrix_kind); L is the local
  // stiffness (gradient) matrix; Lb is the boundary-flux matrix. Not yet
  // used by any assembly code.
  const Eigen::Matrix2d M;
  const Eigen::Matrix2d L;
  const Eigen::Matrix2d Lb;
};

/// @brief Which of a pair of axes is more significant (outer) in a flat
/// index: the x/spatial axis or the E/energy axis. Used for CornerValues'
/// (x-side, E-side) corner order. Field's own (group, cell) storage order is
/// NOT configurable via this -- it's fixed group-major, since the solver
/// always solves one full energy group at a time (energy coupling is
/// one-directional, high energy to low, never solved simultaneously across
/// groups) and depends on that layout.
enum class AxisOrder { XMajor, EMajor };

/// @brief One of the four corners of a linear-discontinuous (cell, group)
/// element: low/high in x crossed with low/high in energy.
enum class Corner { LeftDown, LeftUp, RightDown, RightUp };

/// @brief Maps a named Corner to its physical slot (0-3) within a 4-value
/// block, given which axis is more significant. This is the single place
/// that mapping is computed -- both CornerValues' accessors and
/// Field::index() go through it, so they can never disagree.
inline int cornerSlot(Corner corner, AxisOrder order) {
  const int x_bit = (corner == Corner::RightDown || corner == Corner::RightUp) ? 1 : 0;
  const int e_bit = (corner == Corner::LeftUp || corner == Corner::RightUp) ? 1 : 0;
  return order == AxisOrder::XMajor ? x_bit * 2 + e_bit : e_bit * 2 + x_bit;
}

/// @class CornerValues
/// @brief References 4 values in a (Field) Array using Eigen::Map
/// @details Non-constant references to four contiguous values in a Field; can be used to directly
/// modify those values with minimal overhead. Which physical slot each named
/// corner maps to depends on the AxisOrder it's constructed with (see
/// cornerSlot()), so it stays consistent with whatever order the owning
/// Field was configured for.
class CornerValues {
public:
  CornerValues(Eigen::Map<Eigen::Vector4d> values, AxisOrder corner_order)
      : values_(values), corner_order_(corner_order) {}

  double& leftDown() { return values_(cornerSlot(Corner::LeftDown, corner_order_)); }
  double& leftUp() { return values_(cornerSlot(Corner::LeftUp, corner_order_)); }
  double& rightDown() { return values_(cornerSlot(Corner::RightDown, corner_order_)); }
  double& rightUp() { return values_(cornerSlot(Corner::RightUp, corner_order_)); }

  double leftDown() const { return values_(cornerSlot(Corner::LeftDown, corner_order_)); }
  double leftUp() const { return values_(cornerSlot(Corner::LeftUp, corner_order_)); }
  double rightDown() const { return values_(cornerSlot(Corner::RightDown, corner_order_)); }
  double rightUp() const { return values_(cornerSlot(Corner::RightUp, corner_order_)); }

private:
  Eigen::Map<Eigen::Vector4d> values_;
  AxisOrder corner_order_;
};

/// @class Field
/// @brief Scalar field: holds 4 corner values in a n_x by G rectangular grid
/// @details Storage is always group-major: (group * n_x + cell) * 4 + corner -- fixed, not
/// configurable, since the solver always solves one full energy group at a time. Corner order
/// (which of x-side/E-side is more significant) is still configurable via AxisOrder.
/// field[group][cell] returns a CornerValues view into the Eigen array; index() returns the same
/// location as a plain integer, e.g. for building a stiffness matrix.
class Field {
public:
  // Constructs a Field over n_x cells and G groups, zero-initialized. Both
  // must be positive. corner_order controls whether x-side or E-side is the
  // more significant bit within each 4-value corner block.
  Field(int n_x, int G, AxisOrder corner_order = AxisOrder::XMajor);

  int numCells() const { return n_x_; }
  int numGroups() const { return G_; }
  AxisOrder cornerOrder() const { return corner_order_; }

  // The flat-storage index of one corner of one (group, cell) element --
  // the same location field[group][cell]'s named accessor for `corner`
  // would read/write, but as a plain int rather than a reference. Intended
  // for building a global stiffness matrix (e.g. as Eigen::Triplet
  // row/column indices) without needing a live view into the field. Only
  // meaningful for a matrix spanning this whole Field's (group, cell) range
  // -- a single-group system (e.g. one transport sweep's bilinear form) has
  // its own, group-independent local indexing; don't use this for that.
  // group must be in [0, numGroups()), cell in [0, numCells()).
  int index(int group, int cell, Corner corner) const;

  // Non-owning accessor for one group's row of cells, returned by
  // Field::operator[] so that `field[group][cell]` works.
  class Row {
  public:
    CornerValues operator[](int cell);

  private:
    Row(double* data, int group, int n_x, AxisOrder corner_order)
        : data_(data), group_(group), n_x_(n_x), corner_order_(corner_order) {}
    double* data_;
    int group_;
    int n_x_;
    AxisOrder corner_order_;
    friend class Field;
  };

  Row operator[](int group);

  // Underlying flat storage, laid out for direct use in a linear solve.
  Eigen::VectorXd& values() { return values_; }
  const Eigen::VectorXd& values() const { return values_; }

private:
  // The flat-storage offset of the first (corner 0) value of (group, cell).
  // The single formula index() and operator[] both build on.
  int blockOffset(int group, int cell) const;

  const int n_x_;
  const int G_;
  const AxisOrder corner_order_;
  Eigen::VectorXd values_;
};

#endif
