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

namespace {

int checkedFieldSize(int n_x, int G) {
  if (n_x <= 0) {
    throw std::invalid_argument("Field: n_x must be positive");
  }
  if (G <= 0) {
    throw std::invalid_argument("Field: G must be positive");
  }
  return 4 * n_x * G;
}

} // namespace

Field::Field(int n_x, int G, AxisOrder axis_order, AxisOrder corner_order)
    : n_x_(n_x), G_(G), axis_order_(axis_order), corner_order_(corner_order),
      values_(Eigen::VectorXd::Zero(checkedFieldSize(n_x, G))) {}

int Field::blockOffset(int group, int cell) const {
  return axis_order_ == AxisOrder::EMajor ? (group * n_x_ + cell) * 4 : (cell * G_ + group) * 4;
}

// Bounds are checked with a throw rather than assert so out-of-range access
// fails the same way in Release as in Debug. If this ever shows up as a
// hot-path cost, it can be swapped for `assert` (compiled out in Release)
// at the expense of UB on misuse in Release builds.
Field::Row Field::operator[](int group) {
  if (group < 0 || group >= G_) {
    throw std::out_of_range("Field: group index out of range");
  }
  const int base_offset = blockOffset(group, 0);
  const int stride_cell = axis_order_ == AxisOrder::EMajor ? 4 : G_ * 4;
  return Row(values_.data(), base_offset, stride_cell, n_x_, corner_order_);
}

CornerValues Field::Row::operator[](int cell) {
  if (cell < 0 || cell >= n_x_) {
    throw std::out_of_range("Field: cell index out of range");
  }
  const int offset = base_offset_ + cell * stride_cell_;
  return CornerValues(Eigen::Map<Eigen::Vector4d>(data_ + offset), corner_order_);
}

int Field::index(int group, int cell, Corner corner) const {
  if (group < 0 || group >= G_) {
    throw std::out_of_range("Field: group index out of range");
  }
  if (cell < 0 || cell >= n_x_) {
    throw std::out_of_range("Field: cell index out of range");
  }
  return blockOffset(group, cell) + cornerSlot(corner, corner_order_);
}
