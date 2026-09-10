// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "fe_space.h"

#include <stdexcept>

namespace {

Mat2 massMatrix(MassMatrixKind kind) {
  switch (kind) {
  case MassMatrixKind::Consistent:
    return Mat2{{2.0 / 6.0, 1.0 / 6.0}, {1.0 / 6.0, 2.0 / 6.0}};
  case MassMatrixKind::Lumped:
    return Mat2{{0.5, 0.0}, {0.0, 0.5}};
  }
  return Mat2{};
}

} // namespace

FESpace::FESpace(MassMatrixKind mass_matrix_kind)
    : mass_matrix_kind(mass_matrix_kind), M(massMatrix(mass_matrix_kind)),
      L(Mat2{{0.5, 0.5}, {-0.5, -0.5}}), Lb(Mat2{{-1.0, 0.0}, {0.0, 1.0}}) {}

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

Field::Field(int n_x, int G, AxisOrder corner_order)
    : n_x_(n_x), G_(G), corner_order_(corner_order),
      values_(Eigen::VectorXd::Zero(checkedFieldSize(n_x, G))) {}

int Field::blockOffset(int group, int cell) const { return (group * n_x_ + cell) * 4; }

// Bounds are checked with a throw rather than assert so out-of-range access
// fails the same way in Release as in Debug. If this ever shows up as a
// hot-path cost, it can be swapped for `assert` (compiled out in Release)
// at the expense of UB on misuse in Release builds.
Field::Row Field::operator[](int group) {
  if (group < 0 || group >= G_) {
    throw std::out_of_range("Field: group index out of range");
  }
  return Row(values_.data(), group, n_x_, corner_order_);
}

CornerValues Field::Row::operator[](int cell) {
  if (cell < 0 || cell >= n_x_) {
    throw std::out_of_range("Field: cell index out of range");
  }
  const int offset = (group_ * n_x_ + cell) * 4;
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
