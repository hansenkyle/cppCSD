// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "fe_space.h"

#include <stdexcept>

namespace {

Eigen::Matrix2d massMatrix(MassMatrixKind kind) {
  Eigen::Matrix2d m;
  switch (kind) {
  case MassMatrixKind::Consistent:
    m << 2.0, 1.0, 1.0, 2.0;
    break;
  case MassMatrixKind::Lumped:
    m << 3.0, 0.0, 0.0, 3.0;
    break;
  }
  return m / 6.0;
}

} // namespace

FESpace::FESpace(int spatial_degree, int energy_degree, MassMatrixKind mass_matrix_kind)
    : spatial_degree(spatial_degree), energy_degree(energy_degree),
      mass_matrix_kind(mass_matrix_kind), M(massMatrix(mass_matrix_kind)),
      L((Eigen::Matrix2d() << 0.5, 0.5, -0.5, -0.5).finished()),
      Lb((Eigen::Matrix2d() << -1.0, 0.0, 0.0, 1.0).finished()) {
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

Field::ConstRow Field::operator[](int group) const {
  if (group < 0 || group >= G_) {
    throw std::out_of_range("Field: group index out of range");
  }
  return ConstRow(values_.data(), group, n_x_, corner_order_);
}

ConstCornerValues Field::ConstRow::operator[](int cell) const {
  if (cell < 0 || cell >= n_x_) {
    throw std::out_of_range("Field: cell index out of range");
  }
  const int offset = (group_ * n_x_ + cell) * 4;
  return ConstCornerValues(Eigen::Map<const Eigen::Vector4d>(data_ + offset), corner_order_);
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
