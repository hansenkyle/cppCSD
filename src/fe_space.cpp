// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "fe_space.h"

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
