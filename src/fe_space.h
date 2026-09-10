// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#ifndef FE_SPACE_H
#define FE_SPACE_H

// Which local mass matrix a linear-discontinuous FESpace uses: the standard
// (consistent) Galerkin mass matrix, or its row-summed (lumped) variant.
enum class MassMatrixKind { Consistent, Lumped };

// A 2x2 local matrix on the reference element's two nodes, accessed by node
// name (left/right) rather than numeric index -- so callers never need to
// know or agree on an axis ordering to read an entry.
struct Row2 {
  double left;
  double right;
};
struct Mat2 {
  Row2 left;
  Row2 right;
};

/// @class FESpace
/// @brief The linear-discontinuous reference-element matrices, shared by the x and E directions
class FESpace {
public:
  // Constructs an FESpace. mass_matrix_kind selects the local mass matrix
  // M; L and Lb (the local stiffness and boundary-flux matrices) are fixed
  // regardless of it.
  explicit FESpace(MassMatrixKind mass_matrix_kind = MassMatrixKind::Consistent);

  MassMatrixKind mass_matrix_kind;

  // Local 1D linear-basis matrices on the reference element [-1, 1], shared
  // by both the x and E directions (both use the same linear-discontinuous
  // basis). M is the mass matrix (per mass_matrix_kind); L is the local
  // stiffness (gradient) matrix; Lb is the boundary-flux matrix.
  const Mat2 M;
  const Mat2 L;
  const Mat2 Lb;
};

#endif
