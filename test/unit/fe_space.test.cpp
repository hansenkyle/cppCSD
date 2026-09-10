// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "fe_space.h"

#include <initializer_list>

#include <doctest.h>

TEST_SUITE("FESpace") {
  TEST_CASE("defaults to a consistent (unlumped) mass matrix") {
    const FESpace space;

    CHECK(space.mass_matrix_kind == MassMatrixKind::Consistent);
    CHECK(space.M.left.left == doctest::Approx(1.0 / 3.0));
    CHECK(space.M.left.right == doctest::Approx(1.0 / 6.0));
    CHECK(space.M.right.left == doctest::Approx(1.0 / 6.0));
    CHECK(space.M.right.right == doctest::Approx(1.0 / 3.0));
  }

  TEST_CASE("uses the lumped mass matrix when requested") {
    const FESpace space(MassMatrixKind::Lumped);

    CHECK(space.mass_matrix_kind == MassMatrixKind::Lumped);
    CHECK(space.M.left.left == doctest::Approx(0.5));
    CHECK(space.M.left.right == doctest::Approx(0.0));
    CHECK(space.M.right.left == doctest::Approx(0.0));
    CHECK(space.M.right.right == doctest::Approx(0.5));
  }

  TEST_CASE("L and Lb are fixed regardless of mass_matrix_kind") {
    const FESpace consistent(MassMatrixKind::Consistent);
    const FESpace lumped(MassMatrixKind::Lumped);

    for (const FESpace& space : {consistent, lumped}) {
      CHECK(space.L.left.left == doctest::Approx(0.5));
      CHECK(space.L.left.right == doctest::Approx(0.5));
      CHECK(space.L.right.left == doctest::Approx(-0.5));
      CHECK(space.L.right.right == doctest::Approx(-0.5));

      CHECK(space.Lb.left.left == doctest::Approx(-1.0));
      CHECK(space.Lb.left.right == doctest::Approx(0.0));
      CHECK(space.Lb.right.left == doctest::Approx(0.0));
      CHECK(space.Lb.right.right == doctest::Approx(1.0));
    }
  }
}
