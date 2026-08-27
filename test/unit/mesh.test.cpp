// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "mesh.h"

#include <stdexcept>

#include <doctest.h>

TEST_SUITE("Mesh") {
  TEST_CASE("derives sizes, widths, and centers from boundary arrays") {
    const Mesh mesh({0.0, 1.0, 3.0}, {5.0, 2.0, 0.0});

    CHECK(mesh.n_x == 2);
    CHECK(mesh.G == 2);

    CHECK(mesh.x_boundary == std::vector<double>{0.0, 1.0, 3.0});
    CHECK(mesh.dx == std::vector<double>{1.0, 2.0});
    CHECK(mesh.x_center == std::vector<double>{0.5, 2.0});

    CHECK(mesh.E_boundary == std::vector<double>{5.0, 2.0, 0.0});
    CHECK(mesh.dE == std::vector<double>{3.0, 2.0});
    CHECK(mesh.E_center == std::vector<double>{3.5, 1.0});
  }

  TEST_CASE("rejects a non-ascending x_boundary") {
    CHECK_THROWS_AS(Mesh({0.0, 2.0, 1.0}, {5.0, 0.0}), std::invalid_argument);
  }

  TEST_CASE("rejects a non-descending E_boundary") {
    CHECK_THROWS_AS(Mesh({0.0, 1.0}, {2.0, 5.0, 0.0}), std::invalid_argument);
  }

  TEST_CASE("rejects a negative E_boundary entry") {
    CHECK_THROWS_AS(Mesh({0.0, 1.0}, {5.0, -1.0, 0.0}), std::invalid_argument);
  }

  TEST_CASE("rejects an E_boundary that doesn't end at 0") {
    CHECK_THROWS_AS(Mesh({0.0, 1.0}, {5.0, 1.0}), std::invalid_argument);
  }
}
