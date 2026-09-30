// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "smm.h"

#include <filesystem>

#include <doctest.h>
#include <highfive/eigen.hpp>
#include <highfive/highfive.hpp>

TEST_CASE("SecondMoment::writeH5 round-trips the solution with the documented shapes") {
  InputDeck deck;
  REQUIRE(deck.read(std::filesystem::path(TEST_DATA_DIR) / "sample_input.yaml") == 0);
  SecondMoment sm(deck);
  sm.solve(1e-10);

  const std::filesystem::path path = std::filesystem::temp_directory_path() / "ldcsd_h5_test.h5";
  sm.writeH5(path, "2026-01-01 00:00:00");

  const HighFive::File file(path.string(), HighFive::File::ReadOnly);
  const std::size_t nx = deck.mesh.n_x, G = deck.energy.G, M = deck.angle.M;

  CHECK(file.getAttribute("execution_datetime").read<std::string>() == "2026-01-01 00:00:00");
  CHECK(file.getAttribute("n_groups").read<int>() == static_cast<int>(G));

  const auto scalar = file.getDataSet("/solution/scalar_flux").read<Eigen::MatrixXd>();
  CHECK(scalar == sm.solution.scalar_flux);

  const auto angular = file.getDataSet("/solution/angular_flux");
  CHECK(angular.getDimensions() == std::vector<std::size_t>{G, 4 * nx, M});
  const auto psi = angular.read<std::vector<std::vector<std::vector<double>>>>();
  CHECK(psi[G - 1][4 * nx - 1][M - 1] == sm.solution.angular_flux[G - 1](4 * nx - 1, M - 1));

  CHECK(file.getDataSet("/residuals/second_moment").getDimensions() ==
        std::vector<std::size_t>{G, 8 * nx});
  CHECK(file.getDataSet("/solution/closures/K+").getDimensions() ==
        std::vector<std::size_t>{G, 4 * nx});
  CHECK(file.getDataSet("/input/mesh/x_boundary").getDimensions() ==
        std::vector<std::size_t>{nx + 1});
  CHECK(!file.getDataSet("/convergence/g001/delta_l2").read<std::vector<double>>().empty());
}
