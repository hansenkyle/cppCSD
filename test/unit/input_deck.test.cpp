// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "input_deck.h"

#include <fstream>

#include <doctest.h>

namespace {

std::filesystem::path sampleInputPath() {
  return std::filesystem::path(TEST_DATA_DIR) / "sample_input.yaml";
}

// Writes `contents` to a temp file and returns its path; used for cases that
// need to deviate from the checked-in sample input.
std::filesystem::path writeTempYaml(const std::string& contents) {
  const std::filesystem::path path =
      std::filesystem::temp_directory_path() / "ldcsd_input_deck_test_input.yaml";
  std::ofstream out(path);
  out << contents;
  return path;
}

} // namespace

TEST_SUITE("InputDeck") {
  TEST_CASE("reads a well-formed input deck") {
    InputDeck deck;
    CHECK(deck.read(sampleInputPath()) == 0);

    REQUIRE(deck.mesh.has_value());
    CHECK(deck.mesh->x_boundary == std::vector<double>{0.0, 1.0, 2.0, 3.0});
    CHECK(deck.region_materials == std::vector<std::string>{"water", "water", "lead"});
    CHECK(deck.mesh->E_boundary == std::vector<double>{5.0, 1.0, 0.5, 0.0});

    REQUIRE(deck.materials.contains("water"));
    const MaterialData& water = deck.materials.at("water");
    CHECK(water.sigma_t == std::vector<double>{1.2, 1.0, 0.8});
    REQUIRE(water.scattering.size() == 5);
    CHECK(water.scattering[0].from == 0);
    CHECK(water.scattering[0].to == 0);
    CHECK(water.scattering[0].value == doctest::Approx(1.1));
    CHECK(water.scattering[1].from == 0);
    CHECK(water.scattering[1].to == 1);
    CHECK(water.scattering[1].value == doctest::Approx(0.05));
    CHECK(water.stopping_power_average == std::vector<double>{2.0, 1.8, 1.5});
    CHECK(water.stopping_power_boundary == std::vector<double>{2.2, 1.9, 1.6, 1.3});

    REQUIRE(deck.materials.contains("lead"));

    REQUIRE(deck.xs.has_value());
    CHECK(deck.xs->material == std::vector<std::string>{"water", "water", "lead"});
    CHECK(deck.xs->total[0] == std::vector<double>{1.2, 1.2, 3.2});
    CHECK(deck.xs->total[2] == std::vector<double>{0.8, 0.8, 2.8});
    // Cells 0/1 are water, cell 2 is lead -- scattering is indexed [cell].
    CHECK(deck.xs->scattering[0] == deck.materials.at("water").scattering);
    CHECK(deck.xs->scattering[1] == deck.materials.at("water").scattering);
    CHECK(deck.xs->scattering[2] == deck.materials.at("lead").scattering);
    CHECK(deck.xs->stop_power[1] == std::vector<double>{1.8, 1.8, 4.8});
    CHECK(deck.xs->stop_power_boundary[3] == std::vector<double>{1.3, 1.3, 4.3});

    CHECK(deck.convergence.max_iters == 200);
    CHECK(deck.convergence.epsilon == doctest::Approx(1.0e-8));

    REQUIRE(deck.angular_quadrature.has_value());
    CHECK(deck.angular_quadrature->mu == std::vector<double>{-0.9, -0.3, 0.3, 0.9});
    CHECK(deck.angular_quadrature->w == std::vector<double>{0.5, 0.5, 0.5, 0.5});

    // boundary_conditions.left/right follow down = 10*m + g (left) or
    // 100*m + g (right), up = down + 0.5 -- see sample_input.yaml.
    REQUIRE(deck.boundary_conditions.has_value());
    const BoundaryConditions& bc = *deck.boundary_conditions;
    REQUIRE(bc.left.size() == 4);
    REQUIRE(bc.left[0].size() == 3);
    for (int m = 0; m < 4; ++m) {
      for (int g = 0; g < 3; ++g) {
        CHECK(bc.left[m][g].down == doctest::Approx(10 * m + g));
        CHECK(bc.left[m][g].up == doctest::Approx(10 * m + g + 0.5));
        CHECK(bc.right[m][g].down == doctest::Approx(100 * m + g));
        CHECK(bc.right[m][g].up == doctest::Approx(100 * m + g + 0.5));
      }
    }
  }

  TEST_CASE("rejects a non-ascending spatial mesh") {
    const auto path = writeTempYaml(R"(
spatial_mesh: [0.0, 2.0, 1.0]
regions:
  materials: [water, water]
energy_mesh: [1.0, 0.0]
materials:
  water:
    sigma_t: [1.0]
    scattering: [{from: 0, to: 0, value: 1.0}]
    stopping_power:
      group_average: [1.0]
      group_boundary: [1.0, 1.0]
angular_quadrature:
  mu: [-0.5, 0.5]
  w: [1.0, 1.0]
boundary_conditions:
  left:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
  right:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
convergence:
  max_iters: 1
  epsilon: 1.0
)");
    InputDeck deck;
    CHECK(deck.read(path) == 1);
  }

  TEST_CASE("rejects a region count that does not match the spatial mesh") {
    const auto path = writeTempYaml(R"(
spatial_mesh: [0.0, 1.0, 2.0]
regions:
  materials: [water]
energy_mesh: [1.0, 0.0]
materials:
  water:
    sigma_t: [1.0]
    scattering: [{from: 0, to: 0, value: 1.0}]
    stopping_power:
      group_average: [1.0]
      group_boundary: [1.0, 1.0]
angular_quadrature:
  mu: [-0.5, 0.5]
  w: [1.0, 1.0]
boundary_conditions:
  left:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
  right:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
convergence:
  max_iters: 1
  epsilon: 1.0
)");
    InputDeck deck;
    CHECK(deck.read(path) == 1);
  }

  TEST_CASE("rejects a region referencing an undefined material") {
    const auto path = writeTempYaml(R"(
spatial_mesh: [0.0, 1.0]
regions:
  materials: [copper]
energy_mesh: [1.0, 0.0]
materials:
  water:
    sigma_t: [1.0]
    scattering: [{from: 0, to: 0, value: 1.0}]
    stopping_power:
      group_average: [1.0]
      group_boundary: [1.0, 1.0]
angular_quadrature:
  mu: [-0.5, 0.5]
  w: [1.0, 1.0]
boundary_conditions:
  left:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
  right:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
convergence:
  max_iters: 1
  epsilon: 1.0
)");
    InputDeck deck;
    CHECK(deck.read(path) == 1);
  }

  TEST_CASE("rejects cross sections with the wrong number of groups") {
    const auto path = writeTempYaml(R"(
spatial_mesh: [0.0, 1.0]
regions:
  materials: [water]
energy_mesh: [1.0, 0.5, 0.0]
materials:
  water:
    sigma_t: [1.0]
    scattering: [{from: 0, to: 0, value: 1.0}]
    stopping_power:
      group_average: [1.0]
      group_boundary: [1.0, 1.0]
angular_quadrature:
  mu: [-0.5, 0.5]
  w: [1.0, 1.0]
boundary_conditions:
  left:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
  right:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
convergence:
  max_iters: 1
  epsilon: 1.0
)");
    InputDeck deck;
    CHECK(deck.read(path) == 1);
  }

  TEST_CASE("rejects an energy mesh that isn't strictly descending") {
    const auto path = writeTempYaml(R"(
spatial_mesh: [0.0, 1.0]
regions:
  materials: [water]
energy_mesh: [0.5, 1.0, 0.0]
materials:
  water:
    sigma_t: [1.0, 1.0]
    scattering: [{from: 0, to: 0, value: 1.0}, {from: 1, to: 1, value: 1.0}]
    stopping_power:
      group_average: [1.0, 1.0]
      group_boundary: [1.0, 1.0, 1.0]
angular_quadrature:
  mu: [-0.5, 0.5]
  w: [1.0, 1.0]
boundary_conditions:
  left:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
  right:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
convergence:
  max_iters: 1
  epsilon: 1.0
)");
    InputDeck deck;
    CHECK(deck.read(path) == 1);
  }

  TEST_CASE("rejects an energy mesh that doesn't end at 0") {
    const auto path = writeTempYaml(R"(
spatial_mesh: [0.0, 1.0]
regions:
  materials: [water]
energy_mesh: [1.0, 0.5]
materials:
  water:
    sigma_t: [1.0]
    scattering: [{from: 0, to: 0, value: 1.0}]
    stopping_power:
      group_average: [1.0]
      group_boundary: [1.0, 1.0]
angular_quadrature:
  mu: [-0.5, 0.5]
  w: [1.0, 1.0]
boundary_conditions:
  left:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
  right:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
convergence:
  max_iters: 1
  epsilon: 1.0
)");
    InputDeck deck;
    CHECK(deck.read(path) == 1);
  }

  TEST_CASE("rejects a negative cross section") {
    const auto path = writeTempYaml(R"(
spatial_mesh: [0.0, 1.0]
regions:
  materials: [water]
energy_mesh: [1.0, 0.0]
materials:
  water:
    sigma_t: [-1.0]
    scattering: [{from: 0, to: 0, value: 1.0}]
    stopping_power:
      group_average: [1.0]
      group_boundary: [1.0, 1.0]
angular_quadrature:
  mu: [-0.5, 0.5]
  w: [1.0, 1.0]
boundary_conditions:
  left:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
  right:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
convergence:
  max_iters: 1
  epsilon: 1.0
)");
    InputDeck deck;
    CHECK(deck.read(path) == 1);
  }

  TEST_CASE("rejects an input file missing a required key") {
    const auto path = writeTempYaml(R"(
spatial_mesh: [0.0, 1.0]
regions:
  materials: [water]
materials:
  water:
    sigma_t: [1.0]
    scattering: [{from: 0, to: 0, value: 1.0}]
    stopping_power:
      group_average: [1.0]
      group_boundary: [1.0, 1.0]
angular_quadrature:
  mu: [-0.5, 0.5]
  w: [1.0, 1.0]
boundary_conditions:
  left:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
  right:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
convergence:
  max_iters: 1
  epsilon: 1.0
)");
    InputDeck deck;
    CHECK(deck.read(path) == 1);
  }

  TEST_CASE("rejects a non-ascending angular_quadrature.mu") {
    const auto path = writeTempYaml(R"(
spatial_mesh: [0.0, 1.0]
regions:
  materials: [water]
energy_mesh: [1.0, 0.0]
materials:
  water:
    sigma_t: [1.0]
    scattering: [{from: 0, to: 0, value: 1.0}]
    stopping_power:
      group_average: [1.0]
      group_boundary: [1.0, 1.0]
angular_quadrature:
  mu: [0.5, -0.5]
  w: [1.0, 1.0]
boundary_conditions:
  left:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
  right:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
convergence:
  max_iters: 1
  epsilon: 1.0
)");
    InputDeck deck;
    CHECK(deck.read(path) == 1);
  }

  TEST_CASE("rejects angular_quadrature.w with the wrong number of entries") {
    const auto path = writeTempYaml(R"(
spatial_mesh: [0.0, 1.0]
regions:
  materials: [water]
energy_mesh: [1.0, 0.0]
materials:
  water:
    sigma_t: [1.0]
    scattering: [{from: 0, to: 0, value: 1.0}]
    stopping_power:
      group_average: [1.0]
      group_boundary: [1.0, 1.0]
angular_quadrature:
  mu: [-0.5, 0.5]
  w: [1.0, 1.0, 1.0]
boundary_conditions:
  left:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
  right:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
convergence:
  max_iters: 1
  epsilon: 1.0
)");
    InputDeck deck;
    CHECK(deck.read(path) == 1);
  }

  TEST_CASE("normalizes angular_quadrature.w to sum to 2") {
    const auto path = writeTempYaml(R"(
spatial_mesh: [0.0, 1.0]
regions:
  materials: [water]
energy_mesh: [1.0, 0.0]
materials:
  water:
    sigma_t: [1.0]
    scattering: [{from: 0, to: 0, value: 1.0}]
    stopping_power:
      group_average: [1.0]
      group_boundary: [1.0, 1.0]
angular_quadrature:
  mu: [-0.5, 0.5]
  w: [1.0, 3.0]
boundary_conditions:
  left:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
  right:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
convergence:
  max_iters: 1
  epsilon: 1.0
)");
    InputDeck deck;
    CHECK(deck.read(path) == 0);

    REQUIRE(deck.angular_quadrature.has_value());
    CHECK(deck.angular_quadrature->w[0] == doctest::Approx(0.5));
    CHECK(deck.angular_quadrature->w[1] == doctest::Approx(1.5));

    REQUIRE(deck.boundary_conditions.has_value());
    CHECK(deck.boundary_conditions->left.size() == 2);
    CHECK(deck.boundary_conditions->left[0].size() == 1);
    CHECK(deck.boundary_conditions->left[0][0].down == doctest::Approx(0.0));
    CHECK(deck.boundary_conditions->left[0][0].up == doctest::Approx(0.0));
  }

  TEST_CASE("rejects boundary_conditions.left.down with the wrong number of ordinates") {
    const auto path = writeTempYaml(R"(
spatial_mesh: [0.0, 1.0]
regions:
  materials: [water]
energy_mesh: [1.0, 0.0]
materials:
  water:
    sigma_t: [1.0]
    scattering: [{from: 0, to: 0, value: 1.0}]
    stopping_power:
      group_average: [1.0]
      group_boundary: [1.0, 1.0]
angular_quadrature:
  mu: [-0.5, 0.5]
  w: [1.0, 1.0]
boundary_conditions:
  left:
    down: [[0.0]]
    up: [[0.0], [0.0]]
  right:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
convergence:
  max_iters: 1
  epsilon: 1.0
)");
    InputDeck deck;
    CHECK(deck.read(path) == 1);
  }

  TEST_CASE("rejects boundary_conditions.right.up with the wrong number of groups") {
    const auto path = writeTempYaml(R"(
spatial_mesh: [0.0, 1.0]
regions:
  materials: [water]
energy_mesh: [1.0, 0.0]
materials:
  water:
    sigma_t: [1.0]
    scattering: [{from: 0, to: 0, value: 1.0}]
    stopping_power:
      group_average: [1.0]
      group_boundary: [1.0, 1.0]
angular_quadrature:
  mu: [-0.5, 0.5]
  w: [1.0, 1.0]
boundary_conditions:
  left:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
  right:
    down: [[0.0], [0.0]]
    up: [[0.0, 0.0], [0.0, 0.0]]
convergence:
  max_iters: 1
  epsilon: 1.0
)");
    InputDeck deck;
    CHECK(deck.read(path) == 1);
  }

  TEST_CASE("rejects a material scattering entry with an out-of-range group") {
    const auto path = writeTempYaml(R"(
spatial_mesh: [0.0, 1.0]
regions:
  materials: [water]
energy_mesh: [1.0, 0.0]
materials:
  water:
    sigma_t: [1.0]
    scattering: [{from: 0, to: 1, value: 1.0}]
    stopping_power:
      group_average: [1.0]
      group_boundary: [1.0, 1.0]
angular_quadrature:
  mu: [-0.5, 0.5]
  w: [1.0, 1.0]
boundary_conditions:
  left:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
  right:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
convergence:
  max_iters: 1
  epsilon: 1.0
)");
    InputDeck deck;
    CHECK(deck.read(path) == 1);
  }

  TEST_CASE("rejects a duplicate material scattering entry") {
    const auto path = writeTempYaml(R"(
spatial_mesh: [0.0, 1.0]
regions:
  materials: [water]
energy_mesh: [1.0, 0.0]
materials:
  water:
    sigma_t: [1.0]
    scattering: [{from: 0, to: 0, value: 1.0}, {from: 0, to: 0, value: 2.0}]
    stopping_power:
      group_average: [1.0]
      group_boundary: [1.0, 1.0]
angular_quadrature:
  mu: [-0.5, 0.5]
  w: [1.0, 1.0]
boundary_conditions:
  left:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
  right:
    down: [[0.0], [0.0]]
    up: [[0.0], [0.0]]
convergence:
  max_iters: 1
  epsilon: 1.0
)");
    InputDeck deck;
    CHECK(deck.read(path) == 1);
  }
}

TEST_SUITE("AngularQuadrature") {
  TEST_CASE("rejects mismatched mu/w sizes") {
    CHECK_THROWS_AS(AngularQuadrature({-0.5, 0.5}, {1.0}), std::invalid_argument);
  }

  TEST_CASE("rejects a non-ascending mu") {
    CHECK_THROWS_AS(AngularQuadrature({0.5, -0.5}, {1.0, 1.0}), std::invalid_argument);
  }

  TEST_CASE("normalizes w to sum to 2") {
    const AngularQuadrature quadrature({-0.5, 0.5}, {1.0, 3.0});

    CHECK(quadrature.w[0] == doctest::Approx(0.5));
    CHECK(quadrature.w[1] == doctest::Approx(1.5));
  }
}

TEST_SUITE("BoundaryConditions") {
  TEST_CASE("rejects left with the wrong number of ordinates") {
    CHECK_THROWS_AS(BoundaryConditions({{{0.0, 0.0}}}, {{{0.0, 0.0}}, {{0.0, 0.0}}}, 2, 1),
                    std::invalid_argument);
  }

  TEST_CASE("rejects right with the wrong number of groups") {
    // left: 2 ordinates x 1 group each (matches expected). right: 2
    // ordinates x 2 groups each (expected only 1 group).
    CHECK_THROWS_AS(BoundaryConditions({{{0.0, 0.0}}, {{0.0, 0.0}}},
                                       {{{0.0, 0.0}, {0.0, 0.0}}, {{0.0, 0.0}, {0.0, 0.0}}}, 2, 1),
                    std::invalid_argument);
  }

  TEST_CASE("accepts left/right matching num_ordinates x num_groups") {
    const BoundaryConditions bc({{{1.0, 2.0}}, {{3.0, 4.0}}}, {{{5.0, 6.0}}, {{7.0, 8.0}}}, 2, 1);

    CHECK(bc.left[0][0].down == 1.0);
    CHECK(bc.left[0][0].up == 2.0);
    CHECK(bc.right[1][0].down == 7.0);
    CHECK(bc.right[1][0].up == 8.0);
  }
}

TEST_SUITE("InputDeck setters") {
  TEST_CASE("setMesh replaces mesh and leaves xs unset if it was never set") {
    InputDeck deck;
    deck.setMesh(Mesh({0.0, 1.0}, {1.0, 0.0}));

    REQUIRE(deck.mesh.has_value());
    CHECK(deck.mesh->n_x == 1);
    CHECK_FALSE(deck.xs.has_value());
  }

  TEST_CASE("setMesh clears an existing xs") {
    InputDeck deck;
    REQUIRE(deck.read(sampleInputPath()) == 0);
    REQUIRE(deck.xs.has_value());

    deck.setMesh(Mesh({0.0, 1.0, 2.0, 3.0}, {5.0, 1.0, 0.5, 0.0}));

    CHECK_FALSE(deck.xs.has_value());
  }

  TEST_CASE("setAngularQuadrature replaces angular_quadrature") {
    InputDeck deck;
    deck.setAngularQuadrature(AngularQuadrature({-0.5, 0.5}, {1.0, 1.0}));

    REQUIRE(deck.angular_quadrature.has_value());
    CHECK(deck.angular_quadrature->mu == std::vector<double>{-0.5, 0.5});
  }

  TEST_CASE("setBoundaryConditions replaces boundary_conditions") {
    InputDeck deck;
    deck.setBoundaryConditions(BoundaryConditions({{{1.0, 2.0}}}, {{{3.0, 4.0}}}, 1, 1));

    REQUIRE(deck.boundary_conditions.has_value());
    CHECK(deck.boundary_conditions->left[0][0].down == 1.0);
  }

  TEST_CASE("setConvergence replaces convergence") {
    InputDeck deck;
    deck.setConvergence(ConvergenceCriteria{42, 1.0e-6});

    CHECK(deck.convergence.max_iters == 42);
    CHECK(deck.convergence.epsilon == doctest::Approx(1.0e-6));
  }
}
