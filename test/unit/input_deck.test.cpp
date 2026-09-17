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
// A minimally-valid deck: n_x = 2 cells, G = 2 groups, M = 2 ordinates.
InputDeck makeValidDeck() {
  InputDeck deck;
  deck.mesh.n_x = 2;
  deck.mesh.x_boundary = Eigen::Vector3d(0.0, 1.0, 2.0);
  deck.energy.G = 2;
  deck.energy.E_boundary = Eigen::Vector3d(2.0, 1.0, 0.0);
  deck.angle.M = 2;
  deck.angle.mu = Eigen::Vector2d(-0.5, 0.5);
  deck.angle.w = Eigen::Vector2d(1.0, 1.0);
  deck.xs.total = Eigen::MatrixXd::Constant(2, 2, 1.0);
  deck.xs.scatter = {{{0, 0, 0.5}}, {{0, 0, 0.5}}}; // one in-group entry per cell
  deck.xs.S = Eigen::MatrixXd::Constant(2, 2, 1.0);
  deck.xs.S_bound = Eigen::MatrixXd::Constant(3, 2, 1.0);
  deck.bc.values = Eigen::MatrixXd::Constant(4, 2, 0.0);
  deck.source.values = {Eigen::MatrixXd::Constant(8, 2, 0.0), Eigen::MatrixXd::Constant(8, 2, 0.0)};
  return deck;
}

// Writes a minimal one-cell, one-group deck YAML with the given
// `scattering:` list substituted in, and reads it back.
int readWithScattering(const std::string& scattering_yaml) {
  const std::filesystem::path path =
      std::filesystem::temp_directory_path() / "ldcsd_scattering_test.yaml";
  {
    std::ofstream out(path);
    out << R"(
spatial_mesh: [0.0, 1.0]
regions:
  materials: [water]
energy_mesh: [1.0, 0.0]
materials:
  water:
    sigma_t: [1.0]
    scattering: )"
        << scattering_yaml << R"(
    stopping_power:
      group_average: [1.0]
      group_boundary: [1.0, 1.0]
angular_quadrature:
  mu: [-0.5, 0.5]
  w: [1.0, 1.0]
boundary_conditions:
  down: [[0.0, 0.0]]
  up: [[0.0, 0.0]]
source:
  - [[{up_left: 0.0, up_right: 0.0, down_left: 0.0, down_right: 0.0}],
     [{up_left: 0.0, up_right: 0.0, down_left: 0.0, down_right: 0.0}]]
)";
  }
  InputDeck deck;
  const int result = deck.read(path);
  std::filesystem::remove(path);
  return result;
}
} // namespace

TEST_CASE("InputDeck::read parses the sample deck") {
  InputDeck deck;

  REQUIRE(deck.read(std::filesystem::path(TEST_DATA_DIR) / "sample_input.yaml") == 0);

  CHECK(deck.mesh.n_x == 3);
  CHECK(deck.mesh.x_boundary[3] == doctest::Approx(3.0));
  CHECK(deck.mesh.dx[0] == doctest::Approx(1.0));

  CHECK(deck.energy.G == 3);
  CHECK(deck.energy.E_boundary[0] == doctest::Approx(5.0));
  CHECK(deck.energy.dE[0] == doctest::Approx(4.0));

  CHECK(deck.angle.M == 4);
  CHECK(deck.angle.mu[0] == doctest::Approx(-0.9));
  CHECK(deck.angle.w.sum() == doctest::Approx(2.0));

  // Cell 0 and 1 are water, cell 2 is lead.
  CHECK(deck.xs.total(0, 0) == doctest::Approx(1.2));
  CHECK(deck.xs.total(0, 2) == doctest::Approx(3.2));
  CHECK(deck.xs.S(2, 0) == doctest::Approx(1.5));
  CHECK(deck.xs.S_bound(3, 2) == doctest::Approx(4.3));

  // Every cell using a material gets a copy of that material's sparse
  // scattering entries, in file order; entry 3 is the 1->2 downscatter term.
  REQUIRE(deck.xs.scatter[1].size() == 5); // water
  CHECK(deck.xs.scatter[1][3].from == 1);
  CHECK(deck.xs.scatter[1][3].to == 2);
  CHECK(deck.xs.scatter[1][3].value == doctest::Approx(0.03));

  REQUIRE(deck.xs.scatter[2].size() == 5); // lead
  CHECK(deck.xs.scatter[2][3].from == 1);
  CHECK(deck.xs.scatter[2][3].to == 2);
  CHECK(deck.xs.scatter[2][3].value == doctest::Approx(0.08));

  // bc.values row 2*g is group g's up row, row 2*g + 1 is down; down[g][m] =
  // 10*m + g, up = down + 0.5.
  CHECK(deck.bc.values(0, 1) == doctest::Approx(10.5)); // group 0, up, m=1
  CHECK(deck.bc.values(1, 1) == doctest::Approx(10.0)); // group 0, down, m=1
  CHECK(deck.bc.values(4, 3) == doctest::Approx(32.5)); // group 2, up, m=3

  // source.values[g](4*c + corner, m) = 100*m + 10*g + c + corner offset
  // (0.1/0.2/0.3/0.4 for up_left/up_right/down_left/down_right).
  REQUIRE(deck.source.values.size() == 3);
  CHECK(deck.source.values[0].rows() == 4 * 3);                         // 4 * n_x
  CHECK(deck.source.values[0].cols() == 4);                             // M
  CHECK(deck.source.values[0](0, 0) == doctest::Approx(0.1));           // g=0, m=0, c=0, up_left
  CHECK(deck.source.values[0](1, 0) == doctest::Approx(0.2));           // g=0, m=0, c=0, up_right
  CHECK(deck.source.values[0](2, 0) == doctest::Approx(0.3));           // g=0, m=0, c=0, down_left
  CHECK(deck.source.values[0](3, 0) == doctest::Approx(0.4));           // g=0, m=0, c=0, down_right
  CHECK(deck.source.values[1](4 * 1, 2) == doctest::Approx(211.1));     // g=1, m=2, c=1, up_left
  CHECK(deck.source.values[2](4 * 2 + 3, 3) == doctest::Approx(322.4)); // g=2, m=3, c=2, down_right
}

TEST_CASE("InputDeck::read overrides an already-populated deck") {
  InputDeck deck;
  deck.mesh.n_x = 99;
  deck.mesh.x_boundary = Eigen::VectorXd::LinSpaced(100, 0.0, 99.0);
  deck.energy.G = 42;

  REQUIRE(deck.read(std::filesystem::path(TEST_DATA_DIR) / "sample_input.yaml") == 0);

  CHECK(deck.mesh.n_x == 3);
  CHECK(deck.energy.G == 3);
}

TEST_CASE("InputDeck::read returns 1 and leaves no crash on a missing file") {
  InputDeck deck;

  CHECK(deck.read(std::filesystem::path(TEST_DATA_DIR) / "does_not_exist.yaml") == 1);
}

TEST_CASE("InputDeck::read returns 1 on a missing required key") {
  const std::filesystem::path path =
      std::filesystem::temp_directory_path() / "ldcsd_missing_key.yaml";
  {
    std::ofstream out(path);
    out << "spatial_mesh: [0.0, 1.0]\n";
  }
  InputDeck deck;

  CHECK(deck.read(path) == 1);

  std::filesystem::remove(path);
}

TEST_CASE("InputDeck::read returns 1 when a region references an undefined material") {
  const std::filesystem::path path =
      std::filesystem::temp_directory_path() / "ldcsd_bad_region.yaml";
  {
    std::ofstream out(path);
    out << R"(
spatial_mesh: [0.0, 1.0]
regions:
  materials: [unobtainium]
energy_mesh: [1.0, 0.0]
materials:
  water:
    sigma_t: [1.0]
    scattering: [{from: 1, to: 1, value: 0.5}]
    stopping_power:
      group_average: [1.0]
      group_boundary: [1.0, 1.0]
angular_quadrature:
  mu: [-0.5, 0.5]
  w: [1.0, 1.0]
boundary_conditions:
  down: [[0.0, 0.0]]
  up: [[0.0, 0.0]]
)";
  }
  InputDeck deck;

  CHECK(deck.read(path) == 1);

  std::filesystem::remove(path);
}

TEST_CASE("InputDeck::read accepts 1-indexed scattering group indices") {
  CHECK(readWithScattering("[{from: 1, to: 1, value: 0.5}]") == 0);
}

TEST_CASE("InputDeck::read returns 1 when scattering 'from' is 0 (not 1-indexed)") {
  CHECK(readWithScattering("[{from: 0, to: 1, value: 0.5}]") == 1);
}

TEST_CASE("InputDeck::read returns 1 when scattering 'to' exceeds num_groups") {
  CHECK(readWithScattering("[{from: 1, to: 2, value: 0.5}]") == 1); // G = 1
}

TEST_CASE("InputDeck::read returns 1 on a duplicate scattering (from, to) pair") {
  CHECK(readWithScattering("[{from: 1, to: 1, value: 0.5}, {from: 1, to: 1, value: 0.1}]") == 1);
}

TEST_CASE("InputDeck::read accepts a dense scattering matrix") {
  CHECK(readWithScattering("[[0.5]]") == 0);
}

TEST_CASE("InputDeck::read returns 1 when a dense scattering matrix has the wrong row count") {
  CHECK(readWithScattering("[[0.5], [0.1]]") == 1); // G = 1, but two rows given
}

TEST_CASE("InputDeck::read returns 1 when a dense scattering matrix row has the wrong length") {
  CHECK(readWithScattering("[[0.5, 0.1]]") == 1); // G = 1, but row has two columns
}

TEST_CASE("InputDeck::validate accepts a consistent deck") {
  InputDeck deck = makeValidDeck();

  CHECK_NOTHROW(deck.validate());
}

TEST_CASE("InputDeck::validate rejects xs shaped against the wrong number of groups") {
  InputDeck deck = makeValidDeck();
  deck.xs.total = Eigen::MatrixXd::Constant(3, 2, 1.0);

  CHECK_THROWS_AS(deck.validate(), std::runtime_error);
}

TEST_CASE("InputDeck::validate rejects xs shaped against the wrong number of cells") {
  InputDeck deck = makeValidDeck();
  deck.xs.scatter.resize(3); // mesh.n_x is 2

  CHECK_THROWS_AS(deck.validate(), std::runtime_error);
}

TEST_CASE("InputDeck::validate rejects S_bound with the wrong number of rows") {
  InputDeck deck = makeValidDeck();
  deck.xs.S_bound = Eigen::MatrixXd::Constant(2, 2, 1.0);

  CHECK_THROWS_AS(deck.validate(), std::runtime_error);
}

TEST_CASE("InputDeck::validate rejects bc shaped against the wrong number of groups") {
  InputDeck deck = makeValidDeck();
  deck.bc.values = Eigen::MatrixXd::Constant(2, 2, 0.0);

  CHECK_THROWS_AS(deck.validate(), std::runtime_error);
}

TEST_CASE("InputDeck::validate rejects bc shaped against the wrong number of ordinates") {
  InputDeck deck = makeValidDeck();
  deck.bc.values = Eigen::MatrixXd::Constant(4, 3, 0.0);

  CHECK_THROWS_AS(deck.validate(), std::runtime_error);
}

TEST_CASE("InputDeck::validate rejects source with the wrong number of groups") {
  InputDeck deck = makeValidDeck();
  deck.source.values = {Eigen::MatrixXd::Constant(8, 2, 0.0)}; // energy.G is 2

  CHECK_THROWS_AS(deck.validate(), std::runtime_error);
}

TEST_CASE("InputDeck::validate rejects a source group shaped against the wrong number of "
          "cells") {
  InputDeck deck = makeValidDeck();
  deck.source.values[0] = Eigen::MatrixXd::Constant(4, 2, 0.0); // mesh.n_x is 2, expects 8 rows

  CHECK_THROWS_AS(deck.validate(), std::runtime_error);
}

TEST_CASE("InputDeck::validate rejects a source group shaped against the wrong number of "
          "ordinates") {
  InputDeck deck = makeValidDeck();
  deck.source.values[0] = Eigen::MatrixXd::Constant(8, 3, 0.0); // angle.M is 2

  CHECK_THROWS_AS(deck.validate(), std::runtime_error);
}

TEST_CASE("InputDeck::Source::validate accepts non-negative values") {
  InputDeck::Source source;
  source.values = {Eigen::MatrixXd::Constant(4, 2, 1.0)};

  CHECK_NOTHROW(source.validate());
}

TEST_CASE("InputDeck::Source::validate rejects a negative entry") {
  InputDeck::Source source;
  source.values = {Eigen::MatrixXd::Constant(4, 2, 1.0)};
  source.values[0](0, 0) = -1.0;

  CHECK_THROWS_AS(source.validate(), std::runtime_error);
}

TEST_CASE("InputDeck::Xs::validate accepts non-negative tables") {
  InputDeck::Xs xs;
  xs.total = Eigen::MatrixXd::Constant(2, 3, 1.0); // G=2, n_x=3
  xs.scatter = {{{0, 1, 0.2}}, {}, {{1, 0, 0.1}}};
  xs.S = Eigen::MatrixXd::Constant(2, 3, 1.0);
  xs.S_bound = Eigen::MatrixXd::Constant(3, 3, 1.0); // G+1=3

  CHECK_NOTHROW(xs.validate());
}

TEST_CASE("InputDeck::Xs::validate rejects a negative total entry") {
  InputDeck::Xs xs;
  xs.total = Eigen::MatrixXd::Constant(2, 3, 1.0);
  xs.total(0, 0) = -1.0;
  xs.scatter = {{}, {}, {}};
  xs.S = Eigen::MatrixXd::Constant(2, 3, 1.0);
  xs.S_bound = Eigen::MatrixXd::Constant(3, 3, 1.0);

  CHECK_THROWS_AS(xs.validate(), std::runtime_error);
}

TEST_CASE("InputDeck::Xs::validate rejects a negative S_bound entry") {
  InputDeck::Xs xs;
  xs.total = Eigen::MatrixXd::Constant(2, 3, 1.0);
  xs.scatter = {{}, {}, {}};
  xs.S = Eigen::MatrixXd::Constant(2, 3, 1.0);
  xs.S_bound = Eigen::MatrixXd::Constant(3, 3, 1.0);
  xs.S_bound(2, 1) = -1.0;

  CHECK_THROWS_AS(xs.validate(), std::runtime_error);
}

TEST_CASE("InputDeck::Xs::validate rejects a negative scatter entry") {
  InputDeck::Xs xs;
  xs.total = Eigen::MatrixXd::Constant(2, 3, 1.0);
  xs.scatter = {{{0, 1, -0.2}}, {}, {}};
  xs.S = Eigen::MatrixXd::Constant(2, 3, 1.0);
  xs.S_bound = Eigen::MatrixXd::Constant(3, 3, 1.0);

  CHECK_THROWS_AS(xs.validate(), std::runtime_error);
}

TEST_CASE("InputDeck::BoundaryConditions::operator[] views a group's up/down rows") {
  InputDeck::BoundaryConditions bc;
  bc.values = Eigen::MatrixXd(4, 2); // G = 2, M = 2
  bc.values << 1.0, 2.0,             // group 0, up
      3.0, 4.0,                      // group 0, down
      5.0, 6.0,                      // group 1, up
      7.0, 8.0;                      // group 1, down

  const Eigen::MatrixXd group1 = bc[1];

  REQUIRE(group1.rows() == 2);
  REQUIRE(group1.cols() == 2);
  CHECK(group1(0, 0) == doctest::Approx(5.0)); // up
  CHECK(group1(1, 1) == doctest::Approx(8.0)); // down
}

TEST_CASE("InputDeck::Energy::validate derives dE from E_boundary") {
  InputDeck::Energy energy;
  energy.G = 3;
  energy.E_boundary = Eigen::Vector4d(6.0, 3.0, 1.0, 0.0);

  energy.validate();

  REQUIRE(energy.dE.size() == 3);
  CHECK(energy.dE[0] == doctest::Approx(3.0));
  CHECK(energy.dE[1] == doctest::Approx(2.0));
  CHECK(energy.dE[2] == doctest::Approx(1.0));
}

TEST_CASE("InputDeck::Energy::validate rejects non-positive G") {
  InputDeck::Energy energy;
  energy.G = 0;
  energy.E_boundary = Eigen::VectorXd(1);
  energy.E_boundary[0] = 0.0;

  CHECK_THROWS_AS(energy.validate(), std::runtime_error);
}

TEST_CASE("InputDeck::Energy::validate rejects mismatched E_boundary size") {
  InputDeck::Energy energy;
  energy.G = 3;
  energy.E_boundary = Eigen::Vector2d(1.0, 0.0);

  CHECK_THROWS_AS(energy.validate(), std::runtime_error);
}

TEST_CASE("InputDeck::Energy::validate rejects non-descending E_boundary") {
  InputDeck::Energy energy;
  energy.G = 2;
  energy.E_boundary = Eigen::Vector3d(1.0, 1.0, 0.0);

  CHECK_THROWS_AS(energy.validate(), std::runtime_error);
}

TEST_CASE("InputDeck::Angle::validate accepts an exact quadrature unchanged") {
  InputDeck::Angle angle;
  angle.M = 2;
  angle.mu = Eigen::Vector2d(-0.5773502692, 0.5773502692);
  angle.w = Eigen::Vector2d(1.0, 1.0);

  angle.validate();

  CHECK(angle.w[0] == doctest::Approx(1.0));
  CHECK(angle.w[1] == doctest::Approx(1.0));
}

TEST_CASE("InputDeck::Angle::validate renormalizes w within tolerance") {
  InputDeck::Angle angle;
  angle.M = 2;
  angle.mu = Eigen::Vector2d(-0.5773502692, 0.5773502692);
  angle.w = Eigen::Vector2d(0.99995, 0.99995); // sums to 1.9999, rel diff 5e-5 < 1e-4

  angle.validate();

  CHECK(angle.w.sum() == doctest::Approx(2.0));
}

TEST_CASE("InputDeck::Angle::validate rejects w outside tolerance") {
  InputDeck::Angle angle;
  angle.M = 2;
  angle.mu = Eigen::Vector2d(-0.5773502692, 0.5773502692);
  angle.w = Eigen::Vector2d(0.9, 0.9); // sums to 1.8, rel diff 0.1 > 1e-4

  CHECK_THROWS_AS(angle.validate(), std::runtime_error);
}

TEST_CASE("InputDeck::Angle::validate rejects non-ascending mu") {
  InputDeck::Angle angle;
  angle.M = 2;
  angle.mu = Eigen::Vector2d(0.5, 0.5);
  angle.w = Eigen::Vector2d(1.0, 1.0);

  CHECK_THROWS_AS(angle.validate(), std::runtime_error);
}

TEST_CASE("InputDeck::Mesh::validate derives dx from x_boundary") {
  InputDeck::Mesh mesh;
  mesh.n_x = 3;
  mesh.x_boundary = Eigen::Vector4d(0.0, 1.0, 3.0, 6.0);

  mesh.validate();

  REQUIRE(mesh.dx.size() == 3);
  CHECK(mesh.dx[0] == doctest::Approx(1.0));
  CHECK(mesh.dx[1] == doctest::Approx(2.0));
  CHECK(mesh.dx[2] == doctest::Approx(3.0));
}

TEST_CASE("InputDeck::Mesh::validate rejects non-positive n_x") {
  InputDeck::Mesh mesh;
  mesh.n_x = 0;
  mesh.x_boundary = Eigen::VectorXd(1);
  mesh.x_boundary[0] = 0.0;

  CHECK_THROWS_AS(mesh.validate(), std::runtime_error);
}

TEST_CASE("InputDeck::Mesh::validate rejects mismatched x_boundary size") {
  InputDeck::Mesh mesh;
  mesh.n_x = 3;
  mesh.x_boundary = Eigen::Vector2d(0.0, 1.0);

  CHECK_THROWS_AS(mesh.validate(), std::runtime_error);
}

TEST_CASE("InputDeck::Mesh::validate rejects non-ascending x_boundary") {
  InputDeck::Mesh mesh;
  mesh.n_x = 2;
  mesh.x_boundary = Eigen::Vector3d(0.0, 1.0, 1.0);

  CHECK_THROWS_AS(mesh.validate(), std::runtime_error);
}

TEST_CASE("InputDeck::Angle::validate rejects mismatched sizes") {
  InputDeck::Angle angle;
  angle.M = 2;
  angle.mu = Eigen::Vector2d(-0.5, 0.5);
  angle.w = Eigen::VectorXd::Constant(3, 1.0);

  CHECK_THROWS_AS(angle.validate(), std::runtime_error);
}

TEST_SUITE("InputDeck::echo") {
  TEST_CASE("renders every section in order, with no material-wise section") {
    InputDeck deck = makeValidDeck();
    deck.validate(); // populates dx/dE, which echo() reads

    const std::string result = deck.echo();

    const std::vector<std::string> expected_headers = {"Problem Size",
                                                       "Spatial Discretization",
                                                       "Energy Discretization",
                                                       "Angular Quadrature",
                                                       "Cross Sections - Group 0",
                                                       "Cross Sections - Group 1",
                                                       "Scattering Matrix - Cell 0",
                                                       "Scattering Matrix - Cell 1",
                                                       "Boundary Conditions",
                                                       "External Source - Group 0, Ordinate 0",
                                                       "External Source - Group 0, Ordinate 1",
                                                       "External Source - Group 1, Ordinate 0",
                                                       "External Source - Group 1, Ordinate 1"};

    std::size_t last_pos = 0;
    for (const std::string& header : expected_headers) {
      const std::size_t pos = result.find(header);
      CHECK(pos != std::string::npos);
      CHECK(pos >= last_pos);
      last_pos = pos;
    }

    // InputDeck doesn't retain material names/grouping, so there's nothing
    // to echo a material-wise section from.
    CHECK(result.find("material") == std::string::npos);
    CHECK(result.find("Material") == std::string::npos);
  }

  TEST_CASE("reports problem size and formats cross sections at 6-digit scientific precision") {
    InputDeck deck = makeValidDeck();
    deck.validate();

    const std::string result = deck.echo();

    CHECK(result.find("Spatial cells") != std::string::npos);
    CHECK(result.find("1.000000e+00") != std::string::npos); // xs.total / xs.S entries
    CHECK(result.find("5.000000e-01") != std::string::npos); // within-group scatter, 0.5
  }

  TEST_CASE("renders external source corner values") {
    InputDeck deck = makeValidDeck();
    deck.source.values[1](4, 0) = 7.5; // group 1, cell 1, up_left, ordinate 0
    deck.validate();

    const std::string result = deck.echo();

    CHECK(result.find("up_left") != std::string::npos);
    CHECK(result.find("up_right") != std::string::npos);
    CHECK(result.find("down_left") != std::string::npos);
    CHECK(result.find("down_right") != std::string::npos);
    CHECK(result.find("7.500000e+00") != std::string::npos);
  }
}
