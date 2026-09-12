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
  deck.xs.scatter = Eigen::MatrixXd::Constant(2, 2, 0.5);
  deck.xs.S = Eigen::MatrixXd::Constant(2, 2, 1.0);
  deck.xs.S_bound = Eigen::MatrixXd::Constant(3, 2, 1.0);
  deck.bc.values = Eigen::MatrixXd::Constant(4, 2, 0.0);
  return deck;
}
} // namespace

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
  deck.xs.scatter = Eigen::MatrixXd::Constant(2, 3, 0.5);

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

TEST_CASE("InputDeck::Xs::validate accepts non-negative tables") {
  InputDeck::Xs xs;
  xs.total = Eigen::MatrixXd::Constant(2, 3, 1.0); // G=2, n_x=3
  xs.scatter = Eigen::MatrixXd::Constant(2, 3, 0.0);
  xs.S = Eigen::MatrixXd::Constant(2, 3, 1.0);
  xs.S_bound = Eigen::MatrixXd::Constant(3, 3, 1.0); // G+1=3

  CHECK_NOTHROW(xs.validate());
}

TEST_CASE("InputDeck::Xs::validate rejects a negative total entry") {
  InputDeck::Xs xs;
  xs.total = Eigen::MatrixXd::Constant(2, 3, 1.0);
  xs.total(0, 0) = -1.0;
  xs.scatter = Eigen::MatrixXd::Constant(2, 3, 0.0);
  xs.S = Eigen::MatrixXd::Constant(2, 3, 1.0);
  xs.S_bound = Eigen::MatrixXd::Constant(3, 3, 1.0);

  CHECK_THROWS_AS(xs.validate(), std::runtime_error);
}

TEST_CASE("InputDeck::Xs::validate rejects a negative S_bound entry") {
  InputDeck::Xs xs;
  xs.total = Eigen::MatrixXd::Constant(2, 3, 1.0);
  xs.scatter = Eigen::MatrixXd::Constant(2, 3, 0.0);
  xs.S = Eigen::MatrixXd::Constant(2, 3, 1.0);
  xs.S_bound = Eigen::MatrixXd::Constant(3, 3, 1.0);
  xs.S_bound(2, 1) = -1.0;

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

TEST_CASE("InputDeck::BoundaryConditions::validate accepts non-negative values") {
  InputDeck::BoundaryConditions bc;
  bc.values = Eigen::MatrixXd::Constant(4, 2, 1.0);

  CHECK_NOTHROW(bc.validate());
}

TEST_CASE("InputDeck::BoundaryConditions::validate rejects negative values") {
  InputDeck::BoundaryConditions bc;
  bc.values = Eigen::MatrixXd::Constant(4, 2, 1.0);
  bc.values(0, 0) = -1.0;

  CHECK_THROWS_AS(bc.validate(), std::runtime_error);
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
