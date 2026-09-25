// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "method.h"

#include <doctest.h>

#include <cmath>

namespace {
// Method is abstract with a protected constructor; this stub exposes it so the non-virtual
// helpers can be exercised directly.
class StubMethod : public Method {
public:
  explicit StubMethod(InputDeck deck) : Method("stub", deck) {}
  void solve(double, int) override {}
};

// Two cells of unequal width (dx = 0.5, 1.5), so a norm that ignores or misindexes dx gives a
// different answer. The norms only read the mesh, so the rest of the deck is left empty.
InputDeck makeNormDeck() {
  InputDeck deck;
  deck.mesh.n_x = 2;
  deck.mesh.x_boundary = Eigen::Vector3d(0.0, 0.5, 2.0);
  deck.mesh.validate();
  deck.energy.G = 1;
  return deck;
}

// Corner values for both cells, laid out per cell as [up_L, up_R, down_L, down_R].
Eigen::VectorXd makeCorners(std::initializer_list<double> values) {
  Eigen::VectorXd v(values.size());
  Eigen::Index k = 0;
  for (double x : values) {
    v(k++) = x;
  }
  return v;
}
} // namespace

TEST_CASE("Method::l2norm is zero for a zero vector") {
  const StubMethod method(makeNormDeck());
  CHECK(method.l2norm(Eigen::VectorXd::Zero(8)) == 0.0);
}

TEST_CASE("Method::l2norm integrates a constant field exactly") {
  StubMethod method(makeNormDeck());
  // int_0^2 int_dE 3^2 dE dx = 9 * 2 * 0.5 = 9
  Eigen::VectorXd v = Eigen::VectorXd::Constant(8, 3.0);
  CHECK(method.l2norm(v, 0.5) == doctest::Approx(3.0));
}

TEST_CASE("Method::l2norm integrates a field linear in x exactly") {
  StubMethod method(makeNormDeck());
  // f(x) = x across both cells: int_0^2 x^2 dx = 8/3
  Eigen::VectorXd v = makeCorners({0, 0.5, 0, 0.5, 0.5, 2, 0.5, 2});
  CHECK(method.l2norm(v) == doctest::Approx(std::sqrt(8.0 / 3.0)));
}

TEST_CASE("Method::l2norm integrates a field linear in energy exactly") {
  StubMethod method(makeNormDeck());
  // f = 1 at the upper group bound, 0 at the lower: int_dE f^2 dE = dE / 3, times L = 2
  Eigen::VectorXd v = makeCorners({1, 1, 0, 0, 1, 1, 0, 0});
  CHECK(method.l2norm(v, 0.6) == doctest::Approx(std::sqrt(2.0 * 0.6 / 3.0)));
}

TEST_CASE("Method::l2norm weights each cell by its width") {
  StubMethod method(makeNormDeck());
  // cell 0: 0.5 * (240 / 36) = 10/3, cell 1: 1.5 * 1 = 3/2
  Eigen::VectorXd v = makeCorners({1, 2, 3, 4, 1, 1, 1, 1});
  CHECK(method.l2norm(v) == doctest::Approx(std::sqrt(29.0 / 6.0)));
}

TEST_CASE("Method::l2norm reads each cell's own corners") {
  StubMethod method(makeNormDeck());
  // A single unit corner integrates to dx / 9 in the cell it lives in.
  CHECK(method.l2norm(makeCorners({1, 0, 0, 0, 0, 0, 0, 0})) ==
        doctest::Approx(std::sqrt(0.5 / 9.0)));
  CHECK(method.l2norm(makeCorners({0, 0, 0, 0, 0, 0, 0, 1})) ==
        doctest::Approx(std::sqrt(1.5 / 9.0)));
}

TEST_CASE("Method::l2norm scales with sqrt(dE)") {
  StubMethod method(makeNormDeck());
  Eigen::VectorXd v = makeCorners({1, 2, 3, 4, 1, 1, 1, 1});
  CHECK(method.l2norm(v, 4.0) == doctest::Approx(2.0 * method.l2norm(v)));
  CHECK(method.l2norm(v, 1.0) == doctest::Approx(method.l2norm(v)));
}

TEST_CASE("Method::l2norm is independent of sign") {
  StubMethod method(makeNormDeck());
  Eigen::VectorXd v = makeCorners({1, -2, 3, -4, 0.5, 0, -1, 2});
  CHECK(method.l2norm(-v) == doctest::Approx(method.l2norm(v)));
}

TEST_CASE("Method::linfnorm is zero for a zero vector") {
  CHECK(Method::linfnorm(Eigen::VectorXd::Zero(8)) == 0.0);
}

TEST_CASE("Method::linfnorm returns the largest magnitude, including negative entries") {
  CHECK(Method::linfnorm(makeCorners({1, 2, 3, 4, 1, 1, 1, 1})) == 4.0);
  CHECK(Method::linfnorm(makeCorners({1, -5, 3, 4, 1, 1, 1, 1})) == 5.0);
}
