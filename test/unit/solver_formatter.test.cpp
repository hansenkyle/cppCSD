// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "solver_formatter.h"

#include <doctest.h>

namespace {

// A minimally-valid, 1-cell/1-group/1-ordinate deck, small enough that
// formatted numbers can be checked by literal string.
InputDeck makeDeck() {
  InputDeck deck;
  deck.mesh.n_x = 1;
  deck.mesh.x_boundary = Eigen::Vector2d(0.0, 1.0);
  deck.energy.G = 1;
  deck.energy.E_boundary = Eigen::Vector2d(1.0, 0.0);
  deck.angle.M = 1;
  deck.angle.mu = Eigen::VectorXd::Constant(1, 0.5);
  deck.angle.w = Eigen::VectorXd::Constant(1, 2.0);
  deck.xs.total = Eigen::MatrixXd::Constant(1, 1, 1.0);
  deck.xs.scatter = {Eigen::SparseMatrix<double>(1, 1)};
  deck.xs.S = Eigen::MatrixXd::Constant(1, 1, 1.0);
  deck.xs.S_bound = Eigen::MatrixXd::Constant(2, 1, 1.0);
  deck.bc.values = Eigen::MatrixXd::Constant(2, 1, 0.0);
  deck.source.values = {Eigen::MatrixXd::Constant(4, 1, 0.0)};
  deck.validate();
  return deck;
}

} // namespace

TEST_SUITE("SolverFormatter") {
  TEST_CASE("formatRunMetadata renders a run time") {
    const std::string result = SolverFormatter::formatRunMetadata();

    CHECK(result.find("Run Metadata") != std::string::npos);
    CHECK(result.find("Run time") != std::string::npos);
  }

  TEST_CASE("formatResults renders scalar flux as a 2x2 corner cluster per cell") {
    Eigen::MatrixXd scalar_flux(4, 1);
    scalar_flux << 1.0, 2.0, 3.0, 4.0; // up_left, up_right, down_left, down_right

    const std::string result = SolverFormatter::formatResults(scalar_flux, {}, makeDeck());

    CHECK(result.find("Scalar Flux") != std::string::npos);
    CHECK(result.find("1.000000e+00") != std::string::npos);
    CHECK(result.find("2.000000e+00") != std::string::npos);
    CHECK(result.find("3.000000e+00") != std::string::npos);
    CHECK(result.find("4.000000e+00") != std::string::npos);
  }

  TEST_CASE("formatResults renders one angular flux table per group, with each ordinate's mu and "
            "weight") {
    Eigen::MatrixXd group_flux(4, 1); // 1 ordinate
    group_flux << 5.0, 6.0, 7.0, 8.0;

    const std::string result =
        SolverFormatter::formatResults(Eigen::MatrixXd::Zero(4, 1), {group_flux}, makeDeck());

    CHECK(result.find("Angular Flux - Group 0") != std::string::npos);
    CHECK(result.find("mu=5.000000e-01") != std::string::npos);
    CHECK(result.find("w=2.000000e+00") != std::string::npos);
    CHECK(result.find("5.000000e+00") != std::string::npos);
    CHECK(result.find("8.000000e+00") != std::string::npos);
  }

  TEST_CASE("formatResiduals renders one table per group") {
    Eigen::MatrixXd residual(4, 1);
    residual << 0.1, 0.2, 0.3, 0.4;

    const std::string result = SolverFormatter::formatResiduals({residual}, makeDeck());

    CHECK(result.find("Angular Residual - Group 0") != std::string::npos);
    CHECK(result.find("1.000000e-01") != std::string::npos);
    CHECK(result.find("4.000000e-01") != std::string::npos);
  }

  TEST_CASE("formatResiduals renders nothing for an empty residual list") {
    CHECK(SolverFormatter::formatResiduals({}, makeDeck()).empty());
  }
}
