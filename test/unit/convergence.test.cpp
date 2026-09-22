// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include <doctest.h>

#include "convergence.h"

TEST_SUITE("convergence") {

  TEST_CASE("an empty history is trivially converged and reports nothing") {
    ConvergenceHistory history;
    CHECK(history.empty());
    CHECK(history.totalIterations() == 0);
    CHECK(history.totalSeconds() == doctest::Approx(0.0));
    CHECK(history.allConverged());
    CHECK(history.unconvergedGroups().empty());
    CHECK(history.worstResidual() == doctest::Approx(-1.0));
  }

  TEST_CASE("record() derives the relative difference from the norms it's given") {
    ConvergenceHistory history;
    history.record(IterationRecord{0, 1, 4.0, 1.0, 0.0, 1e-15});

    REQUIRE(history.records().size() == 1);
    CHECK(history.records()[0].rel_diff == doctest::Approx(0.25));
  }

  TEST_CASE("a zero flux norm gives a zero relative difference, not a division by zero") {
    ConvergenceHistory history;
    history.record(IterationRecord{0, 1, 0.0, 0.0, 0.0, 0.0});

    CHECK(history.records()[0].rel_diff == doctest::Approx(0.0));
  }

  TEST_CASE("the group roll-up tracks the latest iteration but the worst residual") {
    ConvergenceHistory history;
    history.record(IterationRecord{0, 1, 10.0, 5.0, 0.0, 1e-9});
    history.record(IterationRecord{0, 2, 10.0, 1.0, 0.0, 1e-4}); // residual spike
    history.record(IterationRecord{0, 3, 10.0, 1e-9, 0.0, 1e-12});
    history.finishGroup(0, true, 0.5);

    REQUIRE(history.groups().size() == 1);
    const GroupSummary& summary = history.groups()[0];
    CHECK(summary.iterations == 3);
    CHECK(summary.abs_diff == doctest::Approx(1e-9));
    CHECK(summary.max_residual == doctest::Approx(1e-4));
    CHECK(summary.seconds == doctest::Approx(0.5));
    CHECK(summary.converged);
  }

  TEST_CASE("groups appear in the order they were first recorded") {
    ConvergenceHistory history;
    history.record(IterationRecord{0, 1, 1.0, 1.0, 0.0, 0.0});
    history.record(IterationRecord{1, 1, 1.0, 1.0, 0.0, 0.0});
    history.record(IterationRecord{2, 1, 1.0, 1.0, 0.0, 0.0});
    history.record(IterationRecord{1, 2, 1.0, 0.5, 0.0, 0.0}); // back to group 1

    REQUIRE(history.groups().size() == 3);
    CHECK(history.groups()[0].group == 0);
    CHECK(history.groups()[1].group == 1);
    CHECK(history.groups()[2].group == 2);
    CHECK(history.groups()[1].iterations == 2);
    CHECK(history.totalIterations() == 4);
  }

  TEST_CASE("unconverged groups are reported by index") {
    ConvergenceHistory history;
    history.record(IterationRecord{0, 1, 1.0, 1e-12, 0.0, 0.0});
    history.finishGroup(0, true, 0.1);
    history.record(IterationRecord{1, 1000, 1.0, 1.0, 0.0, 0.0});
    history.finishGroup(1, false, 9.0);

    CHECK_FALSE(history.allConverged());
    CHECK(history.unconvergedGroups() == std::vector<int>{1});
    CHECK(history.totalSeconds() == doctest::Approx(9.1));
  }
}
