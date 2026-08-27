// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "output.h"

#include <stdexcept>

#include <doctest.h>

TEST_SUITE("OutputBlock1d") {
  TEST_CASE("rejects a non-positive n_rows") {
    CHECK_THROWS_AS(OutputBlock1d("TEST", 0), std::invalid_argument);
    CHECK_THROWS_AS(OutputBlock1d("TEST", -1), std::invalid_argument);
  }

  TEST_CASE("addColumn rejects a column with the wrong number of values") {
    OutputBlock1d block("TEST", 2);

    CHECK_THROWS_AS(
        block.addColumn("x", {OutputBlock::Notation::Fixed, 4}, std::vector<double>{1.0}),
        std::invalid_argument);
  }

  TEST_CASE("txt() renders a title, right-justified header, and right-justified rows") {
    OutputBlock1d block("TEST", 2);
    block.addColumn("i", {OutputBlock::Notation::Integer, 0}, {0.0, 1.0});
    block.addColumn("x", {OutputBlock::Notation::Fixed, 4}, {0.5, 1.5});

    const std::string expected = "--- TEST ---\n"
                                 "  i       x\n"
                                 "  0  0.5000\n"
                                 "  1  1.5000\n";
    CHECK(block.txt() == expected);
  }

  TEST_CASE("csv() renders a title row, comma-separated header, and comma-separated rows") {
    OutputBlock1d block("TEST", 2);
    block.addColumn("i", {OutputBlock::Notation::Integer, 0}, {0.0, 1.0});
    block.addColumn("x", {OutputBlock::Notation::Fixed, 4}, {0.5, 1.5});

    const std::string expected = "TEST\n"
                                 "i,x\n"
                                 "0,0.5000\n"
                                 "1,1.5000\n";
    CHECK(block.csv() == expected);
  }

  TEST_CASE("column width grows to fit the widest value, not just the header") {
    OutputBlock1d block("TEST", 2);
    block.addColumn("phi", {OutputBlock::Notation::Scientific, 2}, {1.0, -20.0});

    // "-2.00e+01" (9 chars) is wider than the header "phi" (3 chars), so the
    // column must widen to fit it, not truncate or misalign.
    const std::string expected = "--- TEST ---\n"
                                 "        phi\n"
                                 "   1.00e+00\n"
                                 "  -2.00e+01\n";
    CHECK(block.txt() == expected);
  }
}

TEST_SUITE("OutputBlock2d") {
  TEST_CASE("rejects a non-positive n_x or G") {
    CHECK_THROWS_AS(OutputBlock2d("TEST", 0, 2), std::invalid_argument);
    CHECK_THROWS_AS(OutputBlock2d("TEST", 2, 0), std::invalid_argument);
  }

  TEST_CASE("addColumn rejects a value with the wrong number of groups") {
    OutputBlock2d block("TEST", 2, 2);

    CHECK_THROWS_AS(block.addColumn("phi", {OutputBlock::Notation::Scientific, 2},
                                    {std::vector<double>{1.0, 2.0}}),
                    std::invalid_argument);
  }

  TEST_CASE("addColumn rejects a group with the wrong number of cells") {
    OutputBlock2d block("TEST", 2, 2);

    CHECK_THROWS_AS(
        block.addColumn("phi", {OutputBlock::Notation::Scientific, 2}, {{1.0, 2.0}, {3.0}}),
        std::invalid_argument);
  }

  TEST_CASE("txt() auto-generates the i column with no columns added") {
    OutputBlock2d block("TEST", 2, 1);

    const std::string expected = "--- TEST ---\n"
                                 "g=0\n"
                                 "  i\n"
                                 "  0\n"
                                 "  1\n";
    CHECK(block.txt() == expected);
  }

  TEST_CASE("txt() renders one subtitled sub-table per group, blank-line separated") {
    OutputBlock2d block("FLUX", 2, 2);
    block.addColumn("phi", {OutputBlock::Notation::Scientific, 2}, {{1.0, 2.0}, {3.0, 4.0}});

    const std::string expected = "--- FLUX ---\n"
                                 "g=0\n"
                                 "  i       phi\n"
                                 "  0  1.00e+00\n"
                                 "  1  2.00e+00\n"
                                 "\n"
                                 "g=1\n"
                                 "  i       phi\n"
                                 "  0  3.00e+00\n"
                                 "  1  4.00e+00\n";
    CHECK(block.txt() == expected);
  }

  TEST_CASE("csv() renders one subtitle row and table per group, blank-line separated") {
    OutputBlock2d block("FLUX", 2, 2);
    block.addColumn("phi", {OutputBlock::Notation::Scientific, 2}, {{1.0, 2.0}, {3.0, 4.0}});

    const std::string expected = "FLUX\n"
                                 "g=0\n"
                                 "i,phi\n"
                                 "0,1.00e+00\n"
                                 "1,2.00e+00\n"
                                 "\n"
                                 "g=1\n"
                                 "i,phi\n"
                                 "0,3.00e+00\n"
                                 "1,4.00e+00\n";
    CHECK(block.csv() == expected);
  }

  TEST_CASE("a column added once applies to every group's sub-table") {
    OutputBlock2d block("TEST", 1, 3);
    block.addColumn("phi", {OutputBlock::Notation::Fixed, 1}, {{1.0}, {2.0}, {3.0}});

    const std::string result = block.txt();
    CHECK(result.find("g=0") != std::string::npos);
    CHECK(result.find("g=1") != std::string::npos);
    CHECK(result.find("g=2") != std::string::npos);
    CHECK(result.find("1.0") != std::string::npos);
    CHECK(result.find("2.0") != std::string::npos);
    CHECK(result.find("3.0") != std::string::npos);
  }
}

TEST_SUITE("OutputBlockMetadata") {
  TEST_CASE("txt() left-justifies keys and right-justifies values, colon-separated") {
    OutputBlockMetadata block("META");
    block.addEntry("A", "1");
    block.addEntry("BB", "22");

    const std::string expected = "--- META ---\n"
                                 "A  :  1\n"
                                 "BB : 22\n";
    CHECK(block.txt() == expected);
  }

  TEST_CASE("csv() renders plain, unaligned key,value rows") {
    OutputBlockMetadata block("META");
    block.addEntry("A", "1");
    block.addEntry("BB", "22");

    const std::string expected = "META\n"
                                 "A,1\n"
                                 "BB,22\n";
    CHECK(block.csv() == expected);
  }

  TEST_CASE("entries render in insertion order") {
    OutputBlockMetadata block("META");
    block.addEntry("Date", "8-24-2026");
    block.addEntry("Time", "13:01:32");
    block.addEntry("Energy groups", "13");

    const std::string result = block.txt();
    CHECK(result.find("Date") < result.find("Time"));
    CHECK(result.find("Time") < result.find("Energy groups"));
  }

  TEST_CASE("a block with no entries still renders its title") {
    OutputBlockMetadata block("EMPTY");

    CHECK(block.txt() == "--- EMPTY ---\n");
    CHECK(block.csv() == "EMPTY\n");
  }

  TEST_CASE("polymorphic access through OutputBlock dispatches to the override") {
    OutputBlockMetadata block("META");
    block.addEntry("A", "1");

    const OutputBlock& base = block;
    CHECK(base.txt() == block.txt());
  }
}
