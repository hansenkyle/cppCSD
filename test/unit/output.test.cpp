// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "output.h"

#include <stdexcept>

#include <doctest.h>

TEST_SUITE("OutputTable") {
  TEST_CASE("rejects a non-positive n_rows") {
    CHECK_THROWS_AS(OutputTable("TEST", 0), std::invalid_argument);
    CHECK_THROWS_AS(OutputTable("TEST", -1), std::invalid_argument);
  }

  TEST_CASE("addColumn rejects a column with the wrong number of values") {
    OutputTable table("TEST", 2);

    CHECK_THROWS_AS(
        table.addColumn("x", {OutputTable::Notation::Fixed, 4}, std::vector<double>{1.0}),
        std::invalid_argument);
  }

  TEST_CASE("txt() renders a title, right-justified header, and right-justified rows") {
    OutputTable table("TEST", 2);
    table.addColumn("i", {OutputTable::Notation::Integer, 0}, {0.0, 1.0});
    table.addColumn("x", {OutputTable::Notation::Fixed, 4}, {0.5, 1.5});

    const std::string expected = "--- TEST ---\n"
                                 "  i       x\n"
                                 "  0  0.5000\n"
                                 "  1  1.5000\n";
    CHECK(table.txt() == expected);
  }

  TEST_CASE("csv() renders a title row, comma-separated header, and comma-separated rows") {
    OutputTable table("TEST", 2);
    table.addColumn("i", {OutputTable::Notation::Integer, 0}, {0.0, 1.0});
    table.addColumn("x", {OutputTable::Notation::Fixed, 4}, {0.5, 1.5});

    const std::string expected = "TEST\n"
                                 "i,x\n"
                                 "0,0.5000\n"
                                 "1,1.5000\n";
    CHECK(table.csv() == expected);
  }

  TEST_CASE("column width grows to fit the widest value, not just the header") {
    OutputTable table("TEST", 2);
    table.addColumn("phi", {OutputTable::Notation::Scientific, 2}, {1.0, -20.0});

    // "-2.00e+01" (9 chars) is wider than the header "phi" (3 chars), so the
    // column must widen to fit it, not truncate or misalign.
    const std::string expected = "--- TEST ---\n"
                                 "        phi\n"
                                 "   1.00e+00\n"
                                 "  -2.00e+01\n";
    CHECK(table.txt() == expected);
  }

  TEST_CASE("addRow rejects a row with the wrong number of values") {
    OutputTable table("TEST", 2);

    CHECK_THROWS_AS(table.addRow("x", {OutputTable::Notation::Fixed, 4}, std::vector<double>{1.0}),
                    std::invalid_argument);
  }

  TEST_CASE("addColumn and addRow cannot be mixed on the same table") {
    OutputTable column_first("TEST", 2);
    column_first.addColumn("x", {OutputTable::Notation::Fixed, 1}, {0.5, 1.5});
    CHECK_THROWS_AS(column_first.addRow("y", {OutputTable::Notation::Fixed, 1}, {0.5, 1.5}),
                    std::invalid_argument);

    OutputTable row_first("TEST", 2);
    row_first.addRow("x", {OutputTable::Notation::Fixed, 1}, {0.5, 1.5});
    CHECK_THROWS_AS(row_first.addColumn("y", {OutputTable::Notation::Fixed, 1}, {0.5, 1.5}),
                    std::invalid_argument);
  }

  TEST_CASE("row-major txt() numbers columns in the header and labels each row") {
    OutputTable table("TEST", 2);
    table.addRow("x", {OutputTable::Notation::Fixed, 1}, {0.5, 1.5});
    table.addRow("y", {OutputTable::Notation::Integer, 0}, {2.0, 3.0});

    const std::string expected = "--- TEST ---\n"
                                 "       0    1\n"
                                 "x    0.5  1.5\n"
                                 "y      2    3\n";
    CHECK(table.txt() == expected);
  }

  TEST_CASE("row-major column width grows to fit the widest value in that column") {
    OutputTable table("TEST", 2);
    table.addRow("phi", {OutputTable::Notation::Scientific, 2}, {1.0, -20.0});

    const std::string expected = "--- TEST ---\n"
                                 "              0          1\n"
                                 "phi    1.00e+00  -2.00e+01\n";
    CHECK(table.txt() == expected);
  }

  TEST_CASE("row-major csv() renders a blank corner, numbered header, and labeled rows") {
    OutputTable table("TEST", 2);
    table.addRow("x", {OutputTable::Notation::Fixed, 1}, {0.5, 1.5});
    table.addRow("y", {OutputTable::Notation::Integer, 0}, {2.0, 3.0});

    const std::string expected = "TEST\n"
                                 ",0,1\n"
                                 "x,0.5,1.5\n"
                                 "y,2,3\n";
    CHECK(table.csv() == expected);
  }
}

TEST_SUITE("OutputMetadata") {
  TEST_CASE("txt() left-justifies keys and right-justifies values, colon-separated") {
    OutputMetadata metadata("META");
    metadata.addEntry("A", "1");
    metadata.addEntry("BB", "22");

    const std::string expected = "--- META ---\n"
                                 "A  :  1\n"
                                 "BB : 22\n";
    CHECK(metadata.txt() == expected);
  }

  TEST_CASE("csv() renders plain, unaligned key,value rows") {
    OutputMetadata metadata("META");
    metadata.addEntry("A", "1");
    metadata.addEntry("BB", "22");

    const std::string expected = "META\n"
                                 "A,1\n"
                                 "BB,22\n";
    CHECK(metadata.csv() == expected);
  }

  TEST_CASE("entries render in insertion order") {
    OutputMetadata metadata("META");
    metadata.addEntry("Date", "8-24-2026");
    metadata.addEntry("Time", "13:01:32");
    metadata.addEntry("Energy groups", "13");

    const std::string result = metadata.txt();
    CHECK(result.find("Date") < result.find("Time"));
    CHECK(result.find("Time") < result.find("Energy groups"));
  }

  TEST_CASE("metadata with no entries still renders its title") {
    OutputMetadata metadata("EMPTY");

    CHECK(metadata.txt() == "--- EMPTY ---\n");
    CHECK(metadata.csv() == "EMPTY\n");
  }
}
