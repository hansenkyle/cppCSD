// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "output_block.h"

#include <Eigen/Core>
#include <doctest.h>

TEST_SUITE("KeyValueOutput") {
  TEST_CASE("render_txt brackets the title and pads the dot runs to a common value column") {
    KeyValueOutput meta("META");
    meta.add("A", "1");
    meta.add("BB", "22");

    const std::string expected = "[META]\n"
                                 "\n"
                                 "A  ....  1\n"
                                 "BB  ...  22\n";
    CHECK(meta.render_txt() == expected);
  }

  TEST_CASE("integer keys and values are converted to strings") {
    KeyValueOutput meta("META");
    meta.add("groups", 13);
    meta.add(7, "seven");
    meta.add(1, 2);

    const std::string expected = "[META]\n"
                                 "\n"
                                 "groups  ...  13\n"
                                 "7  ........  seven\n"
                                 "1  ........  2\n";
    CHECK(meta.render_txt() == expected);
  }

  TEST_CASE("entries render in insertion order") {
    KeyValueOutput meta("META");
    meta.add("Date", "8-24-2026");
    meta.add("Time", "13:01:32");
    meta.add("Energy groups", 13);

    const std::string result = meta.render_txt();
    CHECK(result.find("Date") < result.find("Time"));
    CHECK(result.find("Time") < result.find("Energy groups"));
  }

  TEST_CASE("an empty block still renders its title") {
    CHECK(KeyValueOutput("EMPTY").render_txt() == "[EMPTY]\n\n");
  }

  TEST_CASE("an untitled block renders only its entries") {
    KeyValueOutput meta;
    meta.add("A", "1");
    meta.add("BB", "22");

    const std::string expected = "A  ....  1\n"
                                 "BB  ...  22\n";
    CHECK(meta.render_txt() == expected);
  }

  TEST_CASE("tabs indent every line of the block by four spaces each") {
    KeyValueOutput meta("META");
    meta.add("A", "1");
    meta.add("BB", "22");

    const std::string expected = "        [META]\n"
                                 "\n"
                                 "        A  ....  1\n"
                                 "        BB  ...  22\n";
    CHECK(meta.render_txt(2) == expected);
  }
}

TEST_SUITE("VerticalTable") {
  TEST_CASE("render_txt puts the column names on top, one row per index") {
    VerticalTable table("TEST");
    table.add_column("i", {0.0, 1.0});
    table.add_column("x", {0.5, 1.5});

    const std::string expected = "[TEST]\n"
                                 "\n"
                                 "  i    x\n"
                                 "0.0  0.5\n"
                                 "1.0  1.5\n";
    CHECK(table.render_txt("{:.1f}") == expected);
  }

  TEST_CASE("the default format is 4-decimal e-notation") {
    VerticalTable table("TEST");
    table.add_column("scalar_flux", {1.2345e12});

    const std::string expected = "[TEST]\n"
                                 "\n"
                                 "scalar_flux\n"
                                 " 1.2345e+12\n";
    CHECK(table.render_txt() == expected);
  }

  TEST_CASE("columns widen to fit whichever is longer, the name or a value") {
    VerticalTable table("TEST");
    table.add_column("phi", {1.0, -20.0});
    table.add_column("n", {3.0, 4.0});

    const std::string expected = "[TEST]\n"
                                 "\n"
                                 "      phi         n\n"
                                 " 1.00e+00  3.00e+00\n"
                                 "-2.00e+01  4.00e+00\n";
    CHECK(table.render_txt("{:.2e}") == expected);
  }

  TEST_CASE("a negative value widens its column so the digits still line up") {
    VerticalTable table("TEST");
    table.add_column("a", {1.23, -1.23});
    table.add_column("b", {-4.5, 6.5});

    // The minus sign hangs into the column's extra space rather than
    // shifting 1.23 over: both "1" glyphs sit in the same text column.
    const std::string expected = "[TEST]\n"
                                 "\n"
                                 "    a      b\n"
                                 " 1.23  -4.50\n"
                                 "-1.23   6.50\n";
    CHECK(table.render_txt("{:.2f}") == expected);
  }

  TEST_CASE("a short column runs out partway down and leaves the rest of its cells blank") {
    VerticalTable table("TEST");
    table.add_column("x", {0.5, 1.5, 2.5});
    table.add_column("y", {1.0});

    const std::string expected = "[TEST]\n"
                                 "\n"
                                 "  x    y\n"
                                 "0.5  1.0\n"
                                 "1.5\n"
                                 "2.5\n";
    CHECK(table.render_txt("{:.1f}") == expected);
  }

  TEST_CASE("a gap in the middle of a row keeps the columns after it aligned") {
    VerticalTable table("TEST");
    table.add_column("x", {0.5, 1.5});
    table.add_column("shortish", {1.0});
    table.add_column("z", {2.0, 3.0});

    const std::string expected = "[TEST]\n"
                                 "\n"
                                 "  x  shortish    z\n"
                                 "0.5       1.0  2.0\n"
                                 "1.5            3.0\n";
    CHECK(table.render_txt("{:.1f}") == expected);
  }

  TEST_CASE("add_column accepts an Eigen vector") {
    Eigen::VectorXd x(2);
    x << 0.5, 1.5;

    VerticalTable table("TEST");
    table.add_column("i", {0.0, 1.0});
    table.add_column("x", x);

    const std::string expected = "[TEST]\n"
                                 "\n"
                                 "  i    x\n"
                                 "0.0  0.5\n"
                                 "1.0  1.5\n";
    CHECK(table.render_txt("{:.1f}") == expected);
  }

  TEST_CASE("add_column accepts a matrix row, which is not a contiguous vector") {
    Eigen::MatrixXd m(2, 2);
    m << 0.5, 1.5, 2.5, 3.5;

    VerticalTable table;
    table.add_column("row0", m.row(0));
    table.add_column("col1", m.col(1));

    const std::string expected = "row0  col1\n"
                                 " 0.5   1.5\n"
                                 " 1.5   3.5\n";
    CHECK(table.render_txt("{:.1f}") == expected);
  }

  TEST_CASE("an empty table still renders its title") {
    CHECK(VerticalTable("EMPTY").render_txt() == "[EMPTY]\n\n");
  }

  TEST_CASE("an untitled table renders only its rows") {
    VerticalTable table;
    table.add_column("i", {0.0, 1.0});
    table.add_column("x", {0.5, 1.5});

    const std::string expected = "  i    x\n"
                                 "0.0  0.5\n"
                                 "1.0  1.5\n";
    CHECK(table.render_txt("{:.1f}") == expected);
  }

  TEST_CASE("tabs indent every line of the table by four spaces each") {
    VerticalTable table("TEST");
    table.add_column("i", {0.0, 1.0});
    table.add_column("x", {0.5, 1.5});

    const std::string expected = "    [TEST]\n"
                                 "\n"
                                 "      i    x\n"
                                 "    0.0  0.5\n"
                                 "    1.0  1.5\n";
    CHECK(table.render_txt("{:.1f}", 1) == expected);
  }
}

TEST_SUITE("HorizontalTable") {
  TEST_CASE("render_txt puts the row names down the left, one row per series") {
    HorizontalTable table("TEST");
    table.add_row("x", {0.5, 1.5});
    table.add_row("yy", {2.0, 3.0});

    const std::string expected = "[TEST]\n"
                                 "\n"
                                 "x   0.5  1.5\n"
                                 "yy  2.0  3.0\n";
    CHECK(table.render_txt("{:.1f}") == expected);
  }

  TEST_CASE("the default format is 4-decimal e-notation") {
    HorizontalTable table("TEST");
    table.add_row("scalar_flux", {1.2345e12, 1.5e12});

    const std::string expected = "[TEST]\n"
                                 "\n"
                                 "scalar_flux  1.2345e+12  1.5000e+12\n";
    CHECK(table.render_txt() == expected);
  }

  TEST_CASE("a negative value widens its column so the digits still line up") {
    HorizontalTable table("TEST");
    table.add_row("a", {1.23, -1.23});
    table.add_row("bb", {-4.5, 6.5});

    // Same as the vertical case, per value column: 1.23 and -4.50 line up
    // on their digits, with the minus sign in the leading pad.
    const std::string expected = "[TEST]\n"
                                 "\n"
                                 "a    1.23  -1.23\n"
                                 "bb  -4.50   6.50\n";
    CHECK(table.render_txt("{:.2f}") == expected);
  }

  TEST_CASE("rows of different lengths each just end where their data does") {
    HorizontalTable table("TEST");
    table.add_row("x", {0.5, 1.5, 2.5});
    table.add_row("yy", {2.0});

    const std::string expected = "[TEST]\n"
                                 "\n"
                                 "x   0.5  1.5  2.5\n"
                                 "yy  2.0\n";
    CHECK(table.render_txt("{:.1f}") == expected);
  }

  TEST_CASE("a longer row added after a short one still widens the shared columns") {
    HorizontalTable table;
    table.add_row("x", {0.5});
    table.add_row("yy", {-20.0, 3.0});

    const std::string expected = "x     0.5\n"
                                 "yy  -20.0  3.0\n";
    CHECK(table.render_txt("{:.1f}") == expected);
  }

  TEST_CASE("add_row accepts an Eigen vector") {
    Eigen::VectorXd values(2);
    values << 1.2345e12, 1.5e12;

    HorizontalTable table("TEST");
    table.add_row("scalar_flux", values);

    const std::string expected = "[TEST]\n"
                                 "\n"
                                 "scalar_flux  1.2345e+12  1.5000e+12\n";
    CHECK(table.render_txt() == expected);
  }

  TEST_CASE("add_row accepts an Eigen expression") {
    Eigen::VectorXd values(2);
    values << 1.0, 2.0;

    HorizontalTable table;
    table.add_row("doubled", 2.0 * values);

    CHECK(table.render_txt("{:.1f}") == "doubled  2.0  4.0\n");
  }

  TEST_CASE("an empty table still renders its title") {
    CHECK(HorizontalTable("EMPTY").render_txt() == "[EMPTY]\n\n");
  }

  TEST_CASE("an untitled table renders only its rows") {
    HorizontalTable table;
    table.add_row("x", {0.5, 1.5});
    table.add_row("yy", {2.0, 3.0});

    const std::string expected = "x   0.5  1.5\n"
                                 "yy  2.0  3.0\n";
    CHECK(table.render_txt("{:.1f}") == expected);
  }

  TEST_CASE("tabs indent every line of the table by four spaces each") {
    HorizontalTable table("TEST");
    table.add_row("x", {0.5, 1.5});
    table.add_row("yy", {2.0, 3.0});

    const std::string expected = "    [TEST]\n"
                                 "\n"
                                 "    x   0.5  1.5\n"
                                 "    yy  2.0  3.0\n";
    CHECK(table.render_txt("{:.1f}", 1) == expected);
  }
}

TEST_SUITE("OutputUnit") {
  TEST_CASE("set_title replaces the title used by render_txt") {
    KeyValueOutput meta;
    meta.set_title("LATER");

    CHECK(meta.title() == "LATER");
    CHECK(meta.render_txt() == "[LATER]\n\n");
  }

  TEST_CASE("set_title(\"\") drops the title line from a block that had one") {
    KeyValueOutput meta("META");
    meta.add("A", "1");
    REQUIRE(meta.render_txt() == "[META]\n\nA  ...  1\n");

    meta.set_title("");

    CHECK(meta.render_txt() == "A  ...  1\n");
  }

  TEST_CASE("an untitled unit with no data renders nothing at all") {
    CHECK(KeyValueOutput().render_txt() == "");
    CHECK(VerticalTable().render_txt() == "");
    CHECK(HorizontalTable().render_txt() == "");
  }

  TEST_CASE("tabs indent an untitled block's body, with no leading blank line") {
    KeyValueOutput meta;
    meta.add("A", "1");

    CHECK(meta.render_txt(1) == "    A  ...  1\n");
  }

  TEST_CASE("the blank line under the title stays blank rather than becoming whitespace") {
    KeyValueOutput meta("META");
    meta.add("A", "1");

    CHECK(meta.render_txt(1) == "    [META]\n\n    A  ...  1\n");
  }

  TEST_CASE("a non-positive tab count leaves the block unindented") {
    KeyValueOutput meta("META");
    meta.add("A", "1");

    CHECK(meta.render_txt(0) == meta.render_txt());
    CHECK(meta.render_txt(-1) == meta.render_txt());
  }
}
