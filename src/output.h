// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#ifndef OUTPUT_H
#define OUTPUT_H

#include <string>
#include <vector>

// A titled table of named, formatted columns, one row per index -- for
// solver output files (e.g. solution.txt / solution.csv), not run metadata;
// see docs/output-layout.md and OutputMetadata for that. txt()/csv() render
// the table to a self-contained string that the caller appends to whatever
// file they're assembling; OutputTable never touches the filesystem itself.
class OutputTable {
public:
  // How a column's values are rendered: as an integer, fixed-point, or
  // scientific notation with `precision` digits after the decimal point
  // (ignored for Integer).
  enum class Notation { Integer, Fixed, Scientific };
  struct Format {
    Notation notation;
    int precision;
  };

  // Constructs a table with `n_rows` rows. n_rows must be positive; every
  // column added afterward must have exactly n_rows values.
  OutputTable(std::string title, int n_rows);

  // Appends a named column. `values` must have exactly n_rows entries.
  void addColumn(std::string name, Format format, std::vector<double> values);

  // Renders as human-readable, whitespace-aligned text: a title line, then
  // a right-justified header row and data rows. Column widths are computed
  // from content, not hard-coded.
  std::string txt() const;

  // Renders as CSV: the title as a plain (uncommented) row, then
  // comma-separated header and data rows.
  std::string csv() const;

private:
  struct Column {
    std::string name;
    Format format;
    std::vector<double> values;
  };

  std::string title_;
  int n_rows_;
  std::vector<Column> columns_;
};

// A titled list of key-value pairs (e.g. an output file header), in
// insertion order.
class OutputMetadata {
public:
  explicit OutputMetadata(std::string title);

  // Appends one key: value pair.
  void addEntry(std::string key, std::string value);

  // Renders as "key : value" lines, with keys left-justified and values
  // right-justified to widths computed from content.
  std::string txt() const;

  // Renders as plain "key,value" rows, unaligned.
  std::string csv() const;

private:
  struct Entry {
    std::string key;
    std::string value;
  };

  std::string title_;
  std::vector<Entry> entries_;
};

#endif
