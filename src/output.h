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

// A titled table of named, formatted values -- for solver output files
// (e.g. solution.txt / solution.csv), not run metadata; see
// docs/output-layout.md and OutputMetadata for that. txt()/csv() render the
// table to a self-contained string that the caller appends to whatever file
// they're assembling; OutputTable never touches the filesystem itself.
//
// A table is either column-major (addColumn: a named column of values, one
// per row; rows are numbered 0..n_entries-1 implicitly) or row-major
// (addRow: a named row of values, one per column; columns are numbered
// 0..n_entries-1 in the header, each row prefixed with its name) --
// row-major suits a "one row per variable, one column per index" layout
// (e.g. cross sections tabulated one column per spatial cell). Orientation
// is decided by whichever of addColumn()/addRow() is called first; mixing
// the two on one table is an error.
class OutputTable {
public:
  // How a value is rendered: as an integer, fixed-point, or scientific
  // notation with `precision` digits after the decimal point (ignored for
  // Integer).
  enum class Notation { Integer, Fixed, Scientific };
  struct Format {
    Notation notation;
    int precision;
  };

  // Constructs a table with `n_entries` values along its indexed axis --
  // rows for a column-major table, columns for a row-major one. n_entries
  // must be positive.
  OutputTable(std::string title, int n_entries);

  // Appends a named column of n_entries values, one per row. Column-major;
  // cannot be mixed with addRow() on the same table.
  void addColumn(std::string name, Format format, std::vector<double> values);

  // Appends a named row of n_entries values, one per column. Row-major;
  // cannot be mixed with addColumn() on the same table.
  void addRow(std::string name, Format format, std::vector<double> values);

  // Renders as human-readable, whitespace-aligned text: a title line, then
  // a right-justified header row and data rows (column-major), or a header
  // row of numbered columns and left-labeled data rows (row-major). Column
  // widths are computed from content, not hard-coded.
  std::string txt() const;

  // Renders as CSV: the title as a plain (uncommented) row, then
  // comma-separated header and data rows.
  std::string csv() const;

private:
  struct Entry {
    std::string name;
    Format format;
    std::vector<double> values;
  };
  enum class Orientation { Unset, ColumnMajor, RowMajor };

  std::string txtColumnMajor() const;
  std::string txtRowMajor() const;

  std::string title_;
  int n_entries_;
  Orientation orientation_ = Orientation::Unset;
  std::vector<Entry> entries_;
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
