// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#ifndef OUTPUT_BLOCK_H
#define OUTPUT_BLOCK_H

#include <string>
#include <vector>

// A titled, formatted table of numerical results for solver output files
// (e.g. solution.txt / solution.csv) -- not run metadata; see
// docs/output-layout.md and OutputManager for that. txt()/csv() render the
// block to a self-contained string that the caller appends to whatever file
// they're assembling; OutputBlock never touches the filesystem itself.
//
// Storage and rendering live entirely in this base class, built around one
// shared unit: a "sub-table" (an optional subtitle plus its own named,
// formatted columns). A 1D block is just one sub-table with no subtitle; a
// 2D block is one sub-table per second-axis value (e.g. energy group), each
// subtitled accordingly. Derived classes differ only in how they populate
// subtables_ via addColumn.

/// @class OutputBlock
/// @brief Formatted, titled table with member functions to prepare data for output files
/// @details Base class handles data storage and rendering, the derived classes handle different
/// methods of populating the data.
class OutputBlock {
public:
  // How a column's values are rendered: as an integer, fixed-point, or
  // scientific notation with `precision` digits after the decimal point
  // (ignored for Integer).
  enum class Notation { Integer, Fixed, Scientific };
  struct Format {
    Notation notation;
    int precision;
  };

  explicit OutputBlock(std::string title);
  virtual ~OutputBlock() = default;

  // Renders the block as human-readable, whitespace-aligned text: a title
  // line, then each sub-table's optional subtitle, header row, and data
  // rows, blank-line separated. Column widths are computed from content, not
  // hard-coded. The default implementation renders subtables_; overridden by
  // block types (e.g. OutputBlockMetadata) whose data doesn't fit that shape.
  virtual std::string txt() const;

  // Renders the block as CSV: the title and each sub-table's subtitle as
  // plain (uncommented) rows, comma-separated header and data rows, with a
  // blank line between sub-tables. See txt() re: overriding.
  virtual std::string csv() const;

  struct Column {
    std::string name;
    Format format;
    std::vector<double> values;
  };

  struct SubTable {
    std::string subtitle; // empty => no subtitle line rendered
    std::vector<Column> columns;
  };

protected:
  std::string title_;
  std::vector<SubTable> subtables_;
};

/// @class OutputBlock1d
/// @brief 1D result table: one row per index, 1+ named columns
/// @details A 1D result table: one row per index, one or more named columns, each with its own
/// Format (e.g. an integer index column alongside fixed- and scientific-notation columns).
class OutputBlock1d : public OutputBlock {
public:
  // Constructs a block with `n_rows` rows. n_rows must be positive; every
  // column added afterward must have exactly n_rows values.
  OutputBlock1d(std::string title, int n_rows);

  // Appends a named column. `values` must have exactly n_rows entries.
  void addColumn(std::string name, Format format, std::vector<double> values);

private:
  int n_rows_;
};

/// @class OutputBlock2d
/// @brief 2d result table: one sub-table per energy group (subtitle is hard-coded)
/// @details A 2D, group-resolved result table: one sub-table per energy group (subtitled
/// "g=<group>"), each with an auto-generated integer "i" column (0..n_x-1) and whatever named
/// columns are added. A column added once applies across all groups, matching the [group][cell]
/// layout used elsewhere in this codebase (CrossSection).
class OutputBlock2d : public OutputBlock {
public:
  // Constructs a block over n_x cells and G groups. Both must be positive.
  OutputBlock2d(std::string title, int n_x, int G);

  // Appends a named column shared across all groups. `values` must be
  // [group][cell]-shaped: G outer entries, each with n_x values.
  void addColumn(std::string name, Format format, const std::vector<std::vector<double>>& values);

private:
  int n_x_;
  int G_;
};

/// @class OutputBlockMetadata
/// @brief A titled list of key-value pairs (ex. output file header)
/// @details A titled list of key: value pairs (run metadata -- date, time, dimensions, etc.), in
/// insertion order. Unlike OutputBlock1d/2d, values are caller-supplied strings rather than
/// formatted doubles -- this block's job is alignment, not number formatting -- so it overrides
/// txt()/csv() rather than populating subtables_.
class OutputBlockMetadata : public OutputBlock {
public:
  explicit OutputBlockMetadata(std::string title);

  // Appends one key: value pair.
  void addEntry(std::string key, std::string value);

  // Renders as "key : value" lines, with keys left-justified and values
  // right-justified to widths computed from content.
  std::string txt() const override;

  // Renders as plain "key,value" rows, unaligned.
  std::string csv() const override;

private:
  struct Entry {
    std::string key;
    std::string value;
  };

  std::vector<Entry> entries_;
};

#endif
