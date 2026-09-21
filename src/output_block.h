// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#ifndef OUTPUT_BLOCK_H
#define OUTPUT_BLOCK_H

#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <Eigen/Core>

// A block of output. Holds the title and builds the "[title]" + blank
// line preamble that every render_txt() starts with; the derived classes
// own the data and the body. The title is optional: an empty one renders
// no preamble at all, so the block is just its body -- for a unit nested
// under a heading something else already wrote.
class OutputUnit {
public:
  explicit OutputUnit(std::string title = "") : title_(std::move(title)) {}

  void set_title(std::string title) { title_ = std::move(title); }
  const std::string& title() const { return title_; }

protected:
  std::string header() const { return title_.empty() ? "" : "[" + title_ + "]\n\n"; }

  std::string title_;
};

// Key/value block, for metadata. Values are strings by the time they land
// here -- there's no formatting to defer, unlike a table.
class KeyValueOutput : public OutputUnit {
public:
  using OutputUnit::OutputUnit;

  void add(std::string key, std::string value);
  void add(std::string key, int value);
  void add(int key, std::string value);
  void add(int key, int value);

  // "key  ....  value" per entry, in insertion order, with the dot runs
  // sized so every value starts in the same column. `tabs` indents the
  // whole block by that many four-space runs.
  std::string render_txt(int tabs = 0) const;

private:
  std::vector<std::pair<std::string, std::string>> entries_;
};

// Storage shared by the two table orientations: named series of doubles,
// kept numeric until render time. Series need not all be the same length;
// a short one simply runs out, and render time leaves the cells past its
// end blank.
class Table : public OutputUnit {
public:
  using OutputUnit::OutputUnit;

protected:
  // Appends a named series.
  void add_series(std::string name, std::vector<double> values);

  // Same, from an Eigen vector. Ref<const> also binds to expressions that
  // aren't contiguous VectorXd -- a matrix row, a block, a coefficient-wise
  // expression -- by evaluating them into a temporary first.
  void add_series(std::string name, const Eigen::Ref<const Eigen::VectorXd>& values);

  std::vector<std::string> names_;
  std::vector<std::vector<double>> series_;
};

// Table with its names across the top and one row per index:
//
//   x           scalar_flux  sigma_total
//   1.0000e+00  1.2345e+12   1.5000e+12
//
// A column shorter than the longest one runs out partway down; its cells
// below that are left blank.
class VerticalTable : public Table {
public:
  using Table::Table;

  void add_column(std::string name, std::vector<double> values) {
    add_series(std::move(name), std::move(values));
  }

  void add_column(std::string name, const Eigen::Ref<const Eigen::VectorXd>& values) {
    add_series(std::move(name), values);
  }

  // `format` is the std::format spec applied to every value; `tabs`
  // indents the whole block by that many four-space runs.
  std::string render_txt(std::string_view format = "{:.4e}", int tabs = 0) const;
};

// The same table rotated: names down the left, one row per series. Rows
// may differ in length; a short row just ends early.
class HorizontalTable : public Table {
public:
  using Table::Table;

  void add_row(std::string name, std::vector<double> values) {
    add_series(std::move(name), std::move(values));
  }

  void add_row(std::string name, const Eigen::Ref<const Eigen::VectorXd>& values) {
    add_series(std::move(name), values);
  }

  // `format` is the std::format spec applied to every value; `tabs`
  // indents the whole block by that many four-space runs.
  std::string render_txt(std::string_view format = "{:.4e}", int tabs = 0) const;
};

#endif
