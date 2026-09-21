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
#include <type_traits>
#include <utility>
#include <variant>
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

// Storage shared by the two table orientations: named series of cells,
// kept in insertion order. Series need not all be the same length; a short
// one simply runs out, and render time leaves the cells past its end
// blank.
class Table : public OutputUnit {
public:
  using OutputUnit::OutputUnit;

  // A cell is either numeric -- held as a double so a format spec can be
  // applied at render time -- or literal text. Strings and integers have
  // nothing left to format, so they become text when they're added. Cells
  // are per-value, not per-series, so one series can mix the two.
  using Cell = std::variant<double, std::string>;

protected:
  // Appends a named series. A braced list resolves here, so mixed literals
  // -- {1.0, "n/a", 2.5} -- work directly.
  void add_series(std::string name, std::vector<Cell> values);

  // Same, from a vector of anything a Cell holds. Integers are stringified
  // now; everything else converts to a Cell as-is.
  template <typename T> void add_series(std::string name, const std::vector<T>& values) {
    std::vector<Cell> cells;
    cells.reserve(values.size());
    for (const T& value : values) {
      if constexpr (std::is_integral_v<T>) {
        cells.emplace_back(std::to_string(value));
      } else {
        cells.emplace_back(value);
      }
    }
    add_series(std::move(name), std::move(cells));
  }

  // Same, from an Eigen vector. Ref<const> also binds to expressions that
  // aren't contiguous VectorXd -- a matrix row, a block, a coefficient-wise
  // expression -- by evaluating them into a temporary first.
  void add_series(std::string name, const Eigen::Ref<const Eigen::VectorXd>& values);

  std::vector<std::string> names_;
  std::vector<std::vector<Cell>> series_;
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

  void add_column(std::string name, std::vector<Cell> values) {
    add_series(std::move(name), std::move(values));
  }

  template <typename T> void add_column(std::string name, const std::vector<T>& values) {
    add_series(std::move(name), values);
  }

  void add_column(std::string name, const Eigen::Ref<const Eigen::VectorXd>& values) {
    add_series(std::move(name), values);
  }

  // `format` is the std::format spec applied to every numeric cell; text
  // cells are already rendered. `tabs` indents the whole block by that
  // many four-space runs.
  std::string render_txt(std::string_view format = "{:.4e}", int tabs = 0) const;
};

// The same table rotated: names down the left, one row per series. Rows
// may differ in length; a short row just ends early.
class HorizontalTable : public Table {
public:
  using Table::Table;

  void add_row(std::string name, std::vector<Cell> values) {
    add_series(std::move(name), std::move(values));
  }

  template <typename T> void add_row(std::string name, const std::vector<T>& values) {
    add_series(std::move(name), values);
  }

  void add_row(std::string name, const Eigen::Ref<const Eigen::VectorXd>& values) {
    add_series(std::move(name), values);
  }

  // `format` is the std::format spec applied to every numeric cell; text
  // cells are already rendered. `tabs` indents the whole block by that
  // many four-space runs.
  std::string render_txt(std::string_view format = "{:.4e}", int tabs = 0) const;
};

// One or more rendered units stacked into a single block, separated by
// blank lines and indented one level under the group's own title. A group
// is itself a unit, so groups nest -- each layer indents its contents one
// step further than the layer above.
class UnitGroup : public OutputUnit {
public:
  using OutputUnit::OutputUnit;

  // Renders `unit` now and keeps the result. Extra arguments are passed
  // through to the unit's own render_txt(), so a table can carry its
  // format spec (e.g. add(table, "{:.2e}")). Anything with a render_txt()
  // works, another UnitGroup included.
  template <typename Unit, typename... Args> void add(const Unit& unit, Args&&... args) {
    bodies_.push_back(unit.render_txt(std::forward<Args>(args)...));
  }

  // `tabs` indents the whole group by that many four-space runs, on top of
  // the one level its contents already sit at.
  std::string render_txt(int tabs = 0) const;

private:
  std::vector<std::string> bodies_;
};

#endif
