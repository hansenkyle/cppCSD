// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "output_block.h"

#include <algorithm>
#include <format>
#include <utility>

namespace {

constexpr std::size_t kColumnGap = 2;
constexpr std::size_t kMinDots = 3;
constexpr std::size_t kTabWidth = 4;

// Prefixes every line of `text` with `tabs` four-space runs. Blank lines
// are left as-is, so indenting a block never gives it a line made only of
// whitespace.
std::string indent(const std::string& text, int tabs) {
  if (tabs <= 0) {
    return text;
  }

  const std::string pad(kTabWidth * static_cast<std::size_t>(tabs), ' ');
  std::string out;
  for (std::size_t start = 0; start < text.size();) {
    const std::size_t newline = text.find('\n', start);
    const std::size_t stop = (newline == std::string::npos) ? text.size() : newline + 1;
    if (text[start] != '\n') {
      out += pad;
    }
    out.append(text, start, stop - start);
    start = stop;
  }
  return out;
}

// Renders a grid of already-stringified cells as whitespace-aligned
// columns. The first `n_left` columns are row labels and are
// left-justified; the rest hold values and are right-justified, so digits
// line up under their header regardless of sign or width. Rows may be of
// different lengths -- each column is sized from whichever rows reach it.
// Nothing is padded past the last cell on a line, so no line carries
// trailing whitespace.
std::string render_grid(const std::vector<std::vector<std::string>>& rows, std::size_t n_left) {
  if (rows.empty()) {
    return "";
  }

  std::vector<std::size_t> widths;
  for (const std::vector<std::string>& row : rows) {
    widths.resize(std::max(widths.size(), row.size()), 0);
    for (std::size_t c = 0; c < row.size(); ++c) {
      widths[c] = std::max(widths[c], row[c].size());
    }
  }

  std::string out;
  for (const std::vector<std::string>& row : rows) {
    for (std::size_t c = 0; c < row.size(); ++c) {
      const std::string& cell = row[c];
      const std::string fill(widths[c] - cell.size(), ' ');
      const bool last = c + 1 == row.size();
      if (c < n_left) {
        out += last ? cell : cell + fill;
      } else {
        out += fill + cell;
      }
      if (!last) {
        out += std::string(kColumnGap, ' ');
      }
    }
    out += '\n';
  }
  return out;
}

// One cell as it displays: a numeric cell gets `format` applied now, a
// text cell is already what it renders as.
std::string render_cell(const Table::Cell& cell, std::string_view format) {
  if (const std::string* text = std::get_if<std::string>(&cell)) {
    return *text;
  }

  double value = std::get<double>(cell);
  return std::vformat(format, std::make_format_args(value));
}

} // namespace

void KeyValueOutput::add(std::string key, std::string value) {
  entries_.emplace_back(std::move(key), std::move(value));
}

void KeyValueOutput::add(std::string key, int value) { add(std::move(key), std::to_string(value)); }

void KeyValueOutput::add(int key, std::string value) { add(std::to_string(key), std::move(value)); }

void KeyValueOutput::add(int key, int value) { add(std::to_string(key), std::to_string(value)); }

std::string KeyValueOutput::render_txt(int tabs) const {
  std::size_t key_width = 0;
  for (const auto& [key, value] : entries_) {
    key_width = std::max(key_width, key.size());
  }

  std::string out = header();
  for (const auto& [key, value] : entries_) {
    const std::size_t dots = key_width - key.size() + kMinDots;
    out += key + std::string(kColumnGap, ' ') + std::string(dots, '.') +
           std::string(kColumnGap, ' ') + value + '\n';
  }
  return indent(out, tabs);
}

void Table::add_series(std::string name, std::vector<Cell> values) {
  names_.push_back(std::move(name));
  series_.push_back(std::move(values));
}

void Table::add_series(std::string name, const Eigen::Ref<const Eigen::VectorXd>& values) {
  add_series(std::move(name), std::vector<Cell>(values.begin(), values.end()));
}

std::string VerticalTable::render_txt(std::string_view format, int tabs) const {
  if (series_.empty()) {
    return indent(header(), tabs);
  }

  std::size_t n_rows = 0;
  for (const std::vector<Cell>& column : series_) {
    n_rows = std::max(n_rows, column.size());
  }

  std::vector<std::vector<std::string>> rows;
  rows.push_back(names_);
  for (std::size_t i = 0; i < n_rows; ++i) {
    std::vector<std::string> row;
    std::size_t filled = 0;
    for (const std::vector<Cell>& column : series_) {
      if (i < column.size()) {
        row.push_back(render_cell(column[i], format));
        filled = row.size();
      } else {
        row.emplace_back();
      }
    }
    // A column that has run out contributes a blank cell, but trailing
    // blanks would render as trailing whitespace -- drop them instead.
    row.resize(filled);
    rows.push_back(std::move(row));
  }
  return indent(header() + render_grid(rows, 0), tabs);
}

std::string HorizontalTable::render_txt(std::string_view format, int tabs) const {
  std::vector<std::vector<std::string>> rows;
  for (std::size_t i = 0; i < series_.size(); ++i) {
    std::vector<std::string> row{names_[i]};
    for (const Cell& cell : series_[i]) {
      row.push_back(render_cell(cell, format));
    }
    rows.push_back(std::move(row));
  }
  return indent(header() + render_grid(rows, 1), tabs);
}

void MatrixTable::set_data(const std::vector<std::vector<Cell>>& data,
                           std::vector<std::string> column_labels,
                           std::vector<std::string> row_labels) {
  names_.clear();
  series_.clear();
  column_levels_.clear();
  row_levels_.clear();
  if (!column_labels.empty()) {
    column_levels_.push_back(std::move(column_labels));
  }
  if (!row_labels.empty()) {
    row_levels_.push_back(std::move(row_labels));
  }

  for (const std::vector<Cell>& row : data) {
    add_series(std::string(), row);
  }
}

void MatrixTable::set_data(const Eigen::Ref<const Eigen::MatrixXd>& data,
                           std::vector<std::string> column_labels,
                           std::vector<std::string> row_labels) {
  std::vector<std::vector<Cell>> rows(static_cast<std::size_t>(data.rows()));
  for (Eigen::Index r = 0; r < data.rows(); ++r) {
    std::vector<Cell>& row = rows[static_cast<std::size_t>(r)];
    row.reserve(static_cast<std::size_t>(data.cols()));
    for (Eigen::Index c = 0; c < data.cols(); ++c) {
      row.emplace_back(data(r, c));
    }
  }
  set_data(rows, std::move(column_labels), std::move(row_labels));
}

void MatrixTable::add_column_label(std::vector<std::string> labels) {
  column_levels_.push_back(std::move(labels));
}

void MatrixTable::add_row_label(std::vector<std::string> labels) {
  row_levels_.push_back(std::move(labels));
}

std::string MatrixTable::render_txt(std::string_view format, int tabs) const {
  if (series_.empty()) {
    return indent(header(), tabs);
  }

  const std::size_t n_label_cols = row_levels_.size();
  std::size_t n_head_rows = column_levels_.size();
  // Unlabelled columns leave no header row for the corner text to sit in.
  // Give it one, since the row labels under it still want naming.
  if (n_head_rows == 0 && n_label_cols > 0 && !corner_.empty()) {
    n_head_rows = 1;
  }

  std::vector<std::vector<std::string>> rows;

  for (std::size_t level = 0; level < n_head_rows; ++level) {
    std::vector<std::string> head;
    // The corner text goes in the last header row, so it sits directly
    // above the row labels no matter how many column levels there are.
    const bool last_head_row = level + 1 == n_head_rows;
    for (std::size_t col = 0; col < n_label_cols; ++col) {
      head.push_back(col == 0 && last_head_row ? corner_ : std::string());
    }
    if (level < column_levels_.size()) {
      const std::vector<std::string>& names = column_levels_[level];
      head.insert(head.end(), names.begin(), names.end());
    }
    while (!head.empty() && head.back().empty()) {
      head.pop_back();
    }
    // A header row with nothing left in it -- a level that names no
    // columns, or no column labels at all -- would render as a line of
    // nothing, so drop it rather than opening the table with an empty
    // line. Dropping the trailing blanks first also keeps a level that
    // names only the leftmost columns from trailing whitespace across the
    // rest.
    if (!head.empty()) {
      rows.push_back(std::move(head));
    }
  }

  for (std::size_t r = 0; r < series_.size(); ++r) {
    std::vector<std::string> row;
    for (std::size_t col = 0; col < n_label_cols; ++col) {
      const bool labelled = r < row_levels_[col].size();
      row.push_back(labelled ? row_levels_[col][r] : std::string());
    }
    for (const Cell& cell : series_[r]) {
      row.push_back(render_cell(cell, format));
    }
    rows.push_back(std::move(row));
  }

  return indent(header() + render_grid(rows, n_label_cols), tabs);
}

std::string UnitGroup::render_txt(int tabs) const {
  std::string out = header();
  for (std::size_t i = 0; i < bodies_.size(); ++i) {
    if (i > 0) {
      out += '\n';
    }
    out += indent(bodies_[i], 1);
  }
  return indent(out, tabs);
}
