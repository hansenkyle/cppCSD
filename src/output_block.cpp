// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "output_block.h"

#include <algorithm>
#include <format>
#include <stdexcept>

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

// Renders a rectangular grid of already-stringified cells as whitespace-
// aligned columns. The first `n_left` columns are row labels and are
// left-justified; the rest hold values and are right-justified, so digits
// line up under their header regardless of sign or width. Nothing is
// padded past the last cell on a line, so no line carries trailing
// whitespace.
std::string render_grid(const std::vector<std::vector<std::string>>& rows, std::size_t n_left) {
  if (rows.empty()) {
    return "";
  }

  std::vector<std::size_t> widths(rows.front().size(), 0);
  for (const std::vector<std::string>& row : rows) {
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

void Table::add_series(std::string name, std::vector<double> values) {
  if (!series_.empty() && values.size() != series_.front().size()) {
    throw std::invalid_argument("Table: '" + name + "' has " + std::to_string(values.size()) +
                                " values, expected " + std::to_string(series_.front().size()));
  }
  names_.push_back(std::move(name));
  series_.push_back(std::move(values));
}

std::string VerticalTable::render_txt(std::string_view format, int tabs) const {
  if (series_.empty()) {
    return indent(header(), tabs);
  }

  std::vector<std::vector<std::string>> rows;
  rows.push_back(names_);
  for (std::size_t i = 0; i < series_.front().size(); ++i) {
    std::vector<std::string> row;
    for (const std::vector<double>& column : series_) {
      double value = column[i];
      row.push_back(std::vformat(format, std::make_format_args(value)));
    }
    rows.push_back(std::move(row));
  }
  return indent(header() + render_grid(rows, 0), tabs);
}

std::string HorizontalTable::render_txt(std::string_view format, int tabs) const {
  std::vector<std::vector<std::string>> rows;
  for (std::size_t i = 0; i < series_.size(); ++i) {
    std::vector<std::string> row{names_[i]};
    for (double value : series_[i]) {
      row.push_back(std::vformat(format, std::make_format_args(value)));
    }
    rows.push_back(std::move(row));
  }
  return indent(header() + render_grid(rows, 1), tabs);
}
