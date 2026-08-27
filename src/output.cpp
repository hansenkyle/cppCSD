// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "output.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace {

constexpr std::size_t kColumnPadding = 2;

std::string formatValue(double value, const OutputBlock::Format& format) {
  std::ostringstream oss;
  switch (format.notation) {
  case OutputBlock::Notation::Integer:
    oss << std::llround(value);
    break;
  case OutputBlock::Notation::Fixed:
    oss << std::fixed << std::setprecision(format.precision) << value;
    break;
  case OutputBlock::Notation::Scientific:
    oss << std::scientific << std::setprecision(format.precision) << value;
    break;
  }
  return oss.str();
}

std::size_t numRows(const OutputBlock::SubTable& subtable) {
  return subtable.columns.empty() ? 0 : subtable.columns[0].values.size();
}

void appendTxtSubTable(std::ostringstream& out, const OutputBlock::SubTable& subtable) {
  const std::size_t num_columns = subtable.columns.size();
  std::vector<std::vector<std::string>> formatted(num_columns);
  std::vector<std::size_t> widths(num_columns);

  for (std::size_t c = 0; c < num_columns; ++c) {
    const OutputBlock::Column& column = subtable.columns[c];
    formatted[c].reserve(column.values.size());
    std::size_t width = column.name.size();
    for (double value : column.values) {
      std::string text = formatValue(value, column.format);
      width = std::max(width, text.size());
      formatted[c].push_back(std::move(text));
    }
    widths[c] = width + kColumnPadding;
  }

  for (std::size_t c = 0; c < num_columns; ++c) {
    out << std::setw(static_cast<int>(widths[c])) << subtable.columns[c].name;
  }
  out << "\n";

  const std::size_t num_rows = numRows(subtable);
  for (std::size_t r = 0; r < num_rows; ++r) {
    for (std::size_t c = 0; c < num_columns; ++c) {
      out << std::setw(static_cast<int>(widths[c])) << formatted[c][r];
    }
    out << "\n";
  }
}

void appendCsvSubTable(std::ostringstream& out, const OutputBlock::SubTable& subtable) {
  const std::size_t num_columns = subtable.columns.size();
  for (std::size_t c = 0; c < num_columns; ++c) {
    if (c > 0) {
      out << ",";
    }
    out << subtable.columns[c].name;
  }
  out << "\n";

  const std::size_t num_rows = numRows(subtable);
  for (std::size_t r = 0; r < num_rows; ++r) {
    for (std::size_t c = 0; c < num_columns; ++c) {
      if (c > 0) {
        out << ",";
      }
      out << formatValue(subtable.columns[c].values[r], subtable.columns[c].format);
    }
    out << "\n";
  }
}

} // namespace

OutputBlock::OutputBlock(std::string title) : title_(std::move(title)) {}

std::string OutputBlock::txt() const {
  std::ostringstream out;
  out << "--- " << title_ << " ---\n";
  for (std::size_t s = 0; s < subtables_.size(); ++s) {
    if (s > 0) {
      out << "\n";
    }
    const SubTable& subtable = subtables_[s];
    if (!subtable.subtitle.empty()) {
      out << subtable.subtitle << "\n";
    }
    appendTxtSubTable(out, subtable);
  }
  return out.str();
}

std::string OutputBlock::csv() const {
  std::ostringstream out;
  out << title_ << "\n";
  for (std::size_t s = 0; s < subtables_.size(); ++s) {
    if (s > 0) {
      out << "\n";
    }
    const SubTable& subtable = subtables_[s];
    if (!subtable.subtitle.empty()) {
      out << subtable.subtitle << "\n";
    }
    appendCsvSubTable(out, subtable);
  }
  return out.str();
}

OutputBlock1d::OutputBlock1d(std::string title, int n_rows)
    : OutputBlock(std::move(title)), n_rows_(n_rows) {
  if (n_rows <= 0) {
    throw std::invalid_argument("OutputBlock1d: n_rows must be positive");
  }
  subtables_.push_back(SubTable{"", {}});
}

void OutputBlock1d::addColumn(std::string name, Format format, std::vector<double> values) {
  if (static_cast<int>(values.size()) != n_rows_) {
    throw std::invalid_argument("OutputBlock1d::addColumn: '" + name + "' has " +
                                std::to_string(values.size()) + " values, expected " +
                                std::to_string(n_rows_));
  }
  subtables_[0].columns.push_back(Column{std::move(name), format, std::move(values)});
}

OutputBlock2d::OutputBlock2d(std::string title, int n_x, int G)
    : OutputBlock(std::move(title)), n_x_(n_x), G_(G) {
  if (n_x <= 0) {
    throw std::invalid_argument("OutputBlock2d: n_x must be positive");
  }
  if (G <= 0) {
    throw std::invalid_argument("OutputBlock2d: G must be positive");
  }

  for (int g = 0; g < G; ++g) {
    std::vector<double> index(n_x);
    std::iota(index.begin(), index.end(), 0.0);
    subtables_.push_back(SubTable{"g=" + std::to_string(g),
                                  {Column{"i", Format{Notation::Integer, 0}, std::move(index)}}});
  }
}

void OutputBlock2d::addColumn(std::string name, Format format,
                              const std::vector<std::vector<double>>& values) {
  if (static_cast<int>(values.size()) != G_) {
    throw std::invalid_argument("OutputBlock2d::addColumn: '" + name + "' has " +
                                std::to_string(values.size()) + " groups, expected " +
                                std::to_string(G_));
  }
  for (int g = 0; g < G_; ++g) {
    if (static_cast<int>(values[g].size()) != n_x_) {
      throw std::invalid_argument("OutputBlock2d::addColumn: '" + name + "' group " +
                                  std::to_string(g) + " has " + std::to_string(values[g].size()) +
                                  " values, expected " + std::to_string(n_x_));
    }
  }
  for (int g = 0; g < G_; ++g) {
    subtables_[g].columns.push_back(Column{name, format, values[g]});
  }
}

OutputBlockMetadata::OutputBlockMetadata(std::string title) : OutputBlock(std::move(title)) {}

void OutputBlockMetadata::addEntry(std::string key, std::string value) {
  entries_.push_back(Entry{std::move(key), std::move(value)});
}

std::string OutputBlockMetadata::txt() const {
  std::size_t key_width = 0;
  std::size_t value_width = 0;
  for (const Entry& entry : entries_) {
    key_width = std::max(key_width, entry.key.size());
    value_width = std::max(value_width, entry.value.size());
  }

  std::ostringstream out;
  out << "--- " << title_ << " ---\n";
  for (const Entry& entry : entries_) {
    out << std::left << std::setw(static_cast<int>(key_width)) << entry.key << " : " << std::right
        << std::setw(static_cast<int>(value_width)) << entry.value << "\n";
  }
  return out.str();
}

std::string OutputBlockMetadata::csv() const {
  std::ostringstream out;
  out << title_ << "\n";
  for (const Entry& entry : entries_) {
    out << entry.key << "," << entry.value << "\n";
  }
  return out.str();
}
