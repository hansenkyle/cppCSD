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
#include <sstream>
#include <stdexcept>
#include <utility>

namespace {

constexpr std::size_t kColumnPadding = 2;

std::string formatValue(double value, const OutputTable::Format& format) {
  std::ostringstream oss;
  switch (format.notation) {
  case OutputTable::Notation::Integer:
    oss << std::llround(value);
    break;
  case OutputTable::Notation::Fixed:
    oss << std::fixed << std::setprecision(format.precision) << value;
    break;
  case OutputTable::Notation::Scientific:
    oss << std::scientific << std::setprecision(format.precision) << value;
    break;
  }
  return oss.str();
}

} // namespace

OutputTable::OutputTable(std::string title, int n_rows)
    : title_(std::move(title)), n_rows_(n_rows) {
  if (n_rows <= 0) {
    throw std::invalid_argument("OutputTable: n_rows must be positive");
  }
}

void OutputTable::addColumn(std::string name, Format format, std::vector<double> values) {
  if (static_cast<int>(values.size()) != n_rows_) {
    throw std::invalid_argument("OutputTable::addColumn: '" + name + "' has " +
                                std::to_string(values.size()) + " values, expected " +
                                std::to_string(n_rows_));
  }
  columns_.push_back(Column{std::move(name), format, std::move(values)});
}

std::string OutputTable::txt() const {
  const std::size_t num_columns = columns_.size();
  std::vector<std::vector<std::string>> formatted(num_columns);
  std::vector<std::size_t> widths(num_columns);

  for (std::size_t c = 0; c < num_columns; ++c) {
    const Column& column = columns_[c];
    formatted[c].reserve(column.values.size());
    std::size_t width = column.name.size();
    for (double value : column.values) {
      std::string text = formatValue(value, column.format);
      width = std::max(width, text.size());
      formatted[c].push_back(std::move(text));
    }
    widths[c] = width + kColumnPadding;
  }

  std::ostringstream out;
  out << "--- " << title_ << " ---\n";
  for (std::size_t c = 0; c < num_columns; ++c) {
    out << std::setw(static_cast<int>(widths[c])) << columns_[c].name;
  }
  out << "\n";

  for (int r = 0; r < n_rows_; ++r) {
    for (std::size_t c = 0; c < num_columns; ++c) {
      out << std::setw(static_cast<int>(widths[c])) << formatted[c][r];
    }
    out << "\n";
  }
  return out.str();
}

std::string OutputTable::csv() const {
  std::ostringstream out;
  out << title_ << "\n";

  const std::size_t num_columns = columns_.size();
  for (std::size_t c = 0; c < num_columns; ++c) {
    if (c > 0) {
      out << ",";
    }
    out << columns_[c].name;
  }
  out << "\n";

  for (int r = 0; r < n_rows_; ++r) {
    for (std::size_t c = 0; c < num_columns; ++c) {
      if (c > 0) {
        out << ",";
      }
      out << formatValue(columns_[c].values[r], columns_[c].format);
    }
    out << "\n";
  }
  return out.str();
}

OutputMetadata::OutputMetadata(std::string title) : title_(std::move(title)) {}

void OutputMetadata::addEntry(std::string key, std::string value) {
  entries_.push_back(Entry{std::move(key), std::move(value)});
}

std::string OutputMetadata::txt() const {
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

std::string OutputMetadata::csv() const {
  std::ostringstream out;
  out << title_ << "\n";
  for (const Entry& entry : entries_) {
    out << entry.key << "," << entry.value << "\n";
  }
  return out.str();
}
