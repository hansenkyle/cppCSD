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

OutputTable::OutputTable(std::string title, int n_entries)
    : title_(std::move(title)), n_entries_(n_entries) {
  if (n_entries <= 0) {
    throw std::invalid_argument("OutputTable: n_entries must be positive");
  }
}

void OutputTable::addColumn(std::string name, Format format, std::vector<double> values) {
  if (orientation_ == Orientation::RowMajor) {
    throw std::invalid_argument("OutputTable::addColumn: table already uses addRow (row-major); "
                                "the two can't be mixed on one table");
  }
  orientation_ = Orientation::ColumnMajor;
  if (static_cast<int>(values.size()) != n_entries_) {
    throw std::invalid_argument("OutputTable::addColumn: '" + name + "' has " +
                                std::to_string(values.size()) + " values, expected " +
                                std::to_string(n_entries_));
  }
  entries_.push_back(Entry{std::move(name), format, std::move(values)});
}

void OutputTable::addRow(std::string name, Format format, std::vector<double> values) {
  if (orientation_ == Orientation::ColumnMajor) {
    throw std::invalid_argument("OutputTable::addRow: table already uses addColumn (column-major); "
                                "the two can't be mixed on one table");
  }
  orientation_ = Orientation::RowMajor;
  if (static_cast<int>(values.size()) != n_entries_) {
    throw std::invalid_argument("OutputTable::addRow: '" + name + "' has " +
                                std::to_string(values.size()) + " values, expected " +
                                std::to_string(n_entries_));
  }
  entries_.push_back(Entry{std::move(name), format, std::move(values)});
}

std::string OutputTable::txtColumnMajor() const {
  const std::size_t num_columns = entries_.size();
  std::vector<std::vector<std::string>> formatted(num_columns);
  std::vector<std::size_t> widths(num_columns);

  for (std::size_t c = 0; c < num_columns; ++c) {
    const Entry& column = entries_[c];
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
    out << std::setw(static_cast<int>(widths[c])) << entries_[c].name;
  }
  out << "\n";

  for (int r = 0; r < n_entries_; ++r) {
    for (std::size_t c = 0; c < num_columns; ++c) {
      out << std::setw(static_cast<int>(widths[c])) << formatted[c][r];
    }
    out << "\n";
  }
  return out.str();
}

std::string OutputTable::txtRowMajor() const {
  std::size_t label_width = 0;
  for (const Entry& row : entries_) {
    label_width = std::max(label_width, row.name.size());
  }
  label_width += kColumnPadding;

  std::vector<std::vector<std::string>> formatted(entries_.size());
  std::vector<std::size_t> widths(n_entries_);
  for (int c = 0; c < n_entries_; ++c) {
    widths[c] = std::to_string(c).size();
  }
  for (std::size_t r = 0; r < entries_.size(); ++r) {
    const Entry& row = entries_[r];
    formatted[r].reserve(n_entries_);
    for (int c = 0; c < n_entries_; ++c) {
      std::string text = formatValue(row.values[c], row.format);
      widths[c] = std::max(widths[c], text.size());
      formatted[r].push_back(std::move(text));
    }
  }
  for (std::size_t c = 0; c < widths.size(); ++c) {
    widths[c] += kColumnPadding;
  }

  std::ostringstream out;
  out << "--- " << title_ << " ---\n";
  out << std::setw(static_cast<int>(label_width)) << "";
  for (int c = 0; c < n_entries_; ++c) {
    out << std::setw(static_cast<int>(widths[c])) << c;
  }
  out << "\n";

  for (std::size_t r = 0; r < entries_.size(); ++r) {
    out << std::left << std::setw(static_cast<int>(label_width)) << entries_[r].name << std::right;
    for (int c = 0; c < n_entries_; ++c) {
      out << std::setw(static_cast<int>(widths[c])) << formatted[r][c];
    }
    out << "\n";
  }
  return out.str();
}

std::string OutputTable::txt() const {
  if (orientation_ == Orientation::RowMajor) {
    return txtRowMajor();
  }
  return txtColumnMajor();
}

std::string OutputTable::csv() const {
  std::ostringstream out;
  out << title_ << "\n";

  if (orientation_ == Orientation::RowMajor) {
    for (int c = 0; c < n_entries_; ++c) {
      out << "," << c;
    }
    out << "\n";
    for (const Entry& row : entries_) {
      out << row.name;
      for (int c = 0; c < n_entries_; ++c) {
        out << "," << formatValue(row.values[c], row.format);
      }
      out << "\n";
    }
    return out.str();
  }

  const std::size_t num_columns = entries_.size();
  for (std::size_t c = 0; c < num_columns; ++c) {
    if (c > 0) {
      out << ",";
    }
    out << entries_[c].name;
  }
  out << "\n";

  for (int r = 0; r < n_entries_; ++r) {
    for (std::size_t c = 0; c < num_columns; ++c) {
      if (c > 0) {
        out << ",";
      }
      out << formatValue(entries_[c].values[r], entries_[c].format);
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
