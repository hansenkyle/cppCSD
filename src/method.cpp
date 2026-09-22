// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "method.h"
#include "output_block.h"
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

namespace {
std::string timestamp() {
  const std::time_t now = std::time(nullptr);
  const std::tm* tm = std::localtime(&now);
  std::ostringstream oss;
  oss << std::put_time(tm, "%Y-%m-%d %H:%M:%S");
  return oss.str();
}
} // namespace

void Method::appendToFile(const std::filesystem::path& file_path, const std::string& text) const {
  std::ofstream out(file_path, std::ios::app);
  if (!out.is_open()) {
    throw std::runtime_error("Solver: failed to open '" + file_path.string() + "' for writing");
  }
  out << text;
}
void Method::writeMetadata(const std::filesystem::path& file_path) const {
  KeyValueOutput metadata("run info");
  metadata.add("execution date/time", timestamp());
  metadata.add("n_groups", input_deck.energy.G);
  metadata.add("n_cells", input_deck.mesh.n_x);
  metadata.add("n_angles", input_deck.angle.M);
  metadata.add("method", name);
  appendToFile(file_path, metadata.render_txt() + "\n");
}