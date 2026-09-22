// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "method.h"
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

void Method::appendToFile(const std::filesystem::path& file_path, const std::string& text) const {
  std::ofstream out(file_path, std::ios::app);
  if (!out.is_open()) {
    throw std::runtime_error("Solver: failed to open '" + file_path.string() + "' for writing");
  }
  out << text;
}