// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#ifndef OUTPUT_MANAGER_H
#define OUTPUT_MANAGER_H

#include <filesystem>

/// @class OutputManager
/// @brief Contains the locations of all output files, per @ref docs/output-layout.md
/// @details Bookkeeping only. OutputManager contains no functionality for modifying or reading the
/// files (yet). Automatically assigns run ID based on time and date
class OutputManager {
public:
  // Determines the next run ID under <deck_dir>/runs, creates that run
  // directory (and runs/ itself, if needed), and repoints runs/latest at it.
  explicit OutputManager(const std::filesystem::path& deck_dir);

  std::filesystem::path run_dir;
  std::filesystem::path info_path;
  std::filesystem::path log_path;
  std::filesystem::path out_path;
  std::filesystem::path results_path;
  std::filesystem::path residuals_path;
};

#endif
