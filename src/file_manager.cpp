// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "file_manager.h"

#include <algorithm>
#include <charconv>
#include <iomanip>
#include <optional>
#include <sstream>
#include <string>

#include "logger.h"

namespace fs = std::filesystem;

namespace {

constexpr char kRunPrefix[] = "run_";
constexpr int kRunIdWidth = 4;

// Parses a run directory name of the form "run_NNNN", returning the numeric
// ID, or nullopt if the name doesn't match that pattern.
std::optional<int> parseRunId(const std::string& dirname) {
  if (dirname.rfind(kRunPrefix, 0) != 0) {
    return std::nullopt;
  }
  const std::string suffix = dirname.substr(std::string(kRunPrefix).size());
  if (suffix.empty()) {
    return std::nullopt;
  }

  int id = 0;
  const auto* begin = suffix.data();
  const auto* end = suffix.data() + suffix.size();
  const auto [ptr, ec] = std::from_chars(begin, end, id);
  if (ec != std::errc() || ptr != end) {
    return std::nullopt;
  }
  return id;
}

std::string formatRunId(int id) {
  std::ostringstream oss;
  oss << kRunPrefix << std::setfill('0') << std::setw(kRunIdWidth) << id;
  return oss.str();
}

// Scans runs_dir for existing run_NNNN entries and returns one past the
// largest ID found (or 1 if none exist / runs_dir doesn't exist yet).
int nextRunId(const fs::path& runs_dir) {
  int max_id = 0;
  if (fs::exists(runs_dir)) {
    for (const fs::directory_entry& entry : fs::directory_iterator(runs_dir)) {
      if (!entry.is_directory()) {
        continue;
      }
      if (const std::optional<int> id = parseRunId(entry.path().filename().string())) {
        max_id = std::max(max_id, *id);
      }
    }
  } else {
    LDCSD_LOG_DEBUG("runs directory '" + runs_dir.string() + "' does not exist yet");
  }
  LDCSD_LOG_DEBUG("next run id will be " + std::to_string(max_id + 1));
  return max_id + 1;
}

// (Re)points runs_dir/latest at run_dir, replacing any existing symlink.
void updateLatestSymlink(const fs::path& runs_dir, const fs::path& run_dir) {
  const fs::path latest = runs_dir / "latest";
  if (fs::is_symlink(latest) || fs::exists(latest)) {
    LDCSD_LOG_DEBUG("replacing existing 'latest' symlink at '" + latest.string() + "'");
    fs::remove(latest);
  }
  fs::create_directory_symlink(run_dir.filename(), latest);
}

} // namespace

OutputManager::OutputManager(const fs::path& deck_dir) {
  try {
    const fs::path runs_dir = deck_dir / "runs";
    fs::create_directories(runs_dir);

    run_dir = runs_dir / formatRunId(nextRunId(runs_dir));
    fs::create_directory(run_dir);

    updateLatestSymlink(runs_dir, run_dir);
  } catch (const fs::filesystem_error& e) {
    LDCSD_LOG_ERROR(std::string("failed to set up output directory under '") + deck_dir.string() +
                    "': " + e.what());
    throw;
  }

  LDCSD_LOG_INFO("created run directory '" + run_dir.string() + "'");

  info_path = run_dir / "info.txt";
  log_path = run_dir / "log.txt";
  out_path = run_dir / "out.txt";
  results_path = run_dir / "results.txt";
  residuals_path = run_dir / "residuals.txt";
}
