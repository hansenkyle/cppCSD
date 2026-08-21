#ifndef OUTPUT_MANAGER_H
#define OUTPUT_MANAGER_H

#include <filesystem>

// Owns the on-disk layout for a single solver run's output, per
// docs/output-layout.md: <deck_dir>/runs/run_NNNN/{info,log,out}, with
// run_NNNN chosen automatically as the next unused monotonic counter and
// <deck_dir>/runs/latest kept pointing at it. Agnostic of what ends up
// written into those files -- this is path bookkeeping only.
class OutputManager {
public:
  // Determines the next run ID under <deck_dir>/runs, creates that run
  // directory (and runs/ itself, if needed), and repoints runs/latest at it.
  explicit OutputManager(const std::filesystem::path &deck_dir);

  std::filesystem::path run_dir;
  std::filesystem::path info_path;
  std::filesystem::path log_path;
  std::filesystem::path out_path;
};

#endif
