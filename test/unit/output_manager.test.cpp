// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "output_manager.h"

#include <fstream>

#include <doctest.h>

namespace {

// Creates a fresh, empty temp directory to act as a deck directory and
// returns its path; each call gets a distinct directory so tests don't
// interfere with each other's run counters.
std::filesystem::path makeTempDeckDir(const std::string& name) {
  const std::filesystem::path dir =
      std::filesystem::temp_directory_path() / "ldcsd_output_manager_test" / name;
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);
  return dir;
}

} // namespace

TEST_SUITE("OutputManager") {
  TEST_CASE("creates run_0001 and its files' paths on a fresh deck directory") {
    const std::filesystem::path deck_dir = makeTempDeckDir("fresh");

    const OutputManager manager(deck_dir);

    CHECK(manager.run_dir == deck_dir / "runs" / "run_0001");
    CHECK(manager.info_path == manager.run_dir / "info");
    CHECK(manager.log_path == manager.run_dir / "log");
    CHECK(manager.out_path == manager.run_dir / "out");
    CHECK(std::filesystem::is_directory(manager.run_dir));
  }

  TEST_CASE("picks the next unused monotonic run ID") {
    const std::filesystem::path deck_dir = makeTempDeckDir("increment");

    const OutputManager first(deck_dir);
    const OutputManager second(deck_dir);
    const OutputManager third(deck_dir);

    CHECK(first.run_dir.filename() == "run_0001");
    CHECK(second.run_dir.filename() == "run_0002");
    CHECK(third.run_dir.filename() == "run_0003");
  }

  TEST_CASE("resumes numbering above the highest existing run directory") {
    const std::filesystem::path deck_dir = makeTempDeckDir("resume");
    std::filesystem::create_directories(deck_dir / "runs" / "run_0005");
    std::filesystem::create_directories(deck_dir / "runs" / "run_0002");

    const OutputManager manager(deck_dir);

    CHECK(manager.run_dir.filename() == "run_0006");
  }

  TEST_CASE("ignores non-conforming entries under runs/") {
    const std::filesystem::path deck_dir = makeTempDeckDir("ignore-junk");
    std::filesystem::create_directories(deck_dir / "runs" / "notes");
    std::ofstream(deck_dir / "runs" / "run_0003") << "not a directory";

    const OutputManager manager(deck_dir);

    CHECK(manager.run_dir.filename() == "run_0001");
  }

  TEST_CASE("points runs/latest at the newest run directory") {
    const std::filesystem::path deck_dir = makeTempDeckDir("latest");

    const OutputManager first(deck_dir);
    const OutputManager second(deck_dir);

    const std::filesystem::path latest = deck_dir / "runs" / "latest";
    REQUIRE(std::filesystem::is_symlink(latest));
    CHECK(std::filesystem::read_symlink(latest) == second.run_dir.filename());
    CHECK_FALSE(std::filesystem::equivalent(latest, first.run_dir));
    CHECK(std::filesystem::equivalent(latest, second.run_dir));
  }
}
