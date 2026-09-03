// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "cli.h"

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include <doctest.h>

namespace {

// Builds a char** argv view over storage, whose elements must outlive the
// parseArgs() call using the returned pointers.
std::vector<char*> makeArgv(std::vector<std::string>& storage) {
  std::vector<char*> argv;
  argv.reserve(storage.size());
  for (std::string& arg : storage) {
    argv.push_back(arg.data());
  }
  return argv;
}

} // namespace

TEST_SUITE("parseArgs") {
  TEST_CASE("parses a single positional argument as the yaml path") {
    std::vector<std::string> storage = {"ldcsd", "deck.yaml"};
    std::vector<char*> argv = makeArgv(storage);

    int exit_code = -1;
    const std::optional<std::filesystem::path> yaml_path =
        parseArgs(static_cast<int>(argv.size()), argv.data(), exit_code);

    REQUIRE(yaml_path.has_value());
    CHECK(*yaml_path == std::filesystem::path("deck.yaml"));
    CHECK(exit_code == -1); // left untouched on success
  }

  TEST_CASE("accepts a nested path") {
    std::vector<std::string> storage = {"ldcsd", "some/nested/deck.yaml"};
    std::vector<char*> argv = makeArgv(storage);

    int exit_code = -1;
    const std::optional<std::filesystem::path> yaml_path =
        parseArgs(static_cast<int>(argv.size()), argv.data(), exit_code);

    REQUIRE(yaml_path.has_value());
    CHECK(*yaml_path == std::filesystem::path("some/nested/deck.yaml"));
  }

  TEST_CASE("returns nullopt and a nonzero exit code when the yaml path is missing") {
    std::vector<std::string> storage = {"ldcsd"};
    std::vector<char*> argv = makeArgv(storage);

    int exit_code = -1;
    const std::optional<std::filesystem::path> yaml_path =
        parseArgs(static_cast<int>(argv.size()), argv.data(), exit_code);

    CHECK_FALSE(yaml_path.has_value());
    CHECK(exit_code != 0);
  }

  TEST_CASE("returns nullopt and a nonzero exit code with extra positional arguments") {
    std::vector<std::string> storage = {"ldcsd", "deck.yaml", "extra.yaml"};
    std::vector<char*> argv = makeArgv(storage);

    int exit_code = -1;
    const std::optional<std::filesystem::path> yaml_path =
        parseArgs(static_cast<int>(argv.size()), argv.data(), exit_code);

    CHECK_FALSE(yaml_path.has_value());
    CHECK(exit_code != 0);
  }

  TEST_CASE("--help returns nullopt with exit code 0") {
    std::vector<std::string> storage = {"ldcsd", "--help"};
    std::vector<char*> argv = makeArgv(storage);

    int exit_code = -1;
    const std::optional<std::filesystem::path> yaml_path =
        parseArgs(static_cast<int>(argv.size()), argv.data(), exit_code);

    CHECK_FALSE(yaml_path.has_value());
    CHECK(exit_code == 0);
  }
}
