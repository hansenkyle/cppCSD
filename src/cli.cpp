// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "cli.h"

#include <CLI/CLI.hpp>

std::optional<std::filesystem::path> parseArgs(int argc, char** argv, int& exit_code) {
  CLI::App app{"ldcsd: solves the Continuous Slowing-Down equation via Lewis and Miller's Second "
               "Moment Method"};

  std::string yaml_path;
  app.add_option("yaml_path", yaml_path, "Path to the input deck YAML file")->required();

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    exit_code = app.exit(e);
    return std::nullopt;
  }

  return std::filesystem::path(yaml_path);
}
