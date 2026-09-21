// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include <optional>

#include "cli.h"
#include "file_manager.h"
#include "input_deck.h"
#include "logger.h"
#include "solver.h"
#include "terminal.h"

int main(int argc, char** argv) {
  int exit_code = 0;

  const std::optional<std::filesystem::path> yaml_path = parseArgs(argc, argv, exit_code);
  if (!yaml_path.has_value()) {
    return exit_code;
  }

  const std::filesystem::path deck_dir = std::filesystem::absolute(*yaml_path).parent_path();
  const OutputManager output(deck_dir);
  Logger::configure(output.log_path);
  Terminal::configure(output.out_path);

  LDCSD_LOG(Channel::General, Level::Info, "General/info message");
  LDCSD_LOG_INFO("ldcsd starting, input deck: " + yaml_path->string());

  InputDeck deck;
  if (deck.read(*yaml_path) != 0) {
    // read() has already logged the specific failure.
    return 1;
  }

  Solver solver(deck);
  LDCSD_LOG_INFO("constructed Solver");

  solver.writeMetadata(output.info_path);
  solver.writeInputDeckEcho(output.info_path);

  const Eigen::MatrixXd phi = solver.sourceIterate(1e-8);

  solver.writeConvergence(output.results_path);
  solver.writeResults(output.results_path, phi, {});

  return solver.convergence().allConverged() ? 0 : 1;
}
