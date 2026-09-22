// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#ifndef SOLVER_FORMATTER_H
#define SOLVER_FORMATTER_H

#include <filesystem>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "convergence.h"
#include "input_deck.h"

// Renders Solver's output blocks (run metadata, results, residuals) to
// text. Pure formatting: takes the data it needs as arguments, returns a
// string, and never touches the filesystem or knows when/how often it's
// called -- that's Solver's job. Kept separate so Solver stays responsible
// for *what* to write and *when* (it owns the run's data and decides the
// write schedule), not *how* it's laid out as text.
namespace SolverFormatter {

std::string formatInputEcho(const InputDeck& deck);

// "Scalar Flux" block, followed by one "Angular Flux - Group N" block per
// entry in angular_flux (indexed like InputDeck::Source: one entry per
// group, one ordinate per column).
std::string formatResults(const Eigen::MatrixXd& scalar_flux,
                          const std::vector<Eigen::MatrixXd>& angular_flux, const InputDeck& deck);

// One "Angular Residual - Group N" block per entry in residuals (indexed
// like InputDeck::Source).
std::string formatResiduals(const std::vector<Eigen::MatrixXd>& residuals, const InputDeck& deck);

// "Convergence Summary" (one row per energy group) followed by
// "Iteration History" (one row per iteration of every group). Returns a
// short note instead of tables if the history is empty.
std::string formatConvergence(const ConvergenceHistory& history);

} // namespace SolverFormatter

#endif
