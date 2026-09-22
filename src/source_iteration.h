// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#ifndef SOLVER_H
#define SOLVER_H

#include <filesystem>
#include <vector>

#include <Eigen/Dense>

#include <boost/multiprecision/eigen.hpp>
#include <boost/multiprecision/float128.hpp>

#include "convergence.h"
#include "input_deck.h"
#include "method.h"
#include "transport_operator.h"

using HighPrecision = boost::multiprecision::float128;

class SourceIteration : public Method {
public:
  // constructor from input deck (copy)
  SourceIteration(InputDeck input_deck)
      : Method("source iteration", input_deck), transport_operator(input_deck) {}
  InputDeck input_deck;

  /// @brief Compute scalar flux for all space, all energy groups using Source Iteration.
  ///
  /// Every iteration of every group is recorded into convergence(), which
  /// survives the call and is what gets written to the results file.
  ///
  /// @param epsilon Convergence criterion. Transport iteration stops when |phi_old - phi_new|_2 >
  /// |phi_new|_2*epsilon
  /// @param max_iterations Per-group iteration cap. A group that hits it is marked unconverged
  /// (logged as a warning) and the solve moves on to the next group rather than spinning forever.
  /// @return Eigen::MatrixXd. Scalar flux in all energy groups. [4nx by G]
  void solve(double epsilon, int max_iterations = kDefaultMaxIterations);

  // When true, every iteration re-evaluates the largest residual in the
  // group term-by-term and logs the breakdown. Extremely verbose (it was
  // unconditionally on before this became a switch) -- for hunting a
  // discretization bug, not for normal runs.
  bool log_residual_terms = false;

  static constexpr int kDefaultMaxIterations = 1000;

  // Appends the "run info" block -- run time, deck path, problem
  // dimensions -- to the results file. deck_path is only recorded, never
  // read, so it can be whatever spelling the caller wants shown.
  void write_metadata(const std::filesystem::path& results_path,
                      const std::filesystem::path& deck_path) const;

  // Appends an echo of input_deck (InputDeck::echo()).
  void writeInputDeckEcho(const std::filesystem::path& results_path) const;

  // Appends the results block: scalar flux, then angular flux. angular_flux
  // is indexed like InputDeck::Source (values[g], 4*n_x rows x M cols).
  void writeResults(const std::filesystem::path& file_path, const Eigen::MatrixXd& scalar_flux,
                    const std::vector<Eigen::MatrixXd>& angular_flux) const;

  // // Appends the residuals block (one table per group, indexed like
  // // InputDeck::Source).
  // void writeResiduals(const std::filesystem::path& file_path,
  //                     const std::vector<Eigen::MatrixXd>& residuals) const;

  // // Appends the convergence blocks (per-group summary, then the full
  // // per-iteration history) for the most recent solve.
  // void writeConvergence(const std::filesystem::path& file_path) const;

protected:
  TransportOperator transport_operator;
  ConvergenceHistory convergence_;
};

#endif
