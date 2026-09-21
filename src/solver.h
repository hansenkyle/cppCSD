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
#include "transport_operator.h"

using HighPrecision = boost::multiprecision::float128;

class Solver {
public:
  // constructor from input deck (copy)
  Solver(InputDeck input_deck) : input_deck(input_deck), transport_operator(input_deck) {}
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
  Eigen::MatrixXd sourceIterate(double epsilon, int max_iterations = kDefaultMaxIterations);

  // Per-iteration convergence record from the most recent sourceIterate().
  // Empty before the first solve.
  const ConvergenceHistory& convergence() const { return convergence_; }

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

  class Kernel {
    // contains mass matrices, etc.
    // functions include:
    // solveBLD()
    Eigen::Matrix2d M;
    Eigen::Matrix2d L;
    Eigen::Matrix2d Lb;

    Eigen::Matrix4d A;
    Eigen::Vector4d b;

  public:
    Kernel();

    /// @brief Solve the high-order transport equation in a single cell. Forms 4x4 system, solved
    /// with Eigen direct solver. Problem is rotated internally, forces mu>0
    /// @param cosine Angle cosine, "mu"
    /// @param dx Cell width
    /// @param dE Energy group width
    /// @param xs Total cross section in this cell and energy group
    /// @param S  Group-average stopping power in this cell and energy group
    /// @param S_up Stopping power at upper energy boundary
    /// @param S_down Stopping power at lower energy boundary
    /// @param psi_in_E L/R pair; "D" moment of incoming-in-E flux
    /// @param psi_in_x_down Incoming flux, "D" moment
    /// @param psi_in_x_up Incoming flux, "U" moment
    /// @param q_up External source, "U" moment, L/R pair
    /// @param q_down External source, "D" moment, L/R pair
    /// @param sigma_sdEprime Sigma_s(g' -> g) times dE(g') for all g'
    /// @param phi_gprime_up Scalar flux, L/R pair, "U" moment, all energy groups
    /// @param phi_gprime_down Scalar flux, L/R pair, "D" moment, all energy groups
    /// @param check_condition Print matrix's condition number to LOG_INFO, default false
    /// @return Angular flux in this cell: [U_L, U_R, D_L, D_R]^T
    Eigen::Vector4d solveDirect(double cosine, double dx, double dE, double xs, double S,
                                double S_up, double S_down, Eigen::Vector2d psi_in_E,
                                double psi_in_x_down, double psi_in_x_up, Eigen::Vector2d q_up,
                                Eigen::Vector2d q_down, const Eigen::VectorXd& sigma_sdEprime,
                                const Eigen::MatrixXd& phi_gprime_up,
                                const Eigen::MatrixXd& phi_gprime_down,
                                bool check_condition = false);
  };

protected:
  TransportOperator transport_operator;
  Kernel kernel;
  ConvergenceHistory convergence_;
};

#endif
