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
#include <Eigen/Sparse>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/eigen.hpp>

#include "input_deck.h"

using HighPrecision = boost::multiprecision::cpp_bin_float_50;

class Solver {
public:
  // constructor from input deck (copy)
  Solver(InputDeck input_deck) : input_deck(input_deck) {}
  InputDeck input_deck;

  /// @brief Solve high-order transport equation for all angles in a single energy group
  /// @param g Group index. Needed for slicing sigma_s(g' -> g)
  /// @param psi_in_E Angular flux, next-highest energy group. [4nx by M]
  /// @param scalar_flux Scalar flux, all groups. [4nx by G]
  /// @return Angular flux. [4nx by M]
  Eigen::MatrixXd transportSweep(int g, Eigen::MatrixXd psi_in_E, Eigen::MatrixXd scalar_flux);

  Eigen::VectorXd integrateAngle(Eigen::MatrixXd psi);

  /// @brief Compute scalar flux for all space, all energy groups using Source Iteration.
  ///
  /// @param epsilon Convergence criterion. Transport iteration stops when |phi_old - phi_new|_2 >
  /// |phi_new|_2*epsilon
  /// @return Eigen::MatrixXd. Scalar flux in all energy groups. [4nx by G]
  Eigen::MatrixXd sourceIterate(double epsilon);

  /// @brief Calculate residuals (Ax-b) given a solution and all coefficients; checks that transport
  /// equation was solved correctly in a single cell. Fully independent from solver methods.
  /// @param mu Cos(theta)
  /// @param dx Spatial cell width
  /// @param dE Energy cell width
  /// @param sigma_t Total cross section (cm-1)
  /// @param S_bar Group average stopping power
  /// @param S_Eg Stopping power at lower enegy bound
  /// @param S_Egm1 Stopping power at higher energy bound
  /// @param psi_gm1_d Angular flux, "down" for next-highest energy group
  /// @param psi_b_u Flux AT boundary-- use upwinding conditions
  /// @param psi_b_d Flux AT boundary-- use upwinding conditions
  /// @param q_u External source, L/R values -- upper energy moment
  /// @param q_d External source, L/R values -- higher energy moment
  /// @param sigma_sdEprime sigma_s(g' -> g) * dE_g' for all g'
  /// @param phi_gprime_u scalar flux for all g' -- 'U' moment
  /// @param phi_gprime_d scalar flux for all g' -- 'D' moment
  /// @param psi_up Angular flux solution -- 'U' moment
  /// @param psi_down Angular flux solution -- 'U' moment
  /// @param verbose Prints all residual components for this cell to LOG_INFO
  /// @return Eigen::Vector<HighPrecision, 4>: 4 residual values at extended precision
  Eigen::Vector<HighPrecision, 4>
  cellResidual(double mu, double dx, double dE, double sigma_t, double S_bar, double S_Eg,
               double S_Egm1, const Eigen::Vector2d& psi_gm1_d, const Eigen::Vector2d& psi_b_u,
               const Eigen::Vector2d& psi_b_d, const Eigen::Vector2d& q_u,
               const Eigen::Vector2d& q_d, const Eigen::VectorXd& sigma_sdEprime,
               const Eigen::MatrixXd& phi_gprime_u, const Eigen::MatrixXd& phi_gprime_d,
               const Eigen::Vector2d& psi_up, const Eigen::Vector2d& psi_down,
               bool verbose = false) const;

  /// @brief Calculate residuals for each of 4 equations, all space, one energy group.
  /// @param g Group index, needed to slice scattering matrix
  /// @param angular Psi: angular flux for this energy group. [4nx by M]
  /// @param psi_gm1 Angular flux in next-highest energy group. [4nx by M]
  /// @param scalar Scalar flux, all energy groups. [4nx by G]
  /// @param debug_max cellResidual set to verbose for every cell, default false
  /// @return Eigen::MatrixXd [4nx by M]. Values demoted to double-precision
  Eigen::MatrixXd calculateResiduals(int g, const Eigen::MatrixXd& angular,
                                     const Eigen::MatrixXd& psi_gm1, const Eigen::MatrixXd& scalar,
                                     bool debug_max = false);

  // Appends the run metadata block.
  void writeMetadata(const std::filesystem::path& file_path) const;

  // Appends an echo of input_deck (InputDeck::echo()).
  void writeInputDeckEcho(const std::filesystem::path& file_path) const;

  // Appends the results block: scalar flux, then angular flux. angular_flux
  // is indexed like InputDeck::Source (values[g], 4*n_x rows x M cols).
  void writeResults(const std::filesystem::path& file_path, const Eigen::MatrixXd& scalar_flux,
                    const std::vector<Eigen::MatrixXd>& angular_flux) const;

  // Appends the residuals block (one table per group, indexed like
  // InputDeck::Source).
  void writeResiduals(const std::filesystem::path& file_path,
                      const std::vector<Eigen::MatrixXd>& residuals) const;

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
  Kernel kernel;
};

#endif
