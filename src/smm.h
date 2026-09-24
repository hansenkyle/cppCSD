// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#ifndef SMM_H
#define SMM_H

#include "convergence.h"
#include "input_deck.h"
#include "method.h"
#include "transport_operator.h"

#include <Eigen/Dense>
#include <boost/multiprecision/float128.hpp>
#include <filesystem>
#include <string>

using HighPrecision = boost::multiprecision::float128;

struct SMMResult : public MethodResult {
  Eigen::MatrixXd reconstructed_scalar;
  Eigen::MatrixXd current;

  Eigen::MatrixXd cell_average_current() const;
};

class SecondMoment : public Method {
public:
  SecondMoment(InputDeck input_deck);

  static constexpr int kDefaultMaxIterations = 1000;

  void solve(double epsilon, int max_iterations = kDefaultMaxIterations);

  void writeResults(const std::filesystem::path& file_path) const;
  void writeConvergence(const std::filesystem::path& file_path) const;
  void writeResiduals(const std::filesystem::path& file_path, std::string timestamp) const;
  SMMResult solution;
  Residuals residuals;

  // Residuals of the SM equations for group g over the whole mesh: rows 8i..8i+7 hold cell i's
  // cellResidual. scalar and current are the LO solution [4nx x G], filled for every group up to
  // and including g (g-1 feeds the CSD source, g' <= g the scattering source). psi is group g's
  // angular flux [4nx x M], which the closures F, K and T are computed from. debug_max
  // re-evaluates the cell with the largest residual term-by-term and logs the breakdown.
  Eigen::VectorXd calculateResiduals(int g, const Eigen::MatrixXd& scalar,
                                     const Eigen::MatrixXd& current, const Eigen::MatrixXd& psi,
                                     bool debug_max = false);

  // Residuals (LHS - RHS) of the eight SM equations (53a-h) in one cell, in that order. Pairs
  // are (L, R) at corners and (x_{i-1}, x_i) on faces; _u/_d name the energy edge.
  //   dx, dE                 : cell widths
  //   sigma_t                : total xs
  //   S_bar, S_Eg, S_Egm1    : group-average S, S(E_g), S(E_{g-1})
  //   phi_gm1_d, J_gm1_d     : group g-1's LO solution, d edge (zero when g = 0)
  //   phi_b_*, J_b_*         : face values from the interface conditions (57), (59)-(62)
  //   F_b_*, F_*             : SM closure on faces (55) and at corners (54)
  //   q0_*, q1_*             : zeroth and first angular moments of the external source
  //   sigma_sdEprime         : sigma_s0(g' -> g) * dE_g' for all g'      : [Gx1]
  //   phi_gprime_u/d         : scalar flux in all groups                 : [2xG]
  //   phi_*, J_*             : this group's LO unknowns
  Eigen::Vector<HighPrecision, 8>
  cellResidual(double dx, double dE, double sigma_t, double S_bar, double S_Eg, double S_Egm1,
               const Eigen::Vector2d& phi_gm1_d, const Eigen::Vector2d& J_gm1_d,
               const Eigen::Vector2d& phi_b_u, const Eigen::Vector2d& phi_b_d,
               const Eigen::Vector2d& J_b_u, const Eigen::Vector2d& J_b_d,
               const Eigen::Vector2d& F_b_u, const Eigen::Vector2d& F_b_d,
               const Eigen::Vector2d& F_u, const Eigen::Vector2d& F_d, const Eigen::Vector2d& q0_u,
               const Eigen::Vector2d& q0_d, const Eigen::Vector2d& q1_u,
               const Eigen::Vector2d& q1_d, const Eigen::VectorXd& sigma_sdEprime,
               const Eigen::MatrixXd& phi_gprime_u, const Eigen::MatrixXd& phi_gprime_d,
               const Eigen::Vector2d& phi_u, const Eigen::Vector2d& phi_d,
               const Eigen::Vector2d& J_u, const Eigen::Vector2d& J_d, bool verbose = false) const;

private:
  TransportOperator transport_operator;
  ConvergenceHistory convergence_;

  std::vector<Eigen::MatrixXd> J_in_positive;
  std::vector<Eigen::MatrixXd> J_in_negative;
  std::vector<Eigen::MatrixXd> phi_in_positive;
  std::vector<Eigen::MatrixXd> phi_in_negative;

  Eigen::VectorXd calculateK(Eigen::MatrixXd psi_slice, int sign) const;
  Eigen::VectorXd calculateT(Eigen::MatrixXd psi_slice, int sign) const;
  Eigen::VectorXd solveSM(Eigen::VectorXd Kpos, Eigen::VectorXd Kneg, Eigen::VectorXd Tpos,
                          Eigen::VectorXd Tneg) const;
};

#endif