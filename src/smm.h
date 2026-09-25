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
#include <Eigen/Sparse>
#include <boost/multiprecision/float128.hpp>
#include <filesystem>
#include <string>
#include <utility>

using HighPrecision = boost::multiprecision::float128;

struct SMMResult : public MethodResult {
  Eigen::MatrixXd reconstructed_scalar;
  Eigen::MatrixXd current;

  Eigen::MatrixXd cell_average_current() const;
};

// Closures of the SM equations for one group, evaluated from its angular flux at every corner
// [4nx], in the same layout as scalar_flux. The +/- terms are kept separate: a face value takes
// the + term from the R corner left of the face and the - term from the L corner right of it.
struct SMClosures {
  Eigen::VectorXd F;     // (54)
  Eigen::VectorXd F_pos; // (56a)
  Eigen::VectorXd F_neg; // (56b)
  Eigen::VectorXd K_pos; // (58b)
  Eigen::VectorXd K_neg; // (58a)
  Eigen::VectorXd T_pos; // (58d)
  Eigen::VectorXd T_neg; // (58c)
};

class SecondMoment : public Method {
public:
  SecondMoment(InputDeck input_deck);

  static constexpr int kDefaultMaxIterations = 1000;

  void solve(double epsilon, int max_iterations = kDefaultMaxIterations);

  // Appends the results block (see Method::solutionBlock), then every group's closures.
  void writeResults(const std::filesystem::path& file_path) const;
  void writeResiduals(const std::filesystem::path& file_path, std::string timestamp);
  SMMResult solution;
  std::vector<SMClosures> closures;
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
  //   sigma_sdEprime         : sigma_s0(g' -> g) * dE_g' / dE_g, all g' : [Gx1]
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

  // Closures F, F+-, K+- and T+- at every corner, from one group's angular flux psi [4nx x M].
  SMClosures computeClosures(const Eigen::MatrixXd& psi) const;

  // Group g's LO matrix [8nx x 8nx]: cell i's unknowns [phi_i; J_i] are columns 8i..8i+7 and its
  // equations (53a-h) rows 8i..8i+7. Block tridiagonal; only depends on the problem data, never
  // on the closures, so one factorization serves every iteration within the group.
  Eigen::SparseMatrix<double> buildGroupMatrix(int g);

  // Group g's LO right-hand side [8nx]: every closure, boundary, CSD, scattering and external
  // source term of (53). scalar and current are the LO solution [4nx x G], filled for every group
  // before g; column g itself is not read (within-group scattering is in the matrix).
  Eigen::VectorXd buildGroupRHS(int g, const SMClosures& closures, const Eigen::MatrixXd& scalar,
                                const Eigen::MatrixXd& current);

  // Builds and factorizes group g's LO matrix. Every solveGroup call reuses this factorization
  // until factorizeGroup is called again.
  void factorizeGroup(int g);

  // Solves the factorized group's LO system for rhs, returning (phi_g, J_g), each [4nx].
  std::pair<Eigen::VectorXd, Eigen::VectorXd> solveGroup(const Eigen::VectorXd& rhs) const;

private:
  TransportOperator transport_operator;
  Eigen::SparseLU<Eigen::SparseMatrix<double>> lo_solver_;
  bool lo_factorized_ = false;

  // Half-range moments of group g's incoming boundary flux, (63)-(64), indexed (u, d). The same
  // bc[g] feeds both boundaries: mu > 0 enters at x_0, mu < 0 at x_I. F_in is not in the notes,
  // but the boundary faces' F^b needs it just like (55) does on interior faces.
  struct IncomingMoments {
    Eigen::Vector2d J_pos, J_neg;
    Eigen::Vector2d phi_pos, phi_neg;
    Eigen::Vector2d F_pos, F_neg;
  };
  IncomingMoments incomingMoments(int g) const;

  // Values on every face x_0..x_I [nx+1 x 2], columns (u, d): pos taken at the R corner of the
  // cell left of the face, neg at the L corner of the cell right of it. On an outer face the
  // missing side is the incoming boundary moment in_pos (at x_0) or in_neg (at x_I) instead.
  Eigen::MatrixXd assembleFaces(const Eigen::VectorXd& pos, const Eigen::VectorXd& neg,
                                const Eigen::Vector2d& in_pos, const Eigen::Vector2d& in_neg) const;

  Eigen::VectorXd calculateK(Eigen::MatrixXd psi_slice, int sign) const;
  Eigen::VectorXd calculateT(Eigen::MatrixXd psi_slice, int sign) const;
};

#endif