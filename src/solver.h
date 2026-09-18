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

  Eigen::MatrixXd transportSweep(int g, Eigen::MatrixXd psi_in_E, Eigen::MatrixXd scalar_flux);

  Eigen::VectorXd integrateAngle(Eigen::MatrixXd psi);

  Eigen::MatrixXd sourceIterate(double epsilon);

  // Independently re-derives one cell's local system (mass/streaming
  // matrices, cross sections, sources) in extended precision and evaluates
  // the discretized BTE's residual (eq. 41a for rows 0,1 [up_left,
  // up_right]; eq. 41b for rows 2,3 [down_left, down_right]) against
  // candidate, the corner values solveDirect produced for this cell.
  // Deliberately rebuilds M/L/Lb here rather than reusing Kernel's, so a
  // bug in that assembly can't hide from its own residual check -- same
  // reasoning as residual.h.
  //
  // Unlike Kernel::solveDirect, mu keeps its sign here (no abs(), no L/R
  // swap): eq. 41c's mu-dependent upwind selection is resolved by the
  // caller (the phase-space loop) before calling this, so psi_b_u/psi_b_d
  // arrive already correct for whatever sign mu has.
  //
  // Symbol correspondence with the writeup:
  //   mu                    : mu_m
  //   dx                    : Delta x_i
  //   dE                    : Delta E_g
  //   sigma_t               : sigma_{t,g,i}
  //   S_bar                 : bar S_{g,i}            (group-average stopping power)
  //   S_Eg                  : S_i(E_g)                (stopping power at this group's lower edge)
  //   S_Egm1                : S_i(E_{g-1})             (stopping power at this group's upper edge)
  //   psi_gm1_d             : Psi_{m,g-1,d,i}          (previous group's "d"-edge flux, CSD source)
  //   psi_b_u, psi_b_d      : Psi^b_{m,g,i,u}, Psi^b_{m,g,i,d}   (eq. 41c, already resolved)
  //   q_u, q_d              : q_{m,g,i,u}, q_{m,g,i,d}
  //   sigma_sdEprime        : Delta E_g' * sigma_{s0,g'->g,i}, all g'
  //   phi_gprime_u/d        : Phi_{g',u,i}, Phi_{g',d,i}, all g'  (each [L,R] x g' matrix)
  //   candidate             : [up_left, up_right, down_left, down_right] = [Psi_u,L Psi_u,R Psi_d,L
  //   Psi_d,R]
  // verbose: if true, logs every named intermediate term (streaming,
  // absorption+CSD-loss, CSD source, scattering source, external source,
  // and the final residual) for both eq. 41a and eq. 41b, one LDCSD_LOG_INFO
  // line each -- for tracking down which term disagrees with solveDirect
  // at a specific cell. Off by default since it's very noisy.
  Eigen::Vector<HighPrecision, 4>
  cellResidual(double mu, double dx, double dE, double sigma_t, double S_bar, double S_Eg,
               double S_Egm1, const Eigen::Vector2d& psi_gm1_d, const Eigen::Vector2d& psi_b_u,
               const Eigen::Vector2d& psi_b_d, const Eigen::Vector2d& q_u,
               const Eigen::Vector2d& q_d, const Eigen::VectorXd& sigma_sdEprime,
               const Eigen::MatrixXd& phi_gprime_u, const Eigen::MatrixXd& phi_gprime_d,
               const Eigen::Vector2d& psi_up, const Eigen::Vector2d& psi_down,
               bool verbose = false) const;

  // Residuals for one energy group g. angular is that group's own psi
  // (4*n_x x M); psi_gm1 is the previous group's converged psi, same
  // shape (needed for the CSD source -- pass a zero matrix for g==0), the
  // same value sourceIterate already tracks as psi_up. scalar is still
  // all groups (4*n_x x G), since the scattering source sums over g'.
  //
  // debug_max: if true, after computing the full grid, finds the cell/
  // ordinate with the largest |residual|, logs which one it picked, and
  // re-evaluates cellResidual there with verbose=true.
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
    // check_condition: if true, logs A's condition number (via JacobiSVD,
    // largest/smallest singular value) before solving -- a diagnostic for
    // whether this cell's system is too ill-conditioned for a plain double
    // partialPivLu solve to be trusted. Off by default since it's not free.
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
