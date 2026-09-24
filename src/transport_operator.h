// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#ifndef TRANSPORT_OPERATOR_H
#define TRANSPORT_OPERATOR_H

#include "input_deck.h"
#include <Eigen/Dense>
#include <boost/multiprecision/float128.hpp>

using HighPrecision = boost::multiprecision::float128;

class TransportOperator {
public:
  TransportOperator();
  TransportOperator(InputDeck input_deck);
  Eigen::MatrixXd sweep(int g, Eigen::MatrixXd psi_in_E, Eigen::MatrixXd scalar_flux);
  Eigen::VectorXd integrateAngle(Eigen::MatrixXd psi);
  Eigen::VectorXd integrateAngle(Eigen::MatrixXd psi, Eigen::VectorXd weight);

  Eigen::MatrixXd calculateResiduals(int g, const Eigen::MatrixXd& angular,
                                     const Eigen::MatrixXd& psi_gm1, const Eigen::MatrixXd& scalar,
                                     bool debug_max = false);

  Eigen::Vector<HighPrecision, 4>
  cellResidual(double mu, double dx, double dE, double sigma_t, double S_bar, double S_Eg,
               double S_Egm1, const Eigen::Vector2d& psi_gm1_d, const Eigen::Vector2d& psi_b_u,
               const Eigen::Vector2d& psi_b_d, const Eigen::Vector2d& q_u,
               const Eigen::Vector2d& q_d, const Eigen::VectorXd& sigma_sdEprime,
               const Eigen::MatrixXd& phi_gprime_u, const Eigen::MatrixXd& phi_gprime_d,
               const Eigen::Vector2d& psi_up, const Eigen::Vector2d& psi_down,
               bool verbose = false) const;
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

private:
  InputDeck input_deck;
  Kernel kernel;
};

#endif