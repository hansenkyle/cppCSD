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

private:
  TransportOperator transport_operator;
  ConvergenceHistory convergence_;

  std::vector<Eigen::MatrixXd> J_in_positive;
  std::vector<Eigen::MatrixXd> J_in_negative;
  std::vector<Eigen::MatrixXd> phi_in_positive;
  std::vector<Eigen::MatrixXd> phi_in_negative;
};

#endif