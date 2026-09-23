// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "smm.h"

#include <Eigen/Dense>

Eigen::MatrixXd SMMResult::cell_average_current() const {
  using Eigen::seqN;
  using Eigen::placeholders::all;
  int I = scalar_flux.rows() / 4;
  int G = scalar_flux.cols();

  Eigen::MatrixXd leftsum = (current(seqN(0, I, 4), all) + current(seqN(2, I, 4), all));
  Eigen::MatrixXd rightsum = (current(seqN(1, I, 4), all) + current(seqN(3, I, 4), all));
  Eigen::MatrixXd result = (leftsum + rightsum) / 4;
  return result;
}

SecondMoment::SecondMoment(InputDeck input_deck)
    : Method("second moment method", input_deck), transport_operator(input_deck),
      convergence_(input_deck.energy.G) {

  J_in_positive = std::vector<Eigen::MatrixXd>(input_deck.energy.G);
  J_in_negative = J_in_positive;
  phi_in_positive = J_in_positive;
  phi_in_negative = J_in_positive;

  auto& bc = input_deck.bc.values;
  auto& w = input_deck.angle.w;
  auto& mu = input_deck.angle.mu;
  auto& w_pos = input_deck.angle.w_positive;
  auto& w_neg = input_deck.angle.w_negative;
  for (int g = 0; g < input_deck.energy.G; g++) {
    phi_in_positive[g] = transport_operator.integrateAngle(bc, w_pos);
    phi_in_negative[g] = transport_operator.integrateAngle(bc, w_neg);
    J_in_positive[g] = transport_operator.integrateAngle(bc, mu.cwiseProduct(w_pos));
    J_in_negative[g] = transport_operator.integrateAngle(bc, mu.cwiseProduct(w_neg));
  }
}
