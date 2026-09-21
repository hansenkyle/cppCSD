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

class TransportOperator {
public:
  TransportOperator(InputDeck input_deck);
  Eigen::MatrixXd sweepTransport();
  Eigen::VectorXd integrateAngle();

private:
  InputDeck input_deck;
};

#endif