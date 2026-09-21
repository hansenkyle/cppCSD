// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#ifndef SOURCE_ITERATION_H
#define SOURCE_ITERATION_H

#include "input_deck.h"
#include "method.h"
#include "transport_operator.h"

class SourceIteration : public Method {
public:
  SourceIteration(InputDeck input_deck);

private:
  InputDeck input_deck;
  TransportOperator transport_operator;
};

#endif