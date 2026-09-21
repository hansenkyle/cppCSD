// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#ifndef METHOD_H
#define METHOD_H

#include <Eigen/Dense>
#include <string>
#include <vector>

class Method {
public:
  std::string name;
  virtual void solve() = 0;
};

struct MethodResult {
  Eigen::MatrixXd scalar_flux;
  std::vector<Eigen::MatrixXd> angular_flux;
};

#endif