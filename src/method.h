// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#ifndef METHOD_H
#define METHOD_H

#include "input_deck.h"
#include <Eigen/Dense>
#include <filesystem>
#include <string>
#include <vector>

struct MethodResult {
  Eigen::MatrixXd scalar_flux;
  std::vector<Eigen::MatrixXd> angular_flux;
  Eigen::MatrixXd spectrum() const;
  Eigen::MatrixXd multigroup() const;
  Eigen::MatrixXd cell_average_scalar() const;
  std::vector<Eigen::MatrixXd> cell_average_angular() const;
};

struct Residuals {
  std::vector<Eigen::MatrixXd> high_order;
};

class Method {
public:
  std::string name;

  virtual void solve(double epsilon, int max_iterations) = 0;

  void writeMetadata(const std::filesystem::path& file_path, std::string timestamp) const;
  void writeInputEcho(const std::filesystem::path& file_path) const;

protected:
  Method(std::string name, InputDeck input_deck) : name(name), input_deck(input_deck) {}

  InputDeck input_deck;
  void appendToFile(const std::filesystem::path& file_path, const std::string& text) const;
};

#endif