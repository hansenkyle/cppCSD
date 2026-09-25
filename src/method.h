// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#ifndef METHOD_H
#define METHOD_H

#include "convergence.h"
#include "input_deck.h"
#include "output_block.h"
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
  // Averages a corner-valued field [4nx x G] over each space-energy cell -> [nx x G].
  static Eigen::MatrixXd cell_average(const Eigen::MatrixXd& corners);
  std::vector<Eigen::MatrixXd> cell_average_angular() const;
};

struct Residuals {
  std::vector<Eigen::MatrixXd> high_order;
  // Low-order (SM equation) residuals per group, each [8nx], as returned by
  // SecondMoment::calculateResiduals.
  std::vector<Eigen::VectorXd> low_order;
};

class Method {
public:
  std::string name;

  virtual void solve(double epsilon, int max_iterations) = 0;

  void writeMetadata(const std::filesystem::path& file_path, std::string timestamp) const;
  void writeInputEcho(const std::filesystem::path& file_path) const;
  void writeConvergence(const std::filesystem::path& file_path) const;

  // integrates a vector of corner values over x (and dE if provided)
  double l2norm(const Eigen::VectorXd vector, double dE = 1);
  double linfnorm(const Eigen::VectorXd vector);

protected:
  Method(std::string name, InputDeck input_deck)
      : name(name), input_deck(input_deck), convergence_(input_deck.energy.G) {}

  InputDeck input_deck;
  ConvergenceHistory convergence_;
  // The results block shared by every method: cell-average scalar flux, cell-average angular
  // flux, then the scalar flux slices (multigroup, spectrum). Returned unrendered so a method can
  // append its own units before writing it.
  UnitGroup solutionBlock(const MethodResult& solution) const;
  // Table of a cell-averaged field [nx x G], laid out and labeled like the cell-average scalar
  // flux.
  MatrixTable cellAverageTable(const std::string& title, const Eigen::MatrixXd& cell_average) const;
  void appendToFile(const std::filesystem::path& file_path, const std::string& text) const;
};

#endif