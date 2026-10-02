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
#include <utility>
#include <vector>

namespace HighFive {
class File;
}

struct MethodResult {
  Eigen::MatrixXd scalar_flux;
  std::vector<Eigen::MatrixXd> angular_flux;
  Eigen::MatrixXd spectrum() const;
  Eigen::MatrixXd multigroup() const;
  Eigen::MatrixXd cell_average_scalar() const;
  // Averages a corner-valued field [4nx x ..] over each space-energy cell -> [nx x ..].
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

  // L2 norm of a corner-valued field [4nx] for one group: the exact integral of its square over
  // x (and over the group's energy width dE, if provided) under the bilinear corner basis.
  double l2norm(const Eigen::VectorXd& vector, double dE = 1) const;
  static double linfnorm(const Eigen::VectorXd& vector);

protected:
  // Validates the deck's cross-struct checks (and derived fields) before anything reads it, so a
  // deck built setter by setter can't reach a solver half-consistent.
  Method(std::string name, InputDeck deck)
      : name(name), input_deck(validated(std::move(deck))), convergence_(input_deck.energy.G) {}

  static InputDeck validated(InputDeck deck) {
    deck.validate();
    return deck;
  }

  InputDeck input_deck;
  ConvergenceHistory convergence_;
  // The results block shared by every method: cell-average scalar flux, cell-average angular
  // flux, then the scalar flux slices (multigroup, spectrum). Returned unrendered so a method can
  // append its own units before writing it.
  UnitGroup solutionBlock(const MethodResult& solution) const;
  // Table of a cell-averaged field [nx x G], laid out and labeled like the cell-average scalar
  // flux.
  MatrixTable cellAverageTable(const std::string& title, const Eigen::MatrixXd& cell_average) const;
  // Corner values; 1 horizontal table per group, 1 row in tables per corner value
  UnitGroup cornerValueTable(const std::string& title, const Eigen::MatrixXd& corner_values,
                             std::string outer_idx = "g") const;

  void appendToFile(const std::filesystem::path& file_path, const std::string& text) const;

  // Row layout of every corner-valued field [4nx], stored as an attribute on its .h5 dataset.
  static constexpr const char* kCornerOrder =
      "row 4i+k is cell i, k = (up L, up R, down L, down R)";
  // Writes the .h5 contents every method shares: run metadata as root attributes, the input deck
  // under /input, and the convergence history under /convergence.
  void writeH5Common(HighFive::File& file, const std::string& timestamp) const;
  // Stacks per-group matrices [R x C] into one [G x R x C] array, and per-group vectors [n] into
  // one [G x n] matrix, for writing as a single .h5 dataset.
  static std::vector<std::vector<std::vector<double>>>
  stackGroups(const std::vector<Eigen::MatrixXd>& per_group);
  static Eigen::MatrixXd stackGroups(const std::vector<Eigen::VectorXd>& per_group);
};

#endif