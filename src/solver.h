#ifndef SOLVER_H
#define SOLVER_H

#include <filesystem>
#include <vector>

#include <Eigen/Sparse>

#include "input_deck.h"

class Solver {
public:
  // constructor from input deck (copy)
  Solver(InputDeck input_deck) : input_deck(input_deck) {}
  InputDeck input_deck;

protected:
  class Kernel {
    // contains mass matrices, etc.
    // functions include:
    // solveBLD()
  };

  Kernel kernel();
};

class SourceIterationSolver : public Solver {
public:
  using Solver::Solver;

  // Runs source iteration to convergence (or until params.max_iters), per
  // params. Takes the concrete params type directly -- no base
  // IterationParameters type, no runtime check that it's "the right"
  // derived type: callers already have a SourceIterationSolver, so they
  // already know which params type matches it.
  void solve(const SourceIterationParams& params);

  // Writes this solve's results to path. Not virtual: each solver's
  // results are shaped differently, so there's no shared interface for
  // this to satisfy.
  void writeResults(const std::filesystem::path& path) const;
};

#endif
