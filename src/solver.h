#ifndef SOLVER_H
#define SOLVER_H

#include <Eigen/Sparse>

#include "cross_section.h"
#include "input_deck.h"
#include "mesh.h"

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
};

#endif
