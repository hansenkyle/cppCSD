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

#endif
