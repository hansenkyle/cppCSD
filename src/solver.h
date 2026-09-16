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
    Kernel();
    // contains mass matrices, etc.
    // functions include:
    // solveBLD()
    Eigen::Matrix2d M;
    Eigen::Matrix2d L;
    Eigen::Matrix2d Lb;

    Eigen::Matrix4d A;
    Eigen::Vector4d b;
  };

  Kernel kernel();
};

#endif
