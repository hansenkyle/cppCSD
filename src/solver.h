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

    public:
    Eigen::Vector4d solveDirect(double cosine, double dx, double dE, double xs,
                                            double S, double S_up, double S_down,
                                            Eigen::Vector2d psi_in_E, double psi_in_x_down,
                                            double psi_in_x_up, Eigen::Vector2d q_up,
                                            Eigen::Vector2d q_down, Eigen::VectorXd sigma_sdEprime, Eigen::MatrixXd phi_gprime_up,
                                            Eigen::MatrixXd phi_gprime_down);
  };

  Kernel kernel();
};

#endif
