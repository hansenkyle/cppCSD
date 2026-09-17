#ifndef SOLVER_H
#define SOLVER_H

#include <filesystem>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "input_deck.h"

class Solver {
public:
  // constructor from input deck (copy)
  Solver(InputDeck input_deck) : input_deck(input_deck) {}
  InputDeck input_deck;

  Eigen::VectorXd transportSweep(int g, Eigen::MatrixXd psi_in_E, Eigen::MatrixXd scalar_flux);

  // Appends the run metadata block.
  void writeMetadata(const std::filesystem::path& file_path) const;

  // Appends an echo of input_deck (InputDeck::echo()).
  void writeInputDeckEcho(const std::filesystem::path& file_path) const;

  // Appends the results block: scalar flux, then angular flux. angular_flux
  // is indexed like InputDeck::Source (values[g], 4*n_x rows x M cols).
  void writeResults(const std::filesystem::path& file_path, const Eigen::MatrixXd& scalar_flux,
                    const std::vector<Eigen::MatrixXd>& angular_flux) const;

  // Appends the residuals block (one table per group, indexed like
  // InputDeck::Source).
  void writeResiduals(const std::filesystem::path& file_path,
                      const std::vector<Eigen::MatrixXd>& residuals) const;

  class Kernel {
    // contains mass matrices, etc.
    // functions include:
    // solveBLD()
    Eigen::Matrix2d M;
    Eigen::Matrix2d L;
    Eigen::Matrix2d Lb;

    Eigen::Matrix4d A;
    Eigen::Vector4d b;

  public:
    Kernel();
    Eigen::Vector4d solveDirect(double cosine, double dx, double dE, double xs, double S,
                                double S_up, double S_down, Eigen::Vector2d psi_in_E,
                                double psi_in_x_down, double psi_in_x_up, Eigen::Vector2d q_up,
                                Eigen::Vector2d q_down, const Eigen::VectorXd& sigma_sdEprime,
                                const Eigen::MatrixXd& phi_gprime_up,
                                const Eigen::MatrixXd& phi_gprime_down);
  };

protected:
  Kernel kernel;
};

#endif
