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

  // Flux solution fields, as produced by the solver. scalar_flux and each
  // angular_flux entry are (4 * input_deck.mesh.n_x) rows by
  // input_deck.energy.G columns: rows are cell-major, and within a cell
  // ordered up_left, up_right, down_left, down_right. angular_flux has one
  // entry per ordinate (input_deck.angle.M), in ordinate order.
  struct Results {
    Eigen::MatrixXd scalar_flux;
    std::vector<Eigen::MatrixXd> angular_flux;
  };

  // The four writeX() methods below each append one block to file_path
  // (creating it if it doesn't exist) rather than truncating -- so calling
  // them in sequence over the course of a run, as each block's data
  // becomes available, builds up the same file a single call would have,
  // but leaves a self-describing partial file behind if the run crashes
  // partway through. Callers are responsible for calling them in the
  // intended order: metadata and the input deck echo first (before/at the
  // start of the solve), then results and residuals once the solve
  // produces them. Each formats its block via SolverFormatter
  // (solver_formatter.h) -- Solver supplies the data, not the layout.

  // Appends the run metadata block.
  void writeMetadata(const std::filesystem::path& file_path) const;

  // Appends an echo of input_deck (InputDeck::echo()).
  void writeInputDeckEcho(const std::filesystem::path& file_path) const;

  // Appends the results block (scalar and angular flux).
  void writeResults(const std::filesystem::path& file_path, const Results& results) const;

  // Appends the residuals block (one table per ordinate).
  void writeResiduals(const std::filesystem::path& file_path,
                      const std::vector<Eigen::MatrixXd>& residuals) const;

protected:
  class Kernel {
    // contains mass matrices, etc.
    // functions include:
    // solveBLD()
  };

  Kernel kernel();
};

#endif
