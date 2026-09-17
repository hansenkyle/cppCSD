#include "solver.h"
#include "solver_formatter.h"

#include <array>
#include <cstdlib>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <vector>

Solver::Kernel::Kernel() {
  M << 2.0, 1.0, 1.0, 2.0;
  M *= (1.0 / 6);
  L << 0.5, 0.5, -0.5, -0.5;
  Lb << -1, 0, 0, -1;

  A = Eigen::Matrix4d::Zero();
  b = Eigen::Vector4d::Zero();
}

std::vector<Eigen::VectorXd> Solver::transportSweep(int g){

  // Compute angular flux in a single energy group given a known source (and known scalar flux in all groups).
  
  // references: delete these
  // M = input_deck.angle.M
  // mu = input_deck.angle.mu
  // sigma_total = input_deck.xs.total // Eigen matrix(G, x)
  // boundary = input_deck.bc[group](up/down, m)
  // dx = input_deck.mesh.dx / .mesh.n_x   // Eigen::VectorXd
  // dE = input_deck.energy.dE  // Eigen::VectorXd



  // initialize guess
  std::vector<Eigen::VectorXd> psi(input_deck.angle.M, Eigen::VectorXd::Zero(4*input_deck.mesh.n_x));

  // prepare data (dE*sigma_s)
  auto& dx = input_deck.mesh.dx;
  Eigen::VectorXd& sigma_t = input_deck.xs.total(g, Eigen::placeholds::all);
  Eigen::
  double dE = input_deck.energy.dE(g);


  for (int m=0; m<input_deck.angle.M; m++){
    auto mu = input_deck.angle.mu[m];
    switch (mu > 0) {
      case true: // left-to-right
      // solve leftmost cell using boundary conditions
      auto dx = input_deck.mesh.dx(0);
      auto dE = input_deck.energy.dE(g);
      
      psi[m]({0,1,2,3}) = kernel.solveDirect(mu,)
      // loop through all other cells
      break;
      case false: // right-to-left
      // solve rightmost cell using boundary conditions

      // loop through all other cells
      break;
    }
    // switch (mu > 0)
    // case true
    // case false
  }




  // return
  return psi;
}

Eigen::Vector4d Solver::Kernel::solveDirect(
    double cosine, double dx, double dE, double xs, double S, double S_up, double S_down,
    Eigen::Vector2d psi_in_E, double psi_in_x_down, double psi_in_x_up, Eigen::Vector2d q_up,
    Eigen::Vector2d q_down, const Eigen::VectorXd& sigma_sdEprime,
    const Eigen::MatrixXd& phi_gprime_up, const Eigen::MatrixXd& phi_gprime_down) {
  A = Eigen::Matrix4d::Zero();
  b = Eigen::Vector4d::Zero();

  double mu = std::abs(cosine);

  // The element matrices below are built assuming flow travels L->R; for
  // cosine < 0 the L/R labeling of every spatially-structured input must be
  // swapped to match before assembly (scalar psi_in_x_up/down are already
  // direction-relative "upwind" values, so they're left alone).
  if (cosine < 0) {
    psi_in_E.reverseInPlace();
    q_up.reverseInPlace();
    q_down.reverseInPlace();
  }

  /*
  matrix/vector are energy-major, space-minor:
  up:    L
         R

  down:  L
         R
  */

  // Parameters:
  //   cosine             : angle consine (mu)
  //   dx, dE             : cell widhts
  //   xs                 : total xs
  //   S                  : group-average stopping power
  //   S_up, S_down       : S(g-1), S(g)
  //   psi_in_E           : flux at next-higher energy group, L/R    : [2x1]
  //   psi_in_x_up        : upwind flux in same energy group, up     : scalar
  //   psi_in_x_down      : "                              ", down   : scalar
  //   q_up               : external source, up (L/R)                : [2x1]
  //   q_down             : "             ", down (L/R)              : [2x1]
  //   sigma_sdEprime     : sigma_s(g' -> g) * dE_g' for all g'      : [Gx1]
  //   phi_gprime_up/down : scalar flux in all groups                : [2xG]

  // "Up" LHS
  // streaming (L)
  A({0, 1}, {0, 1}) += (mu / 6) * 2 * L;
  A({0, 1}, {2, 3}) += (mu / 6) * L;
  // streaming (Lb)
  A(1, 1) += (mu / 6) * 2;
  A(1, 3) += (mu / 6);
  // absorption + CSD loss
  A({0, 1}, {0, 1}) += dx * (xs / 3 + S / (2 * dE)) * M;
  A({0, 1}, {2, 3}) += dx * (xs / 6 + S / (2 * dE)) * M;

  // "Up" RHS
  // streaming source
  b(0) += (mu / 6) * (2 * psi_in_x_up + psi_in_x_down);
  // CSD source
  b({0, 1}) += (dx / dE) * S_up * M * psi_in_E;
  // Scattering source
  b({0, 1}) += (0.125) * M * (phi_gprime_down + phi_gprime_up) * sigma_sdEprime;
  // External source
  b({0, 1}) += (dx / 6) * M * (2 * q_up + q_down);

  // "Down" LHS
  // streaming (L)
  A({2, 3}, {0, 1}) += (mu / 6) * L;
  A({2, 3}, {2, 3}) += (mu / 6) * 2 * L;
  // Streaming (Lb)
  A(3, 1) += (mu / 6);
  A(3, 3) += (mu / 6) * 2;
  // absorption + CSD loss
  A({2, 3}, {0, 1}) += dx * (xs / 6 - S / (2 * dE)) * M;
  A({2, 3}, {2, 3}) += dx * (xs / 3 + (S_down - S / 2) / dE) * M;

  // "Down" RHS
  // streaming source
  b(2) += (mu / 6) * (psi_in_x_up + 2 * psi_in_x_down);
  // Scattering source
  b({2, 3}) += (0.125) * M * (phi_gprime_down + phi_gprime_up) * sigma_sdEprime;
  // External source
  b({2, 3}) += (dx / 6) * M * (q_up + 2 * q_down);

  Eigen::Vector4d x = A.partialPivLu().solve(b);

  if (cosine < 0) {
    return x({1, 0, 3, 2});
  }

  return x;

  return Eigen::Vector4d::Zero();
}

namespace {

void appendToFile(const std::filesystem::path& file_path, const std::string& text) {
  std::ofstream out(file_path, std::ios::app);
  if (!out.is_open()) {
    throw std::runtime_error("Solver: failed to open '" + file_path.string() + "' for writing");
  }
  out << text;
}

} // namespace

void Solver::writeMetadata(const std::filesystem::path& file_path) const {
  appendToFile(file_path, SolverFormatter::formatRunMetadata());
}

void Solver::writeInputDeckEcho(const std::filesystem::path& file_path) const {
  appendToFile(file_path, input_deck.echo());
}

void Solver::writeResults(const std::filesystem::path& file_path, const Results& results) const {
  appendToFile(file_path, SolverFormatter::formatResults(results, input_deck));
}

void Solver::writeResiduals(const std::filesystem::path& file_path,
                            const std::vector<Eigen::MatrixXd>& residuals) const {
  appendToFile(file_path, SolverFormatter::formatResiduals(residuals, input_deck));
}
