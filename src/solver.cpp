#include "solver.h"

#include <array>
#include <cstdlib>
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

Eigen::Vector4d Solver::Kernel::solveDirect(double cosine, double dx, double dE, double xs,
                                            double S, double S_up, double S_down,
                                            Eigen::Vector2d psi_in_E, double psi_in_x_down,
                                            double psi_in_x_up, Eigen::Vector2d q_up,
                                            Eigen::Vector2d q_down, Eigen::VectorXd sigma_sdEprime, Eigen::MatrixXd phi_gprime_up,
                                            Eigen::MatrixXd phi_gprime_down) {
  A = Eigen::Matrix4d::Zero();
  b = Eigen::Vector4d::Zero();

  double mu = std::abs(cosine);

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
  b({0, 1}) += (0.125) *M * (phi_gprime_down + phi_gprime_up)*sigma_sdEprime;
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
  b({2,3}) += (0.125) *M * (phi_gprime_down + phi_gprime_up)*sigma_sdEprime;
  // External source
  b({2, 3}) += (dx / 6) * M * (q_up + 2 * q_down);

  Eigen::Vector4d x = A.partialPivLu().solve(b);

  if (cosine < 0) {
    return x({1, 0, 3, 2});
  }

  return x;

  return Eigen::Vector4d::Zero();
}