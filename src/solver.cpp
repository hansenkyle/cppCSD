#include "solver.h"
#include "logger.h"
#include "solver_formatter.h"

#include <array>
#include <cstdlib>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

Solver::Kernel::Kernel() {
  M << 2.0, 1.0, 1.0, 2.0;
  M *= (1.0 / 6);
  L << 0.5, 0.5, -0.5, -0.5;
  Lb << -1, 0, 0, 1;

  A = Eigen::Matrix4d::Zero();
  b = Eigen::Vector4d::Zero();
}

Eigen::MatrixXd Solver::transportSweep(int g, Eigen::MatrixXd psi_in_E,
                                       Eigen::MatrixXd scalar_flux) {
  // Compute angular flux in a single energy group given a known source (and known scalar flux in
  // all groups).

  // Parameters:
  //   g                  : energy group. Needed for indexing xs and dE.
  //   psi_in_E           : angular flux at higher energy group  : [4nx x M]
  //   scalar_flux        : scalar flux, all energy groups       : [4nx x G]

  using Eigen::seqN;
  using Eigen::placeholders::all;

  // initialize guess
  Eigen::MatrixXd psi = Eigen::MatrixXd::Zero(4 * input_deck.mesh.n_x, input_deck.angle.M);

  // prepare data (dE*sigma_s)
  double dE = input_deck.energy.dE(g);
  auto dx = input_deck.mesh.dx;
  auto sigma_t = input_deck.xs.total(g, all);
  auto S = input_deck.xs.S(g, all);
  auto S_up = input_deck.xs.S_bound(g, all);
  auto S_down = input_deck.xs.S_bound(g + 1, all);
  double bc_up, bc_down;
  Eigen::Vector2d q_up, q_down;
  Eigen::VectorXd sigmaSdEprime;
  Eigen::MatrixXd phi;

  // loop over all angles
  for (int m = 0; m < input_deck.angle.M; m++) {
    bool print_condition_number=false;
    auto mu = input_deck.angle.mu[m];
    auto q = input_deck.source.values[g](all, m);

    int i = 0;

    switch (mu > 0) {

    case true: // left-to-right
      // solve leftmost cell using boundary conditions
      i = 0;
      // slice data
      q_up = q(seqN(i * 4, 2));
      q_down = q(seqN(i * 4 + 2, 2));
      phi = scalar_flux(seqN(i * 4, 4), all);
      sigmaSdEprime = input_deck.xs.scatter[i].col(g).cwiseProduct(input_deck.energy.dE);
      bc_up = input_deck.bc[g](0, m);
      bc_down = input_deck.bc[g](1, m);

      psi(seqN(i * 4, 4), m) = kernel.solveDirect(
          mu, dx[i], dE, sigma_t[i], S[i], S_up[i], S_down[i], psi_in_E(seqN(i * 4 + 2, 2), m),
          bc_down, bc_up, q_up, q_down, sigmaSdEprime, phi({0, 1}, all), phi({2, 3}, all));
      // loop through all other cells
      for (i = 1; i < input_deck.mesh.n_x; i++) {
        // slice data
        q_up = q(seqN(i * 4, 2));
        q_down = q(seqN(i * 4 + 2, 2));
        phi = scalar_flux(seqN(i * 4, 4), all);
        sigmaSdEprime = input_deck.xs.scatter[i].col(g).cwiseProduct(input_deck.energy.dE);
        bc_up = psi((i - 1) * 4 + 1, m);
        bc_down = psi((i - 1) * 4 + 3, m);

        psi(seqN(i * 4, 4), m) = kernel.solveDirect(
            mu, dx[i], dE, sigma_t[i], S[i], S_up[i], S_down[i], psi_in_E(seqN(i * 4 + 2, 2), m),
            bc_down, bc_up, q_up, q_down, sigmaSdEprime, phi({0, 1}, all), phi({2, 3}, all));
      }
      break;    // left-to-right
    case false: // right-to-left
      // solve rightmost cell using boundary conditions

      i = input_deck.mesh.n_x - 1;
      // slice data
      q_up = q(seqN(i * 4, 2));
      q_down = q(seqN(i * 4 + 2, 2));
      phi = scalar_flux(seqN(i * 4, 4), all);
      sigmaSdEprime = input_deck.xs.scatter[i].col(g).cwiseProduct(input_deck.energy.dE);
      bc_up = input_deck.bc[g](0, m);
      bc_down = input_deck.bc[g](1, m);
      psi(seqN(i * 4, 4), m) = kernel.solveDirect(
          mu, dx[i], dE, sigma_t[i], S[i], S_up[i], S_down[i], psi_in_E(seqN(i * 4 + 2, 2), m),
          bc_down, bc_up, q_up, q_down, sigmaSdEprime, phi({0, 1}, all), phi({2, 3}, all),
          print_condition_number);

      // loop through all other cells
      for (i = input_deck.mesh.n_x - 2; i > -1; i--) {
        // slice data
        q_up = q(seqN(i * 4, 2));
        q_down = q(seqN(i * 4 + 2, 2));
        phi = scalar_flux(seqN(i * 4, 4), all);
        sigmaSdEprime = input_deck.xs.scatter[i].col(g).cwiseProduct(input_deck.energy.dE);
        bc_up = psi((i + 1) * 4, m);
        bc_down = psi((i + 1) * 4 + 2, m);

        psi(seqN(i * 4, 4), m) = kernel.solveDirect(
            mu, dx[i], dE, sigma_t[i], S[i], S_up[i], S_down[i], psi_in_E(seqN(i * 4 + 2, 2), m),
            bc_down, bc_up, q_up, q_down, sigmaSdEprime, phi({0, 1}, all), phi({2, 3}, all));
      }
      break;
    }
  }
  // return
  return psi;
}

Eigen::VectorXd Solver::integrateAngle(Eigen::MatrixXd psi) {
  // psi: [4nx by M]

  return psi * input_deck.angle.w;
}

Eigen::MatrixXd Solver::sourceIterate(double epsilon) {
  // solve transport equqation in all groups via source iteration

  LDCSD_LOG_INFO("Begin source iteration", true);

  // initial guess (maybe provided)
  Eigen::MatrixXd phi = Eigen::MatrixXd::Zero(4 * input_deck.mesh.n_x, input_deck.energy.G);
  Eigen::VectorXd phi_g = Eigen::VectorXd::Zero(4 * input_deck.mesh.n_x);

  Eigen::MatrixXd psi = Eigen::MatrixXd::Zero(4 * input_deck.mesh.n_x, input_deck.angle.M);
  Eigen::MatrixXd psi_up = Eigen::MatrixXd::Zero(4 * input_deck.mesh.n_x, input_deck.angle.M);

  // for each E:
  for (int g = 0; g < input_deck.energy.G; g++) {
    LDCSD_LOG_INFO("Beginning group " + std::to_string(g));
    // while not converged:
    // while norm(phi_latest - phi_old) > norm(phi_latest)epsilon
    auto i = 0;
    double abs_diff_norm = 1e10;
    while (abs_diff_norm > (phi_g.norm() * epsilon)) { //  TODO calculate these
      i++;
      phi.col(g) = phi_g;
      // solve transport using known phi
      psi = transportSweep(g, psi_up, phi);

      // Independent check: re-derive the discretized equation the sweep
      // just solved, in extended precision, using the exact psi/phi it was
      // solved with (not the updated phi_g below) -- should be near-zero
      // regardless of whether the source itself is physically correct,
      // since this only checks that solveDirect's assembly matches the
      // governing equations, not that the equations model the right
      // problem.
      Eigen::MatrixXd residuals = calculateResiduals(g, psi, psi_up, phi, true);
      Eigen::Index max_row, max_col, min_row, min_col;
      double max_residual = residuals.cwiseAbs().maxCoeff(&max_row, &max_col);
      double min_residual = residuals.cwiseAbs().minCoeff(&min_row, &min_col);

      int max_cell = int(max_row) / 4;
      int cv_index = max_row - (4 * max_cell);
      LDCSD_LOG_INFO("Group " + std::to_string(g) + ", iteration " + std::to_string(i) +
                     " complete.");
      LDCSD_LOG_INFO("\t\tResiduals:");
      LDCSD_LOG_INFO("\t\t-------------------");
      LDCSD_LOG_INFO("\t\t\tmax = " + std::format("{:.4e}", max_residual) + " at (" +
                     std::to_string(max_row) + ", " + std::to_string(max_col) +
                     ") (eqn, ordinate)");
      LDCSD_LOG_INFO("\t\t\t\tCell " + std::to_string(max_cell) + ", Corner value " +
                     std::to_string(cv_index));
      LDCSD_LOG_INFO("\t\t\tmin = " + std::format("{:.4e}", min_residual) + " at (" +
                     std::to_string(min_row) + ", " + std::to_string(min_col) + ")\n");

      // compute new phi
      phi_g = integrateAngle(psi);

      abs_diff_norm = (phi_g - phi.col(g)).norm();
    }
    phi.col(g) = phi_g;
    LDCSD_LOG_INFO("Converged with abs. norm = " + std::format("{:.4e}", abs_diff_norm) + "in " +
                       std::to_string(i) + " iterations",
                   true);
    psi_up = psi;
  }
  return phi;
}
Eigen::Vector<HighPrecision, 4>
Solver::cellResidual(double mu, double dx, double dE, double sigma_t, double S_bar, double S_Eg,
                     double S_Egm1, const Eigen::Vector2d& psi_gm1_d,
                     const Eigen::Vector2d& psi_b_u, const Eigen::Vector2d& psi_b_d,
                     const Eigen::Vector2d& q_u, const Eigen::Vector2d& q_d,
                     const Eigen::VectorXd& sigma_sdEprime, const Eigen::MatrixXd& phi_gprime_u,
                     const Eigen::MatrixXd& phi_gprime_d, const Eigen::Vector2d& psi_up,
                     const Eigen::Vector2d& psi_down, bool verbose) const {
  using Matrix2hp = Eigen::Matrix<HighPrecision, 2, 2>;
  using Vector2hp = Eigen::Vector<HighPrecision, 2>;
  using Vector4hp = Eigen::Vector<HighPrecision, 4>;
  using VectorXhp = Eigen::Matrix<HighPrecision, Eigen::Dynamic, 1>;
  using MatrixXhp = Eigen::Matrix<HighPrecision, Eigen::Dynamic, Eigen::Dynamic>;

  // M, L, Lb exactly as written in the method derivation -- re-derived
  // independently of Kernel's (double-precision) copies; see the
  // class-declaration comment for why.
  Matrix2hp M_hp;
  M_hp << 2, 1, 1, 2;
  M_hp *= HighPrecision(1) / 6;

  Matrix2hp L_hp;
  L_hp << 1, 1, -1, -1;
  L_hp *= HighPrecision(1) / 2;

  Matrix2hp Lb_hp;
  Lb_hp << -1, 0, 0, 1;

  // Promote every scalar input. mu keeps its sign (see the declaration
  // comment) -- no std::abs, no L/R swap.
  const HighPrecision mu_hp = mu;
  const HighPrecision dx_hp = dx;
  const HighPrecision dE_hp = dE;
  const HighPrecision sigma_t_hp = sigma_t;
  const HighPrecision S_bar_hp = S_bar;
  const HighPrecision S_Eg_hp = S_Eg;
  const HighPrecision S_Egm1_hp = S_Egm1;

  // Promote every vector/matrix input.
  const Vector2hp psi_gm1_d_hp = psi_gm1_d.cast<HighPrecision>();
  const Vector2hp psi_b_u_hp = psi_b_u.cast<HighPrecision>();
  const Vector2hp psi_b_d_hp = psi_b_d.cast<HighPrecision>();
  const Vector2hp q_u_hp = q_u.cast<HighPrecision>();
  const Vector2hp q_d_hp = q_d.cast<HighPrecision>();
  const VectorXhp sigma_sdEprime_hp = sigma_sdEprime.cast<HighPrecision>();
  const MatrixXhp phi_gprime_u_hp = phi_gprime_u.cast<HighPrecision>();
  const MatrixXhp phi_gprime_d_hp = phi_gprime_d.cast<HighPrecision>();

  // candidate is laid out like everywhere else in Solver: [up_left,
  // up_right, down_left, down_right] = [Psi_u,L Psi_u,R Psi_d,L Psi_d,R].
  const Vector2hp Psi_u_hp = psi_up.cast<HighPrecision>();
  const Vector2hp Psi_d_hp = psi_down.cast<HighPrecision>();

  // Scattering source is the same expression in both equations, just with
  // opposite external-source weighting below -- computed once.
  const Vector2hp scatter_source =
      (dx_hp / 8) * (M_hp * ((phi_gprime_d_hp + phi_gprime_u_hp) * sigma_sdEprime_hp));

  // ---- eq. 41a ("u" edge; candidate rows 0,1) ----
  // No dx on the streaming term (writeup typo).
  const Vector2hp streaming_u =
      (mu_hp / 6) * (Lb_hp * (psi_b_d_hp + 2 * psi_b_u_hp) + L_hp * (Psi_d_hp + 2 * Psi_u_hp));
  const Vector2hp absorption_csd_loss_u =
      dx_hp * (sigma_t_hp / 6 + S_bar_hp / (2 * dE_hp)) * (M_hp * Psi_d_hp) +
      dx_hp * (sigma_t_hp / 3 + S_bar_hp / (2 * dE_hp)) * (M_hp * Psi_u_hp);
  const Vector2hp csd_source_u = (dx_hp / dE_hp) * S_Egm1_hp * (M_hp * psi_gm1_d_hp);
  const Vector2hp external_source_u = (dx_hp / 6) * (M_hp * (q_d_hp + 2 * q_u_hp));

  const Vector2hp residual_u =
      streaming_u + absorption_csd_loss_u - csd_source_u - scatter_source - external_source_u;

  // ---- eq. 41b ("d" edge; candidate rows 2,3) ----
  // No dx on the streaming term (writeup typo); sigma_t's 1/6 and 1/3
  // swapped from how the writeup prints them (writeup typo) -- the "u"
  // coefficient is 1/3 and the "d" coefficient is 1/6, same pattern as
  // eq. 41a, matching Kernel::solveDirect's existing xs/3, xs/6 split.
  const Vector2hp streaming_d =
      (mu_hp / 6) * (Lb_hp * (2 * psi_b_d_hp + psi_b_u_hp) + L_hp * (2 * Psi_d_hp + Psi_u_hp));
  const Vector2hp absorption_csd_loss_d =
      dx_hp * (sigma_t_hp / 3 + (S_Eg_hp - S_bar_hp / 2) / dE_hp) * (M_hp * Psi_d_hp) +
      dx_hp * (sigma_t_hp / 6 - S_bar_hp / (2 * dE_hp)) * (M_hp * Psi_u_hp);
  // No CSD source here -- unlike eq. 41a, eq. 41b has no dependence on the
  // previous group.
  const Vector2hp external_source_d = (dx_hp / 6) * (M_hp * (2 * q_d_hp + q_u_hp));

  const Vector2hp residual_d =
      streaming_d + absorption_csd_loss_d - scatter_source - external_source_d;

  if (verbose) {
    auto logTerm = [](const std::string& name, const Vector2hp& v) {
      LDCSD_LOG_INFO("\t\t\t\t" + name + " = [" + std::format("{:.6e}", static_cast<double>(v(0))) +
                     ", " + std::format("{:.6e}", static_cast<double>(v(1))) + "]");
    };
    LDCSD_LOG_INFO("\t\t\teq. 41a (\"u\" edge):");
    logTerm("streaming_u          ", streaming_u);
    logTerm("absorption_csd_loss_u", absorption_csd_loss_u);
    logTerm("csd_source_u         ", csd_source_u);
    logTerm("scatter_source       ", scatter_source);
    logTerm("external_source_u    ", external_source_u);
    logTerm("residual_u           ", residual_u);
    LDCSD_LOG_INFO("\t\t\teq. 41b (\"d\" edge):");
    logTerm("streaming_d          ", streaming_d);
    logTerm("absorption_csd_loss_d", absorption_csd_loss_d);
    logTerm("scatter_source       ", scatter_source);
    logTerm("external_source_d    ", external_source_d);
    logTerm("residual_d           ", residual_d);
  }

  Vector4hp result;
  result << residual_u, residual_d;
  return result;
}

Eigen::MatrixXd Solver::calculateResiduals(int g, const Eigen::MatrixXd& angular,
                                           const Eigen::MatrixXd& psi_gm1,
                                           const Eigen::MatrixXd& scalar, bool debug_max) {
  // angular, psi_gm1: (4nx, M), this group's psi and the previous group's
  // (zero matrix for g==0). scalar: (4nx, G), all groups.
  int L = 0;
  int R = 1;

  int nx = input_deck.mesh.n_x;
  int M = input_deck.angle.M;

  auto mu = input_deck.angle.mu;
  auto dx = input_deck.mesh.dx;
  auto dE = input_deck.energy.dE;
  auto xs = input_deck.xs.total;
  auto S = input_deck.xs.S;
  auto S_bound = input_deck.xs.S_bound;

  Eigen::MatrixXd residuals = Eigen::MatrixXd::Zero(4 * nx, M);

  Eigen::Vector2d psi_b_up, psi_b_down, psi_in_E;

  // Gathers cell (i, m)'s inputs and evaluates cellResidual there -- shared
  // by the main grid pass below and the debug_max re-evaluation, so the two
  // can't drift apart from each other.
  auto cellResidualAt = [&](int i, int m, bool verbose) -> Eigen::Vector4d {
    auto psi_up = angular(Eigen::seqN(4 * i, 2), m);
    auto psi_down = angular(Eigen::seqN(4 * i + 2, 2), m);

    auto q_up = input_deck.source.values[g](Eigen::seqN(4 * i, 2), m);
    auto q_down = input_deck.source.values[g](Eigen::seqN(4 * i + 2, 2), m);

    if (g == 0) {
      psi_in_E = Eigen::Vector2d::Zero();
    } else {
      psi_in_E = psi_gm1(Eigen::seqN(4 * i + 2, 2), m);
    }

    auto phi_gprime_up = scalar(Eigen::seqN(4 * i, 2), Eigen::placeholders::all);
    auto phi_gprime_down = scalar(Eigen::seqN(4 * i + 2, 2), Eigen::placeholders::all);

    auto sigma_s = input_deck.xs.scatter[i].col(g);

    // construct appropraite psi^b
    switch (mu[m] > 0) {
    case true:
      psi_b_up(R) = psi_up(R);
      psi_b_down(R) = psi_down(R);

      if (i == 0) {
        psi_b_up(L) = input_deck.bc[g](0, m);
        psi_b_down(L) = input_deck.bc[g](1, m);
      } else {
        psi_b_up(L) = angular((4 * (i - 1) + 1), m);
        psi_b_down(L) = angular((4 * (i - 1) + 3), m);
      }
      break;
    case false:
      psi_b_up(L) = psi_up(L);
      psi_b_down(L) = psi_down(L);

      if (i == nx - 1) {
        psi_b_up(R) = input_deck.bc[g](0, m);
        psi_b_down(R) = input_deck.bc[g](1, m);
      } else {
        psi_b_up(R) = angular((4 * (i + 1)), m);
        psi_b_down(R) = angular((4 * (i + 1) + 2), m);
      }
      break;
    }

    return cellResidual(mu[m], dx[i], dE[g], xs(g, i), S(g, i), S_bound(g + 1, i), S_bound(g, i),
                        psi_in_E, psi_b_up, psi_b_down, q_up, q_down, dE.cwiseProduct(sigma_s),
                        phi_gprime_up, phi_gprime_down, psi_up, psi_down, verbose)
        .cast<double>();
  };

  for (int m = 0; m < M; m++) {
    for (int i = 0; i < nx; i++) {
      residuals(Eigen::seqN(4 * i, 4), m) = cellResidualAt(i, m, false);
    }
  }

  if (debug_max) {
    Eigen::Index max_row, max_col;
    residuals.cwiseAbs().maxCoeff(&max_row, &max_col);
    int max_cell = static_cast<int>(max_row) / 4;
    int corner = static_cast<int>(max_row) - 4 * max_cell;
    int max_m = static_cast<int>(max_col);
    LDCSD_LOG_INFO("\t\tdebug_max: largest |residual| at group " + std::to_string(g) + ", cell " +
                   std::to_string(max_cell) + ", corner " + std::to_string(corner) +
                   ", ordinate " + std::to_string(max_m) + " (mu=" +
                   std::format("{:.4e}", mu[max_m]) + ") -- recomputing term-by-term:");
    cellResidualAt(max_cell, max_m, true);
  }

  return residuals;
}

Eigen::Vector4d Solver::Kernel::solveDirect(
    double cosine, double dx, double dE, double xs, double S, double S_up, double S_down,
    Eigen::Vector2d psi_in_E, double psi_in_x_down, double psi_in_x_up, Eigen::Vector2d q_up,
    Eigen::Vector2d q_down, const Eigen::VectorXd& sigma_sdEprime,
    const Eigen::MatrixXd& phi_gprime_up, const Eigen::MatrixXd& phi_gprime_down,
    bool check_condition) {
  A = Eigen::Matrix4d::Zero();
  b = Eigen::Vector4d::Zero();

  double mu = std::abs(cosine);

  // The element matrices below are built assuming flow travels L->R; for
  // cosine < 0 the L/R labeling of every spatially-structured input must be
  // swapped to match before assembly (scalar psi_in_x_up/down are already
  // direction-relative "upwind" values, so they're left alone).
  Eigen::MatrixXd phi_gprime_up_local = phi_gprime_up;
  Eigen::MatrixXd phi_gprime_down_local = phi_gprime_down;
  if (cosine < 0) {
    psi_in_E.reverseInPlace();
    q_up.reverseInPlace();
    q_down.reverseInPlace();
    phi_gprime_up_local = phi_gprime_up_local.colwise().reverse().eval();
    phi_gprime_down_local = phi_gprime_down_local.colwise().reverse().eval();
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
  b({0, 1}) += (dx / 8) * M * (phi_gprime_down_local + phi_gprime_up_local) * sigma_sdEprime;
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
  b({2, 3}) += (dx / 8) * M * (phi_gprime_down_local + phi_gprime_up_local) * sigma_sdEprime;
  // External source
  b({2, 3}) += (dx / 6) * M * (q_up + 2 * q_down);

  if (check_condition) {
    Eigen::JacobiSVD<Eigen::Matrix4d> svd(A);
    const Eigen::Vector4d& singular_values = svd.singularValues();
    double condition_number = singular_values(0) / singular_values(singular_values.size() - 1);
    LDCSD_LOG_INFO("solveDirect: condition number = " + std::format("{:.4e}", condition_number));
  }

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

void Solver::writeResults(const std::filesystem::path& file_path,
                          const Eigen::MatrixXd& scalar_flux,
                          const std::vector<Eigen::MatrixXd>& angular_flux) const {
  appendToFile(file_path, SolverFormatter::formatResults(scalar_flux, angular_flux, input_deck));
}

void Solver::writeResiduals(const std::filesystem::path& file_path,
                            const std::vector<Eigen::MatrixXd>& residuals) const {
  appendToFile(file_path, SolverFormatter::formatResiduals(residuals, input_deck));
}
