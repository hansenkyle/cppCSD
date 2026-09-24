// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "smm.h"

#include "logger.h"

#include <Eigen/Dense>
#include <format>
#include <stdexcept>

Eigen::MatrixXd SMMResult::cell_average_current() const {
  using Eigen::seqN;
  using Eigen::placeholders::all;
  int I = scalar_flux.rows() / 4;
  int G = scalar_flux.cols();

  Eigen::MatrixXd leftsum = (current(seqN(0, I, 4), all) + current(seqN(2, I, 4), all));
  Eigen::MatrixXd rightsum = (current(seqN(1, I, 4), all) + current(seqN(3, I, 4), all));
  Eigen::MatrixXd result = (leftsum + rightsum) / 4;
  return result;
}

SecondMoment::SecondMoment(InputDeck input_deck)
    : Method("second moment method", input_deck), transport_operator(input_deck),
      convergence_(input_deck.energy.G) {

  J_in_positive = std::vector<Eigen::MatrixXd>(input_deck.energy.G);
  J_in_negative = J_in_positive;
  phi_in_positive = J_in_positive;
  phi_in_negative = J_in_positive;

  auto& bc = input_deck.bc.values;
  auto& mu = input_deck.angle.mu;
  auto& w_pos = input_deck.angle.w_positive;
  auto& w_neg = input_deck.angle.w_negative;
  for (int g = 0; g < input_deck.energy.G; g++) {
    phi_in_positive[g] = transport_operator.integrateAngle(bc, w_pos);
    phi_in_negative[g] = transport_operator.integrateAngle(bc, w_neg);
    J_in_positive[g] = transport_operator.integrateAngle(bc, mu.cwiseProduct(w_pos));
    J_in_negative[g] = transport_operator.integrateAngle(bc, mu.cwiseProduct(w_neg));
  }
}

Eigen::VectorXd SecondMoment::calculateK(Eigen::MatrixXd psi_slice, int sign) const {
  // calculates
  // K = 0.25 sum_m (w_m (1-2|mu_m|) psi)
  // psi_slice must have M columns

  auto& mu = input_deck.angle.mu;
  return sign * 0.25 * psi_slice *
         input_deck.angle.w.cwiseProduct((2 * mu.cwiseAbs()) - Eigen::VectorXd::Ones(mu.size()));
}

Eigen::VectorXd SecondMoment::calculateT(Eigen::MatrixXd psi_slice, int sign) const {
  // calculates
  // if sign = -1:
  // T = 0.5 sum_m (w_m ( (3mu_m / 2) - (mu_m / abs(mu_m))) psi_slice)
  // if sign = +1:
  // T = 0.5 sum_m (w_m (-(3mu_m / 2) + (mu_m / abs(mu_m))) psi_slice)

  auto& mu = input_deck.angle.mu;
  auto& w = input_deck.angle.w;
  // (sign) * (mu/abs(mu) - 3mu/2)
  auto mu_term = sign * ((mu.cwiseProduct(mu.cwiseAbs().cwiseInverse())) - (1.5 * mu));

  return 0.5 * psi_slice * mu_term.cwiseProduct(w);
}

Eigen::VectorXd SecondMoment::solveSM(Eigen::VectorXd Kpos, Eigen::VectorXd Kneg,
                                      Eigen::VectorXd Tpos, Eigen::VectorXd Tneg) const {

  // solves for scalar flux and current given closure terms K and T.
  int I = input_deck.mesh.n_x;
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(8 * I, 8 * I);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(8 * I);

  // balance, LU (53a)

  // balance, RU (53b)

  // balance, LD (53c)

  // balance, RD (53d)

  // first-moment, LU (53e)

  // first-moment, RU (53f)

  // first-moment, LD (53g)

  // first-moment, RD (53h)
}

void SecondMoment::solve(double /*epsilon*/, int /*max_iterations*/) {
  // Placeholder so SecondMoment is constructible (e.g. by the residual unit tests).
  throw std::logic_error("SecondMoment::solve is not implemented yet");
}

Eigen::VectorXd SecondMoment::calculateResiduals(int g, const Eigen::MatrixXd& scalar,
                                                 const Eigen::MatrixXd& current,
                                                 const Eigen::MatrixXd& psi, bool debug_max) {
  using Eigen::seqN;
  using Eigen::placeholders::all;

  int u = 0;
  int d = 1;

  int nx = input_deck.mesh.n_x;

  auto mu = input_deck.angle.mu;
  auto w = input_deck.angle.w;
  auto dx = input_deck.mesh.dx;
  auto dE = input_deck.energy.dE;
  InputDeck::Xs& xs = input_deck.xs;

  Eigen::VectorXd residuals = Eigen::VectorXd::Zero(8 * nx);

  Eigen::VectorXd phi = scalar.col(g);
  Eigen::VectorXd J = current.col(g);

  // quadrature weights restricted to each half-range
  Eigen::VectorXd w_pos = (mu.array() > 0).select(w, 0.0);
  Eigen::VectorXd w_neg = (mu.array() < 0).select(w, 0.0);
  Eigen::VectorXd third_minus_mu2 = (1.0 / 3) - mu.array().square();

  // closures at every corner of group g : [4nx]
  Eigen::VectorXd F = psi * w.cwiseProduct(third_minus_mu2);         // (54)
  Eigen::VectorXd F_pos = psi * w_pos.cwiseProduct(third_minus_mu2); // (56a)
  // (56b) typo corrected: PDF sums m = M/2+1..M (mu > 0); F^- is the mu < 0 half-range.
  Eigen::VectorXd F_neg = psi * w_neg.cwiseProduct(third_minus_mu2);
  Eigen::VectorXd K_pos = calculateK(psi, +1); // (58b)
  Eigen::VectorXd K_neg = calculateK(psi, -1); // (58a)
  Eigen::VectorXd T_pos = calculateT(psi, +1); // (58d)
  Eigen::VectorXd T_neg = calculateT(psi, -1); // (58c)

  // Half-range moments of the incoming boundary flux, (63)-(64), indexed (u, d). The same bc[g]
  // feeds both boundaries: mu > 0 enters at x_0, mu < 0 at x_I. F_in is not in the notes, but
  // the boundary faces' F^b needs it just like (55) does on interior faces.
  Eigen::MatrixXd bc = input_deck.bc[g];
  Eigen::Vector2d J_in_pos = bc * mu.cwiseProduct(w_pos);
  Eigen::Vector2d J_in_neg = bc * mu.cwiseProduct(w_neg);
  Eigen::Vector2d phi_in_pos = bc * w_pos;
  Eigen::Vector2d phi_in_neg = bc * w_neg;
  Eigen::Vector2d F_in_pos = bc * w_pos.cwiseProduct(third_minus_mu2);
  Eigen::Vector2d F_in_neg = bc * w_neg.cwiseProduct(third_minus_mu2);

  // construct phi^b, J^b, F^b on every face x_0..x_I : [nx+1 x 2], columns (u, d).
  // Each is the mu > 0 part carried out of the cell left of the face (its R corner) plus the
  // mu < 0 part carried out of the cell right of it (its L corner): (55) and (57). On an outer
  // face the missing side is the incoming boundary moment instead, giving (59)-(62).
  Eigen::MatrixXd phi_b(nx + 1, 2), J_b(nx + 1, 2), F_b(nx + 1, 2);
  for (int k = 0; k <= nx; k++) {
    for (int l : {u, d}) {
      // mu > 0 part
      if (k == 0) {
        J_b(k, l) = J_in_pos(l);
        phi_b(k, l) = phi_in_pos(l);
        F_b(k, l) = F_in_pos(l);
      } else {
        int R = 4 * (k - 1) + 2 * l + 1; // R corner of the cell left of x_k
        J_b(k, l) = 0.5 * J(R) + 0.25 * phi(R) + K_pos(R);
        phi_b(k, l) = 0.75 * J(R) + 0.5 * phi(R) + T_pos(R);
        F_b(k, l) = F_pos(R);
      }
      // mu < 0 part
      if (k == nx) {
        J_b(k, l) += J_in_neg(l);
        phi_b(k, l) += phi_in_neg(l);
        F_b(k, l) += F_in_neg(l);
      } else {
        int L = 4 * k + 2 * l; // L corner of the cell right of x_k
        J_b(k, l) += 0.5 * J(L) - 0.25 * phi(L) + K_neg(L);
        phi_b(k, l) += -0.75 * J(L) + 0.5 * phi(L) + T_neg(L);
        F_b(k, l) += F_neg(L);
      }
    }
  }

  Eigen::Vector2d phi_gm1_d, J_gm1_d;

  // Gathers cell i's inputs and evaluates cellResidual there -- shared by the main pass below
  // and the debug_max re-evaluation, so the two can't drift apart from each other.
  auto cellResidualAt = [&](int i, bool verbose) -> Eigen::Vector<double, 8> {
    auto phi_up = phi(seqN(4 * i, 2));
    auto phi_down = phi(seqN(4 * i + 2, 2));
    auto J_up = J(seqN(4 * i, 2));
    auto J_down = J(seqN(4 * i + 2, 2));
    auto F_up = F(seqN(4 * i, 2));
    auto F_down = F(seqN(4 * i + 2, 2));

    // cell i's faces x_{i-1}, x_i are rows i, i+1 of the face arrays
    auto phi_b_up = phi_b(seqN(i, 2), u);
    auto phi_b_down = phi_b(seqN(i, 2), d);
    auto J_b_up = J_b(seqN(i, 2), u);
    auto J_b_down = J_b(seqN(i, 2), d);
    auto F_b_up = F_b(seqN(i, 2), u);
    auto F_b_down = F_b(seqN(i, 2), d);

    auto q0_up = input_deck.source.q0[g](seqN(4 * i, 2));
    auto q0_down = input_deck.source.q0[g](seqN(4 * i + 2, 2));
    auto q1_up = input_deck.source.q1[g](seqN(4 * i, 2));
    auto q1_down = input_deck.source.q1[g](seqN(4 * i + 2, 2));

    // (48): nothing enters the highest-energy group from above
    if (g == 0) {
      phi_gm1_d = Eigen::Vector2d::Zero();
      J_gm1_d = Eigen::Vector2d::Zero();
    } else {
      phi_gm1_d = scalar(seqN(4 * i + 2, 2), g - 1);
      J_gm1_d = current(seqN(4 * i + 2, 2), g - 1);
    }

    auto phi_gprime_up = scalar(seqN(4 * i, 2), all);
    auto phi_gprime_down = scalar(seqN(4 * i + 2, 2), all);

    auto sigma_s = xs.scatter(i).col(g);

    return cellResidual(dx[i], dE[g], xs.total(g, i), xs.S(g, i), xs.S_down(g, i), xs.S_up(g, i),
                        phi_gm1_d, J_gm1_d, phi_b_up, phi_b_down, J_b_up, J_b_down, F_b_up,
                        F_b_down, F_up, F_down, q0_up, q0_down, q1_up, q1_down,
                        dE.cwiseProduct(sigma_s), phi_gprime_up, phi_gprime_down, phi_up, phi_down,
                        J_up, J_down, verbose)
        .cast<double>();
  };

  for (int i = 0; i < nx; i++) {
    residuals(seqN(8 * i, 8)) = cellResidualAt(i, false);
  }

  if (debug_max) {
    Eigen::Index max_row;
    residuals.cwiseAbs().maxCoeff(&max_row);
    int max_cell = static_cast<int>(max_row) / 8;
    int equation = static_cast<int>(max_row) - 8 * max_cell;
    LDCSD_LOG_INFO("\t\tdebug_max: largest |SM residual| at group " + std::to_string(g) +
                   ", cell " + std::to_string(max_cell) + ", eq. 53" +
                   static_cast<char>('a' + equation) + " -- recomputing term-by-term:");
    cellResidualAt(max_cell, true);
  }

  return residuals;
}

Eigen::Vector<HighPrecision, 8> SecondMoment::cellResidual(
    double dx, double dE, double sigma_t, double S_bar, double S_Eg, double S_Egm1,
    const Eigen::Vector2d& phi_gm1_d, const Eigen::Vector2d& J_gm1_d,
    const Eigen::Vector2d& phi_b_u, const Eigen::Vector2d& phi_b_d, const Eigen::Vector2d& J_b_u,
    const Eigen::Vector2d& J_b_d, const Eigen::Vector2d& F_b_u, const Eigen::Vector2d& F_b_d,
    const Eigen::Vector2d& F_u, const Eigen::Vector2d& F_d, const Eigen::Vector2d& q0_u,
    const Eigen::Vector2d& q0_d, const Eigen::Vector2d& q1_u, const Eigen::Vector2d& q1_d,
    const Eigen::VectorXd& sigma_sdEprime, const Eigen::MatrixXd& phi_gprime_u,
    const Eigen::MatrixXd& phi_gprime_d, const Eigen::Vector2d& phi_u, const Eigen::Vector2d& phi_d,
    const Eigen::Vector2d& J_u, const Eigen::Vector2d& J_d, bool verbose) const {
  using Matrix2hp = Eigen::Matrix<HighPrecision, 2, 2>;
  using Vector2hp = Eigen::Vector<HighPrecision, 2>;
  using Vector8hp = Eigen::Vector<HighPrecision, 8>;
  using VectorXhp = Eigen::Matrix<HighPrecision, Eigen::Dynamic, 1>;
  using MatrixXhp = Eigen::Matrix<HighPrecision, Eigen::Dynamic, Eigen::Dynamic>;

  // Each Vector2hp below is a (b_L, b_R)-weighted pair of equations, i.e. rows (53a, 53b),
  // (53c, 53d), (53e, 53f) or (53g, 53h). Face pairs are (x_{i-1}, x_i), corner pairs (L, R).

  Matrix2hp M_hp;
  M_hp << 2, 1, 1, 2;
  M_hp *= HighPrecision(1) / 6;

  Matrix2hp L_hp;
  L_hp << 1, 1, -1, -1;
  L_hp *= HighPrecision(1) / 2;

  Matrix2hp Lb_hp;
  Lb_hp << -1, 0, 0, 1;

  // Promote every scalar input
  const HighPrecision dx_hp = dx;
  const HighPrecision dE_hp = dE;
  const HighPrecision sigma_t_hp = sigma_t;
  const HighPrecision S_bar_hp = S_bar;
  const HighPrecision S_Eg_hp = S_Eg;
  const HighPrecision S_Egm1_hp = S_Egm1;

  // Promote every vector/matrix input
  const Vector2hp phi_gm1_d_hp = phi_gm1_d.cast<HighPrecision>();
  const Vector2hp J_gm1_d_hp = J_gm1_d.cast<HighPrecision>();
  const Vector2hp phi_b_u_hp = phi_b_u.cast<HighPrecision>();
  const Vector2hp phi_b_d_hp = phi_b_d.cast<HighPrecision>();
  const Vector2hp J_b_u_hp = J_b_u.cast<HighPrecision>();
  const Vector2hp J_b_d_hp = J_b_d.cast<HighPrecision>();
  const Vector2hp F_b_u_hp = F_b_u.cast<HighPrecision>();
  const Vector2hp F_b_d_hp = F_b_d.cast<HighPrecision>();
  const Vector2hp F_u_hp = F_u.cast<HighPrecision>();
  const Vector2hp F_d_hp = F_d.cast<HighPrecision>();
  const Vector2hp q0_u_hp = q0_u.cast<HighPrecision>();
  const Vector2hp q0_d_hp = q0_d.cast<HighPrecision>();
  const Vector2hp q1_u_hp = q1_u.cast<HighPrecision>();
  const Vector2hp q1_d_hp = q1_d.cast<HighPrecision>();
  const VectorXhp sigma_sdEprime_hp = sigma_sdEprime.cast<HighPrecision>();
  const MatrixXhp phi_gprime_u_hp = phi_gprime_u.cast<HighPrecision>();
  const MatrixXhp phi_gprime_d_hp = phi_gprime_d.cast<HighPrecision>();

  const Vector2hp Phi_u_hp = phi_u.cast<HighPrecision>();
  const Vector2hp Phi_d_hp = phi_d.cast<HighPrecision>();
  const Vector2hp J_u_hp = J_u.cast<HighPrecision>();
  const Vector2hp J_d_hp = J_d.cast<HighPrecision>();

  // absorption + CSD loss coefficients of (42) and (46); c_ud multiplies the d unknown in the
  // u-weighted equation, etc.
  const HighPrecision c_uu = sigma_t_hp / 3 + S_bar_hp / (2 * dE_hp);
  const HighPrecision c_ud = sigma_t_hp / 6 + S_bar_hp / (2 * dE_hp);
  const HighPrecision c_du = sigma_t_hp / 6 - S_bar_hp / (2 * dE_hp);
  const HighPrecision c_dd = sigma_t_hp / 3 + (S_Eg_hp - S_bar_hp / 2) / dE_hp;

  // Scattering source is the same expression in both balance equations -- computed once. The
  // first-moment equations have none: ldcsd's scattering is isotropic, so sigma_{s1} = 0.
  const Vector2hp scatter_source =
      (dx_hp / 4) * (M_hp * ((phi_gprime_d_hp + phi_gprime_u_hp) * sigma_sdEprime_hp));

  // ---- Balance equations ----

  // UP equations (53a, 53b) -- balance (42a) weighted by b_L, b_R and integrated.
  // typo: the notes' text before 53b says (42b).
  // typo: the notes' interior streaming term prints J_{g,i,d,u,R}; it is J_{g,i,u,R}.
  const Vector2hp streaming_balance_u =
      (Lb_hp * (J_b_d_hp + 2 * J_b_u_hp) + L_hp * (J_d_hp + 2 * J_u_hp)) / 6;
  const Vector2hp absorption_csd_loss_balance_u =
      dx_hp * c_ud * (M_hp * Phi_d_hp) + dx_hp * c_uu * (M_hp * Phi_u_hp);
  const Vector2hp csd_source_balance_u = (dx_hp / dE_hp) * S_Egm1_hp * (M_hp * phi_gm1_d_hp);
  const Vector2hp external_source_balance_u = (dx_hp / 6) * (M_hp * (q0_d_hp + 2 * q0_u_hp));

  const Vector2hp residual_balance_u = streaming_balance_u + absorption_csd_loss_balance_u -
                                       csd_source_balance_u - scatter_source -
                                       external_source_balance_u;

  // DOWN equations (53c, 53d) -- balance (42b) weighted by b_L, b_R and integrated.
  // typo: the notes' interior streaming term prints J_{g,i,d,u,R}; it is J_{g,i,u,R}.
  const Vector2hp streaming_balance_d =
      (Lb_hp * (2 * J_b_d_hp + J_b_u_hp) + L_hp * (2 * J_d_hp + J_u_hp)) / 6;
  const Vector2hp absorption_csd_loss_balance_d =
      dx_hp * c_dd * (M_hp * Phi_d_hp) + dx_hp * c_du * (M_hp * Phi_u_hp);
  const Vector2hp external_source_balance_d = (dx_hp / 6) * (M_hp * (2 * q0_d_hp + q0_u_hp));

  const Vector2hp residual_balance_d = streaming_balance_d + absorption_csd_loss_balance_d -
                                       scatter_source - external_source_balance_d;

  // ---- First-moment equations ----
  // typos common to 53e-h:
  //   - the interior streaming term repeats phi_d; the second pair is phi_u.
  //   - the mass terms print (2 J_{d,L} + J_{d,L}) etc.; the second corner is R, i.e. M * J.
  //   - the source prefactor is printed 1/6 with last term q^1_{u,L}; it is dx/36 with q^1_{u,R},
  //     i.e. (dx/6) * M * (...), exactly as in the balance equations.

  // UP equations (53e, 53f) -- first moment (46a) weighted by b_L, b_R and integrated.
  const Vector2hp streaming_first_moment_u =
      (Lb_hp * (phi_b_d_hp + 2 * phi_b_u_hp) + L_hp * (Phi_d_hp + 2 * Phi_u_hp)) / 18;
  const Vector2hp absorption_csd_loss_first_moment_u =
      dx_hp * c_ud * (M_hp * J_d_hp) + dx_hp * c_uu * (M_hp * J_u_hp);
  const Vector2hp closure_source_first_moment_u =
      (Lb_hp * (F_b_d_hp + 2 * F_b_u_hp) + L_hp * (F_d_hp + 2 * F_u_hp)) / 6;
  const Vector2hp csd_source_first_moment_u = (dx_hp / dE_hp) * S_Egm1_hp * (M_hp * J_gm1_d_hp);
  const Vector2hp external_source_first_moment_u = (dx_hp / 6) * (M_hp * (q1_d_hp + 2 * q1_u_hp));

  const Vector2hp residual_first_moment_u =
      streaming_first_moment_u + absorption_csd_loss_first_moment_u -
      closure_source_first_moment_u - csd_source_first_moment_u - external_source_first_moment_u;

  // DOWN equations (53g, 53h) -- first moment (46b) weighted by b_L, b_R and integrated.
  // typo: 53h prints F^b(x_{i-1}); the b_R row takes F^b(x_i), as Lb_hp does here.
  const Vector2hp streaming_first_moment_d =
      (Lb_hp * (2 * phi_b_d_hp + phi_b_u_hp) + L_hp * (2 * Phi_d_hp + Phi_u_hp)) / 18;
  const Vector2hp absorption_csd_loss_first_moment_d =
      dx_hp * c_dd * (M_hp * J_d_hp) + dx_hp * c_du * (M_hp * J_u_hp);
  const Vector2hp closure_source_first_moment_d =
      (Lb_hp * (2 * F_b_d_hp + F_b_u_hp) + L_hp * (2 * F_d_hp + F_u_hp)) / 6;
  const Vector2hp external_source_first_moment_d = (dx_hp / 6) * (M_hp * (2 * q1_d_hp + q1_u_hp));

  const Vector2hp residual_first_moment_d =
      streaming_first_moment_d + absorption_csd_loss_first_moment_d -
      closure_source_first_moment_d - external_source_first_moment_d;

  if (verbose) {
    auto logTerm = [](const std::string& name, const Vector2hp& v) {
      LDCSD_LOG_INFO("\t\t\t\t" + name + " = [" + std::format("{:.6e}", static_cast<double>(v(0))) +
                     ", " + std::format("{:.6e}", static_cast<double>(v(1))) + "]");
    };
    LDCSD_LOG_INFO("\t\t\teq. 53a,b (balance, \"u\" edge):");
    logTerm("streaming_balance_u               ", streaming_balance_u);
    logTerm("absorption_csd_loss_balance_u     ", absorption_csd_loss_balance_u);
    logTerm("csd_source_balance_u              ", csd_source_balance_u);
    logTerm("scatter_source                    ", scatter_source);
    logTerm("external_source_balance_u         ", external_source_balance_u);
    logTerm("residual_balance_u                ", residual_balance_u);
    LDCSD_LOG_INFO("\t\t\teq. 53c,d (balance, \"d\" edge):");
    logTerm("streaming_balance_d               ", streaming_balance_d);
    logTerm("absorption_csd_loss_balance_d     ", absorption_csd_loss_balance_d);
    logTerm("scatter_source                    ", scatter_source);
    logTerm("external_source_balance_d         ", external_source_balance_d);
    logTerm("residual_balance_d                ", residual_balance_d);
    LDCSD_LOG_INFO("\t\t\teq. 53e,f (first moment, \"u\" edge):");
    logTerm("streaming_first_moment_u          ", streaming_first_moment_u);
    logTerm("absorption_csd_loss_first_moment_u", absorption_csd_loss_first_moment_u);
    logTerm("closure_source_first_moment_u     ", closure_source_first_moment_u);
    logTerm("csd_source_first_moment_u         ", csd_source_first_moment_u);
    logTerm("external_source_first_moment_u    ", external_source_first_moment_u);
    logTerm("residual_first_moment_u           ", residual_first_moment_u);
    LDCSD_LOG_INFO("\t\t\teq. 53g,h (first moment, \"d\" edge):");
    logTerm("streaming_first_moment_d          ", streaming_first_moment_d);
    logTerm("absorption_csd_loss_first_moment_d", absorption_csd_loss_first_moment_d);
    logTerm("closure_source_first_moment_d     ", closure_source_first_moment_d);
    logTerm("external_source_first_moment_d    ", external_source_first_moment_d);
    logTerm("residual_first_moment_d           ", residual_first_moment_d);
  }

  Vector8hp result;
  result << residual_balance_u, residual_balance_d, residual_first_moment_u,
      residual_first_moment_d;
  return result;
}
