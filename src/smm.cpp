// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "smm.h"
#include "logger.h"
#include "output_block.h"

#include <Eigen/Dense>
#include <array>
#include <chrono>
#include <format>
#include <ranges>
#include <stdexcept>
#include <vector>

namespace {
using Matrix8d = Eigen::Matrix<double, 8, 8>;

// Kronecker product A (x) B. Left factor = energy (u, d) or equation/unknown type, right factor =
// space (L, R) or a 4x4 per-cell block, as in smm.md.
template <int Rows, int Cols>
Eigen::Matrix<double, 2 * Rows, 2 * Cols> kron(const Eigen::Matrix2d& A,
                                               const Eigen::Matrix<double, Rows, Cols>& B) {
  Eigen::Matrix<double, 2 * Rows, 2 * Cols> K;
  for (int r = 0; r < 2; r++) {
    for (int c = 0; c < 2; c++) {
      K.template block<Rows, Cols>(Rows * r, Cols * c) = A(r, c) * B;
    }
  }
  return K;
}

Eigen::Matrix2d mat2(double a, double b, double c, double d) {
  return (Eigen::Matrix2d() << a, b, c, d).finished();
}

// 2x2 building blocks and the data-free 4x4 operators of (53), named as in smm.md
const Eigen::Matrix2d kM = mat2(2, 1, 1, 2) / 6;
const Eigen::Matrix2d kL = mat2(1, 1, -1, -1) / 2;
const Eigen::Matrix2d kLb = mat2(-1, 0, 0, 1);
const Eigen::Matrix2d kN = mat2(0, 1, -1, 0) / 2;
const Eigen::Matrix2d kI2 = Eigen::Matrix2d::Identity();

const Eigen::Matrix4d kS = kron(kM, kL);               // interior streaming
const Eigen::Matrix4d kSb = kron(kM, kLb);             // face streaming
const Eigen::Matrix4d kX = kron(mat2(1, 1, 1, 1), kM); // scattering
const Eigen::Matrix4d kP = kron(mat2(0, 1, 0, 0), kM); // CSD inflow from g - 1
const Eigen::Matrix4d kQ = kron(kM, kM);               // external source

// The data-free part of the LO matrix's blocks: streaming, with the face values' dependence on
// the unknowns (57) substituted in. A0 is a cell's own block (still missing R and within-group
// scattering), A_minus couples it to cell i-1 and A_plus to cell i+1. Outer 2x2 indexing is
// (equation type, unknown type).
const Matrix8d kA0Streaming =
    kron(kI2, kron(kM, kI2)) / 4 + kron(mat2(0, 1, 1.0 / 3, 0), kron(kM, kN));
const Matrix8d kAMinus = kron(mat2(0.25, 0.5, 1.0 / 6, 0.25), kron(kM, mat2(0, -1, 0, 0)));
const Matrix8d kAPlus = kron(mat2(-0.25, 0.5, 1.0 / 6, -0.25), kron(kM, mat2(0, 0, 1, 0)));

// Cell i's 4 face values [u(x_{i-1}), u(x_i), d(x_{i-1}), d(x_i)] from a face array [nx+1 x 2].
Eigen::Vector4d cellFaces(const Eigen::MatrixXd& faces, int i) {
  return Eigen::Vector4d(faces(i, 0), faces(i + 1, 0), faces(i, 1), faces(i + 1, 1));
}
} // namespace

Eigen::MatrixXd SMMResult::cell_average_current() const { return cell_average(current); }

SecondMoment::SecondMoment(InputDeck input_deck)
    : Method("second moment method", input_deck), transport_operator(input_deck),
      closures(input_deck.energy.G) {}

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

SMClosures SecondMoment::computeClosures(const Eigen::MatrixXd& psi) const {
  const Eigen::VectorXd& mu = input_deck.angle.mu;
  const Eigen::VectorXd& w = input_deck.angle.w;

  // quadrature weights restricted to each half-range
  const Eigen::VectorXd w_pos = (mu.array() > 0).select(w, 0.0);
  const Eigen::VectorXd w_neg = (mu.array() < 0).select(w, 0.0);
  const Eigen::VectorXd third_minus_mu2 = (1.0 / 3) - mu.array().square();

  SMClosures closures;
  closures.F = psi * w.cwiseProduct(third_minus_mu2);
  closures.F_pos = psi * w_pos.cwiseProduct(third_minus_mu2);
  // (56b) typo corrected: PDF sums m = M/2+1..M (mu > 0); F^- is the mu < 0 half-range.
  closures.F_neg = psi * w_neg.cwiseProduct(third_minus_mu2);
  closures.K_pos = calculateK(psi, +1);
  closures.K_neg = calculateK(psi, -1);
  closures.T_pos = calculateT(psi, +1);
  closures.T_neg = calculateT(psi, -1);
  return closures;
}

SecondMoment::IncomingMoments SecondMoment::incomingMoments(int g) const {
  const Eigen::VectorXd& mu = input_deck.angle.mu;
  const Eigen::VectorXd& w = input_deck.angle.w;
  const Eigen::VectorXd w_pos = (mu.array() > 0).select(w, 0.0);
  const Eigen::VectorXd w_neg = (mu.array() < 0).select(w, 0.0);
  const Eigen::VectorXd third_minus_mu2 = (1.0 / 3) - mu.array().square();

  const Eigen::MatrixXd bc = input_deck.bc[g];
  IncomingMoments in;
  in.J_pos = bc * mu.cwiseProduct(w_pos);
  in.J_neg = bc * mu.cwiseProduct(w_neg);
  in.phi_pos = bc * w_pos;
  in.phi_neg = bc * w_neg;
  in.F_pos = bc * w_pos.cwiseProduct(third_minus_mu2);
  in.F_neg = bc * w_neg.cwiseProduct(third_minus_mu2);
  return in;
}

Eigen::MatrixXd SecondMoment::assembleFaces(const Eigen::VectorXd& pos, const Eigen::VectorXd& neg,
                                            const Eigen::Vector2d& in_pos,
                                            const Eigen::Vector2d& in_neg) const {
  const int nx = input_deck.mesh.n_x;
  Eigen::MatrixXd faces(nx + 1, 2);
  for (int k = 0; k <= nx; k++) {
    for (int l : {0, 1}) {
      // R corner of the cell left of x_k, and L corner of the cell right of it
      faces(k, l) = (k == 0 ? in_pos(l) : pos(4 * (k - 1) + 2 * l + 1)) +
                    (k == nx ? in_neg(l) : neg(4 * k + 2 * l));
    }
  }
  return faces;
}

Eigen::SparseMatrix<double> SecondMoment::buildGroupMatrix(int g) {
  const int nx = input_deck.mesh.n_x;
  const double dE = input_deck.energy.dE(g);
  InputDeck::Xs& xs = input_deck.xs;

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(96 * nx);
  auto addBlock = [&](int row_cell, int col_cell, const Matrix8d& block) {
    for (int c = 0; c < 8; c++) {
      for (int r = 0; r < 8; r++) {
        if (block(r, c) != 0.0) {
          triplets.emplace_back(8 * row_cell + r, 8 * col_cell + c, block(r, c));
        }
      }
    }
  };

  for (int i = 0; i < nx; i++) {
    const double dx = input_deck.mesh.dx(i);
    const double sigma_t = xs.total(g, i);
    const double S_bar = xs.S(g, i);

    // absorption + CSD loss, the same four coefficients solveDirect uses
    const Eigen::Matrix2d C =
        mat2(sigma_t / 3 + S_bar / (2 * dE), sigma_t / 6 + S_bar / (2 * dE),
             sigma_t / 6 - S_bar / (2 * dE), sigma_t / 3 + (xs.S_down(g, i) - S_bar / 2) / dE);
    const Eigen::Matrix4d R = dx * kron(C, kM);
    // within-group scattering, kept on the LHS; sigma_s1 = 0, so only the balance rows get it
    const double w0 = (dx / 4) * xs.scatter(g, g, i); // dE_g / dE_g = 1

    Matrix8d A0 = kA0Streaming;
    A0.topLeftCorner<4, 4>() += R - w0 * kX;
    A0.bottomRightCorner<4, 4>() += R;

    addBlock(i, i, A0);
    if (i > 0) {
      addBlock(i, i - 1, kAMinus);
    }
    if (i < nx - 1) {
      addBlock(i, i + 1, kAPlus);
    }
  }

  Eigen::SparseMatrix<double> A(8 * nx, 8 * nx);
  A.setFromTriplets(triplets.begin(), triplets.end());
  return A;
}

Eigen::VectorXd SecondMoment::buildGroupRHS(int g, const SMClosures& closures,
                                            const Eigen::MatrixXd& scalar,
                                            const Eigen::MatrixXd& current) {
  const int nx = input_deck.mesh.n_x;
  const double dE = input_deck.energy.dE(g);
  InputDeck::Xs& xs = input_deck.xs;
  const IncomingMoments in = incomingMoments(g);

  // closure (and incoming boundary) parts of the face values, (55) and (57)
  const Eigen::MatrixXd K_b = assembleFaces(closures.K_pos, closures.K_neg, in.J_pos, in.J_neg);
  const Eigen::MatrixXd T_b = assembleFaces(closures.T_pos, closures.T_neg, in.phi_pos, in.phi_neg);
  const Eigen::MatrixXd F_b = assembleFaces(closures.F_pos, closures.F_neg, in.F_pos, in.F_neg);

  Eigen::VectorXd b(8 * nx);
  for (int i = 0; i < nx; i++) {
    const double dx = input_deck.mesh.dx(i);

    // (48): nothing enters the highest-energy group from above
    Eigen::Vector4d phi_gm1 = Eigen::Vector4d::Zero();
    Eigen::Vector4d J_gm1 = Eigen::Vector4d::Zero();
    if (g > 0) {
      phi_gm1 = scalar.col(g - 1).segment<4>(4 * i);
      J_gm1 = current.col(g - 1).segment<4>(4 * i);
    }
    const double csd = (dx / dE) * xs.S_up(g, i);

    // scattering from every other group; g -> g is in the matrix
    Eigen::VectorXd sigma_sdEprime = xs.scatter(i).col(g).cwiseProduct(input_deck.energy.dE) / dE;
    sigma_sdEprime(g) = 0.0;
    const Eigen::Vector4d scatter = scalar.middleRows<4>(4 * i) * sigma_sdEprime;

    const Eigen::Vector4d q0 = input_deck.source.q0[g].segment<4>(4 * i);
    const Eigen::Vector4d q1 = input_deck.source.q1[g].segment<4>(4 * i);
    const Eigen::Vector4d F = closures.F.segment<4>(4 * i);

    // balance (53a-d)
    b.segment<4>(8 * i) =
        csd * kP * phi_gm1 + (dx / 4) * kX * scatter + dx * kQ * q0 - kSb * cellFaces(K_b, i);
    // first moment (53e-h)
    b.segment<4>(8 * i + 4) = kS * F + kSb * cellFaces(F_b, i) + csd * kP * J_gm1 + dx * kQ * q1 -
                              kSb * cellFaces(T_b, i) / 3;
  }
  return b;
}

void SecondMoment::factorizeGroup(int g) {
  lo_solver_.compute(buildGroupMatrix(g));
  if (lo_solver_.info() != Eigen::Success) {
    throw std::runtime_error("SecondMoment: LO matrix factorization failed in group " +
                             std::to_string(g) + ": " + lo_solver_.lastErrorMessage());
  }
  lo_factorized_ = true;
}

std::pair<Eigen::VectorXd, Eigen::VectorXd>
SecondMoment::solveGroup(const Eigen::VectorXd& rhs) const {
  if (!lo_factorized_) {
    throw std::logic_error("SecondMoment::solveGroup called before factorizeGroup");
  }
  const Eigen::VectorXd x = lo_solver_.solve(rhs);

  // x holds [phi_i; J_i] per cell; each column of cells is one cell's 8 unknowns
  const Eigen::Map<const Eigen::Matrix<double, 8, Eigen::Dynamic>> cells(x.data(), 8, x.size() / 8);
  return {cells.topRows<4>().reshaped(), cells.bottomRows<4>().reshaped()};
}

void SecondMoment::solve(double epsilon, int max_iterations) {
  // Solve the transport equation in all groups via the second moment method: identical to
  // source iteration, except each sweep's scalar flux comes from the LO (SM) equations with
  // closures computed from the sweep's angular flux, instead of from integrating it over angle.
  //
  // No iteration over energy groups, assume downscatter only

  int I = input_deck.mesh.n_x;
  int G = input_deck.energy.G;
  int M = input_deck.angle.M;

  LDCSD_LOG_INFO("Begin second moment method");

  Eigen::VectorXd phi_g = Eigen::VectorXd::Zero(4 * I);
  Eigen::MatrixXd psi_up = Eigen::MatrixXd::Zero(4 * I, M);

  // write phi and J directly into this->solution
  solution.scalar_flux = Eigen::MatrixXd::Zero(4 * I, G);
  auto& phi = solution.scalar_flux;
  solution.current = Eigen::MatrixXd::Zero(4 * I, G);
  auto& J = solution.current;
  solution.reconstructed_scalar = Eigen::MatrixXd::Zero(4 * I, G);

  solution.angular_flux =
      std::vector<Eigen::MatrixXd>(input_deck.energy.G, Eigen::MatrixXd::Zero(4 * I, M));
  residuals.high_order = solution.angular_flux;
  residuals.low_order = std::vector<Eigen::VectorXd>(G, Eigen::VectorXd::Zero(8 * I));

  for (int g = 0; g < G; g++) {
    LDCSD_LOG_INFO("Beginning group " + std::to_string(g));
    const auto group_start = std::chrono::steady_clock::now();

    // the LO matrix doesn't depend on the closures: factorize once for the whole group
    factorizeGroup(g);

    // write result directly into this->solution
    auto& psi = solution.angular_flux[g];

    int iteration = 0;
    double delta_phi_l2norm = 0.0;
    Eigen::VectorXd absolute_delta_phi;

    while (iteration < max_iterations) {
      iteration++;
      phi.col(g) = phi_g;
      // solve transport using known phi
      psi = transport_operator.sweep(g, psi_up, phi);
      // compute new phi (and J) from the SM equations
      closures[g] = computeClosures(psi);
      auto [phi_lo, J_lo] = solveGroup(buildGroupRHS(g, closures[g], phi, J));
      phi_g = phi_lo;
      J.col(g) = J_lo;

      absolute_delta_phi = (phi_g - phi.col(g));
      delta_phi_l2norm = l2norm(phi_g - phi.col(g));
      const double phi_norm = l2norm(phi_g);

      convergence_.log_group(g, IterationRecord({l2norm(phi_g), linfnorm(phi_g)},
                                                {delta_phi_l2norm, linfnorm(absolute_delta_phi)}));

      if (delta_phi_l2norm <= phi_norm * epsilon) {
        LDCSD_LOG_INFO("Converged with abs. norm = " + std::format("{:.4e}", delta_phi_l2norm) +
                       " in " + std::to_string(iteration) + " iterations");
        break; // exit while loop
      }
    }

    residuals.high_order[g] = transport_operator.calculateResiduals(g, psi, psi_up, phi);
    phi.col(g) = phi_g;
    residuals.low_order[g] = calculateResiduals(g, phi, J, psi);
    solution.reconstructed_scalar.col(g) = transport_operator.integrateAngle(psi);
    psi_up = psi;

    const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - group_start;
    convergence_.time_group(g, elapsed.count());

    if (iteration == max_iterations) {
      LDCSD_LOG_WARN("group " + std::to_string(g) + " did NOT converge: hit the " +
                     std::to_string(max_iterations) +
                     "-iteration cap with abs. norm = " + std::format("{:.4e}", delta_phi_l2norm));
    }
  }
}

Eigen::VectorXd SecondMoment::calculateResiduals(int g, const Eigen::MatrixXd& scalar,
                                                 const Eigen::MatrixXd& current,
                                                 const Eigen::MatrixXd& psi, bool debug_max) {
  using Eigen::seqN;
  using Eigen::placeholders::all;

  int u = 0;
  int d = 1;

  int nx = input_deck.mesh.n_x;

  auto dx = input_deck.mesh.dx;
  auto dE = input_deck.energy.dE;
  InputDeck::Xs& xs = input_deck.xs;

  Eigen::VectorXd residuals = Eigen::VectorXd::Zero(8 * nx);

  Eigen::VectorXd phi = scalar.col(g);
  Eigen::VectorXd J = current.col(g);

  const SMClosures c = computeClosures(psi);
  const Eigen::VectorXd& F = c.F;
  const IncomingMoments in = incomingMoments(g);

  // phi^b, J^b, F^b on every face x_0..x_I : [nx+1 x 2], columns (u, d), from (55) and (57);
  // (59)-(62) on the outer faces.
  const Eigen::MatrixXd J_b = assembleFaces(0.5 * J + 0.25 * phi + c.K_pos,
                                            0.5 * J - 0.25 * phi + c.K_neg, in.J_pos, in.J_neg);
  const Eigen::MatrixXd phi_b = assembleFaces(
      0.75 * J + 0.5 * phi + c.T_pos, -0.75 * J + 0.5 * phi + c.T_neg, in.phi_pos, in.phi_neg);
  const Eigen::MatrixXd F_b = assembleFaces(c.F_pos, c.F_neg, in.F_pos, in.F_neg);

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
                        dE.cwiseProduct(sigma_s) / dE[g], phi_gprime_up, phi_gprime_down, phi_up,
                        phi_down, J_up, J_down, verbose)
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

void SecondMoment::writeResults(const std::filesystem::path& results_path) const {
  using Eigen::seqN;
  const int I = input_deck.mesh.n_x;

  std::vector<std::string> iseq, x_center;
  for (int i = 0; i < I; i++) {
    iseq.push_back(std::to_string(i + 1));
    x_center.push_back(std::format("{:.4e}", input_deck.mesh.x_center(i)));
  }

  static const std::array<std::pair<const char*, Eigen::VectorXd SMClosures::*>, 7> kClosures = {{
      {"F", &SMClosures::F},
      {"F+", &SMClosures::F_pos},
      {"F-", &SMClosures::F_neg},
      {"K+", &SMClosures::K_pos},
      {"K-", &SMClosures::K_neg},
      {"T+", &SMClosures::T_pos},
      {"T-", &SMClosures::T_neg},
  }};

  UnitGroup close("closures");
  for (int g = 0; g < input_deck.energy.G; g++) {
    UnitGroup group("g = " + std::to_string(g + 1));
    for (const auto& [name, field] : kClosures) {
      const Eigen::VectorXd& c = closures[g].*field;
      HorizontalTable table(name);
      table.add_row("x_i", x_center);
      table.add_row("i", iseq);
      table.add_row("up,left", c(seqN(0, I, 4)));
      table.add_row("up,right", c(seqN(1, I, 4)));
      table.add_row("down,left", c(seqN(2, I, 4)));
      table.add_row("down,right", c(seqN(3, I, 4)));
      group.add(table, "{:.6e}");
    }
    close.add(group);
  }

  UnitGroup block = solutionBlock(solution);
  block.add(close);
  block.add(cellAverageTable("reconstructed cell-average scalar flux",
                             MethodResult::cell_average(solution.reconstructed_scalar)),
            "{:.6e}");
  appendToFile(results_path, block.render_txt());
}

void SecondMoment::writeResiduals(const std::filesystem::path& file_path, std::string timestamp) {
  using Eigen::seqN;
  using Eigen::placeholders::all;
  writeMetadata(file_path, timestamp);

  static constexpr std::array<const char*, 4> kEdgeLabels = {"up,L", "up,R", "down,L", "down,R"};
  static constexpr std::array<const char*, 8> kSMEqLabels = {
      "(balance, up L)",    "(balance, up R)",    "(balance, down L)",    "(balance, down R)",
      "(1st moment, up L)", "(1st moment, up R)", "(1st moment, down L)", "(1st moment, down R)"};

  // Largest |residual| over every group, its (cell, angle, edge), from a per-group Eigen
  // container whose rows are laid out cell-major in blocks of block_size (4 for transport, one
  // edge per row; 8 for the SM equations, one equation per row).
  auto peakResidual = [](const auto& per_group, int block_size) {
    struct Peak {
      double value = -1.0;
      int g = -1, i = -1, sub = -1, m = -1;
    } peak;
    for (int g = 0; g < static_cast<int>(per_group.size()); g++) {
      Eigen::Index row, col;
      const double gmax = per_group[g].cwiseAbs().maxCoeff(&row, &col);
      if (gmax > peak.value) {
        peak = {gmax, g, static_cast<int>(row) / block_size, static_cast<int>(row) % block_size,
                static_cast<int>(col)};
      }
    }
    return peak;
  };

  const auto transport_peak = peakResidual(residuals.high_order, 4);
  const std::string transport_summary =
      std::format("max |residual| = {:.4e} at g= {}, i= {}, angle {}, {}", transport_peak.value,
                  transport_peak.g + 1, transport_peak.i + 1, transport_peak.m + 1,
                  kEdgeLabels[transport_peak.sub]);
  LDCSD_LOG_INFO("transport residuals: " + transport_summary);

  // Re-evaluate the true peak (found above, across all groups) term-by-term.
  const Eigen::MatrixXd zero_psi_gm1 =
      Eigen::MatrixXd::Zero(4 * input_deck.mesh.n_x, input_deck.angle.M);
  transport_operator.calculateResiduals(
      transport_peak.g, solution.angular_flux[transport_peak.g],
      transport_peak.g == 0 ? zero_psi_gm1 : solution.angular_flux[transport_peak.g - 1],
      solution.scalar_flux, /*debug_max=*/true);

  UnitGroup transport("transport residuals", transport_summary);
  int I = input_deck.mesh.n_x;
  std::vector<std::string> x_i, mu_m, x_center;

  for (int i = 0; i < input_deck.mesh.n_x; i++) {
    x_i.push_back(std::to_string(i + 1));
    x_center.push_back(std::format("{:.5e}", input_deck.mesh.x_center(i)));
  }
  for (int m = 0; m < input_deck.angle.M; m++) {
    mu_m.push_back(std::to_string(m + 1));
  }
  for (int g = 0; g < input_deck.energy.G; g++) {
    UnitGroup group("g = " + std::to_string(g + 1));
    for (int m = 0; m < input_deck.angle.M; m++) {
      HorizontalTable angle("m = " + std::to_string(m + 1));
      angle.add_row("x_i", x_center);
      angle.add_row("i", x_i);
      angle.add_row("up,L", residuals.high_order[g](seqN(0, I, 4), m));
      angle.add_row("up,R", residuals.high_order[g](seqN(1, I, 4), m));
      angle.add_row("down,L", residuals.high_order[g](seqN(2, I, 4), m));
      angle.add_row("down,R", residuals.high_order[g](seqN(3, I, 4), m));
      group.add(angle, "{:.6e}");
    }
    transport.add(group);
  }

  appendToFile(file_path, transport.render_txt());

  const auto sm_peak = peakResidual(residuals.low_order, 8);
  const std::string sm_summary =
      std::format("max |residual| = {:.4e} at g= {}, i= {}, equation {}", sm_peak.value,
                  sm_peak.g + 1, sm_peak.i + 1, kSMEqLabels[sm_peak.sub]);
  LDCSD_LOG_INFO("second moment equation residuals: " + sm_summary);

  // Re-evaluate the true peak (found above, across all groups) term-by-term.
  calculateResiduals(sm_peak.g, solution.scalar_flux, solution.current,
                     solution.angular_flux[sm_peak.g], /*debug_max=*/true);

  UnitGroup low_order("second moment equation residuals", sm_summary);
  for (int g = 0; g < input_deck.energy.G; g++) {
    HorizontalTable group("g = " + std::to_string(g + 1));
    group.add_row("x_i", x_center);
    group.add_row("i", x_i);
    group.add_row("(balance, up L)", residuals.low_order[g](seqN(0, I, 8)));
    group.add_row("(balance, up R)", residuals.low_order[g](seqN(1, I, 8)));
    group.add_row("(balance, down L)", residuals.low_order[g](seqN(2, I, 8)));
    group.add_row("(balance, down R)", residuals.low_order[g](seqN(3, I, 8)));
    group.add_row("(1st moment, up L)", residuals.low_order[g](seqN(4, I, 8)));
    group.add_row("(1st moment, up R)", residuals.low_order[g](seqN(5, I, 8)));
    group.add_row("(1st moment, down L)", residuals.low_order[g](seqN(6, I, 8)));
    group.add_row("(1st moment, down R)", residuals.low_order[g](seqN(7, I, 8)));
    low_order.add(group, "{:.6e}");
  }

  appendToFile(file_path, low_order.render_txt());
}
