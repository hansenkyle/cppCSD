// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "smm.h"
#include "source_iteration.h"
#include "transport_operator.h"

#include <doctest.h>

namespace {
// 4 cells with a different material in each, 2 groups, and an S4 Gauss-Legendre quadrature (with
// only mu = +-0.5, 2|mu| - 1 = 0 and the K closure would drop out of every test).
// within_group_scatter = false zeroes the g -> g scattering but keeps the g0 -> g1 downscatter.
InputDeck makeSMDeck(bool within_group_scatter) {
  const int n_x = 4;
  const int G = 2;
  const int M = 4;

  InputDeck deck;
  deck.mesh.n_x = n_x;
  deck.mesh.x_boundary = Eigen::Vector<double, 5>(0.0, 0.6, 1.4, 1.9, 2.7);
  deck.energy.G = G;
  deck.energy.E_boundary = Eigen::Vector3d(2.0, 1.2, 0.5);
  deck.angle.M = M;
  deck.angle.mu = Eigen::Vector4d(-0.8611363115940526, -0.3399810435848563, 0.3399810435848563,
                                  0.8611363115940526);
  deck.angle.w = Eigen::Vector4d(0.3478548451374538, 0.6521451548625461, 0.6521451548625461,
                                 0.3478548451374538);

  std::vector<Material> materials;
  for (int i = 0; i < n_x; ++i) {
    Material mat;
    mat.name = "m" + std::to_string(i);
    mat.total = Eigen::Vector2d(1.3 - 0.2 * i, 0.7 + 0.15 * i);
    mat.S = Eigen::Vector2d(0.6 - 0.05 * i, 0.35 + 0.1 * i);
    mat.S_b = Eigen::Vector3d(0.7 + 0.03 * i, 0.5 - 0.04 * i, 0.3 + 0.05 * i);
    mat.scatter = Eigen::MatrixXd::Zero(G, G);
    mat.scatter(0, 1) = 0.3 + 0.02 * i; // downscatter g0 -> g1
    if (within_group_scatter) {
      mat.scatter(0, 0) = 0.4 + 0.05 * i;
      mat.scatter(1, 1) = 0.5 - 0.03 * i;
    }
    materials.push_back(mat);
  }
  deck.xs.set_materials(materials, {0, 1, 2, 3});

  // rows [u, d] per group, different in every ordinate so every half-range moment is distinct
  deck.bc.values = Eigen::MatrixXd(2 * G, M);
  deck.bc.values << 1.1, 0.9, 1.3, 0.8, //
      0.7, 1.2, 0.6, 1.0,               //
      0.4, 0.3, 0.5, 0.25,              //
      0.2, 0.45, 0.35, 0.15;

  // anisotropic, so q1 != 0
  deck.source.values.clear();
  for (int g = 0; g < G; ++g) {
    Eigen::MatrixXd q(4 * n_x, M);
    for (int r = 0; r < 4 * n_x; ++r) {
      for (int m = 0; m < M; ++m) {
        q(r, m) = (0.3 + 0.07 * ((3 * r + 5 * m) % 11)) * (1.0 + 0.25 * g);
      }
    }
    deck.source.values.push_back(q);
  }

  deck.validate();
  return deck;
}

// The exact LO solution of the discrete transport equations: its zeroth and first moments.
Eigen::VectorXd zerothMoment(const InputDeck& deck, const Eigen::MatrixXd& psi) {
  return psi * deck.angle.w;
}
Eigen::VectorXd firstMoment(const InputDeck& deck, const Eigen::MatrixXd& psi) {
  return psi * deck.angle.w.cwiseProduct(deck.angle.mu);
}

Eigen::Matrix4d kron(const Eigen::Matrix2d& A, const Eigen::Matrix2d& B) {
  Eigen::Matrix4d K;
  for (int r = 0; r < 2; ++r) {
    for (int c = 0; c < 2; ++c) {
      K.block<2, 2>(2 * r, 2 * c) = A(r, c) * B;
    }
  }
  return K;
}

} // namespace

TEST_CASE("SM residuals vanish on the moments of an exact transport solution") {
  // With no within-group scattering, two sweeps give the exact discrete transport solution in
  // both groups, so its moments must satisfy (53a-h) to machine precision -- in every term,
  // since the deck makes each one nonzero (boundaries, CSD coupling, downscatter, q0, q1, F, K,
  // T). The SM equations are moments of the transport equations, not a separate model.
  const InputDeck deck = makeSMDeck(false);
  const int rows = 4 * deck.mesh.n_x;
  TransportOperator transport(deck);
  SecondMoment sm(deck);

  Eigen::MatrixXd scalar = Eigen::MatrixXd::Zero(rows, deck.energy.G);
  Eigen::MatrixXd current = Eigen::MatrixXd::Zero(rows, deck.energy.G);

  const Eigen::MatrixXd psi0 = transport.sweep(0, Eigen::MatrixXd::Zero(rows, 4), scalar);
  scalar.col(0) = zerothMoment(deck, psi0);
  current.col(0) = firstMoment(deck, psi0);

  const Eigen::MatrixXd psi1 = transport.sweep(1, psi0, scalar);
  scalar.col(1) = zerothMoment(deck, psi1);
  current.col(1) = firstMoment(deck, psi1);

  const Eigen::VectorXd r0 = sm.calculateResiduals(0, scalar, current, psi0);
  const Eigen::VectorXd r1 = sm.calculateResiduals(1, scalar, current, psi1);

  REQUIRE(r0.size() == 8 * deck.mesh.n_x);
  CHECK(r0.cwiseAbs().maxCoeff() < 1e-13);
  CHECK(r1.cwiseAbs().maxCoeff() < 1e-13);
}

TEST_CASE("SM residuals vanish on a converged source iteration with within-group scattering") {
  // Within-group scattering puts phi_g on both sides of the equations, so the moments only
  // satisfy (53) once source iteration has converged -- to about its tolerance.
  const InputDeck deck = makeSMDeck(true);
  SourceIteration si(deck);
  si.solve(1e-14, 10000);
  SecondMoment sm(deck);

  const Eigen::MatrixXd& scalar = si.solution.scalar_flux;
  Eigen::MatrixXd current(scalar.rows(), scalar.cols());
  for (int g = 0; g < deck.energy.G; ++g) {
    current.col(g) = firstMoment(deck, si.solution.angular_flux[g]);
  }

  for (int g = 0; g < deck.energy.G; ++g) {
    CAPTURE(g);
    const Eigen::VectorXd r =
        sm.calculateResiduals(g, scalar, current, si.solution.angular_flux[g]);
    CHECK(r.cwiseAbs().maxCoeff() < 1e-11);
  }
}

TEST_CASE("SM residuals are nonzero, and local, when the LO solution is perturbed") {
  // Guards the two tests above against passing vacuously, and checks the face bookkeeping:
  // a corner only enters its own cell and the neighbor across the face it touches.
  const InputDeck deck = makeSMDeck(false);
  const int rows = 4 * deck.mesh.n_x;
  TransportOperator transport(deck);
  SecondMoment sm(deck);

  Eigen::MatrixXd scalar = Eigen::MatrixXd::Zero(rows, deck.energy.G);
  Eigen::MatrixXd current = Eigen::MatrixXd::Zero(rows, deck.energy.G);
  const Eigen::MatrixXd psi0 = transport.sweep(0, Eigen::MatrixXd::Zero(rows, 4), scalar);
  scalar.col(0) = zerothMoment(deck, psi0);
  current.col(0) = firstMoment(deck, psi0);

  auto cellNorm = [](const Eigen::VectorXd& r, int i) { return r.segment(8 * i, 8).norm(); };

  SUBCASE("L corner of cell 1 touches face x_1: cells 0 and 1 change") {
    for (int corner : {0, 2}) { // u_L, d_L
      CAPTURE(corner);
      Eigen::MatrixXd perturbed = current;
      perturbed(4 * 1 + corner, 0) += 1e-3;
      const Eigen::VectorXd r = sm.calculateResiduals(0, scalar, perturbed, psi0);
      CHECK(cellNorm(r, 0) > 1e-6);
      CHECK(cellNorm(r, 1) > 1e-6);
      CHECK(cellNorm(r, 2) < 1e-13);
      CHECK(cellNorm(r, 3) < 1e-13);
    }
  }

  SUBCASE("R corner of cell 1 touches face x_2: cells 1 and 2 change") {
    for (int corner : {1, 3}) { // u_R, d_R
      CAPTURE(corner);
      Eigen::MatrixXd perturbed = scalar;
      perturbed(4 * 1 + corner, 0) += 1e-3;
      const Eigen::VectorXd r = sm.calculateResiduals(0, perturbed, current, psi0);
      CHECK(cellNorm(r, 0) < 1e-13);
      CHECK(cellNorm(r, 1) > 1e-6);
      CHECK(cellNorm(r, 2) > 1e-6);
      CHECK(cellNorm(r, 3) < 1e-13);
    }
  }
}

TEST_CASE("cellResidual matches the compact Kronecker form of (53) in smm.md") {
  // An independent second statement of the same equations, checked on arbitrary (unphysical)
  // data so no term can hide behind a cancellation that happens to hold for real solutions.
  using Eigen::Matrix2d;
  using Eigen::Matrix4d;
  using Eigen::Vector2d;
  using Eigen::Vector4d;

  const InputDeck deck = makeSMDeck(false);
  const SecondMoment sm(deck);
  const int G = deck.energy.G;

  std::srand(20260923);
  const double dx = 0.7;
  const double dE = 0.4;
  const double sigma_t = 1.3;
  const double S_bar = 0.6;
  const double S_Eg = 0.45;
  const double S_Egm1 = 0.8;

  // [u_L, u_R, d_L, d_R] at corners, [u(x_{i-1}), u(x_i), d(x_{i-1}), d(x_i)] on faces
  const Vector4d phi = Vector4d::Random(), J = Vector4d::Random(), F = Vector4d::Random();
  const Vector4d phi_b = Vector4d::Random(), J_b = Vector4d::Random(), F_b = Vector4d::Random();
  const Vector4d phi_gm1 = Vector4d::Random(), J_gm1 = Vector4d::Random();
  const Vector4d q0 = Vector4d::Random(), q1 = Vector4d::Random();
  const Eigen::VectorXd sigma_sdEprime = Eigen::VectorXd::Random(G);
  const Eigen::MatrixXd phi_gprime = Eigen::MatrixXd::Random(4, G);

  auto up = [](const Vector4d& v) -> Vector2d { return v.head<2>(); };
  auto down = [](const Vector4d& v) -> Vector2d { return v.tail<2>(); };

  const Eigen::Vector<double, 8> r =
      sm.cellResidual(dx, dE, sigma_t, S_bar, S_Eg, S_Egm1, down(phi_gm1), down(J_gm1), up(phi_b),
                      down(phi_b), up(J_b), down(J_b), up(F_b), down(F_b), up(F), down(F), up(q0),
                      down(q0), up(q1), down(q1), sigma_sdEprime, phi_gprime.topRows<2>(),
                      phi_gprime.bottomRows<2>(), up(phi), down(phi), up(J), down(J))
          .cast<double>();

  Matrix2d M, L, Lb, ones, upper;
  M << 2, 1, 1, 2;
  M /= 6;
  L << 0.5, 0.5, -0.5, -0.5;
  Lb << -1, 0, 0, 1;
  ones << 1, 1, 1, 1;
  upper << 0, 1, 0, 0;

  Matrix2d C;
  C << sigma_t / 3 + S_bar / (2 * dE), sigma_t / 6 + S_bar / (2 * dE),
      sigma_t / 6 - S_bar / (2 * dE), sigma_t / 3 + (S_Eg - S_bar / 2) / dE;

  const Matrix4d S = kron(M, L);
  const Matrix4d Sb = kron(M, Lb);
  const Matrix4d R = dx * kron(C, M);
  const Matrix4d X = kron(ones, M);
  const Matrix4d P = kron(upper, M);
  const Matrix4d Q = kron(M, M);
  auto D = [&](const Vector4d& v, const Vector4d& vb) -> Vector4d { return Sb * vb + S * v; };

  const Vector4d scatter = phi_gprime * sigma_sdEprime;
  const Vector4d balance =
      D(J, J_b) + R * phi - (dx / dE) * S_Egm1 * P * phi_gm1 - (dx / 4) * X * scatter - dx * Q * q0;
  const Vector4d first_moment =
      D(phi, phi_b) / 3 + R * J - D(F, F_b) - (dx / dE) * S_Egm1 * P * J_gm1 - dx * Q * q1;

  for (int k = 0; k < 4; ++k) {
    CAPTURE(k);
    CHECK(r(k) == doctest::Approx(balance(k)).epsilon(1e-14));
    CHECK(r(4 + k) == doctest::Approx(first_moment(k)).epsilon(1e-14));
  }
}

TEST_CASE("calculateResiduals with debug_max returns the same residuals") {
  const InputDeck deck = makeSMDeck(true);
  const int rows = 4 * deck.mesh.n_x;
  SecondMoment sm(deck);

  const Eigen::MatrixXd scalar = Eigen::MatrixXd::Constant(rows, deck.energy.G, 0.8);
  const Eigen::MatrixXd current = Eigen::MatrixXd::Constant(rows, deck.energy.G, 0.1);
  const Eigen::MatrixXd psi = Eigen::MatrixXd::Constant(rows, deck.angle.M, 0.4);

  const Eigen::VectorXd quiet = sm.calculateResiduals(1, scalar, current, psi);
  const Eigen::VectorXd logged = sm.calculateResiduals(1, scalar, current, psi, true);
  CHECK(quiet.cwiseAbs().maxCoeff() > 0);
  CHECK(logged == quiet);
}

TEST_CASE("LO solve reproduces the moments of an exact transport solution") {
  // Closures from the exact discrete transport solution make the LO system exact, so its solution
  // must be that solution's zeroth and first moments -- group 1 also checks the CSD and
  // downscatter coupling to group 0.
  const InputDeck deck = makeSMDeck(false);
  const int rows = 4 * deck.mesh.n_x;
  TransportOperator transport(deck);
  SecondMoment sm(deck);

  Eigen::MatrixXd scalar = Eigen::MatrixXd::Zero(rows, deck.energy.G);
  Eigen::MatrixXd current = Eigen::MatrixXd::Zero(rows, deck.energy.G);
  Eigen::MatrixXd psi_up = Eigen::MatrixXd::Zero(rows, deck.angle.M);

  for (int g = 0; g < deck.energy.G; ++g) {
    CAPTURE(g);
    const Eigen::MatrixXd psi = transport.sweep(g, psi_up, scalar);

    sm.factorizeGroup(g);
    const auto [phi, J] =
        sm.solveGroup(sm.buildGroupRHS(g, sm.computeClosures(psi), scalar, current));

    const Eigen::VectorXd phi_exact = zerothMoment(deck, psi);
    const Eigen::VectorXd J_exact = firstMoment(deck, psi);
    CHECK((phi - phi_exact).norm() < 1e-13 * phi_exact.norm());
    CHECK((J - J_exact).norm() < 1e-13 * J_exact.norm());

    scalar.col(g) = phi;
    current.col(g) = J;
    psi_up = psi;
  }
}

TEST_CASE("LO solve reproduces the moments of a converged source iteration") {
  // Within-group scattering is on the LHS of the LO system, so this checks the matrix's w0 term.
  const InputDeck deck = makeSMDeck(true);
  SourceIteration si(deck);
  si.solve(1e-14, 10000);
  SecondMoment sm(deck);

  const Eigen::MatrixXd& scalar = si.solution.scalar_flux;
  Eigen::MatrixXd current(scalar.rows(), scalar.cols());
  for (int g = 0; g < deck.energy.G; ++g) {
    current.col(g) = firstMoment(deck, si.solution.angular_flux[g]);
  }

  for (int g = 0; g < deck.energy.G; ++g) {
    CAPTURE(g);
    const Eigen::MatrixXd& psi = si.solution.angular_flux[g];
    sm.factorizeGroup(g);
    const auto [phi, J] =
        sm.solveGroup(sm.buildGroupRHS(g, sm.computeClosures(psi), scalar, current));
    CHECK((phi - scalar.col(g)).norm() < 1e-10 * scalar.col(g).norm());
    CHECK((J - current.col(g)).norm() < 1e-10 * current.col(g).norm());
  }
}

TEST_CASE("LO solution zeroes the SM residuals for arbitrary closures") {
  // psi here is not a transport solution, so the LO solution differs from its moments -- but it
  // must still satisfy (53) exactly as calculateResiduals states it, term for term.
  const InputDeck deck = makeSMDeck(true);
  const int rows = 4 * deck.mesh.n_x;
  SecondMoment sm(deck);

  std::srand(20260923);
  Eigen::MatrixXd scalar = Eigen::MatrixXd::Random(rows, deck.energy.G);
  Eigen::MatrixXd current = Eigen::MatrixXd::Random(rows, deck.energy.G);
  const Eigen::MatrixXd psi = Eigen::MatrixXd::Random(rows, deck.angle.M).array() + 1.0;

  for (int g = 0; g < deck.energy.G; ++g) {
    CAPTURE(g);
    sm.factorizeGroup(g);
    const auto [phi, J] =
        sm.solveGroup(sm.buildGroupRHS(g, sm.computeClosures(psi), scalar, current));
    scalar.col(g) = phi;
    current.col(g) = J;
    const Eigen::VectorXd r = sm.calculateResiduals(g, scalar, current, psi);
    CHECK(r.cwiseAbs().maxCoeff() < 1e-13);
  }
}

TEST_CASE("LO matrix is block tridiagonal with 80I - 32 nonzeros") {
  const InputDeck deck = makeSMDeck(true);
  const int I = deck.mesh.n_x;
  SecondMoment sm(deck);

  for (int g = 0; g < deck.energy.G; ++g) {
    CAPTURE(g);
    const Eigen::SparseMatrix<double> A = sm.buildGroupMatrix(g);
    REQUIRE(A.rows() == 8 * I);
    REQUIRE(A.cols() == 8 * I);
    CHECK(A.nonZeros() == 80 * I - 32);
    for (int k = 0; k < A.outerSize(); ++k) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
        CHECK(std::abs(it.row() / 8 - it.col() / 8) <= 1);
      }
    }
  }
}

TEST_CASE("solveGroup before factorizeGroup throws") {
  const InputDeck deck = makeSMDeck(false);
  const SecondMoment sm(deck);
  CHECK_THROWS_AS(sm.solveGroup(Eigen::VectorXd::Zero(8 * deck.mesh.n_x)), std::logic_error);
}
