// Copyright (c) 2026, Kyle Hansen (khansen3@ncsu.edu)
//
// Funded by CARRE (https://carre-psaapiv.org/)
//
// Licensed under BSD 3-Clause License; Redistribution and use in source and binary forms, with
// or without modification are permitted provided that the terms of the license are met.

#include "source_iteration.h"

#include "smm.h"
#include <doctest.h>

namespace {
// Infinite homogeneous medium: pure within-group scatterer (sigma_s,gg = sigma_t), isotropic
// source q per unit energy, and group widths != 1 so a missing dE_g / dE_g' factor in the
// scattering source can't cancel out. Summing the u and d rows of the LD cell equations, the
// within-group scattering cancels the removal exactly and leaves the slowing-down balance
//   S(E_g) psi_d(E_g) = S(E_{g-1}) psi_d(E_{g-1}) + q dE_g  =>  S(E_g) psi_d(E_g) = q (E_0 - E_g).
// The corner values within each group come from the 2x2 energy-only LD system; the boundary
// conditions are set to them, so the exact discrete solution is also flat in x and isotropic.
struct InfiniteMedium {
  InputDeck deck;
  std::vector<Eigen::Vector2d> psi; // [u, d] per group
  Eigen::Vector4d S_b;              // group-boundary stopping powers
};

InfiniteMedium makeInfiniteMediumDeck() {
  const int n_x = 3;
  const int G = 3;
  const int M = 4;
  const double q = 0.9;
  const double sigma = 1.7;

  InfiniteMedium im;
  InputDeck& deck = im.deck;
  deck.mesh.n_x = n_x;
  deck.mesh.x_boundary = Eigen::Vector4d(0.0, 0.23, 0.61, 0.8);
  deck.energy.G = G;
  deck.energy.E_boundary = Eigen::Vector4d(2.0, 1.37, 0.52, 0.11);
  deck.angle.M = M;
  deck.angle.mu = Eigen::Vector4d(-0.8611363115940526, -0.3399810435848563, 0.3399810435848563,
                                  0.8611363115940526);
  deck.angle.w = Eigen::Vector4d(0.3478548451374538, 0.6521451548625461, 0.6521451548625461,
                                 0.3478548451374538);

  Material mat;
  mat.name = "scatterer";
  mat.total = Eigen::Vector3d::Constant(sigma);
  mat.S = Eigen::Vector3d(0.55, 0.8, 1.2);
  mat.S_b = Eigen::Vector4d(0.5, 0.62, 1.0, 1.4);
  mat.scatter = sigma * Eigen::Matrix3d::Identity();
  deck.xs.set_materials({mat}, {0, 0, 0});
  im.S_b = mat.S_b;

  // Per group, the flat isotropic cell equations (divided by dx / 2):
  //   u: (s/3 + S/(2dE)) u + (s/6 + S/(2dE)) d - (s/4)(u + d) = q/2 + S_up psi_in / dE
  //   d: (s/6 - S/(2dE)) u + (s/3 + (S_down - S/2)/dE) d - (s/4)(u + d) = q/2
  deck.bc.values = Eigen::MatrixXd(2 * G, M);
  double psi_in = 0.0;
  for (int g = 0; g < G; ++g) {
    const double dE = deck.energy.E_boundary(g) - deck.energy.E_boundary(g + 1);
    const double S = mat.S(g);
    Eigen::Matrix2d A;
    A << sigma / 12 + S / (2 * dE), -sigma / 12 + S / (2 * dE), //
        -sigma / 12 - S / (2 * dE), sigma / 12 + (mat.S_b(g + 1) - S / 2) / dE;
    const Eigen::Vector2d b(q / 2 + mat.S_b(g) * psi_in / dE, q / 2);
    const Eigen::Vector2d psi_g = A.partialPivLu().solve(b);
    im.psi.push_back(psi_g);
    deck.bc.values.row(2 * g).setConstant(psi_g(0));
    deck.bc.values.row(2 * g + 1).setConstant(psi_g(1));
    psi_in = psi_g(1);
  }

  deck.source.values.assign(G, Eigen::MatrixXd::Constant(4 * n_x, M, q));
  deck.validate();
  return im;
}

// Solution is flat, isotropic, matches the per-group corner values, and satisfies the analytic
// slowing-down balance S(E_g) psi_d = q (E_0 - E_g).
void checkInfiniteMedium(const InfiniteMedium& im, const std::vector<Eigen::MatrixXd>& psi) {
  const InputDeck& deck = im.deck;
  const double q = 0.9;
  for (int g = 0; g < deck.energy.G; ++g) {
    CAPTURE(g);
    const double S_down = im.S_b(g + 1);
    const double E0_minus_Eg = deck.energy.E_boundary(0) - deck.energy.E_boundary(g + 1);
    CHECK(S_down * im.psi[g](1) == doctest::Approx(q * E0_minus_Eg).epsilon(1e-13));
    for (int i = 0; i < deck.mesh.n_x; ++i) {
      for (int m = 0; m < deck.angle.M; ++m) {
        CAPTURE(i);
        CAPTURE(m);
        CHECK(psi[g](4 * i + 0, m) == doctest::Approx(im.psi[g](0)).epsilon(1e-10));
        CHECK(psi[g](4 * i + 1, m) == doctest::Approx(im.psi[g](0)).epsilon(1e-10));
        CHECK(psi[g](4 * i + 2, m) == doctest::Approx(im.psi[g](1)).epsilon(1e-10));
        CHECK(psi[g](4 * i + 3, m) == doctest::Approx(im.psi[g](1)).epsilon(1e-10));
      }
    }
  }
}
} // namespace

TEST_CASE("source iteration conserves particles in an infinite pure scatterer with dE != 1") {
  const InfiniteMedium im = makeInfiniteMediumDeck();
  SourceIteration si(im.deck);
  si.solve(1e-14, 10000);
  checkInfiniteMedium(im, si.solution.angular_flux);
}

TEST_CASE("SMM conserves particles in an infinite pure scatterer with dE != 1") {
  const InfiniteMedium im = makeInfiniteMediumDeck();
  SecondMoment sm(im.deck);
  sm.solve(1e-14, 10000);
  checkInfiniteMedium(im, sm.solution.angular_flux);
}
