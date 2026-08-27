#ifndef SOLVER_H
#define SOLVER_H

#include <Eigen/Sparse>

#include "cross_section.h"
#include "fe_space.h"
#include "input_deck.h"
#include "mesh.h"

// Base class for iterative transport solve strategies (e.g.
// SourceIterationSolver, SecondMomentSolver). Holds functionality shared by
// every derived solve strategy -- the high-order transport sweep chief among
// them -- so it lives in exactly one place rather than being duplicated
// across methods that all need it. Derived classes differ in what they do
// around the sweep (e.g. whether they also run a low-order solve to
// accelerate it), not in the sweep itself. Owns an InputDeck rather than
// its own separate copies of mesh/cross-section data -- input_deck's own
// setters (setMesh, setAngularQuadrature, etc.) are how a caller reconfigures
// part of the problem after construction, e.g. to re-run the same problem
// under a different angular quadrature. FESpace is discretization, not
// problem data, so it stays a separate member rather than living on
// InputDeck.
class Solver {
public:
  // Constructs a Solver over its own copy of input_deck and fe_space.
  // input_deck.mesh must already be set -- there's no legitimate way to
  // construct a solver without one, unlike input_deck.xs, which can
  // legitimately be transiently unset (e.g. right after a mesh change) and
  // is only checked when something actually needs it. corner_order fixes
  // the corner-ordering convention this Solver's own assembly code uses
  // consistently, matching whatever Field the results are eventually
  // written into.
  Solver(InputDeck input_deck, const FESpace& fe_space, AxisOrder corner_order = AxisOrder::XMajor);
  virtual ~Solver() = default;

  InputDeck input_deck;
  const FESpace fe_space;
  const AxisOrder corner_order;

protected:
  // Performs one high-order transport sweep. Shared by every derived solve
  // strategy; not yet implemented.
  void sweep();

  // Builds the high-order transport equation's LHS bilinear form for a
  // single ordinate mu and a single energy group: streaming (upwinded in x,
  // coupling to the neighboring cell) plus absorption and this group's own
  // slowing-down removal term (both local to each cell, no neighbor
  // coupling). Energy-group coupling never appears here -- groups are
  // solved sequentially from high energy to low, so the previous
  // (already-solved) group's inflow is a known RHS source, not an LHS
  // unknown. A is resized to 4*n_x x 4*n_x (n_x from input_deck.mesh); its
  // previous contents are discarded. group must be in [0, G).
  void constructTransportBilinear(Eigen::SparseMatrix<double>& A, double mu, int group) const;

private:
  // Returns input_deck.xs. If it's unset (e.g. left cleared after a mesh
  // change and never reconfigured), logs a fatal error and terminates the
  // program -- solving without cross sections isn't a recoverable
  // condition, so this isn't a throw a caller is expected to catch.
  const CrossSection& requireCrossSection() const;
};

// Solves the high-order transport equation by source iteration (repeated
// sweeps to convergence).
class SourceIterationSolver : public Solver {
public:
  using Solver::Solver;
};

// Solves the transport equation via the Second Moment Method: alternates
// high-order sweeps with a low-order solve to accelerate convergence.
class SecondMomentSolver : public Solver {
public:
  using Solver::Solver;
};

#endif
