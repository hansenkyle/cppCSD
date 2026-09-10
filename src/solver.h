#ifndef SOLVER_H
#define SOLVER_H

#include <Eigen/Sparse>

#include "cross_section.h"
#include "fe_space.h"
#include "mesh.h"

// Base class for iterative transport solve strategies (e.g.
// SourceIterationSolver, SecondMomentSolver). Holds functionality shared by
// every derived solve strategy -- the high-order transport sweep chief among
// them -- so it lives in exactly one place rather than being duplicated
// across methods that all need it. Derived classes differ in what they do
// around the sweep (e.g. whether they also run a low-order solve to
// accelerate it), not in the sweep itself. Owns its own copies of the Mesh,
// CrossSection, and FESpace it was constructed with, so it doesn't depend on
// the caller's originals outliving it.
class Solver {
public:
  // Constructs a Solver over its own copies of mesh, cross_section, and
  // fe_space. cross_section's own Mesh reference is rebound to this
  // Solver's mesh copy, not the caller's original, so the two stay
  // consistent regardless of what happens to the objects passed in.
  Solver(const Mesh& mesh, const CrossSection& cross_section, const FESpace& fe_space);
  virtual ~Solver() = default;

  const Mesh mesh;
  const CrossSection cross_section;
  const FESpace fe_space;

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
  // unknown. A is resized to 4*mesh.n_x x 4*mesh.n_x; its previous contents
  // are discarded. group must be in [0, mesh.G).
  void constructTransportBilinear(Eigen::SparseMatrix<double>& A, double mu, int group) const;
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
