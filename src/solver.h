#ifndef SOLVER_H
#define SOLVER_H

#include "cross_section.h"
#include "mesh.h"

// Base class for iterative transport solve strategies (e.g.
// SourceIterationSolver, SecondMomentSolver). Holds functionality shared by
// every derived solve strategy -- the high-order transport sweep chief among
// them -- so it lives in exactly one place rather than being duplicated
// across methods that all need it. Derived classes differ in what they do
// around the sweep (e.g. whether they also run a low-order solve to
// accelerate it), not in the sweep itself. Owns its own copies of the Mesh
// and CrossSection it was constructed with, so it doesn't depend on the
// caller's originals outliving it.
class Solver {
public:
  // Constructs a Solver over its own copies of mesh and cross_section.
  // cross_section's own Mesh reference is rebound to this Solver's mesh
  // copy, not the caller's original, so the two stay consistent regardless
  // of what happens to the objects passed in.
  Solver(const Mesh& mesh, const CrossSection& cross_section);
  virtual ~Solver() = default;

  const Mesh mesh;
  const CrossSection cross_section;

protected:
  // Performs one high-order transport sweep. Shared by every derived solve
  // strategy; not yet implemented.
  void sweep();
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
