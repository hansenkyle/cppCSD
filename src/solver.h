#ifndef SOLVER_H
#define SOLVER_H

/// @class Solver
/// @brief Base class for iterative transport solve strategies (e.g. SourceIterationSolver,
/// SecondMomentSolver)
/// @details Holds functionality shared by every derived solve strategy -- the high-order transport
/// sweep chief among them -- so it lives in exactly one place rather than being duplicated across
/// methods that all need it. Derived classes differ in what they do around the sweep (e.g. whether
/// they also run a low-order solve to accelerate it), not in the sweep itself.
class Solver {
public:
  virtual ~Solver() = default;

protected:
  // Performs one high-order transport sweep. Shared by every derived solve
  // strategy; not yet implemented.
  void sweep();
};

/// @class SourceIterationSolver
/// @brief Solves the high-order transport equation by source iteration (repeated sweeps to
/// convergence)
class SourceIterationSolver : public Solver {};

/// @class SecondMomentSolver
/// @brief Solves the transport equation via the Second Moment Method: alternates high-order sweeps
/// with a low-order solve to accelerate convergence
class SecondMomentSolver : public Solver {};

#endif
