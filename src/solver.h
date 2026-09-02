#ifndef SOLVER_H
#define SOLVER_H

#include <filesystem>
#include <vector>

#include <Eigen/Sparse>

#include "cross_section.h"
#include "fe_space.h"
#include "input_deck.h"
#include "mesh.h"

/// @enum LinearSolverKind
/// @brief Select at runtime which Eigen sparse linear solver is used
enum class LinearSolverKind {
  SparseLU,   ///< Sparse LU; Direct solver for genral matrices
  BiCGSTAB,   ///< Stabilized Biconjugate Gradient; Iterative solver for general matrices
  GMRES,      ///< General Minimal Residual; Iterative solver for general matrices
  SweepDirect ///< PLACEHOLDER
};

/// @class Solver
/// @brief Base class for iterative electron transport solver
/// @details Holds functions and members common to multiple derived classes (of different methods)
/// to avoid code duplication. Owns an InputDeck, so that Solver callers can declare solvers with
/// the same reconfigured deck.
class Solver {
public:
  // Constructs a Solver over its own copy of input_deck and fe_space.
  // input_deck.mesh must already be set -- there's no legitimate way to
  // construct a solver without one, unlike input_deck.xs, which can
  // legitimately be transiently unset (e.g. right after a mesh change) and
  // is only checked when something actually needs it. corner_order fixes
  // the corner-ordering convention this Solver's own assembly code uses
  // consistently, matching whatever Field the results are eventually
  // written into. linear_solver_kind selects which Eigen solver sweep()
  // uses for each ordinate/group's linear system.
  Solver(InputDeck input_deck, const FESpace& fe_space, AxisOrder corner_order = AxisOrder::XMajor,
         LinearSolverKind linear_solver_kind = LinearSolverKind::SparseLU);
  virtual ~Solver() = default;

  InputDeck input_deck;
  const FESpace fe_space;
  const AxisOrder corner_order;
  const LinearSolverKind linear_solver_kind;

protected:
  /// @brief Solve the high-order transport equation for a given source b.
  ///
  /// Builds the transport bilinear form matrix in A (over-writes current value), then use linear
  /// solver declared at compile-time to solve the system.
  ///
  /// @param &x Vector to write result to (will over-write)
  /// @param &A Reference to Eigen Sparse matrix at which to build the bilinear form
  /// @param &b source vector; must be pre-defined (i.e. by calling constructTransportLinear)
  /// @param mu Angle cosine for this ordinate; required to implement boundary conditions
  /// @param group Group index, used to index cross sections. Must be [0, G).
  void sweep(Eigen::VectorXd& x, Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b,
             double mu, int group) const;

  // Solve transport equation for 1 group, fixed source. A loop over all angles m in angular
  // quadrature. b is source, must be m by 4n_x. (rectangular). throws if b is the wrong size.
  void solveTransport(std::vector<Eigen::VectorXd>& x, Eigen::SparseMatrix<double>& A,
                      const std::vector<Eigen::VectorXd>& b, int group) const;

  /// @brief Construct LHS matrix for high-order transport equation
  ///
  /// Builds the high-order transport equation's LHS (bilinear form) for a single ordinate mu
  /// and a single energy group. Each cell is coupled to the spatially upwind cell (in the -mu
  /// direction) by the advection operator.
  ///
  /// @param &A Reference to the matrix to write the system to. A will be resized if necessary
  /// to 4*n_x by 4*n_x, where n_x is the number of spatiall cells (from input_deck.mesh).
  /// @param mu Cosine of the direction of travel for this ordinate
  /// @param group Energy group. Group must be in [0, G).
  void constructTransportBilinear(Eigen::SparseMatrix<double>& A, double mu, int group) const;

  /// @brief Construct RHS vector for high-order transport equation
  ///
  /// Computes scattering source from scalar_flux and latest_scalar_flux, computes CSD source from
  /// upwind_angular_flux, then combines scattering + CSD + external_source. Writes to b.
  ///
  /// @param &b Reference to the Eigen vector where the term is to be written. Resizes to 4*n_x if
  /// necessary.
  /// @param mu Cosine of the direction of travel for this ordinate. Used to incorporate incoming
  /// flux boundary conditions (i.e. BC only added for mu > 0 on left face, mu < 0 on right face).
  /// @param group Energy group. CSD source is zero for group=0; used for explicit in-group
  /// scattering treatment (use the most recent flux for this group).
  /// @param ordinate_index this ordinate's position in input_deck.angular_quadrature.mu, needed to
  /// look up its boundary condition.
  /// @param upwind_angular_flux Angular flux for previous group (or pass vector of zeros for group
  /// = 1), used to compute CSD source
  /// @param scalar_flux Scalar flux for all groups (only upwind are used), used to compute
  /// scattering source
  /// @param latest_scalar_flux Most recent value of scalar flux, used to calculate within-group
  /// scattering
  /// @param external_source Fixed source for this mu.
  void constructTransportLinear(Eigen::VectorXd& b, double mu, int group, int ordinate_index,
                                Field::ConstRow upwind_angular_flux, const Field& scalar_flux,
                                Field::ConstRow latest_scalar_flux,
                                const Eigen::VectorXd& external_source) const;

private:
  // Returns input_deck.xs. If it's unset (e.g. left cleared after a mesh
  // change and never reconfigured), logs a fatal error and terminates the
  // program -- solving without cross sections isn't a recoverable
  // condition, so this isn't a throw a caller is expected to catch.
  const CrossSection& requireCrossSection() const;

  // Returns input_deck.boundary_conditions. Same fatal-if-unset treatment
  // as requireCrossSection(), for the same reason -- there's no meaningful
  // way to assemble a boundary-facing RHS without it.
  const BoundaryConditions& requireBoundaryConditions() const;

  // Returns input_deck.angular_quadrature. Same fatal-if-unset treatment as
  // requireCrossSection() -- sweep() can't even know how many ordinates to
  // loop over without it.
  const AngularQuadrature& requireAngularQuadrature() const;

protected:
  // Solves A*x = b using whichever Eigen solver linear_solver_kind selects.
  // Throws std::runtime_error if the solve doesn't succeed (e.g. failure to
  // converge, a singular matrix) -- unlike the require*() methods above,
  // this is a per-solve numerical failure a caller could plausibly retry
  // (e.g. with a different linear_solver_kind), not a "this Solver was
  // never fully configured" condition, so it isn't treated as fatal.
  // protected (not private): derived solvers may reasonably want to call
  // this directly too, e.g. for a low-order solve.
  Eigen::VectorXd solveLinearSystem(const Eigen::SparseMatrix<double>& A,
                                    const Eigen::VectorXd& b) const;

  // Same as solveLinearSystem(A, b), but writes the result into x instead of
  // returning it, for callers (e.g. sweep()) that already have a vector to
  // write into.
  void solveLinearSystem(const Eigen::SparseMatrix<double>& A, Eigen::VectorXd& x,
                         const Eigen::VectorXd& b) const;
};

// Iteration parameters for SourceIterationSolver::solve() -- a single flat
// loop (solve transport, update scattering source, repeat), so this is just
// ConvergenceCriteria's shape. Plain aggregate, like ConvergenceCriteria: no
// validation yet, since nothing parses this from YAML yet either.
struct SourceIterationParams {
  int max_iters = 0;
  double epsilon = 0.0;
};

/// @class SourceIterationSolver
/// @brief Transport solver with "naive" source iteration method; No projection, no acceleration
/// @details A baseline method with no acceleration, preconditioning, or projection. Solves
/// transport equation with intial source guess, calculate scattering source, repeat. Exists to
/// compare other methods to
class SourceIterationSolver : public Solver {
public:
  using Solver::Solver;

  // Runs source iteration to convergence (or until params.max_iters), per
  // params. Takes the concrete params type directly -- no base
  // IterationParameters type, no runtime check that it's "the right"
  // derived type: callers already have a SourceIterationSolver, so they
  // already know which params type matches it.
  void solve(const SourceIterationParams& params);

  // Writes this solve's results to path. Not virtual: each solver's
  // results are shaped differently, so there's no shared interface for
  // this to satisfy.
  void writeResults(const std::filesystem::path& path) const;
};

// Iteration parameters for SecondMomentSolver::solve(). Placeholder shape --
// outer loop over transport sweep + closure, inner loop for the low-order
// SMM solve -- to be refined once that iteration structure is actually
// implemented.
struct SecondMomentParams {
  int outer_max_iters = 0;
  double outer_epsilon = 0.0;
  int inner_max_iters = 0;
  double inner_epsilon = 0.0;
};

/// @class SecondMomentSolver
/// @brief Transport solver using the Second Moment method
/// @details Projective method; solves the coupled equations for angular flux, scalar flux, and
/// current. Iterates between solution of the transport (high-order) equation, then calculates
/// scalar flux and current using closure terms (calculated from angular flux), then uses scalar
/// flux to calculate the source for the next transport sweep.
class SecondMomentSolver : public Solver {
public:
  using Solver::Solver;

  // See SourceIterationSolver::solve() -- same reasoning, different params.
  void solve(const SecondMomentParams& params);

  // See SourceIterationSolver::writeResults().
  void writeResults(const std::filesystem::path& path) const;
};

#endif
