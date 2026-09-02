#include "solver.h"

#include <array>
#include <cstdlib>
#include <optional>
#include <stdexcept>
#include <vector>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <unsupported/Eigen/IterativeSolvers>

#include "logger.h"

Solver::Solver(InputDeck input_deck, const FESpace& fe_space, AxisOrder corner_order,
               LinearSolverKind linear_solver_kind)
    : input_deck(std::move(input_deck)), fe_space(fe_space), corner_order(corner_order),
      linear_solver_kind(linear_solver_kind) {
  if (!this->input_deck.mesh.has_value()) {
    throw std::invalid_argument("Solver: input_deck.mesh must be set");
  }
}

const CrossSection& Solver::requireCrossSection() const {
  if (!input_deck.xs.has_value()) {
    LDCSD_LOG_ERROR("Solver: input_deck.xs is not set (cleared by a mesh change and never "
                    "reconfigured?) -- cannot solve");
    std::exit(1);
  }
  return *input_deck.xs;
}

const BoundaryConditions& Solver::requireBoundaryConditions() const {
  if (!input_deck.boundary_conditions.has_value()) {
    LDCSD_LOG_ERROR("Solver: input_deck.boundary_conditions is not set -- cannot solve");
    std::exit(1);
  }
  return *input_deck.boundary_conditions;
}

const AngularQuadrature& Solver::requireAngularQuadrature() const {
  if (!input_deck.angular_quadrature.has_value()) {
    LDCSD_LOG_ERROR("Solver: input_deck.angular_quadrature is not set -- cannot solve");
    std::exit(1);
  }
  return *input_deck.angular_quadrature;
}

Eigen::VectorXd Solver::solveLinearSystem(const Eigen::SparseMatrix<double>& A,
                                          const Eigen::VectorXd& b) const {
  switch (linear_solver_kind) {
  case LinearSolverKind::SparseLU: {
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
      throw std::runtime_error("Solver::solveLinearSystem: SparseLU factorization failed");
    }
    const Eigen::VectorXd x = solver.solve(b);
    if (solver.info() != Eigen::Success) {
      throw std::runtime_error("Solver::solveLinearSystem: SparseLU solve failed");
    }
    return x;
  }
  case LinearSolverKind::BiCGSTAB: {
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    const Eigen::VectorXd x = solver.solve(b);
    if (solver.info() != Eigen::Success) {
      throw std::runtime_error("Solver::solveLinearSystem: BiCGSTAB failed to converge");
    }
    return x;
  }
  case LinearSolverKind::GMRES: {
    Eigen::GMRES<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    const Eigen::VectorXd x = solver.solve(b);
    if (solver.info() != Eigen::Success) {
      throw std::runtime_error("Solver::solveLinearSystem: GMRES failed to converge");
    }
    return x;
  }
  case LinearSolverKind::SweepDirect:
    throw std::runtime_error("Solver::solveLinearSystem: SweepDirect is not yet implemented");
  }
  throw std::runtime_error("Solver::solveLinearSystem: unknown linear_solver_kind");
}

namespace {

// Row/column order for the local dense blocks below -- independent of any
// Field's corner_order, which only matters once these blocks are scattered
// into the global matrix via localIndex().
constexpr std::array<Corner, 4> kLocalCorners = {Corner::LeftDown, Corner::LeftUp,
                                                 Corner::RightDown, Corner::RightUp};

int localIndex(int cell, Corner corner, AxisOrder corner_order) {
  return cell * 4 + cornerSlot(corner, corner_order);
}

// Scatters a cell's full 4x4 local block (streaming's local part plus
// absorption/slowing-down, both local to this cell) into triplets.
void appendSelfBlock(std::vector<Eigen::Triplet<double>>& triplets, int cell,
                     AxisOrder corner_order, const Eigen::Matrix4d& block) {
  for (int r = 0; r < 4; ++r) {
    const int row = localIndex(cell, kLocalCorners[r], corner_order);
    for (int c = 0; c < 4; ++c) {
      const int col = localIndex(cell, kLocalCorners[c], corner_order);
      triplets.emplace_back(row, col, block(r, c));
    }
  }
}

void zeroRow(Field::Row row, int n_x) {
  for (int i = 0; i < n_x; ++i) {
    CornerValues cv = row[i];
    cv.leftDown() = 0.0;
    cv.leftUp() = 0.0;
    cv.rightDown() = 0.0;
    cv.rightUp() = 0.0;
  }
}

// Overwrites row with local's values (local is a 4*n_x solve result in the
// same (cell, corner) layout Field::index()/localIndex() use).
void assignLocalVectorToRow(Field::Row row, int n_x, AxisOrder corner_order,
                            const Eigen::VectorXd& local) {
  for (int i = 0; i < n_x; ++i) {
    CornerValues cv = row[i];
    cv.leftDown() = local(localIndex(i, Corner::LeftDown, corner_order));
    cv.leftUp() = local(localIndex(i, Corner::LeftUp, corner_order));
    cv.rightDown() = local(localIndex(i, Corner::RightDown, corner_order));
    cv.rightUp() = local(localIndex(i, Corner::RightUp, corner_order));
  }
}

// Adds scale * local into row, e.g. accumulating a quadrature-weighted sum
// of several ordinates' solves into a scalar flux.
void accumulateLocalVectorIntoRow(Field::Row row, int n_x, AxisOrder corner_order,
                                  const Eigen::VectorXd& local, double scale) {
  for (int i = 0; i < n_x; ++i) {
    CornerValues cv = row[i];
    cv.leftDown() += scale * local(localIndex(i, Corner::LeftDown, corner_order));
    cv.leftUp() += scale * local(localIndex(i, Corner::LeftUp, corner_order));
    cv.rightDown() += scale * local(localIndex(i, Corner::RightDown, corner_order));
    cv.rightUp() += scale * local(localIndex(i, Corner::RightUp, corner_order));
  }
}

} // namespace

void Solver::constructTransportBilinear(Eigen::SparseMatrix<double>& A, double mu,
                                        int group) const {
  const Mesh& mesh = *input_deck.mesh;
  const CrossSection& cross_section = requireCrossSection();

  if (group < 0 || group >= mesh.G) {
    throw std::out_of_range("Solver::constructTransportBilinear: group index out of range");
  }

  const int n_x = mesh.n_x;
  const double dE = mesh.dE[group];
  const std::vector<double>& xs_total = cross_section.total[group];
  const std::vector<double>& stop_power = cross_section.stop_power[group];
  // This group's own lower-energy-edge stopping power -- the flux this
  // group loses into the next (lower-energy) group. mesh's E_boundary
  // descends from high energy to 0, so group g's lower edge is boundary
  // row g+1.
  const std::vector<double>& stop_power_bound_down = cross_section.stop_power_boundary[group + 1];

  const double m00 = fe_space.M(0, 0);
  const double m01 = fe_space.M(0, 1);
  const double m10 = fe_space.M(1, 0);
  const double m11 = fe_space.M(1, 1);

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(static_cast<std::size_t>(n_x) * 20);

  for (int i = 0; i < n_x; ++i) {
    const double dx = mesh.dx[i];
    const double v1 = dx * ((1.0 / 6.0) * xs_total[i] + (1.0 / (2.0 * dE)) * stop_power[i]);
    const double v2 = dx * ((1.0 / 3.0) * xs_total[i] + (1.0 / (2.0 * dE)) * stop_power[i]);
    const double v3 = dx * ((1.0 / 3.0) * xs_total[i] +
                            (1.0 / dE) * (-stop_power[i] / 2.0 + stop_power_bound_down[i]));
    const double v4 = dx * ((1.0 / 6.0) * xs_total[i] + (1.0 / dE) * (-stop_power[i] / 2.0));

    // Row/col order within this block: 0=LeftDown, 1=LeftUp, 2=RightDown,
    // 3=RightUp (matches kLocalCorners).
    Eigen::Matrix4d self_block = Eigen::Matrix4d::Zero();

    // streaming, local part
    self_block(1, 0) += (1.0 / 12.0) * mu;
    self_block(1, 2) += (1.0 / 12.0) * mu;
    self_block(1, 1) += (1.0 / 6.0) * mu;
    self_block(1, 3) += (1.0 / 6.0) * mu;

    self_block(3, 0) += (-1.0 / 12.0) * mu;
    self_block(3, 2) += (-1.0 / 12.0) * mu;
    self_block(3, 1) += (-1.0 / 6.0) * mu;
    self_block(3, 3) += (-1.0 / 6.0) * mu;

    self_block(0, 0) += (1.0 / 6.0) * mu;
    self_block(0, 2) += (1.0 / 6.0) * mu;
    self_block(0, 1) += (1.0 / 12.0) * mu;
    self_block(0, 3) += (1.0 / 12.0) * mu;

    self_block(2, 0) += (-1.0 / 6.0) * mu;
    self_block(2, 2) += (-1.0 / 6.0) * mu;
    self_block(2, 1) += (-1.0 / 12.0) * mu;
    self_block(2, 3) += (-1.0 / 12.0) * mu;

    // absorption + slowing-down, local part
    self_block(1, 0) += v1 * m00;
    self_block(1, 2) += v1 * m01;
    self_block(1, 1) += v2 * m00;
    self_block(1, 3) += v2 * m01;

    self_block(3, 0) += v1 * m10;
    self_block(3, 2) += v1 * m11;
    self_block(3, 1) += v2 * m10;
    self_block(3, 3) += v2 * m11;

    self_block(0, 0) += v3 * m00;
    self_block(0, 2) += v3 * m01;
    self_block(0, 1) += v4 * m00;
    self_block(0, 3) += v4 * m01;

    self_block(2, 0) += v3 * m10;
    self_block(2, 2) += v3 * m11;
    self_block(2, 1) += v4 * m10;
    self_block(2, 3) += v4 * m11;

    appendSelfBlock(triplets, i, corner_order, self_block);

    // streaming, upwind coupling part -- 2x2, only when a neighbor exists
    // (at a domain boundary, the missing inflow is a boundary condition on
    // the RHS instead, not modeled here).
    if (mu > 0.0) {
      if (i > 0) {
        Eigen::Matrix2d coupling;
        coupling << (-1.0 / 6.0) * mu, (-1.0 / 3.0) * mu, (-1.0 / 3.0) * mu, (-1.0 / 6.0) * mu;
        const int row_left_up = localIndex(i, Corner::LeftUp, corner_order);
        const int row_left_down = localIndex(i, Corner::LeftDown, corner_order);
        const int col_right_down = localIndex(i - 1, Corner::RightDown, corner_order);
        const int col_right_up = localIndex(i - 1, Corner::RightUp, corner_order);
        triplets.emplace_back(row_left_up, col_right_down, coupling(0, 0));
        triplets.emplace_back(row_left_up, col_right_up, coupling(0, 1));
        triplets.emplace_back(row_left_down, col_right_down, coupling(1, 0));
        triplets.emplace_back(row_left_down, col_right_up, coupling(1, 1));
      }
    } else {
      if (i < n_x - 1) {
        Eigen::Matrix2d coupling;
        coupling << (1.0 / 6.0) * mu, (1.0 / 3.0) * mu, (1.0 / 3.0) * mu, (1.0 / 6.0) * mu;
        const int row_right_up = localIndex(i, Corner::RightUp, corner_order);
        const int row_right_down = localIndex(i, Corner::RightDown, corner_order);
        const int col_left_down = localIndex(i + 1, Corner::LeftDown, corner_order);
        const int col_left_up = localIndex(i + 1, Corner::LeftUp, corner_order);
        triplets.emplace_back(row_right_up, col_left_down, coupling(0, 0));
        triplets.emplace_back(row_right_up, col_left_up, coupling(0, 1));
        triplets.emplace_back(row_right_down, col_left_down, coupling(1, 0));
        triplets.emplace_back(row_right_down, col_left_up, coupling(1, 1));
      }
    }
  }

  const int size = 4 * n_x;
  A = Eigen::SparseMatrix<double>(size, size);
  A.setFromTriplets(triplets.begin(), triplets.end());
}

void Solver::constructTransportLinear(Eigen::VectorXd& b, double mu, int group, int ordinate_index,
                                      Field::ConstRow upwind_angular_flux, const Field& scalar_flux,
                                      Field::ConstRow latest_scalar_flux,
                                      const Eigen::VectorXd& external_source) const {
  const Mesh& mesh = *input_deck.mesh;
  const CrossSection& cross_section = requireCrossSection();
  const BoundaryConditions& bc = requireBoundaryConditions();

  if (group < 0 || group >= mesh.G) {
    throw std::out_of_range("Solver::constructTransportLinear: group index out of range");
  }
  if (ordinate_index < 0 || ordinate_index >= static_cast<int>(bc.left.size())) {
    throw std::out_of_range("Solver::constructTransportLinear: ordinate_index out of range");
  }

  const int n_x = mesh.n_x;
  if (external_source.size() != 4 * n_x) {
    throw std::invalid_argument(
        "Solver::constructTransportLinear: external_source has the wrong size");
  }

  // This group's own upper-energy-edge stopping power -- the flux this
  // group receives from the previous (higher-energy, already-solved)
  // group. Row `group` (not group + 1, which is the LHS's lower-edge
  // value).
  const std::vector<double>& stop_power_bound_up = cross_section.stop_power_boundary[group];

  const double m00 = fe_space.M(0, 0);
  const double m01 = fe_space.M(0, 1);
  const double m10 = fe_space.M(1, 0);
  const double m11 = fe_space.M(1, 1);

  b = Eigen::VectorXd::Zero(4 * n_x);

  for (int i = 0; i < n_x; ++i) {
    const double dx = mesh.dx[i];
    const int idx_left_down = localIndex(i, Corner::LeftDown, corner_order);
    const int idx_left_up = localIndex(i, Corner::LeftUp, corner_order);
    const int idx_right_down = localIndex(i, Corner::RightDown, corner_order);
    const int idx_right_up = localIndex(i, Corner::RightUp, corner_order);

    // external source
    const double q_left_down = external_source(idx_left_down);
    const double q_left_up = external_source(idx_left_up);
    const double q_right_down = external_source(idx_right_down);
    const double q_right_up = external_source(idx_right_up);

    b(idx_left_up) += (dx / 6.0) * (m00 * (q_left_down + 2.0 * q_left_up) +
                                    m01 * (q_right_down + 2.0 * q_right_up));
    b(idx_right_up) += (dx / 6.0) * (m10 * (q_left_down + 2.0 * q_left_up) +
                                     m11 * (q_right_down + 2.0 * q_right_up));
    b(idx_left_down) += (dx / 6.0) * (m00 * (2.0 * q_left_down + q_left_up) +
                                      m01 * (2.0 * q_right_down + q_right_up));
    b(idx_right_down) += (dx / 6.0) * (m10 * (2.0 * q_left_down + q_left_up) +
                                       m11 * (2.0 * q_right_down + q_right_up));

    // CSD source: inflow from the previous (higher-energy) group. For
    // group == 0, upwind_angular_flux is expected to be an all-zero view,
    // so this naturally contributes nothing rather than needing a
    // special case here.
    const double dE = mesh.dE[group];
    const ConstCornerValues upwind = upwind_angular_flux[i];
    b(idx_left_up) +=
        (dx / dE) * stop_power_bound_up[i] * (m00 * upwind.leftDown() + m01 * upwind.rightDown());
    b(idx_right_up) +=
        (dx / dE) * stop_power_bound_up[i] * (m10 * upwind.leftDown() + m11 * upwind.rightDown());

    // scattering source: sum over source groups gp with a nonzero transfer
    // into `group`. Only in-group (gp == group) uses latest_scalar_flux
    // rather than scalar_flux[gp] -- see the doc comment on this function.
    for (const ScatterEntry& entry : cross_section.scattering[i]) {
      if (entry.to != group) {
        continue;
      }
      const int gp = entry.from;
      const ConstCornerValues sc = (gp == group) ? latest_scalar_flux[i] : scalar_flux[gp][i];
      const double sc_left = sc.leftDown() + sc.leftUp();
      const double sc_right = sc.rightDown() + sc.rightUp();

      const double contribution_left =
          (dx * mesh.dE[gp] / 8.0) * entry.value * (m00 * sc_left + m01 * sc_right);
      const double contribution_right =
          (dx * mesh.dE[gp] / 8.0) * entry.value * (m10 * sc_left + m11 * sc_right);

      b(idx_left_down) += contribution_left;
      b(idx_left_up) += contribution_left;
      b(idx_right_down) += contribution_right;
      b(idx_right_up) += contribution_right;
    }
  }

  // boundary condition, at whichever edge mu points away from (the edge mu
  // points into is handled by the LHS's upwind coupling instead).
  if (mu > 0.0) {
    const DownUp& incoming = bc.left[ordinate_index][group];
    const int idx_left_down = localIndex(0, Corner::LeftDown, corner_order);
    const int idx_left_up = localIndex(0, Corner::LeftUp, corner_order);
    b(idx_left_up) += (mu / 6.0) * (incoming.down + 2.0 * incoming.up);
    b(idx_left_down) += (mu / 6.0) * (2.0 * incoming.down + incoming.up);
  } else {
    const DownUp& incoming = bc.right[ordinate_index][group];
    const int idx_right_up = localIndex(n_x - 1, Corner::RightUp, corner_order);
    const int idx_right_down = localIndex(n_x - 1, Corner::RightDown, corner_order);
    b(idx_right_up) += -(mu / 6.0) * (incoming.down + 2.0 * incoming.up);
    b(idx_right_down) += (-mu / 6.0) * (2.0 * incoming.down + incoming.up);
  }
}

void Solver::sweep(Eigen::VectorXd& x,Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b, double mu, int group) {
  // Construct bilinear
  constructTransportBilinear(A, mu, group);

  // solveLinearSystem
  x = solveLinearSystem(A, b);
}

// void Solver::sweep(int group, Field& scalar_flux, std::vector<Field>& angular_flux,
//                    Field::ConstRow latest_scalar_flux,
//                    const std::vector<Eigen::VectorXd>& external_source) const {
//   const Mesh& mesh = *input_deck.mesh;
//   const AngularQuadrature& quadrature = requireAngularQuadrature();
//
//   if (group < 0 || group >= mesh.G) {
//     throw std::out_of_range("Solver::sweep: group index out of range");
//   }
//
//   const int num_ordinates = static_cast<int>(quadrature.mu.size());
//   if (static_cast<int>(angular_flux.size()) != num_ordinates) {
//     throw std::invalid_argument("Solver::sweep: angular_flux has the wrong number of ordinates");
//   }
//   if (static_cast<int>(external_source.size()) != num_ordinates) {
//     throw std::invalid_argument("Solver::sweep: external_source has the wrong number of ordinates");
//   }
//
//   const int n_x = mesh.n_x;
//
//   zeroRow(scalar_flux[group], n_x);
//
//   // group == 0 has no previous group; substitute an all-zero view rather
//   // than requiring the caller to supply one.
//   std::optional<Field> zero_upwind;
//   if (group == 0) {
//     zero_upwind.emplace(n_x, 1);
//   }
//
//   Eigen::SparseMatrix<double> A;
//   Eigen::VectorXd b;
//
//   for (int m = 0; m < num_ordinates; ++m) {
//     const double mu = quadrature.mu[m];
//     const double w = quadrature.w[m];
//
//     const Field& angular_flux_m = angular_flux[m];
//     const Field::ConstRow upwind =
//         (group == 0) ? static_cast<const Field&>(*zero_upwind)[0] : angular_flux_m[group - 1];
//
//     constructTransportBilinear(A, mu, group);
//     constructTransportLinear(b, mu, group, m, upwind, scalar_flux, latest_scalar_flux,
//                              external_source[m]);
//
//     const Eigen::VectorXd x = solveLinearSystem(A, b);
//
//     assignLocalVectorToRow(angular_flux[m][group], n_x, corner_order, x);
//     accumulateLocalVectorIntoRow(scalar_flux[group], n_x, corner_order, x, w);
//   }
// }
