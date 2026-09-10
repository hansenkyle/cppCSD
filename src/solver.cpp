#include "solver.h"

#include <array>
#include <stdexcept>
#include <vector>

Solver::Solver(const Mesh& mesh, const CrossSection& cross_section, const FESpace& fe_space)
    : mesh(mesh), cross_section(this->mesh, cross_section.total, cross_section.scattering,
                                cross_section.stop_power, cross_section.stop_power_boundary,
                                cross_section.material),
      fe_space(fe_space) {}

void Solver::sweep() {}

namespace {

// One of the four corners of a linear-discontinuous (cell, group) element:
// low/high in x crossed with low/high in energy. Local to assembly -- the
// order below is the only one that exists, nothing outside this file needs
// to agree on it.
enum class Corner { LeftDown, LeftUp, RightDown, RightUp };

constexpr std::array<Corner, 4> kLocalCorners = {Corner::LeftDown, Corner::LeftUp,
                                                 Corner::RightDown, Corner::RightUp};

int localIndex(int cell, Corner corner) { return cell * 4 + static_cast<int>(corner); }

// Scatters a cell's full 4x4 local block (streaming's local part plus
// absorption/slowing-down, both local to this cell) into triplets.
void appendSelfBlock(std::vector<Eigen::Triplet<double>>& triplets, int cell,
                     const Eigen::Matrix4d& block) {
  for (int r = 0; r < 4; ++r) {
    const int row = localIndex(cell, kLocalCorners[r]);
    for (int c = 0; c < 4; ++c) {
      const int col = localIndex(cell, kLocalCorners[c]);
      triplets.emplace_back(row, col, block(r, c));
    }
  }
}

} // namespace

void Solver::constructTransportBilinear(Eigen::SparseMatrix<double>& A, double mu,
                                        int group) const {
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

  const double m00 = fe_space.M.left.left;
  const double m01 = fe_space.M.left.right;
  const double m10 = fe_space.M.right.left;
  const double m11 = fe_space.M.right.right;

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

    appendSelfBlock(triplets, i, self_block);

    // streaming, upwind coupling part -- 2x2, only when a neighbor exists
    // (at a domain boundary, the missing inflow is a boundary condition on
    // the RHS instead, not modeled here).
    if (mu > 0.0) {
      if (i > 0) {
        Eigen::Matrix2d coupling;
        coupling << (-1.0 / 6.0) * mu, (-1.0 / 3.0) * mu, (-1.0 / 3.0) * mu, (-1.0 / 6.0) * mu;
        const int row_left_up = localIndex(i, Corner::LeftUp);
        const int row_left_down = localIndex(i, Corner::LeftDown);
        const int col_right_down = localIndex(i - 1, Corner::RightDown);
        const int col_right_up = localIndex(i - 1, Corner::RightUp);
        triplets.emplace_back(row_left_up, col_right_down, coupling(0, 0));
        triplets.emplace_back(row_left_up, col_right_up, coupling(0, 1));
        triplets.emplace_back(row_left_down, col_right_down, coupling(1, 0));
        triplets.emplace_back(row_left_down, col_right_up, coupling(1, 1));
      }
    } else {
      if (i < n_x - 1) {
        Eigen::Matrix2d coupling;
        coupling << (1.0 / 6.0) * mu, (1.0 / 3.0) * mu, (1.0 / 3.0) * mu, (1.0 / 6.0) * mu;
        const int row_right_up = localIndex(i, Corner::RightUp);
        const int row_right_down = localIndex(i, Corner::RightDown);
        const int col_left_down = localIndex(i + 1, Corner::LeftDown);
        const int col_left_up = localIndex(i + 1, Corner::LeftUp);
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
