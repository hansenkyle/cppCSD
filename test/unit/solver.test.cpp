#include "solver.h"

#include <memory>

#include <doctest.h>

namespace {

Mesh makeMesh() { return Mesh({0.0, 1.0, 3.0}, {5.0, 2.0, 0.0}); }

CrossSection makeCrossSection(const Mesh& mesh) {
  return CrossSection(mesh, {{1.0, 1.5}, {2.0, 2.5}}, {{0.1, 0.2}, {0.3, 0.4}},
                      {{3.0, 3.5}, {4.0, 4.5}}, {{5.0, 5.5}, {6.0, 6.5}, {7.0, 7.5}},
                      {"water", "lead"});
}

FESpace makeFESpace() { return FESpace(); }

// Exposes Solver's protected members for testing.
class TestSolver : public Solver {
public:
  using Solver::constructTransportBilinear;
  using Solver::Solver;
};

} // namespace

TEST_SUITE("Solver") {
  TEST_CASE("stores its own copy of mesh, not a reference to the caller's") {
    Mesh mesh = makeMesh();
    const CrossSection xs = makeCrossSection(mesh);
    const FESpace fe_space = makeFESpace();

    const Solver solver(mesh, xs, fe_space);

    CHECK(&solver.mesh != &mesh);
    CHECK(solver.mesh.n_x == mesh.n_x);
    CHECK(solver.mesh.G == mesh.G);
  }

  TEST_CASE("stores its own copy of cross_section's data") {
    const Mesh mesh = makeMesh();
    const CrossSection xs = makeCrossSection(mesh);
    const FESpace fe_space = makeFESpace();

    const Solver solver(mesh, xs, fe_space);

    CHECK(&solver.cross_section != &xs);
    CHECK(solver.cross_section.total == xs.total);
    CHECK(solver.cross_section.scattering == xs.scattering);
    CHECK(solver.cross_section.stop_power == xs.stop_power);
    CHECK(solver.cross_section.stop_power_boundary == xs.stop_power_boundary);
    CHECK(solver.cross_section.material == xs.material);
  }

  TEST_CASE("rebinds cross_section's Mesh reference to its own mesh, not the caller's") {
    Mesh mesh = makeMesh();
    const CrossSection xs = makeCrossSection(mesh);
    const FESpace fe_space = makeFESpace();

    const Solver solver(mesh, xs, fe_space);

    CHECK(&solver.cross_section.mesh == &solver.mesh);
    CHECK(&solver.cross_section.mesh != &mesh);
  }

  TEST_CASE("outlives the Mesh and CrossSection it was constructed from") {
    auto mesh = std::make_unique<Mesh>(makeMesh());
    auto xs = std::make_unique<CrossSection>(makeCrossSection(*mesh));
    const FESpace fe_space = makeFESpace();

    const Solver solver(*mesh, *xs, fe_space);
    xs.reset();
    mesh.reset();

    CHECK(solver.mesh.n_x == 2);
    CHECK(&solver.cross_section.mesh == &solver.mesh);
  }

  TEST_CASE("stores its own copy of fe_space") {
    const Mesh mesh = makeMesh();
    const CrossSection xs = makeCrossSection(mesh);
    const FESpace fe_space = makeFESpace();

    const Solver solver(mesh, xs, fe_space);

    CHECK(&solver.fe_space != &fe_space);
    CHECK(solver.fe_space.mass_matrix_kind == fe_space.mass_matrix_kind);
  }
}

TEST_SUITE("SourceIterationSolver") {
  TEST_CASE("constructs via the inherited Solver constructor") {
    const Mesh mesh = makeMesh();
    const CrossSection xs = makeCrossSection(mesh);
    const FESpace fe_space = makeFESpace();

    const SourceIterationSolver solver(mesh, xs, fe_space);

    CHECK(solver.mesh.n_x == mesh.n_x);
  }
}

TEST_SUITE("SecondMomentSolver") {
  TEST_CASE("constructs via the inherited Solver constructor") {
    const Mesh mesh = makeMesh();
    const CrossSection xs = makeCrossSection(mesh);
    const FESpace fe_space = makeFESpace();

    const SecondMomentSolver solver(mesh, xs, fe_space);

    CHECK(solver.mesh.n_x == mesh.n_x);
  }
}

TEST_SUITE("Solver::constructTransportBilinear") {
  // 2-cell, 1-group mesh, uniform cross sections -- small enough to verify
  // specific matrix entries against the same closed-form coefficients the
  // implementation uses, computed independently in each test.
  Mesh makeBilinearMesh() { return Mesh({0.0, 1.0, 2.0}, {2.0, 0.0}); }

  CrossSection makeBilinearCrossSection(const Mesh& mesh) {
    return CrossSection(mesh, {{1.0, 1.0}}, {{0.0, 0.0}}, {{0.5, 0.5}}, {{0.2, 0.2}, {0.4, 0.4}},
                        {"water", "water"});
  }

  TEST_CASE("produces a 4*n_x x 4*n_x matrix") {
    const Mesh mesh = makeBilinearMesh();
    const CrossSection xs = makeBilinearCrossSection(mesh);
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(mesh, xs, fe_space);

    Eigen::SparseMatrix<double> A;
    solver.constructTransportBilinear(A, 0.5, 0);

    CHECK(A.rows() == 8);
    CHECK(A.cols() == 8);
  }

  TEST_CASE("rejects an out-of-range group") {
    const Mesh mesh = makeBilinearMesh();
    const CrossSection xs = makeBilinearCrossSection(mesh);
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(mesh, xs, fe_space);

    Eigen::SparseMatrix<double> A;
    CHECK_THROWS_AS(solver.constructTransportBilinear(A, 0.5, -1), std::out_of_range);
    CHECK_THROWS_AS(solver.constructTransportBilinear(A, 0.5, 1), std::out_of_range);
  }

  TEST_CASE("self-block entries match the closed-form streaming + absorption coefficients") {
    const Mesh mesh = makeBilinearMesh();
    const CrossSection xs = makeBilinearCrossSection(mesh);
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(mesh, xs, fe_space);
    const double mu = 0.5;

    Eigen::SparseMatrix<double> A;
    solver.constructTransportBilinear(A, mu, 0);

    // Cell 0's own data (uniform across both cells in this fixture).
    const double dx = mesh.dx[0];
    const double dE = mesh.dE[0];
    const double xs_total = xs.total[0][0];
    const double stop_power = xs.stop_power[0][0];
    const double stop_power_bound_down = xs.stop_power_boundary[1][0];
    const double v1 = dx * ((1.0 / 6.0) * xs_total + (1.0 / (2.0 * dE)) * stop_power);
    const double v2 = dx * ((1.0 / 3.0) * xs_total + (1.0 / (2.0 * dE)) * stop_power);
    const double v3 =
        dx * ((1.0 / 3.0) * xs_total + (1.0 / dE) * (-stop_power / 2.0 + stop_power_bound_down));
    const double m00 = fe_space.M.left.left;
    const double m01 = fe_space.M.left.right;

    // Local corner order: 0=LeftDown, 1=LeftUp, 2=RightDown, 3=RightUp.
    const int row_left_up = 0 * 4 + 1;
    const int col_left_down = 0 * 4 + 0;
    const int col_right_down = 0 * 4 + 2;

    CHECK(A.coeff(row_left_up, col_left_down) == doctest::Approx((1.0 / 12.0) * mu + v1 * m00));
    CHECK(A.coeff(row_left_up, col_right_down) == doctest::Approx((1.0 / 12.0) * mu + v1 * m01));

    const int row_left_down = 0 * 4 + 0;
    CHECK(A.coeff(row_left_down, col_left_down) == doctest::Approx((1.0 / 6.0) * mu + v3 * m00));
  }

  TEST_CASE("mu > 0 couples a cell's left corners to its left neighbor's right corners") {
    const Mesh mesh = makeBilinearMesh();
    const CrossSection xs = makeBilinearCrossSection(mesh);
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(mesh, xs, fe_space);
    const double mu = 0.5;

    Eigen::SparseMatrix<double> A;
    solver.constructTransportBilinear(A, mu, 0);

    // Local corner order: 0=LeftDown, 1=LeftUp, 2=RightDown, 3=RightUp.
    const int row_left_up = 1 * 4 + 1;
    const int row_left_down = 1 * 4 + 0;
    const int col_right_down = 0 * 4 + 2;
    const int col_right_up = 0 * 4 + 3;

    CHECK(A.coeff(row_left_up, col_right_down) == doctest::Approx((-1.0 / 6.0) * mu));
    CHECK(A.coeff(row_left_up, col_right_up) == doctest::Approx((-1.0 / 3.0) * mu));
    CHECK(A.coeff(row_left_down, col_right_down) == doctest::Approx((-1.0 / 3.0) * mu));
    CHECK(A.coeff(row_left_down, col_right_up) == doctest::Approx((-1.0 / 6.0) * mu));

    // Cell 0 has no left neighbor, so its LeftUp row gets only the 4
    // self-block entries -- no additional coupling entry.
    const int row0_left_up = 0 * 4 + 1;
    const Eigen::VectorXd row0_left_up_dense = A.row(row0_left_up);
    CHECK((row0_left_up_dense.array() != 0.0).count() == 4);
  }

  TEST_CASE("mu < 0 couples a cell's right corners to its right neighbor's left corners") {
    const Mesh mesh = makeBilinearMesh();
    const CrossSection xs = makeBilinearCrossSection(mesh);
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(mesh, xs, fe_space);
    const double mu = -0.3;

    Eigen::SparseMatrix<double> A;
    solver.constructTransportBilinear(A, mu, 0);

    // Local corner order: 0=LeftDown, 1=LeftUp, 2=RightDown, 3=RightUp.
    const int row_right_up = 0 * 4 + 3;
    const int row_right_down = 0 * 4 + 2;
    const int col_left_down = 1 * 4 + 0;
    const int col_left_up = 1 * 4 + 1;

    CHECK(A.coeff(row_right_up, col_left_down) == doctest::Approx((1.0 / 6.0) * mu));
    CHECK(A.coeff(row_right_up, col_left_up) == doctest::Approx((1.0 / 3.0) * mu));
    CHECK(A.coeff(row_right_down, col_left_down) == doctest::Approx((1.0 / 3.0) * mu));
    CHECK(A.coeff(row_right_down, col_left_up) == doctest::Approx((1.0 / 6.0) * mu));

    // Cell 1 has no right neighbor, so its RightUp row gets only the 4
    // self-block entries -- no additional coupling entry.
    const int row1_right_up = 1 * 4 + 3;
    const Eigen::VectorXd row1_right_up_dense = A.row(row1_right_up);
    CHECK((row1_right_up_dense.array() != 0.0).count() == 4);
  }
}
