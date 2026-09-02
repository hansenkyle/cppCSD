#include "solver.h"

#include <doctest.h>

namespace {

Mesh makeMesh() { return Mesh({0.0, 1.0, 3.0}, {5.0, 2.0, 0.0}); }

CrossSection makeCrossSection(const Mesh& mesh) {
  return CrossSection(mesh, {{1.0, 1.5}, {2.0, 2.5}}, {{{0, 0, 0.1}}, {{1, 1, 0.4}}},
                      {{3.0, 3.5}, {4.0, 4.5}}, {{5.0, 5.5}, {6.0, 6.5}, {7.0, 7.5}},
                      {"water", "lead"});
}

FESpace makeFESpace() { return FESpace(1, 1); }

// An InputDeck with just mesh and xs set -- enough for Solver's constructor
// and constructTransportBilinear, without going through a full YAML read().
InputDeck makeInputDeck() {
  InputDeck deck;
  deck.setMesh(makeMesh());
  deck.xs.emplace(makeCrossSection(*deck.mesh));
  return deck;
}

// Exposes Solver's protected members for testing.
class TestSolver : public Solver {
public:
  using Solver::constructTransportBilinear;
  using Solver::constructTransportLinear;
  using Solver::solveLinearSystem;
  using Solver::Solver;
  using Solver::sweep;
};

} // namespace

TEST_SUITE("Solver") {
  TEST_CASE("stores its own copy of input_deck") {
    InputDeck deck = makeInputDeck();
    const FESpace fe_space = makeFESpace();

    const Solver solver(deck, fe_space);

    CHECK(&solver.input_deck != &deck);
    REQUIRE(solver.input_deck.mesh.has_value());
    CHECK(solver.input_deck.mesh->n_x == deck.mesh->n_x);
    REQUIRE(solver.input_deck.xs.has_value());
    CHECK(solver.input_deck.xs->total == deck.xs->total);
  }

  TEST_CASE("rejects an input_deck with no mesh set") {
    const InputDeck deck; // fresh, mesh unset
    const FESpace fe_space = makeFESpace();

    CHECK_THROWS_AS(Solver(deck, fe_space), std::invalid_argument);
  }

  TEST_CASE("stores its own copy of fe_space and defaults to XMajor corner order") {
    InputDeck deck = makeInputDeck();
    const FESpace fe_space = makeFESpace();

    const Solver solver(deck, fe_space);

    CHECK(&solver.fe_space != &fe_space);
    CHECK(solver.fe_space.spatial_degree == fe_space.spatial_degree);
    CHECK(solver.corner_order == AxisOrder::XMajor);
  }

  TEST_CASE("stores the corner_order it's constructed with") {
    InputDeck deck = makeInputDeck();
    const FESpace fe_space = makeFESpace();

    const Solver solver(deck, fe_space, AxisOrder::EMajor);

    CHECK(solver.corner_order == AxisOrder::EMajor);
  }

  TEST_CASE("defaults to SparseLU as the linear solver") {
    InputDeck deck = makeInputDeck();
    const FESpace fe_space = makeFESpace();

    const Solver solver(deck, fe_space);

    CHECK(solver.linear_solver_kind == LinearSolverKind::SparseLU);
  }

  TEST_CASE("stores the linear_solver_kind it's constructed with") {
    InputDeck deck = makeInputDeck();
    const FESpace fe_space = makeFESpace();

    const Solver solver(deck, fe_space, AxisOrder::XMajor, LinearSolverKind::GMRES);

    CHECK(solver.linear_solver_kind == LinearSolverKind::GMRES);
  }

  TEST_CASE("input_deck can be reconfigured after construction") {
    InputDeck deck = makeInputDeck();
    const FESpace fe_space = makeFESpace();
    Solver solver(deck, fe_space);

    solver.input_deck.setAngularQuadrature(AngularQuadrature({-0.5, 0.5}, {1.0, 1.0}));

    REQUIRE(solver.input_deck.angular_quadrature.has_value());
    CHECK(solver.input_deck.angular_quadrature->mu == std::vector<double>{-0.5, 0.5});
  }
}

TEST_SUITE("SourceIterationSolver") {
  TEST_CASE("constructs via the inherited Solver constructor") {
    InputDeck deck = makeInputDeck();
    const FESpace fe_space = makeFESpace();

    const SourceIterationSolver solver(deck, fe_space);

    CHECK(solver.input_deck.mesh->n_x == deck.mesh->n_x);
  }
}

TEST_SUITE("SecondMomentSolver") {
  TEST_CASE("constructs via the inherited Solver constructor") {
    InputDeck deck = makeInputDeck();
    const FESpace fe_space = makeFESpace();

    const SecondMomentSolver solver(deck, fe_space);

    CHECK(solver.input_deck.mesh->n_x == deck.mesh->n_x);
  }
}

TEST_SUITE("Solver::constructTransportBilinear") {
  // 2-cell, 1-group mesh, uniform cross sections -- small enough to verify
  // specific matrix entries against the same closed-form coefficients the
  // implementation uses, computed independently in each test.
  Mesh makeBilinearMesh() { return Mesh({0.0, 1.0, 2.0}, {2.0, 0.0}); }

  CrossSection makeBilinearCrossSection(const Mesh& mesh) {
    return CrossSection(mesh, {{1.0, 1.0}}, {{}, {}}, {{0.5, 0.5}}, {{0.2, 0.2}, {0.4, 0.4}},
                        {"water", "water"});
  }

  InputDeck makeBilinearInputDeck() {
    InputDeck deck;
    deck.setMesh(makeBilinearMesh());
    deck.xs.emplace(makeBilinearCrossSection(*deck.mesh));
    return deck;
  }

  TEST_CASE("produces a 4*n_x x 4*n_x matrix") {
    InputDeck deck = makeBilinearInputDeck();
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(deck, fe_space);

    Eigen::SparseMatrix<double> A;
    solver.constructTransportBilinear(A, 0.5, 0);

    CHECK(A.rows() == 8);
    CHECK(A.cols() == 8);
  }

  TEST_CASE("rejects an out-of-range group") {
    InputDeck deck = makeBilinearInputDeck();
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(deck, fe_space);

    Eigen::SparseMatrix<double> A;
    CHECK_THROWS_AS(solver.constructTransportBilinear(A, 0.5, -1), std::out_of_range);
    CHECK_THROWS_AS(solver.constructTransportBilinear(A, 0.5, 1), std::out_of_range);
  }

  TEST_CASE("self-block entries match the closed-form streaming + absorption coefficients") {
    InputDeck deck = makeBilinearInputDeck();
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(deck, fe_space);
    const double mu = 0.5;

    Eigen::SparseMatrix<double> A;
    solver.constructTransportBilinear(A, mu, 0);

    // Cell 0's own data (uniform across both cells in this fixture).
    const Mesh& mesh = *deck.mesh;
    const CrossSection& xs = *deck.xs;
    const double dx = mesh.dx[0];
    const double dE = mesh.dE[0];
    const double xs_total = xs.total[0][0];
    const double stop_power = xs.stop_power[0][0];
    const double stop_power_bound_down = xs.stop_power_boundary[1][0];
    const double v1 = dx * ((1.0 / 6.0) * xs_total + (1.0 / (2.0 * dE)) * stop_power);
    const double v2 = dx * ((1.0 / 3.0) * xs_total + (1.0 / (2.0 * dE)) * stop_power);
    const double v3 =
        dx * ((1.0 / 3.0) * xs_total + (1.0 / dE) * (-stop_power / 2.0 + stop_power_bound_down));
    const double m00 = fe_space.M(0, 0);
    const double m01 = fe_space.M(0, 1);

    const int row_left_up = 0 * 4 + cornerSlot(Corner::LeftUp, solver.corner_order);
    const int col_left_down = 0 * 4 + cornerSlot(Corner::LeftDown, solver.corner_order);
    const int col_right_down = 0 * 4 + cornerSlot(Corner::RightDown, solver.corner_order);

    CHECK(A.coeff(row_left_up, col_left_down) == doctest::Approx((1.0 / 12.0) * mu + v1 * m00));
    CHECK(A.coeff(row_left_up, col_right_down) == doctest::Approx((1.0 / 12.0) * mu + v1 * m01));

    const int row_left_down = 0 * 4 + cornerSlot(Corner::LeftDown, solver.corner_order);
    CHECK(A.coeff(row_left_down, col_left_down) == doctest::Approx((1.0 / 6.0) * mu + v3 * m00));
  }

  TEST_CASE("mu > 0 couples a cell's left corners to its left neighbor's right corners") {
    InputDeck deck = makeBilinearInputDeck();
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(deck, fe_space);
    const double mu = 0.5;

    Eigen::SparseMatrix<double> A;
    solver.constructTransportBilinear(A, mu, 0);

    const AxisOrder order = solver.corner_order;
    const int row_left_up = 1 * 4 + cornerSlot(Corner::LeftUp, order);
    const int row_left_down = 1 * 4 + cornerSlot(Corner::LeftDown, order);
    const int col_right_down = 0 * 4 + cornerSlot(Corner::RightDown, order);
    const int col_right_up = 0 * 4 + cornerSlot(Corner::RightUp, order);

    CHECK(A.coeff(row_left_up, col_right_down) == doctest::Approx((-1.0 / 6.0) * mu));
    CHECK(A.coeff(row_left_up, col_right_up) == doctest::Approx((-1.0 / 3.0) * mu));
    CHECK(A.coeff(row_left_down, col_right_down) == doctest::Approx((-1.0 / 3.0) * mu));
    CHECK(A.coeff(row_left_down, col_right_up) == doctest::Approx((-1.0 / 6.0) * mu));

    // Cell 0 has no left neighbor, so its LeftUp row gets only the 4
    // self-block entries -- no additional coupling entry.
    const int row0_left_up = 0 * 4 + cornerSlot(Corner::LeftUp, order);
    const Eigen::VectorXd row0_left_up_dense = A.row(row0_left_up);
    CHECK((row0_left_up_dense.array() != 0.0).count() == 4);
  }

  TEST_CASE("mu < 0 couples a cell's right corners to its right neighbor's left corners") {
    InputDeck deck = makeBilinearInputDeck();
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(deck, fe_space);
    const double mu = -0.3;

    Eigen::SparseMatrix<double> A;
    solver.constructTransportBilinear(A, mu, 0);

    const AxisOrder order = solver.corner_order;
    const int row_right_up = 0 * 4 + cornerSlot(Corner::RightUp, order);
    const int row_right_down = 0 * 4 + cornerSlot(Corner::RightDown, order);
    const int col_left_down = 1 * 4 + cornerSlot(Corner::LeftDown, order);
    const int col_left_up = 1 * 4 + cornerSlot(Corner::LeftUp, order);

    CHECK(A.coeff(row_right_up, col_left_down) == doctest::Approx((1.0 / 6.0) * mu));
    CHECK(A.coeff(row_right_up, col_left_up) == doctest::Approx((1.0 / 3.0) * mu));
    CHECK(A.coeff(row_right_down, col_left_down) == doctest::Approx((1.0 / 3.0) * mu));
    CHECK(A.coeff(row_right_down, col_left_up) == doctest::Approx((1.0 / 6.0) * mu));

    // Cell 1 has no right neighbor, so its RightUp row gets only the 4
    // self-block entries -- no additional coupling entry.
    const int row1_right_up = 1 * 4 + cornerSlot(Corner::RightUp, order);
    const Eigen::VectorXd row1_right_up_dense = A.row(row1_right_up);
    CHECK((row1_right_up_dense.array() != 0.0).count() == 4);
  }
}

TEST_SUITE("Solver::constructTransportLinear") {
  // 2-cell, 2-group mesh -- group 1 is "the group being solved" throughout,
  // with group 0 available as the other source group for scattering tests.
  // Each test isolates one RHS term by zeroing the others (empty
  // scattering, an all-zero upwind view, or a zero boundary condition).
  Mesh makeLinearMesh() { return Mesh({0.0, 1.0, 2.0}, {2.0, 1.0, 0.0}); }

  CrossSection makeLinearCrossSection(const Mesh& mesh, std::vector<ScatterEntry> cell0_scattering,
                                      std::vector<ScatterEntry> cell1_scattering) {
    return CrossSection(
        mesh, {{1.0, 1.0}, {1.0, 1.0}}, {std::move(cell0_scattering), std::move(cell1_scattering)},
        {{0.5, 0.5}, {0.5, 0.5}}, {{0.1, 0.1}, {0.2, 0.2}, {0.3, 0.3}}, {"water", "water"});
  }

  InputDeck makeLinearInputDeck(std::vector<ScatterEntry> cell0_scattering = {},
                                std::vector<ScatterEntry> cell1_scattering = {}) {
    InputDeck deck;
    deck.setMesh(makeLinearMesh());
    deck.xs.emplace(makeLinearCrossSection(*deck.mesh, std::move(cell0_scattering),
                                           std::move(cell1_scattering)));
    // 2 ordinates x 2 groups, all zero unless a test overrides it.
    std::vector<std::vector<DownUp>> zero_bc(2, std::vector<DownUp>(2, DownUp{0.0, 0.0}));
    deck.setBoundaryConditions(BoundaryConditions(zero_bc, zero_bc, 2, 2));
    return deck;
  }

  TEST_CASE("produces a 4*n_x vector") {
    InputDeck deck = makeLinearInputDeck();
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(deck, fe_space);
    const Field zero_flux(2, 1);
    const Field scalar_flux(2, 2);
    const Eigen::VectorXd external_source = Eigen::VectorXd::Zero(8);

    Eigen::VectorXd b;
    solver.constructTransportLinear(b, 0.5, 1, 0, zero_flux[0], scalar_flux, zero_flux[0],
                                    external_source);

    CHECK(b.size() == 8);
  }

  TEST_CASE("rejects an out-of-range group") {
    InputDeck deck = makeLinearInputDeck();
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(deck, fe_space);
    const Field zero_flux(2, 1);
    const Field scalar_flux(2, 2);
    const Eigen::VectorXd external_source = Eigen::VectorXd::Zero(8);

    Eigen::VectorXd b;
    CHECK_THROWS_AS(solver.constructTransportLinear(b, 0.5, -1, 0, zero_flux[0], scalar_flux,
                                                    zero_flux[0], external_source),
                    std::out_of_range);
    CHECK_THROWS_AS(solver.constructTransportLinear(b, 0.5, 2, 0, zero_flux[0], scalar_flux,
                                                    zero_flux[0], external_source),
                    std::out_of_range);
  }

  TEST_CASE("rejects an out-of-range ordinate_index") {
    InputDeck deck = makeLinearInputDeck();
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(deck, fe_space);
    const Field zero_flux(2, 1);
    const Field scalar_flux(2, 2);
    const Eigen::VectorXd external_source = Eigen::VectorXd::Zero(8);

    Eigen::VectorXd b;
    CHECK_THROWS_AS(solver.constructTransportLinear(b, 0.5, 1, -1, zero_flux[0], scalar_flux,
                                                    zero_flux[0], external_source),
                    std::out_of_range);
    CHECK_THROWS_AS(solver.constructTransportLinear(b, 0.5, 1, 2, zero_flux[0], scalar_flux,
                                                    zero_flux[0], external_source),
                    std::out_of_range);
  }

  TEST_CASE("rejects an external_source with the wrong size") {
    InputDeck deck = makeLinearInputDeck();
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(deck, fe_space);
    const Field zero_flux(2, 1);
    const Field scalar_flux(2, 2);
    const Eigen::VectorXd wrong_size_source = Eigen::VectorXd::Zero(7);

    Eigen::VectorXd b;
    CHECK_THROWS_AS(solver.constructTransportLinear(b, 0.5, 1, 0, zero_flux[0], scalar_flux,
                                                    zero_flux[0], wrong_size_source),
                    std::invalid_argument);
  }

  TEST_CASE("external source term matches the closed-form coefficients") {
    InputDeck deck = makeLinearInputDeck();
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(deck, fe_space);
    const Field zero_flux(2, 1);
    const Field scalar_flux(2, 2);
    Eigen::VectorXd external_source(8);
    external_source << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0;

    Eigen::VectorXd b;
    solver.constructTransportLinear(b, 0.5, 1, 0, zero_flux[0], scalar_flux, zero_flux[0],
                                    external_source);

    const double dx = deck.mesh->dx[0];
    const double m00 = fe_space.M(0, 0);
    const double m01 = fe_space.M(0, 1);
    const double m10 = fe_space.M(1, 0);
    const double m11 = fe_space.M(1, 1);
    const double q_ld = external_source(0);
    const double q_lu = external_source(1);
    const double q_rd = external_source(2);
    const double q_ru = external_source(3);

    const double expected_lu = (dx / 6.0) * (m00 * (q_ld + 2.0 * q_lu) + m01 * (q_rd + 2.0 * q_ru));
    const double expected_ru = (dx / 6.0) * (m10 * (q_ld + 2.0 * q_lu) + m11 * (q_rd + 2.0 * q_ru));
    const double expected_ld = (dx / 6.0) * (m00 * (2.0 * q_ld + q_lu) + m01 * (2.0 * q_rd + q_ru));
    const double expected_rd = (dx / 6.0) * (m10 * (2.0 * q_ld + q_lu) + m11 * (2.0 * q_rd + q_ru));

    CHECK(b(1) == doctest::Approx(expected_lu));
    CHECK(b(3) == doctest::Approx(expected_ru));
    CHECK(b(0) == doctest::Approx(expected_ld));
    CHECK(b(2) == doctest::Approx(expected_rd));
  }

  TEST_CASE("CSD term matches the closed-form coefficients") {
    InputDeck deck = makeLinearInputDeck();
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(deck, fe_space);
    Field upwind_mutable(2, 1);
    upwind_mutable[0][0].leftDown() = 0.1;
    upwind_mutable[0][0].leftUp() = 0.2;
    upwind_mutable[0][0].rightDown() = 0.3;
    upwind_mutable[0][0].rightUp() = 0.4;
    const Field& upwind = upwind_mutable;
    const Field scalar_flux(2, 2);
    const Field zero_latest(2, 1);
    const Eigen::VectorXd external_source = Eigen::VectorXd::Zero(8);

    Eigen::VectorXd b;
    solver.constructTransportLinear(b, 0.5, 1, 0, upwind[0], scalar_flux, zero_latest[0],
                                    external_source);

    const double dx = deck.mesh->dx[0];
    const double dE = deck.mesh->dE[1];
    const double stop_power_bound_up = deck.xs->stop_power_boundary[1][0]; // row `group`
    const double m00 = fe_space.M(0, 0);
    const double m01 = fe_space.M(0, 1);
    const double m10 = fe_space.M(1, 0);
    const double m11 = fe_space.M(1, 1);

    const double expected_lu = (dx / dE) * stop_power_bound_up * (m00 * 0.1 + m01 * 0.3);
    const double expected_ru = (dx / dE) * stop_power_bound_up * (m10 * 0.1 + m11 * 0.3);

    CHECK(b(1) == doctest::Approx(expected_lu));
    CHECK(b(3) == doctest::Approx(expected_ru));
    // No CSD contribution to the down rows.
    CHECK(b(0) == doctest::Approx(0.0));
    CHECK(b(2) == doctest::Approx(0.0));
  }

  TEST_CASE("scattering term (cross-group) reads scalar_flux[from]") {
    InputDeck deck = makeLinearInputDeck({{0, 1, 0.05}}, {});
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(deck, fe_space);
    const Field zero_flux(2, 1);
    Field scalar_flux(2, 2);
    scalar_flux[0][0].leftDown() = 1.0;
    scalar_flux[0][0].leftUp() = 2.0;
    scalar_flux[0][0].rightDown() = 3.0;
    scalar_flux[0][0].rightUp() = 4.0;
    const Eigen::VectorXd external_source = Eigen::VectorXd::Zero(8);

    Eigen::VectorXd b;
    solver.constructTransportLinear(b, 0.5, 1, 0, zero_flux[0], scalar_flux, zero_flux[0],
                                    external_source);

    const double dx = deck.mesh->dx[0];
    const double dE_from = deck.mesh->dE[0]; // dE of the source group (from=0)
    const double xs = 0.05;
    const double m00 = fe_space.M(0, 0);
    const double m01 = fe_space.M(0, 1);
    const double m10 = fe_space.M(1, 0);
    const double m11 = fe_space.M(1, 1);
    const double sc_left = 1.0 + 2.0;
    const double sc_right = 3.0 + 4.0;

    const double expected_left = (dx * dE_from / 8.0) * xs * (m00 * sc_left + m01 * sc_right);
    const double expected_right = (dx * dE_from / 8.0) * xs * (m10 * sc_left + m11 * sc_right);

    CHECK(b(0) == doctest::Approx(expected_left));  // LeftDown
    CHECK(b(1) == doctest::Approx(expected_left));  // LeftUp
    CHECK(b(2) == doctest::Approx(expected_right)); // RightDown
    CHECK(b(3) == doctest::Approx(expected_right)); // RightUp
  }

  TEST_CASE("scattering term (in-group) reads latest_scalar_flux, not scalar_flux[group]") {
    InputDeck deck = makeLinearInputDeck({{1, 1, 0.3}}, {});
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(deck, fe_space);
    const Field zero_flux(2, 1);

    // scalar_flux[1] (group 1's stored, stale value) is a decoy: if the
    // implementation mistakenly used it instead of latest_scalar_flux, the
    // check below would fail.
    Field scalar_flux(2, 2);
    scalar_flux[1][0].leftDown() = 999.0;
    scalar_flux[1][0].leftUp() = 999.0;
    scalar_flux[1][0].rightDown() = 999.0;
    scalar_flux[1][0].rightUp() = 999.0;

    Field latest_mutable(2, 1);
    latest_mutable[0][0].leftDown() = 10.0;
    latest_mutable[0][0].leftUp() = 11.0;
    latest_mutable[0][0].rightDown() = 12.0;
    latest_mutable[0][0].rightUp() = 13.0;
    const Field& latest = latest_mutable;
    const Eigen::VectorXd external_source = Eigen::VectorXd::Zero(8);

    Eigen::VectorXd b;
    solver.constructTransportLinear(b, 0.5, 1, 0, zero_flux[0], scalar_flux, latest[0],
                                    external_source);

    const double dx = deck.mesh->dx[0];
    const double dE_from = deck.mesh->dE[1]; // dE of the source group (from=1, in-group)
    const double xs = 0.3;
    const double m00 = fe_space.M(0, 0);
    const double m01 = fe_space.M(0, 1);
    const double sc_left = 10.0 + 11.0;
    const double sc_right = 12.0 + 13.0;

    const double expected_left = (dx * dE_from / 8.0) * xs * (m00 * sc_left + m01 * sc_right);

    CHECK(b(0) == doctest::Approx(expected_left));
    CHECK(b(1) == doctest::Approx(expected_left));
  }

  TEST_CASE("boundary condition term (mu > 0, left)") {
    InputDeck deck = makeLinearInputDeck();
    std::vector<std::vector<DownUp>> zero_bc(2, std::vector<DownUp>(2, DownUp{0.0, 0.0}));
    std::vector<std::vector<DownUp>> left_bc = zero_bc;
    left_bc[0][1] = DownUp{100.0, 200.0};
    deck.setBoundaryConditions(BoundaryConditions(left_bc, zero_bc, 2, 2));

    const FESpace fe_space = makeFESpace();
    const TestSolver solver(deck, fe_space);
    const Field zero_flux(2, 1);
    const Field scalar_flux(2, 2);
    const Eigen::VectorXd external_source = Eigen::VectorXd::Zero(8);
    const double mu = 0.5;

    Eigen::VectorXd b;
    solver.constructTransportLinear(b, mu, 1, 0, zero_flux[0], scalar_flux, zero_flux[0],
                                    external_source);

    const double expected_lu = (mu / 6.0) * (100.0 + 2.0 * 200.0);
    const double expected_ld = (mu / 6.0) * (2.0 * 100.0 + 200.0);
    CHECK(b(1) == doctest::Approx(expected_lu));
    CHECK(b(0) == doctest::Approx(expected_ld));
    // The right boundary (cell 1) is untouched by mu > 0.
    CHECK(b(6) == doctest::Approx(0.0));
    CHECK(b(7) == doctest::Approx(0.0));
  }

  TEST_CASE("boundary condition term (mu < 0, right)") {
    InputDeck deck = makeLinearInputDeck();
    std::vector<std::vector<DownUp>> zero_bc(2, std::vector<DownUp>(2, DownUp{0.0, 0.0}));
    std::vector<std::vector<DownUp>> right_bc = zero_bc;
    right_bc[0][1] = DownUp{100.0, 200.0};
    deck.setBoundaryConditions(BoundaryConditions(zero_bc, right_bc, 2, 2));

    const FESpace fe_space = makeFESpace();
    const TestSolver solver(deck, fe_space);
    const Field zero_flux(2, 1);
    const Field scalar_flux(2, 2);
    const Eigen::VectorXd external_source = Eigen::VectorXd::Zero(8);
    const double mu = -0.5;

    Eigen::VectorXd b;
    solver.constructTransportLinear(b, mu, 1, 0, zero_flux[0], scalar_flux, zero_flux[0],
                                    external_source);

    const double expected_ru = -(mu / 6.0) * (100.0 + 2.0 * 200.0);
    const double expected_rd = (-mu / 6.0) * (2.0 * 100.0 + 200.0);
    // Cell 1 (index n_x - 1 = 1) occupies local indices 4..7; under the
    // default XMajor corner order, RightDown=slot 2, RightUp=slot 3.
    CHECK(b(7) == doctest::Approx(expected_ru));
    CHECK(b(6) == doctest::Approx(expected_rd));
    // The left boundary (cell 0) is untouched by mu < 0.
    CHECK(b(0) == doctest::Approx(0.0));
    CHECK(b(1) == doctest::Approx(0.0));
  }
}

TEST_SUITE("Solver::solveLinearSystem") {
  // 2x + y = 5, x + 3y = 10 -> x = 1, y = 3.
  Eigen::SparseMatrix<double> makeSimpleMatrix() {
    Eigen::SparseMatrix<double> A(2, 2);
    const std::vector<Eigen::Triplet<double>> triplets = {
        {0, 0, 2.0}, {0, 1, 1.0}, {1, 0, 1.0}, {1, 1, 3.0}};
    A.setFromTriplets(triplets.begin(), triplets.end());
    return A;
  }

  TEST_CASE("SparseLU solves a small system correctly") {
    InputDeck deck = makeInputDeck();
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(deck, fe_space, AxisOrder::XMajor, LinearSolverKind::SparseLU);

    const Eigen::SparseMatrix<double> A = makeSimpleMatrix();
    Eigen::VectorXd b(2);
    b << 5.0, 10.0;

    const Eigen::VectorXd x = solver.solveLinearSystem(A, b);
    CHECK(x(0) == doctest::Approx(1.0));
    CHECK(x(1) == doctest::Approx(3.0));
  }

  TEST_CASE("BiCGSTAB solves a small system correctly") {
    InputDeck deck = makeInputDeck();
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(deck, fe_space, AxisOrder::XMajor, LinearSolverKind::BiCGSTAB);

    const Eigen::SparseMatrix<double> A = makeSimpleMatrix();
    Eigen::VectorXd b(2);
    b << 5.0, 10.0;

    const Eigen::VectorXd x = solver.solveLinearSystem(A, b);
    CHECK(x(0) == doctest::Approx(1.0));
    CHECK(x(1) == doctest::Approx(3.0));
  }

  TEST_CASE("GMRES solves a small system correctly") {
    InputDeck deck = makeInputDeck();
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(deck, fe_space, AxisOrder::XMajor, LinearSolverKind::GMRES);

    const Eigen::SparseMatrix<double> A = makeSimpleMatrix();
    Eigen::VectorXd b(2);
    b << 5.0, 10.0;

    const Eigen::VectorXd x = solver.solveLinearSystem(A, b);
    CHECK(x(0) == doctest::Approx(1.0));
    CHECK(x(1) == doctest::Approx(3.0));
  }

  TEST_CASE("SweepDirect throws, since it isn't implemented yet") {
    InputDeck deck = makeInputDeck();
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(deck, fe_space, AxisOrder::XMajor, LinearSolverKind::SweepDirect);

    const Eigen::SparseMatrix<double> A = makeSimpleMatrix();
    Eigen::VectorXd b(2);
    b << 5.0, 10.0;

    CHECK_THROWS_AS(solver.solveLinearSystem(A, b), std::runtime_error);
  }

  TEST_CASE("x-reference overload matches the return-value overload") {
    InputDeck deck = makeInputDeck();
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(deck, fe_space, AxisOrder::XMajor, LinearSolverKind::SparseLU);

    const Eigen::SparseMatrix<double> A = makeSimpleMatrix();
    Eigen::VectorXd b(2);
    b << 5.0, 10.0;

    Eigen::VectorXd x;
    solver.solveLinearSystem(A, x, b);
    CHECK(x(0) == doctest::Approx(1.0));
    CHECK(x(1) == doctest::Approx(3.0));
  }

  TEST_CASE("x-reference overload overwrites a pre-existing value of x") {
    InputDeck deck = makeInputDeck();
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(deck, fe_space, AxisOrder::XMajor, LinearSolverKind::SparseLU);

    const Eigen::SparseMatrix<double> A = makeSimpleMatrix();
    Eigen::VectorXd b(2);
    b << 5.0, 10.0;

    Eigen::VectorXd x(2);
    x << 100.0, -100.0;
    solver.solveLinearSystem(A, x, b);
    CHECK(x(0) == doctest::Approx(1.0));
    CHECK(x(1) == doctest::Approx(3.0));
  }
}

TEST_SUITE("Solver::sweep") {
  TEST_CASE("builds A via constructTransportBilinear for this mu/group and solves A*x=b") {
    InputDeck deck = makeInputDeck();
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(deck, fe_space);

    const double mu = 0.5;
    const int group = 1;

    Eigen::SparseMatrix<double> expected_A;
    solver.constructTransportBilinear(expected_A, mu, group);
    const Eigen::VectorXd b = Eigen::VectorXd::LinSpaced(expected_A.rows(), 1.0, expected_A.rows());
    const Eigen::VectorXd expected_x = solver.solveLinearSystem(expected_A, b);

    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd x;
    solver.sweep(x, A, b, mu, group);

    REQUIRE(A.rows() == expected_A.rows());
    REQUIRE(A.cols() == expected_A.cols());
    CHECK((A - expected_A).norm() == doctest::Approx(0.0));

    REQUIRE(x.size() == expected_x.size());
    for (int i = 0; i < x.size(); ++i) {
      CHECK(x(i) == doctest::Approx(expected_x(i)));
    }
  }

  TEST_CASE("overwrites a pre-existing value of A rather than accumulating into it") {
    InputDeck deck = makeInputDeck();
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(deck, fe_space);

    const double mu = -0.5;
    const int group = 0;

    Eigen::SparseMatrix<double> expected_A;
    solver.constructTransportBilinear(expected_A, mu, group);
    const Eigen::VectorXd b = Eigen::VectorXd::LinSpaced(expected_A.rows(), 1.0, expected_A.rows());

    // Pre-fill A with a stale matrix of a different size to confirm sweep()
    // fully replaces it rather than adding to whatever was already there.
    Eigen::SparseMatrix<double> A(2, 2);
    A.insert(0, 0) = 123.0;

    Eigen::VectorXd x;
    solver.sweep(x, A, b, mu, group);

    REQUIRE(A.rows() == expected_A.rows());
    REQUIRE(A.cols() == expected_A.cols());
    CHECK((A - expected_A).norm() == doctest::Approx(0.0));
  }

  TEST_CASE("propagates an out-of-range group from constructTransportBilinear") {
    InputDeck deck = makeInputDeck();
    const FESpace fe_space = makeFESpace();
    const TestSolver solver(deck, fe_space);

    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd x;
    const Eigen::VectorXd b = Eigen::VectorXd::Zero(12);

    CHECK_THROWS_AS(solver.sweep(x, A, b, 0.5, -1), std::out_of_range);
    CHECK_THROWS_AS(solver.sweep(x, A, b, 0.5, 2), std::out_of_range);
  }
}
