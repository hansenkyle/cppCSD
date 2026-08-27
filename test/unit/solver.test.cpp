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
  using Solver::Solver;
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
