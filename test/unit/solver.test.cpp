#include "solver.h"

#include <fstream>
#include <sstream>

#include <doctest.h>

namespace {

// A minimally-valid, 1-cell/1-group/1-ordinate deck.
InputDeck makeDeck() {
  InputDeck deck;
  deck.mesh.n_x = 1;
  deck.mesh.x_boundary = Eigen::Vector2d(0.0, 1.0);
  deck.energy.G = 1;
  deck.energy.E_boundary = Eigen::Vector2d(1.0, 0.0);
  deck.angle.M = 1;
  deck.angle.mu = Eigen::VectorXd::Constant(1, 0.5);
  deck.angle.w = Eigen::VectorXd::Constant(1, 2.0);
  deck.xs.total = Eigen::MatrixXd::Constant(1, 1, 1.0);
  deck.xs.scatter = {{}};
  deck.xs.S = Eigen::MatrixXd::Constant(1, 1, 1.0);
  deck.xs.S_bound = Eigen::MatrixXd::Constant(2, 1, 1.0);
  deck.bc.values = Eigen::MatrixXd::Constant(2, 1, 0.0);
  deck.validate();
  return deck;
}

// Reads a whole file into a string.
std::string readFile(const std::filesystem::path& path) {
  std::ifstream in(path);
  std::ostringstream oss;
  oss << in.rdbuf();
  return oss.str();
}

// A fresh path with nothing there yet, so each test starts from a clean
// (nonexistent) file regardless of what earlier runs left behind.
std::filesystem::path freshPath(const std::string& name) {
  const std::filesystem::path path = std::filesystem::temp_directory_path() / name;
  std::filesystem::remove(path);
  return path;
}

} // namespace

TEST_SUITE("Solver") {
  TEST_CASE("writeMetadata creates the file if it doesn't exist, and appends on later calls") {
    const std::filesystem::path path = freshPath("ldcsd_solver_test_metadata_append.txt");
    Solver solver(makeDeck());

    solver.writeMetadata(path);
    const std::string after_first = readFile(path);
    CHECK(after_first.find("Run Metadata") != std::string::npos);

    solver.writeMetadata(path);
    const std::string after_second = readFile(path);

    // The second call's block is appended after the first's, not written
    // over it -- so the metadata header now appears twice.
    CHECK(after_second.size() > after_first.size());
    CHECK(after_second.substr(0, after_first.size()) == after_first);
    const std::size_t first_pos = after_second.find("Run Metadata");
    const std::size_t second_pos = after_second.find("Run Metadata", first_pos + 1);
    CHECK(second_pos != std::string::npos);

    std::filesystem::remove(path);
  }

  TEST_CASE("writeInputDeckEcho appends after existing content rather than truncating it") {
    const std::filesystem::path path = freshPath("ldcsd_solver_test_echo_append.txt");
    Solver solver(makeDeck());

    solver.writeMetadata(path);
    const std::string after_metadata = readFile(path);

    solver.writeInputDeckEcho(path);
    const std::string after_echo = readFile(path);

    CHECK(after_echo.substr(0, after_metadata.size()) == after_metadata);
    CHECK(after_echo.find("Problem Size") != std::string::npos);

    std::filesystem::remove(path);
  }

  TEST_CASE("writeResults appends after existing content rather than truncating it") {
    const std::filesystem::path path = freshPath("ldcsd_solver_test_results_append.txt");
    Solver solver(makeDeck());

    solver.writeMetadata(path);
    const std::string before = readFile(path);

    solver.writeResults(path, Solver::Results{Eigen::MatrixXd::Zero(4, 1), {}});
    const std::string after = readFile(path);

    CHECK(after.substr(0, before.size()) == before);
    CHECK(after.find("Scalar Flux") != std::string::npos);

    std::filesystem::remove(path);
  }

  TEST_CASE("writeResiduals appends after existing content rather than truncating it") {
    const std::filesystem::path path = freshPath("ldcsd_solver_test_residuals_append.txt");
    Solver solver(makeDeck());

    solver.writeResults(path, Solver::Results{Eigen::MatrixXd::Zero(4, 1), {}});
    const std::string before = readFile(path);

    Eigen::MatrixXd residual = Eigen::MatrixXd::Constant(4, 1, 1e-8);
    solver.writeResiduals(path, {residual});
    const std::string after = readFile(path);

    CHECK(after.substr(0, before.size()) == before);
    CHECK(after.find("Angular Residual - Ordinate 0") != std::string::npos);

    std::filesystem::remove(path);
  }

  TEST_CASE("calling all four writers in order builds the full report in that order") {
    const std::filesystem::path path = freshPath("ldcsd_solver_test_full_report.txt");
    Solver solver(makeDeck());

    solver.writeMetadata(path);
    solver.writeInputDeckEcho(path);
    solver.writeResults(path, Solver::Results{Eigen::MatrixXd::Zero(4, 1), {}});
    solver.writeResiduals(path, {Eigen::MatrixXd::Zero(4, 1)});

    const std::string result = readFile(path);
    const std::size_t metadata_pos = result.find("Run Metadata");
    const std::size_t echo_pos = result.find("Problem Size");
    const std::size_t results_pos = result.find("Scalar Flux");
    const std::size_t residuals_pos = result.find("Angular Residual");

    REQUIRE(metadata_pos != std::string::npos);
    REQUIRE(echo_pos != std::string::npos);
    REQUIRE(results_pos != std::string::npos);
    REQUIRE(residuals_pos != std::string::npos);
    CHECK(metadata_pos < echo_pos);
    CHECK(echo_pos < results_pos);
    CHECK(results_pos < residuals_pos);

    std::filesystem::remove(path);
  }
}
