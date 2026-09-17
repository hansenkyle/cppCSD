#include "solver.h"

#include <array>
#include <cstdlib>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <vector>

#include "solver_formatter.h"

namespace {

void appendToFile(const std::filesystem::path& file_path, const std::string& text) {
  std::ofstream out(file_path, std::ios::app);
  if (!out.is_open()) {
    throw std::runtime_error("Solver: failed to open '" + file_path.string() + "' for writing");
  }
  out << text;
}

} // namespace

void Solver::writeMetadata(const std::filesystem::path& file_path) const {
  appendToFile(file_path, SolverFormatter::formatRunMetadata());
}

void Solver::writeInputDeckEcho(const std::filesystem::path& file_path) const {
  appendToFile(file_path, input_deck.echo());
}

void Solver::writeResults(const std::filesystem::path& file_path, const Results& results) const {
  appendToFile(file_path, SolverFormatter::formatResults(results, input_deck));
}

void Solver::writeResiduals(const std::filesystem::path& file_path,
                            const std::vector<Eigen::MatrixXd>& residuals) const {
  appendToFile(file_path, SolverFormatter::formatResiduals(residuals, input_deck));
}
