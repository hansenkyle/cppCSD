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

} // namespace

TEST_SUITE("Solver") {
  TEST_CASE("stores its own copy of mesh, not a reference to the caller's") {
    Mesh mesh = makeMesh();
    const CrossSection xs = makeCrossSection(mesh);

    const Solver solver(mesh, xs);

    CHECK(&solver.mesh != &mesh);
    CHECK(solver.mesh.n_x == mesh.n_x);
    CHECK(solver.mesh.G == mesh.G);
  }

  TEST_CASE("stores its own copy of cross_section's data") {
    const Mesh mesh = makeMesh();
    const CrossSection xs = makeCrossSection(mesh);

    const Solver solver(mesh, xs);

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

    const Solver solver(mesh, xs);

    CHECK(&solver.cross_section.mesh == &solver.mesh);
    CHECK(&solver.cross_section.mesh != &mesh);
  }

  TEST_CASE("outlives the Mesh and CrossSection it was constructed from") {
    auto mesh = std::make_unique<Mesh>(makeMesh());
    auto xs = std::make_unique<CrossSection>(makeCrossSection(*mesh));

    const Solver solver(*mesh, *xs);
    xs.reset();
    mesh.reset();

    CHECK(solver.mesh.n_x == 2);
    CHECK(&solver.cross_section.mesh == &solver.mesh);
  }
}

TEST_SUITE("SourceIterationSolver") {
  TEST_CASE("constructs via the inherited Solver constructor") {
    const Mesh mesh = makeMesh();
    const CrossSection xs = makeCrossSection(mesh);

    const SourceIterationSolver solver(mesh, xs);

    CHECK(solver.mesh.n_x == mesh.n_x);
  }
}

TEST_SUITE("SecondMomentSolver") {
  TEST_CASE("constructs via the inherited Solver constructor") {
    const Mesh mesh = makeMesh();
    const CrossSection xs = makeCrossSection(mesh);

    const SecondMomentSolver solver(mesh, xs);

    CHECK(solver.mesh.n_x == mesh.n_x);
  }
}
