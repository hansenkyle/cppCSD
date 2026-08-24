#include "cross_section.h"

#include <stdexcept>

#include <doctest.h>

namespace {

// A 2-cell, 2-group mesh shared by these tests.
Mesh makeMesh() { return Mesh({0.0, 1.0, 3.0}, {5.0, 2.0, 0.0}); }

} // namespace

TEST_SUITE("CrossSection") {
  TEST_CASE("stores per-cell, per-group data indexed [group][cell]") {
    const Mesh mesh = makeMesh();
    const CrossSection xs(mesh, {{1.0, 1.5}, {2.0, 2.5}}, {{0.1, 0.2}, {0.3, 0.4}},
                          {{3.0, 3.5}, {4.0, 4.5}}, {{5.0, 5.5}, {6.0, 6.5}, {7.0, 7.5}},
                          {"water", "lead"});

    CHECK(xs.total[0] == std::vector<double>{1.0, 1.5});
    CHECK(xs.total[1] == std::vector<double>{2.0, 2.5});
    CHECK(xs.scattering[0] == std::vector<double>{0.1, 0.2});
    CHECK(xs.stop_power[1] == std::vector<double>{4.0, 4.5});
    CHECK(xs.stop_power_boundary[2] == std::vector<double>{7.0, 7.5});
    CHECK(xs.material == std::vector<std::string>{"water", "lead"});
  }

  TEST_CASE("rejects total with the wrong number of group rows") {
    const Mesh mesh = makeMesh();
    CHECK_THROWS_AS(CrossSection(mesh, {{1.0, 1.5}}, {{0.1, 0.2}, {0.3, 0.4}},
                                 {{3.0, 3.5}, {4.0, 4.5}}, {{5.0, 5.5}, {6.0, 6.5}, {7.0, 7.5}},
                                 {"water", "lead"}),
                    std::invalid_argument);
  }

  TEST_CASE("rejects a group row with the wrong number of cells") {
    const Mesh mesh = makeMesh();
    CHECK_THROWS_AS(CrossSection(mesh, {{1.0}, {2.0, 2.5}}, {{0.1, 0.2}, {0.3, 0.4}},
                                 {{3.0, 3.5}, {4.0, 4.5}}, {{5.0, 5.5}, {6.0, 6.5}, {7.0, 7.5}},
                                 {"water", "lead"}),
                    std::invalid_argument);
  }

  TEST_CASE("rejects stop_power_boundary with G rows instead of G + 1") {
    const Mesh mesh = makeMesh();
    CHECK_THROWS_AS(CrossSection(mesh, {{1.0, 1.5}, {2.0, 2.5}}, {{0.1, 0.2}, {0.3, 0.4}},
                                 {{3.0, 3.5}, {4.0, 4.5}}, {{5.0, 5.5}, {6.0, 6.5}},
                                 {"water", "lead"}),
                    std::invalid_argument);
  }

  TEST_CASE("rejects a negative value") {
    const Mesh mesh = makeMesh();
    CHECK_THROWS_AS(CrossSection(mesh, {{-1.0, 1.5}, {2.0, 2.5}}, {{0.1, 0.2}, {0.3, 0.4}},
                                 {{3.0, 3.5}, {4.0, 4.5}}, {{5.0, 5.5}, {6.0, 6.5}, {7.0, 7.5}},
                                 {"water", "lead"}),
                    std::invalid_argument);
  }

  TEST_CASE("rejects a material vector with the wrong size") {
    const Mesh mesh = makeMesh();
    CHECK_THROWS_AS(CrossSection(mesh, {{1.0, 1.5}, {2.0, 2.5}}, {{0.1, 0.2}, {0.3, 0.4}},
                                 {{3.0, 3.5}, {4.0, 4.5}}, {{5.0, 5.5}, {6.0, 6.5}, {7.0, 7.5}},
                                 {"water"}),
                    std::invalid_argument);
  }
}
