#include "input_deck.h"

#include <fstream>

#include <doctest.h>

namespace {

std::filesystem::path sampleInputPath() {
  return std::filesystem::path(TEST_DATA_DIR) / "sample_input.yaml";
}

// Writes `contents` to a temp file and returns its path; used for cases that
// need to deviate from the checked-in sample input.
std::filesystem::path writeTempYaml(const std::string& contents) {
  const std::filesystem::path path =
      std::filesystem::temp_directory_path() / "ldcsd_input_deck_test_input.yaml";
  std::ofstream out(path);
  out << contents;
  return path;
}

} // namespace

TEST_SUITE("InputDeck") {
  TEST_CASE("reads a well-formed input deck") {
    InputDeck deck;
    CHECK(deck.read(sampleInputPath()) == 0);

    CHECK(deck.spatial_mesh == std::vector<double>{0.0, 1.0, 2.0, 3.0});
    CHECK(deck.region_materials == std::vector<std::string>{"water", "water", "lead"});
    CHECK(deck.energy_mesh == std::vector<double>{0.1, 0.5, 1.0, 5.0});

    REQUIRE(deck.materials.contains("water"));
    const MaterialData& water = deck.materials.at("water");
    CHECK(water.sigma_t == std::vector<double>{1.2, 1.0, 0.8});
    CHECK(water.sigma_s == std::vector<double>{1.1, 0.9, 0.7});
    CHECK(water.stopping_power_average == std::vector<double>{2.0, 1.8, 1.5});
    CHECK(water.stopping_power_boundary == std::vector<double>{2.2, 1.9, 1.6, 1.3});

    REQUIRE(deck.materials.contains("lead"));

    CHECK(deck.convergence.max_iters == 200);
    CHECK(deck.convergence.epsilon == doctest::Approx(1.0e-8));
  }

  TEST_CASE("rejects a non-ascending spatial mesh") {
    const auto path = writeTempYaml(R"(
spatial_mesh: [0.0, 2.0, 1.0]
regions:
  materials: [water, water]
energy_mesh: [0.1, 1.0]
materials:
  water:
    sigma_t: [1.0]
    sigma_s: [1.0]
    stopping_power:
      group_average: [1.0]
      group_boundary: [1.0, 1.0]
convergence:
  max_iters: 1
  epsilon: 1.0
)");
    InputDeck deck;
    CHECK(deck.read(path) == 1);
  }

  TEST_CASE("rejects a region count that does not match the spatial mesh") {
    const auto path = writeTempYaml(R"(
spatial_mesh: [0.0, 1.0, 2.0]
regions:
  materials: [water]
energy_mesh: [0.1, 1.0]
materials:
  water:
    sigma_t: [1.0]
    sigma_s: [1.0]
    stopping_power:
      group_average: [1.0]
      group_boundary: [1.0, 1.0]
convergence:
  max_iters: 1
  epsilon: 1.0
)");
    InputDeck deck;
    CHECK(deck.read(path) == 1);
  }

  TEST_CASE("rejects a region referencing an undefined material") {
    const auto path = writeTempYaml(R"(
spatial_mesh: [0.0, 1.0]
regions:
  materials: [copper]
energy_mesh: [0.1, 1.0]
materials:
  water:
    sigma_t: [1.0]
    sigma_s: [1.0]
    stopping_power:
      group_average: [1.0]
      group_boundary: [1.0, 1.0]
convergence:
  max_iters: 1
  epsilon: 1.0
)");
    InputDeck deck;
    CHECK(deck.read(path) == 1);
  }

  TEST_CASE("rejects cross sections with the wrong number of groups") {
    const auto path = writeTempYaml(R"(
spatial_mesh: [0.0, 1.0]
regions:
  materials: [water]
energy_mesh: [0.1, 0.5, 1.0]
materials:
  water:
    sigma_t: [1.0]
    sigma_s: [1.0]
    stopping_power:
      group_average: [1.0]
      group_boundary: [1.0, 1.0]
convergence:
  max_iters: 1
  epsilon: 1.0
)");
    InputDeck deck;
    CHECK(deck.read(path) == 1);
  }

  TEST_CASE("rejects an input file missing a required key") {
    const auto path = writeTempYaml(R"(
spatial_mesh: [0.0, 1.0]
regions:
  materials: [water]
materials:
  water:
    sigma_t: [1.0]
    sigma_s: [1.0]
    stopping_power:
      group_average: [1.0]
      group_boundary: [1.0, 1.0]
convergence:
  max_iters: 1
  epsilon: 1.0
)");
    InputDeck deck;
    CHECK(deck.read(path) == 1);
  }
}
