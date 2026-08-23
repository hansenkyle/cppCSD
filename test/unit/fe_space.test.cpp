#include "fe_space.h"

#include <stdexcept>

#include <doctest.h>

TEST_SUITE("FESpace") {
  TEST_CASE("stores the spatial and energy degrees it's constructed with") {
    const FESpace space(1, 2);

    CHECK(space.spatial_degree == 1);
    CHECK(space.energy_degree == 2);
  }

  TEST_CASE("rejects a negative spatial degree") {
    CHECK_THROWS_AS(FESpace(-1, 0), std::invalid_argument);
  }

  TEST_CASE("rejects a negative energy degree") {
    CHECK_THROWS_AS(FESpace(0, -1), std::invalid_argument);
  }
}
