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

TEST_SUITE("Field") {
  TEST_CASE("stores the dimensions it's constructed with") {
    const Field field(3, 2);

    CHECK(field.numCells() == 3);
    CHECK(field.numGroups() == 2);
  }

  TEST_CASE("rejects a non-positive n_x") {
    CHECK_THROWS_AS(Field(0, 2), std::invalid_argument);
    CHECK_THROWS_AS(Field(-1, 2), std::invalid_argument);
  }

  TEST_CASE("rejects a non-positive G") {
    CHECK_THROWS_AS(Field(2, 0), std::invalid_argument);
    CHECK_THROWS_AS(Field(2, -1), std::invalid_argument);
  }

  TEST_CASE("zero-initializes every corner of every (group, cell)") {
    Field field(2, 2);

    for (int g = 0; g < field.numGroups(); ++g) {
      for (int i = 0; i < field.numCells(); ++i) {
        CHECK(field[g][i].leftDown() == 0.0);
        CHECK(field[g][i].leftUp() == 0.0);
        CHECK(field[g][i].rightDown() == 0.0);
        CHECK(field[g][i].rightUp() == 0.0);
      }
    }
  }

  TEST_CASE("named accessors read back what was written through them") {
    Field field(2, 1);

    field[0][1].leftDown() = 1.0;
    field[0][1].leftUp() = 2.0;
    field[0][1].rightDown() = 3.0;
    field[0][1].rightUp() = 4.0;

    CHECK(field[0][1].leftDown() == 1.0);
    CHECK(field[0][1].leftUp() == 2.0);
    CHECK(field[0][1].rightDown() == 3.0);
    CHECK(field[0][1].rightUp() == 4.0);
  }

  TEST_CASE("writes to one (group, cell) don't leak into another") {
    Field field(2, 2);

    field[0][0].leftDown() = 5.0;

    CHECK(field[0][1].leftDown() == 0.0);
    CHECK(field[1][0].leftDown() == 0.0);
    CHECK(field[1][1].leftDown() == 0.0);
  }

  TEST_CASE("CornerValues is a view: writes land directly in the flat storage") {
    Field field(2, 2);
    const int n_x = field.numCells();

    field[1][0].rightUp() = 42.0;

    // group 1, cell 0, corner 3 (rightUp) -> (group * n_x + cell) * 4 + corner
    const int index = (1 * n_x + 0) * 4 + 3;
    CHECK(field.values()(index) == 42.0);

    field.values()(index) = 7.0;
    CHECK(field[1][0].rightUp() == 7.0);
  }

  TEST_CASE("underlying flat vector is sized 4 * n_x * G") {
    const Field field(3, 5);

    CHECK(field.values().size() == 4 * 3 * 5);
  }

  TEST_CASE("const Field exposes its flat storage read-only") {
    Field field(1, 1);
    field[0][0].leftUp() = 9.0;

    const Field& const_field = field;
    CHECK(const_field.values()(1) == 9.0);
  }

  TEST_CASE("rejects an out-of-range group index") {
    Field field(2, 2);

    CHECK_THROWS_AS(field[-1], std::out_of_range);
    CHECK_THROWS_AS(field[2], std::out_of_range);
  }

  TEST_CASE("rejects an out-of-range cell index") {
    Field field(2, 2);

    CHECK_THROWS_AS(field[0][-1], std::out_of_range);
    CHECK_THROWS_AS(field[0][2], std::out_of_range);
  }
}
