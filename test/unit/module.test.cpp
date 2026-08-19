
#include "module.h"
#include <doctest.h> // or "doctest/doctest.h", depending on which include path you chose

TEST_SUITE("example") {
  TEST_CASE("example") { CHECK(1 + 1 == 2); }
}
