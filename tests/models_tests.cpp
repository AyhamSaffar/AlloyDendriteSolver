#include <catch2/catch_test_macros.hpp>

TEST_CASE("2 plus 2 equals 4", "[basics]")
{
    REQUIRE(2+2==4);
    REQUIRE(1+3==5);
}

