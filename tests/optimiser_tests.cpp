#include <tuple>
#include <catch2/catch_test_macros.hpp>
#include "optimiser.h"
#include "differentials.h"


TEST_CASE("Newton Raphson optimser can solve linear basis in one step", "[optimiser]")
{
    double f1{1.0}, f2{2.0}, df1dV{3.0}, df1dR{4.0}, df2dV{5.0}, df2dR{6.0}, df1{}, df2{};
    diff::Jacobian J{df1dV, df1dR, df2dV, df2dR};
    std::tie(df1, df2) = optimisers::newtonRaphson(f1, f2, J);
    REQUIRE(df1 == 1.0) // TODO finish
}
