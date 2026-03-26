#include <catch2/catch_test_macros.hpp>
#include <tuple>
#include <cmath> // for std::exp
#include "differentials.h"
#include "models.h"
#include "alloy.h"

// easily differentiable function that follows models::ModelFunc prototype
std::tuple<double, double> dummyModel(double V, double R, double dT, double C0, const alloy::Alloy& A)
{
    double f1{V*V*std::exp(R)};
    double f2{std::log(V) + 3*R};
    return std::make_tuple(f1, f2);
}

TEST_CASE("diff::calculateGrads works for analytically differentiable expressions", "[differentials]")
{
    double V{3.0}, R{5.0};
    diff::Jacobian J{diff::calculateGrads<dummyModel>(V, R, 1.0, 1.0, alloy::SnAg)};
    REQUIRE(J.df1dV == 2*V*std::exp(R));
    REQUIRE(J.df1dR == V*V*std::exp(R));
    REQUIRE(J.df2dV == 1/V);
    REQUIRE(J.df2dR == 3);
}

TEST_CASE("diff::calculateGrads works for successive analytically differentiable expressions", "[differentials]")
{
    for (double shift{0}; shift < 2; shift += 0.4)
    {
        double V{3.0-shift}, R{5.0+shift};
        diff::Jacobian J{diff::calculateGrads<dummyModel>(V, R, 1.0, 1.0, alloy::SnAg)};
        REQUIRE(J.df1dR == V*V*std::exp(R));
        REQUIRE(J.df1dV == 2*V*std::exp(R));
        REQUIRE(J.df2dR == 3);
        REQUIRE(J.df2dV == 1/V);
    }
}
