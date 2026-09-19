#include <catch2/catch_test_macros.hpp>
#include <tuple>
#include <cmath> // for std::exp
#include "differentials.h"
#include "models.h"
#include "alloys.h"

// easily differentiable function that follows models::ModelFunc prototype
std::tuple<double, double, models::DTs> dummyModel(double V, double R, double dT, double C0, const alloys::Alloy& A)
{
    double f1{V*V*std::exp(R)};
    double f2{std::log(V) + 3*R};
    models::DTs _{1.0, 1.0, 1.0, 1.0}; // must explicitly pass floats so Enzyme knows what data type to expect in DTs 
    return std::make_tuple(f1, f2, _);
}

TEST_CASE("diff::calculateGrads works for analytically differentiable expressions", "[differentials]")
{
    double V{3}, R{5}, dT{1}, C0{1};
    diff::Jacobian J{diff::calculateGrads<dummyModel>(V, R, dT, C0, alloys::SnAg_wtp)};
    REQUIRE(J.df1dV == 2*V*std::exp(R));
    REQUIRE(J.df1dR == V*V*std::exp(R));
    REQUIRE(J.df2dV == 1/V);
    REQUIRE(J.df2dR == 3);
}

TEST_CASE("diff::calculateGrads works for successive analytically differentiable expressions", "[differentials]")
{
    for (double shift{0}; shift < 2; shift += 0.4)
    {
        double V{3.0-shift}, R{5.0+shift}, dT{1}, C0{1};
        diff::Jacobian J{diff::calculateGrads<dummyModel>(V, R, dT, C0, alloys::SnAg_wtp)};
        REQUIRE(J.df1dR == V*V*std::exp(R));
        REQUIRE(J.df1dV == 2*V*std::exp(R));
        REQUIRE(J.df2dR == 3);
        REQUIRE(J.df2dV == 1/V);
    }
}

TEST_CASE("diff::calculateGrads does not modify Alloy objects passed to it", "[differentials]")
{
    // checks enzyme bug (https://github.com/EnzymeAD/Enzyme/issues/3073) is prevented from modifying Alloy objects
    
    // must create an Alloy that can be used with all models
    std::vector<alloys::Fit> fits{ alloys::Fit{{1, 2, 3}}, alloys::Fit{{4, 5, 6}}};
    const alloys::Alloy A{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, fits, fits, fits, 1};
    const alloys::Alloy ACopy{A};
    diff::Jacobian J{};

    J = diff::calculateGrads<models::LGK>(1, 1, 1, 1, A);
    REQUIRE(A==ACopy);

    J = diff::calculateGrads<models::LKT_BCT>(1, 1, 1, 1, A);
    REQUIRE(A==ACopy);

    J = diff::calculateGrads<models::CLW>(1, 1, 1, 1, A);
    REQUIRE(A==ACopy);

    J = diff::calculateGrads<models::GD>(1, 1, 1, 1, A);
    REQUIRE(A==ACopy);

    J = diff::calculateGrads<models::WLCYZ>(1, 1, 1, 1, A);
    REQUIRE(A==ACopy);
}
