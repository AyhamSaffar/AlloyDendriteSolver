#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <tuple>
#include "approximators.h"
#include "models.h"
#include "alloys.h"

TEST_CASE("LGK model roughly agrees with approximation at low undercooling", "[models]")
{
    double dT{0.001}, C0{5.0}, f1{}, f2{}; 
    alloys::Alloy A{alloys::SnAg_wtp};
    models::DTs _{};

    double V{approx::getV(dT, C0, A)}, R{approx::getR(dT, C0, A)};
    std::tie(f1, f2, _) = models::LGK(V, R, dT, C0, A);
    
    // near zero values for f1 and f2 means the models predict the values V and R are correct for the dT & C0 used
    REQUIRE(std::abs(f1) < 0.1);
    REQUIRE(std::abs(f2) < 0.1);
}

TEST_CASE("LKT_BCT model roughly agrees with approximation at low undercooling", "[models]")
{
    double dT{0.001}, C0{5.0}, f1{}, f2{}; 
    alloys::Alloy A{alloys::SnAg_wtp};
    models::DTs _{};

    double V{approx::getV(dT, C0, A)}, R{approx::getR(dT, C0, A)};
    std::tie(f1, f2, _) = models::LKT_BCT(V, R, dT, C0, A);
    
    // near zero values for f1 and f2 means the models predict the values V and R are correct for the dT & C0 used
    REQUIRE(std::abs(f1) < 0.1);
    REQUIRE(std::abs(f2) < 0.1);
}

TEST_CASE("models do not modify Alloy objects passed to them", "[differentials]")
{
    // checks enzyme bug (https://github.com/EnzymeAD/Enzyme/issues/3073) is prevented from modifying Alloy objects
    
    // must create an Alloy that can be used with all models
    std::vector<alloys::Fit> fits{ alloys::Fit{{1, 2, 3}}, alloys::Fit{{4, 5, 6}}};
    const alloys::Alloy A{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, fits, fits, fits, 1};
    const alloys::Alloy ACopy{A};

    double f1{}, f2{};
    models::DTs dTs{};

    std::tie(f1, f2, dTs) = models::LGK(1, 1, 1, 1, A);
    REQUIRE(A==ACopy);

    std::tie(f1, f2, dTs) = models::LKT_BCT(1, 1, 1, 1, A);
    REQUIRE(A==ACopy);

    std::tie(f1, f2, dTs) = models::CLW(1, 1, 1, 1, A);
    REQUIRE(A==ACopy);

    std::tie(f1, f2, dTs) = models::GD(1, 1, 1, 1, A);
    REQUIRE(A==ACopy);

    std::tie(f1, f2, dTs) = models::WLCYZ(1, 1, 1, 1, A);
    REQUIRE(A==ACopy);
}
