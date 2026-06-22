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

    double V{approx::getV(dT, C0, A)}, R{approx::getR(dT, C0, A)};
    std::tie(f1, f2) = models::LGK(V, R, dT, C0, A);
    
    // near zero values for f1 and f2 means the models predict the values V and R are correct for the dT & C0 used
    REQUIRE(std::abs(f1) < 0.1);
    REQUIRE(std::abs(f2) < 0.1);
}

TEST_CASE("LKT_BCT model roughly agrees with approximation at low undercooling", "[models]")
{
    double dT{0.001}, C0{5.0}, f1{}, f2{}; 
    alloys::Alloy A{alloys::SnAg_wtp};

    double V{approx::getV(dT, C0, A)}, R{approx::getR(dT, C0, A)};
    std::tie(f1, f2) = models::LKT_BCT(V, R, dT, C0, A);
    
    // near zero values for f1 and f2 means the models predict the values V and R are correct for the dT & C0 used
    REQUIRE(std::abs(f1) < 0.1);
    REQUIRE(std::abs(f2) < 0.1);
}
