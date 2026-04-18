#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <tuple>
#include "approximators.h"
#include "models.h"
#include "alloy.h"

TEST_CASE("LGK model agrees with approximation at low undercooling", "[models]")
{
    double dT{0.001}, C0{5.0}, f1{}, f2{}; 
    alloy::Alloy A{alloy::SnAg};

    double V{approx::getTipVelocity(dT, C0, A)}, R{approx::getTipRadius(dT, C0, A)};
    std::tie(f1, f2) = models::LGK(V, R, dT, C0, A);
    
    // near zero values for f1 and f2 means the models predict the values V and R are correct for the dT & C0 used
    REQUIRE(std::abs(f1) < 0.01);
    REQUIRE(std::abs(f2) < 0.01);
}
