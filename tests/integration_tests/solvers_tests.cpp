#include <numbers>
#include <cmath>
#include <tuple>
#include <string>
#include <catch2/catch_test_macros.hpp>
#include "alloys.h"
#include "models.h"
#include "solvers.h"


// numericallly derived LGK fit from https://doi.org/10.1007/s10854-025-14979-6
double getPublishedLGKSnAgVFit(double dT, double C0)
{
    return std::pow(dT*std::pow(C0, -0.59)/28.6, 1/0.35);
}

TEST_CASE("LGK model V prediction agrees with published LGK SnAg numerical fit and gives positive R", "[solvers]")
{
    for (double dT{2.5}; dT<=60; dT+=2.5)
        for (double C0{3.0}; C0<=6.0; C0+=1.0)
        {
            INFO("dT = " + std::to_string(dT) + ", and C0 = " + std::to_string(C0));
            solvers::Result result{solvers::newton<models::LGK>(dT, C0, alloys::SnAg_wtp)};
            REQUIRE(result.hasConverged);
            REQUIRE(result.R > 0);

            double VFit{getPublishedLGKSnAgVFit(dT, C0)};
            REQUIRE(std::abs(result.V-VFit)/VFit < 0.20); // maximum of 20% error as numerical fit
        }
}


TEST_CASE("LKT-BCT model V prediction agrees with LGK at low undercooling and gives positive R", "[solvers]")
{
    // low undercoolings = lower V and R = k does not vary from k0, stability functions equal 1, and negligible kinetic
    // undercooling. Therefore LKT_BCT reduces to LGK 
    for (double dT{2.5}; dT<=20; dT+=2.5)
        for (double C0{3.0}; C0<=6.0; C0+=1.0)
        {
            INFO("dT = " + std::to_string(dT) + ", and C0 = " + std::to_string(C0));

            constexpr bool legacy{false}; // ensures LGK form is consistent with LKT_BCT
            solvers::Result LGKResult{solvers::newton<models::LGK<legacy>>(dT, C0, alloys::SnAg_wtp)};
            solvers::Result LKT_BCTResult{solvers::newton<models::LKT_BCT>(dT, C0, alloys::SnAg_wtp)};

            REQUIRE(LKT_BCTResult.hasConverged);
            REQUIRE(LKT_BCTResult.R > 0);
            REQUIRE(std::abs(LKT_BCTResult.V-LGKResult.V)/LGKResult.V < 0.05); // maximum of 5% error
        }
}

// Currently fails. Does slightly better if in using an atp alloy.
TEST_CASE("Linearised Dynamic model agrees with LKT-BCT at pre solute trapping undercoolings", "[Solvers]")
{
    const alloys::Alloy A{alloys::AgCu_wtp}; // LGK-BCT capable Alloy
    const alloys::Alloy ALin{ // dynamic model with fixed D, k0 and linear Tl
        A.L, A.Cp, A.m, A.k0, A.r, A.D, A.a, A.o, A.a0, A.V0, A.Tm, A.D, 0, {A.Tm, A.m}, {A.k0}
    };
    const double C0{15}, dT0{1}, Vd{A.D/A.a0}; // solute trapping starts as V approaches max diffusive Velocity Vd 
    double V0{approx::getV(dT0, C0, A)}, R0{approx::getR(dT0, C0, A)};
    
    for (double dT{1}; dT<500; ++dT)
    {
        INFO("dT = " + std::to_string(dT));
        solvers::Result LKT_BCTResult{solvers::newton<models::LKT_BCT>(dT, C0, A, V0, R0)};
        solvers::Result dynamicResult{solvers::newton<models::dynamic>(dT, C0, ALin, V0, R0)};
        INFO("V LKT-BCT = " + std::to_string(LKT_BCTResult.V) + ", V dynamic = " + std::to_string(dynamicResult.V) + '\n');
        INFO("R LKT-BCT = " + std::to_string(LKT_BCTResult.R) + ", R dynamic = " + std::to_string(dynamicResult.R) + '\n');

        REQUIRE(LKT_BCTResult.hasConverged);
        REQUIRE(dynamicResult.hasConverged);
        REQUIRE((std::abs(LKT_BCTResult.R - dynamicResult.R)/LKT_BCTResult.R) < 2); // tolerance for rounding errors
        REQUIRE((std::abs(LKT_BCTResult.V - dynamicResult.V)/LKT_BCTResult.V) < 2);

        if (dynamicResult.V > (Vd/10))
            break;  // dynamic model not neccessarily equivalent to LKT_BCT when solute trapping starts
        std::tie(V0, R0) = std::tie(dynamicResult.V, dynamicResult.R);
    }

}
