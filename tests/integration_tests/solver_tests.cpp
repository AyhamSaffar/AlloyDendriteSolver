#include <numbers>
#include <cmath>
#include <tuple>
#include <string>
#include <catch2/catch_test_macros.hpp>
#include "alloys.h"
#include "models.h"
#include "solvers.h"


// numericallly derived fit from https://doi.org/10.1007/s10854-025-14979-6
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
            solvers::Result result{solvers::solve<models::LGK>(dT, C0, alloys::SnAg)};
            REQUIRE(std::abs(result.f1) < 1e-12);
            REQUIRE(std::abs(result.f2) < 1e-12);
            REQUIRE(result.R > 0);

            double VFit{getPublishedLGKSnAgVFit(dT, C0)};
            REQUIRE(std::abs(result.V-VFit)/VFit < 0.20); // maximum of 20% error as numerical fit
        }
}

TEST_CASE(
    "LKT-BCT model V prediction agrees with published LGK SnAg numerical fit at low undercooling and gives positive R",
    "[solvers]"
)
{
    for (double dT{2.5}; dT<=20; dT+=2.5)
        for (double C0{3.0}; C0<=6.0; C0+=1.0)
        {
            INFO("dT = " + std::to_string(dT) + ", and C0 = " + std::to_string(C0));
            solvers::Result result{solvers::solve<models::LKT_BCT>(dT, C0, alloys::SnAg)};
            REQUIRE(std::abs(result.f1) < 1e-12);
            REQUIRE(std::abs(result.f2) < 1e-12);
            REQUIRE(result.R > 0);

            double VFit{getPublishedLGKSnAgVFit(dT, C0)};
            REQUIRE(std::abs(result.V-VFit)/VFit < 0.2); // maximum of 20% error as numerical fit
        }
}
