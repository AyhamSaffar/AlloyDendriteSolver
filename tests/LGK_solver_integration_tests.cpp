#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <tuple>
#include <string>
#include "alloy.h"
#include "approximators.h"
#include "differentials.h"
#include "models.h"
#include "optimiser.h"


// returns iteratively calculated V and R
std::tuple<double, double> LGKSolver(double dT, double C0, const alloy::Alloy& A)
{
    double f1{}, f2{}, dV{}, dR{}, V{approx::getTipVelocity(dT, C0, A)}, R{approx::getTipRadius(dT, C0, A)};
    diff::Jacobian J{};

    for (int step{0}; step<100; ++step)
    {
        std::tie(f1, f2) = models::LGK(V, R, dT, C0, A);
        J = diff::calculateGrads<models::LGK>(V, R, dT, C0, A);
        std::tie(dV, dR) = optimisers::newtonRaphson(f1, f2, J);
        V += 0.1 * dV; // smaller steps improve convergence
        R += 0.1 * dR;
    }

    return std::tuple(V, R);
}

// numericallly derived fit from https://doi.org/10.1007/s10854-025-14979-6
double getPublishedSnAgVFit(double dT, double C0)
{
    return std::pow(dT*std::pow(C0, -0.59)/28.6, 1/0.35);
}

TEST_CASE("LGK model V prediction agrees with published numerical fit and gives positive R", "[solver]")
{
    double VPred{}, RPred{};
    for (double dT{10.0}; dT<=60; dT+=10.0)
    {
        for (double C0{3.0}; C0<=6.0; C0+=1.0)
        {
            INFO("dT = " + std::to_string(dT) + ", and C0 = " + std::to_string(C0));
            std::tie(VPred, RPred) = LGKSolver(dT, C0, alloy::SnAg);
            REQUIRE(RPred > 0);

            double VFit{getPublishedSnAgVFit(dT, C0)};
            REQUIRE(std::abs(VPred-VFit)/VFit < 0.20); // maximum of 20% error as simple solver and numerical fit
        }
    }
}
