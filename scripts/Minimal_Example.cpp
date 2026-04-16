// Script used to demonstrate simplest possible workflow

#include <iostream>
#include <tuple>
#include "alloy.h"
#include "approximators.h"
#include "differentials.h"
#include "models.h"
#include "optimiser.h"


int main()
{
    // initialise variables
    const alloy::Alloy A{alloy::SnAg}; // common solder material
    double f1{}, f2{}, dV{}, dR{}, dT{10.0}, C0{5.0}; // only non SI unit is concentration (wt.%)
    double V{approx::getTipVelocity(dT, C0, A)}, R{approx::getTipRadius(dT, C0, A)};
    diff::Jacobian J{};

    // iteratively solve for V and R
    for (int step{0}; step<100; ++step)
    {
        std::tie(f1, f2) = models::LGK(V, R, dT, C0, A);
        J = diff::calculateGrads<models::LGK>(V, R, dT, C0, A);
        std::tie(dV, dR) = optimisers::newtonRaphson(f1, f2, J);
        V += 0.01 * dV; // smaller steps improve convergence
        R += 0.01 * dR;
    }

    // print result
    std::cout << "R = " << R << " m, V = " << V << " m/s\n";
    return 0;
}
