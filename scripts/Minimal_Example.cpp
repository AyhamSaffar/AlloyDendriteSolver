// Script used to demonstrate simplest possible workflow

#include <iostream>
#include <tuple>
#include "alloys.h"
#include "approximators.h"
#include "differentials.h"
#include "models.h"
#include "optimiser.h"


int main()
{
    // initialise variables
    const alloys::Alloy A{alloys::SnAg_wtp}; // common solder material, with concentrations in wt.%
    double f1{}, f2{}, dV{}, dR{}, dT{10.0}, C0{5.0};
    double V{approx::getV(dT, C0, A)}, R{approx::getR(dT, C0, A)};
    diff::Jacobian J{};

    // iteratively solve for V and R
    for (int step{0}; step<100; ++step)
    {
        std::tie(f1, f2) = models::LGK(V, R, dT, C0, A);
        J = diff::calculateGrads<models::LGK>(V, R, dT, C0, A);
        std::tie(dV, dR) = optimisers::newtonRaphson(f1, f2, J);
        V += 0.1 * dV; // smaller steps increase the range of starting V and R that don't diverge
        R += 0.1 * dR;
    }

    // print result
    std::cout << "R = " << R << " m, V = " << V << " m/s\n";
    return 0;
}
