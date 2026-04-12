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
    const alloy::Alloy A{alloy::SnAg};
    double f1{}, f2{}, dV{}, dR{}, dT{0.1}, C0{1.0};
    double V{approx::getTipVelocity(dT, C0, A)}, R{approx::getTipRadius(dT, C0, A)};
    diff::Jacobian J{};

    for (int step{0}; step<100; ++step)
    {
        std::tie(f1, f2) = models::LGK(V, R, dT, C0, A);
        J = diff::calculateGrads<models::LGK>(V, R, dT, C0, A);
        std::tie(dV, dR) = optimisers::newtonRaphson(f1, f2, J);
        V += dV;
        R += dR;
    }

    std::cout << "R = " << R << ", V = " << V << '\n';
    return 0;
}
