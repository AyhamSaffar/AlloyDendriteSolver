#include <iostream>
#include <tuple>
#include "alloy.h"
#include "approximators.h"
#include "differentials.h"
#include "models.h"
#include "optimiser.h"


int main()
{
    double dT{1.0}, C0{5.0}, f1{}, f2{}, dV{}, dR{};
    alloy::Alloy A{alloy::SnAg};
    double V{approx::getTipVelocity(A.D, A.m, A.k0, A.r, dT, C0)};
    double R{approx::getTipRadius(A.r, A.m, A.k0, C0, dT)};
    diff::Jacobian J{};

    for (int i{0}; i<10; ++i)
    {
        std::tie(f1, f2) = models::LGK(V, R, dT, C0, A);
        J = diff::calculateGrads<models::LGK>(V, R, dT, C0, A);
        std::tie(dV, dR) = optimisers::newtonRaphson(f1, f2, J);
        std::cout << "V: " << V << ", dV: " << dV << ", R: " << R << ", dR: " << dR << ", f1: " << f1 << ", f2: " 
            << f2 << '\n';
        V += dV;
        R += dR;
    }

    return 0;
}
