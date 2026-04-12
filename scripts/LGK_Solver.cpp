// Minimal script used to reproduce figure 3 in https://doi.org/10.1016/0025-5416(84)90199-X.

#include <fstream>
#include <tuple>
#include <string>
#include <cmath>
#include "alloy.h"
#include "approximators.h"
#include "differentials.h"
#include "models.h"
#include "optimiser.h"

struct Result{bool hasDiverged{}; double f1{}; double f2{}; double V{}; double R{};};

Result VRSolver(double dT, double C0, const alloy::Alloy& A)
{
    double f1{}, f2{}, V{approx::getTipVelocity(dT, C0, A)}, R{approx::getTipRadius(dT, C0, A)}, dV{}, dR{};
    diff::Jacobian J{};

    for (int step{0}; step<1000; ++step)
    {
        std::tie(f1, f2) = models::LGK(V, R, dT, C0, A);
        J = diff::calculateGrads<models::LGK>(V, R, dT, C0, A);
        std::tie(dV, dR) = optimisers::newtonRaphson(f1, f2, J);
        if (std::isnan(dV) || std::isnan(dR)) // solver diverges
            return Result{true};
        V += 0.1*dV;
        R += 0.1*dR;
        if ((std::abs(f1)<1e-16) && (std::abs(f2)<1e-16)) // solver converged
            break;
    }
    
    return Result{false, f1, f2, V, R};
}

int main()
{
    std::string dataPath{DATA_PATH};
    std::ofstream outf{dataPath + "/data.csv"};
    outf << "diverged,dT,C0,V,R,f1,f2\n" << std::boolalpha;
    for (double dT{5}; dT<10.0; dT+=4)
    {
        for(double C0Molar{0.01}; C0Molar<=1.0; C0Molar+=0.01)
        {  
            double C0wt{C0Molar*7.252e-3};
            Result result{VRSolver(dT, C0wt, alloy::SucAce)};
            outf << result.hasDiverged << ',' << dT << ',' << C0wt << ',' << result.V << ',' << result.R << ',' <<
                result.f1 << ',' << result.f2 << '\n';
        }
    }
    return 0;
}
