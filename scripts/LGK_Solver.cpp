// Minimal script used to reproduce figure 3 in https://doi.org/10.1016/0025-5416(84)90199-X.

#include <optional>
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
    double f1{}, f2{}, dV{}, dR{};
    double V{approx::getTipVelocity(dT, C0, A)};
    double R{approx::getTipRadius(dT, C0, A)};
    diff::Jacobian J{};

    for (int i{0}; i<1000; ++i)
    {
        std::tie(f1, f2) = models::LGK(V, R, dT, C0, A);
        J = diff::calculateGrads<models::LGK>(V, R, dT, C0, A);
        std::tie(dV, dR) = optimisers::newtonRaphson(f1, f2, J);
        if (std::isnan(dV) || std::isnan(dR)) // solver diverges
            return Result{true};
        V += dV * 0.1;
        R += dR * 0.1;
        if ((f1<1e-8) && (f2<1e-8)) // solver converged
            break;
    }
    
    return Result{false, f1, f2, V, R};
}

int main()
{
    std::string dataPath{DATA_PATH};
    std::ofstream outf{dataPath + "/data.csv"};
    outf << "diverged,dT,C0,V,R,f1,f2\n" << std::boolalpha;
    for (double dT{0.5}; dT<1.0; dT+=0.4)
    {
        for(double C0{0.0}; C0<=1.0; C0+=0.01)
        {
            Result result{VRSolver(dT, C0, alloy::SucAce)};
            outf << result.hasDiverged << ',' << dT << ',' << C0 << ',' << result.V << ',' << result.R << ',' <<
                result.f1 << ',' << result.f2 << '\n';
        }
    }
    return 0;
}
