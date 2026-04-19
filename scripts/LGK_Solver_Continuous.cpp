// Script used to reproduce figure 3 in https://doi.org/10.1016/0025-5416(84)90199-X.

#include <fstream>
#include <tuple>
#include <string>
#include <cmath>
#include "alloy.h"
#include "approximators.h"
#include "differentials.h"
#include "models.h"
#include "optimiser.h"

struct Result{bool hasDiverged{}; bool hasConverged{}; int steps{}; double f1{}; double f2{}; double V{}; double R{};};

Result VRSolver(double dT, double C0, const alloy::Alloy& A, double V0, double R0)
{
    double f1{}, f2{}, V{V0}, R{R0}, dV{}, dR{};
    diff::Jacobian J{};
    const int maxSteps{10'000};

    for (int step{0}; step<maxSteps; ++step)
    {
        std::tie(f1, f2) = models::LGK(V, R, dT, C0, A);
        J = diff::calculateGrads<models::LGK>(V, R, dT, C0, A);
        std::tie(dV, dR) = optimisers::newtonRaphson(f1, f2, J);
        if (std::isnan(dV) || std::isnan(dR)) // solver diverges
            return Result{true, false, step};
        V += 0.1 * dV;
        R += 0.1 * dR;
        if ((std::abs(f1)<1e-8) && (std::abs(f2)<1e-12)) // solver converged. f2 criteria lower as R is small
            return Result{false, true, step, f1, f2, V, R};
    }
    
    return Result{false, false, maxSteps, f1, f2, V, R};
}

int main()
{
    std::string dataPath{DATA_PATH};
    std::ofstream outf{dataPath + "/data.csv"};
    outf << "diverged,converged,steps,dT,C0,V,R,f1,f2\n" << std::boolalpha;

    const alloy::Alloy A{alloy::SucAce};
    for (double dT{5}; dT<10.0; dT+=4)
    {
        bool converged{false};
        double VBest{}, RBest{};
        for(double C0Molar{0.01}; C0Molar<=1.0; C0Molar+=0.001)
        {  
            double C0wt{C0Molar*7.252e-3};
            if (!converged)
            {
                VBest = approx::getTipVelocity(dT, C0wt, A);
                RBest = approx::getTipRadius(dT, C0wt, A);
            }
            Result result{VRSolver(dT, C0wt, A, VBest, RBest)};
            if (result.hasConverged)
            {
                std::tie(VBest, RBest) = std::tie(result.V, result.R);
                converged = true;
            }
            else
                converged = false;

            outf << result.hasDiverged << ',' << result.hasConverged << ',' << result.steps << ',' <<  dT << ',' <<
                C0wt << ',' << result.V << ',' << result.R << ',' << result.f1 << ',' << result.f2 << '\n';
        }
    }
    return 0;
}
