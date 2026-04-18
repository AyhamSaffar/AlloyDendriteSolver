// Script used to reproduce figure 11 in https://link.springer.com/article/10.1007/s10854-025-14979-6, adding
// bounds for how well the model converges given a starting guess for dendrite velocity.

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

Result VRSolver(double dT, double C0, const alloy::Alloy& A, double V0)
{
    double f1{}, f2{}, V{V0}, R{approx::getTipRadius(dT, C0, A)}, dV{}, dR{};
    diff::Jacobian J{};
    int nSteps{1000};

    for (int step{0}; step<nSteps; ++step)
    {
        std::tie(f1, f2) = models::LGK(V, R, dT, C0, A);
        J = diff::calculateGrads<models::LGK>(V, R, dT, C0, A);
        std::tie(dV, dR) = optimisers::newtonRaphson(f1, f2, J);
        if (std::isnan(dV) || std::isnan(dR)) // solver diverges
            return Result{true, false, step};
        V += 0.1*dV;
        R += 0.1*dR;
        if ((std::abs(f1)<1e-12) && (std::abs(f2)<1e-12)) // solver converged
            return Result{false, true, step, f1, f2, V, R};
    }
    
    return Result{false, false, nSteps, f1, f2, V, R};
}

int main()
{
    std::string dataPath{DATA_PATH};
    std::ofstream outfSolver{dataPath + "/solver_data.csv"};
    outfSolver << "diverged,converged,steps,dT,C0,V,R,f1,f2,V0\n" << std::boolalpha;
    std::ofstream outfAprrox{dataPath + "/approx_data.csv"};
    outfAprrox << "dT,C0,V,R\n";
    
    const double C0{3.5};
    alloy::Alloy A{alloy::SnAg};
    for (double dT{1.0}; dT<=50.0; dT+=0.5)
    {
        for(double V0Power{-5.0}; V0Power<=1.0; V0Power+=0.1)
        {
            double V0{std::pow(10.0, V0Power)};
            Result result{VRSolver(dT, C0, A, V0)};
            outfSolver << result.hasDiverged << ',' << result.hasConverged << ',' << result.steps << ',' << dT << ',' <<
                C0 << ',' << result.V << ',' << result.R << ',' << result.f1 << ',' << result.f2 << ',' << V0 << '\n';
        }

        outfAprrox << dT << ',' << C0 << ',' << approx::getTipVelocity(dT, C0, A) << ',' <<
            approx::getTipRadius(dT, C0, A) << '\n';
    }
    return 0;
}
