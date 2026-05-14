// Script used to reproduce figure 2 in page 22 of
// https://search.library.uq.edu.au/discovery/fulldisplay/alma991011497109703131/61UQ_INST:61UQ.

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
    const int maxSteps{1000};

    for (int step{0}; step<maxSteps; ++step)
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
    
    return Result{false, false, maxSteps, f1, f2, V, R};
}

int main()
{
    std::string dataPath{DATA_PATH};
    std::ofstream outf{dataPath + "/data.csv"};
    outf << "diverged,converged,steps,dT,C0,V,R,f1,f2\n" << std::boolalpha;

    const alloy::Alloy A{alloy::CuAg};
    double C0{15}, dT{1.0}, dTStep{1.0};
    double V{approx::getTipVelocity(dT, C0, A)}, R{approx::getTipRadius(dT, C0, A)};

    while( (dT<=320) && (dTStep>=std::pow(2, -100)) )
    {
        Result result{VRSolver(dT, C0, A, V, R)};
        if (result.hasConverged)
        {
            outf << result.hasDiverged << ',' << result.hasConverged << ',' << result.steps << ',' <<  dT << ',' <<
                C0 << ',' << result.V << ',' << result.R << ',' << result.f1 << ',' << result.f2 << '\n';
            std::tie(V, R) = std::tie(result.V, result.R);
            dTStep = 1.0;
            dT += dTStep;
        }
        else
        {
            dT -= dTStep;
            dTStep /= 2;
            dT += dTStep;
        }
    }
    return 0;
}
