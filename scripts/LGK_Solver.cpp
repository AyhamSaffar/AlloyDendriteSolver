// script used to solve V and R for range of dT and C0

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

Result VRSolver(double dT, double C0, const alloy::Alloy& A)
{
    double f1{}, f2{}, V{approx::getTipVelocity(dT, C0, A)}, R{approx::getTipRadius(dT, C0, A)}, dV{}, dR{};
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

    alloy::AlloyTDependant A{alloy::SnAgTDependant};
    for (double dT{1.0}; dT<=60; dT+=0.5)
    {
        for(double C0{1.0}; C0<=5.0; C0+=2.0)
        {
            A.updateDiffusivity(dT, C0);
            Result result{VRSolver(dT, C0, A)};
            outf << result.hasDiverged << ',' << result.hasConverged << ',' << result.steps << ',' <<  dT << ',' <<
                C0 << ',' << result.V << ',' << result.R << ',' << result.f1 << ',' << result.f2 << '\n';
        }
    }
    return 0;
}
