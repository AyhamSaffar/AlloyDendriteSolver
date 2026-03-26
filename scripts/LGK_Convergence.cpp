// Minimal script used track how well the optimiser converges when attempting to reproduce fig.3 in
// https://doi.org/10.1016/0025-5416(84)90199-X.

#include <vector>
#include <fstream>
#include <tuple>
#include <string>
#include <cmath>
#include "alloy.h"
#include "approximators.h"
#include "differentials.h"
#include "models.h"
#include "optimiser.h"

struct Result{std::vector<double> f1s; std::vector<double> f2s;};

Result VRSolver(double dT, double C0, const alloy::Alloy& A)
{
    const int nSteps{1000};
    double f1{}, f2{}, V{approx::getTipVelocity(dT, C0, A)}, R{approx::getTipRadius(dT, C0, A)}, dV{}, dR{};
    std::vector<double> f1s(nSteps), f2s(nSteps);
    diff::Jacobian J{};

    for (int step{0}; step<nSteps; ++step)
    {
        std::tie(f1, f2) = models::LGK(V, R, dT, C0, A);
        f1s[step] = f1;
        f2s[step] = f2;
        J = diff::calculateGrads<models::LGK>(V, R, dT, C0, A);
        std::tie(dV, dR) = optimisers::newtonRaphson(f1, f2, J);
        if (std::isnan(dV) || std::isnan(dR)) // solver diverges
            break;
        V += dV*0.1;
        R += dR*0.1;
        if ((std::abs(f1)<1e-8) && (std::abs(f2)<1e-8)) // solver converged
            break;
    }
    
    return Result{f1s, f2s};
}

int main()
{
    std::string dataPath{DATA_PATH};
    std::ofstream outf{dataPath + "/data.csv"};
    outf << "dT,C0,f1s,f2s\n" << std::boolalpha;
    for (double dT{0.5}; dT<1.0; dT+=0.4)
    {
        for(double C0{0.1}; C0<=1.0; C0+=0.1)
        {
            Result result{VRSolver(dT, C0, alloy::SucAce)};
            outf <<  dT << ',' << C0 << ',';
            for (double f1: result.f1s)
                outf << f1 << ' ';
            outf << ',';
            for (double f2: result.f2s)
                outf << f2 << ' ';
            outf << '\n';
        }
    }
    return 0;
}
