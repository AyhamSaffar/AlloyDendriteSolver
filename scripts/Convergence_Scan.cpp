// Script used to reproduce figure 11 in https://link.springer.com/article/10.1007/s10854-025-14979-6, adding
// bounds for how well the model converges given a starting guess for dendrite velocity.

#include <fstream>
#include <string>
#include <cmath>
#include "solver.h"
#include "alloy.h"
#include "models.h"
#include "approximators.h"


int main()
{
    std::string dataPath{DATA_PATH};
    std::ofstream outfSolver{dataPath + "/solver_data.csv"};
    outfSolver << solver::Result::commaSeparatedColumns << '\n';
    std::ofstream outfAprrox{dataPath + "/approx_data.csv"};
    outfAprrox << "dT,C0,V,R" << '\n';
    
    const alloy::Alloy A{alloy::CuAg};
    const double C0{15};
    for (double dT{1.0}; dT<=320; dT+=1.0)
    {
        for(double V0Power{-7.0}; V0Power<=4.0; V0Power+=0.1)
        {
            double V0{std::pow(10.0, V0Power)};
            solver::Result result{solver::solve<models::LGK>(dT, C0, A, V0, approx::getTipRadius(dT, C0, A))};
            outfSolver << result.commaSeparatedValues() << '\n';
        }

        outfAprrox << dT << ',' << C0 << ',' << approx::getTipVelocity(dT, C0, A) << ',' <<
            approx::getTipRadius(dT, C0, A) << '\n';
    }
    return 0;
}
