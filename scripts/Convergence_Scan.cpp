// Script used to test whether a model will converge for a range of starting guesses for dendrite velocity.

#include <fstream>
#include <string>
#include <cmath>
#include "solver.h"
#include "alloys.h"
#include "models.h"
#include "approximators.h"


int main()
{
    std::string dataPath{DATA_PATH};
    std::ofstream outfSolver{dataPath + "/solver_data.csv"};
    outfSolver << solver::Result::commaSeparatedColumns << '\n';
    std::ofstream outfAprrox{dataPath + "/approx_data.csv"};
    outfAprrox << "dT,C0,V,R" << '\n';
    
    const alloys::Alloy A{alloys::CuAg};
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
