// Script used to test whether a model will converge for a range of starting guesses for dendrite velocity.

#include <fstream>
#include <string>
#include <cmath>
#include <array>
#include "solvers.h"
#include "alloys.h"
#include "models.h"
#include "approximators.h"


int main()
{
    std::string dataPath{DATA_PATH};
    std::ofstream outfSolver{dataPath + "/solver_data.csv"};
    outfSolver << solvers::Result::commaSeparatedColumns << ",V0,R0\n";
    std::ofstream outfAprrox{dataPath + "/approx_data.csv"};
    outfAprrox << "dT,C0,V,R" << '\n';
    
    const alloys::Alloy A{alloys::AgCu};
    const double C0{15};
    for (double dT{50}; dT<=150; dT+=50)
    {
        outfAprrox << dT << ',' << C0 << ',' << approx::getTipVelocity(dT, C0, A) << ',' <<
            approx::getTipRadius(dT, C0, A) << '\n';

        for (double V0Power{-7.0}; V0Power<=4.0; V0Power+=0.1)
        {
            double V0{std::pow(10.0, V0Power)};
            for (double R0Power{-8}; R0Power<=-5.2; R0Power+=0.1)
            {
                double R0{std::pow(10.0, R0Power)};
                outfSolver << solvers::solve<models::LKT_BCT>(dT, C0, A, V0, R0).commaSeparatedValues() << ',';
                outfSolver << V0 << ',' << R0 << '\n';
            } 
        }

    }

    return 0;
}
