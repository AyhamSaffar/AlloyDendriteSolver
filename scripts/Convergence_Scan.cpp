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
    
    const alloys::Alloy A{alloys::NiB_atp};
    const double C0{1};
    for (double dT{100}; dT<=300; dT+=100)
    {
        outfAprrox << dT << ',' << C0 << ',' << approx::getV(dT, C0, A) << ',' <<
            approx::getR(dT, C0, A) << '\n';

        for (double V0Power{-3}; V0Power<=3.0; V0Power+=0.1)
        {
            double V0{std::pow(10.0, V0Power)};
            for (double R0Power{-8}; R0Power<=-5; R0Power+=0.1)
            {
                double R0{std::pow(10.0, R0Power)};
                outfSolver << solvers::newton<models::LKT_BCT>(dT, C0, A, V0, R0).commaSeparatedValues() << ',';
                outfSolver << V0 << ',' << R0 << '\n';
            } 
        }

    }

    return 0;
}
