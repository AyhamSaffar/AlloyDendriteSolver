// Script used to reproduce figure 2 in page 22 of
// https://search.library.uq.edu.au/discovery/fulldisplay/alma991011497109703131/61UQ_INST:61UQ.

#include <fstream>
#include <tuple>
#include <string>
#include <cmath>
#include "alloys.h"
#include "models.h"
#include "solver.h"


int main()
{
    std::string dataPath{DATA_PATH};
    std::ofstream outf{dataPath + "/data.csv"};
    outf << solver::Result::commaSeparatedColumns << '\n';

    const alloys::Alloy A{alloys::CuAg};
    double C0{15}, dT{1.0}, dTStep{1.0};
    double V{approx::getTipVelocity(dT, C0, A)}, R{approx::getTipRadius(dT, C0, A)};

    while( (dT<=320) && (dTStep>=std::pow(2, -100)) )
    {
        solver::Result result{solver::solve<models::LGK>(dT, C0, A, V, R)};
        if (result.hasConverged)
        {
            outf << result.commaSeparatedValues() << '\n';
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
