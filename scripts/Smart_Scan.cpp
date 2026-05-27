// Script used to more robustly solve for higher undercoolings by using the solutions to previous undercoolings as a
// starting guess for V and R.

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

    const alloys::Alloy A{alloys::AgCu};
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
