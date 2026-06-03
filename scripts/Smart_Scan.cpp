// Script used to more robustly solve for higher undercoolings by using the solutions to previous undercoolings as a
// starting guess for V and R.

#include <fstream>
#include <tuple>
#include <string>
#include <cmath>
#include "alloys.h"
#include "models.h"
#include "solvers.h"


int main()
{
    std::string dataPath{DATA_PATH};
    std::ofstream outf{dataPath + "/data.csv"};
    outf << solvers::Result::commaSeparatedColumns << '\n';

    const alloys::Alloy A{alloys::AgCu};
    double C0{15}, dT{1.0}, dTStep{1.0};
    double V0{approx::getTipVelocity(dT, C0, A)}, R0{approx::getTipRadius(dT, C0, A)};

    for (double C0{30}; C0<=50; C0+=10)
        while( (dT<=320) && (dTStep>=std::pow(2, -100)) )
        {
            solvers::Result result{solvers::solve<models::LGK>(dT, C0, A, V0, R0)};
            if (result.hasConverged)
            {
                outf << result.commaSeparatedValues() << '\n';
                std::tie(V0, R0) = std::tie(result.V, result.R);
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
