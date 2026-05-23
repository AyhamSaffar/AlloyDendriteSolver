// script used to solve V and R for range of dT and C0

#include <string>
#include <fstream>
#include "solver.h"
#include "alloys.h"
#include "models.h"


int main()
{
    std::string dataPath{DATA_PATH};
    std::ofstream outf{dataPath + "/data.csv"};
    outf << solver::Result::commaSeparatedColumns << '\n';

    for (double dT{0.5}; dT<=50; dT+=0.5)
        for(double C0{3.5}; C0<=5.0; C0+=1.5)
        {
            solver::Result result{solver::solve<models::LGK>(dT, C0, alloys::SnAg)};
            outf << result.commaSeparatedValues() << '\n';
        }

    return 0;
}
