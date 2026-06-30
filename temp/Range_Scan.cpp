// script used to newton V and R for range of dT and C0

#include <string>
#include <fstream>
#include "solvers.h"
#include "alloys.h"
#include "models.h"


int main()
{
    std::string dataPath{DATA_PATH};
    std::ofstream outf{dataPath + "/data.csv"};
    outf << solvers::Result::commaSeparatedColumns << '\n';

    for (double dT{0.5}; dT<=50; dT+=0.5)
        for(double C0{3.5}; C0<=5.0; C0+=1.5)
        {
            solvers::Result result{solvers::newton<models::LGK>(dT, C0, alloys::SnAg_wtp)};
            outf << result.commaSeparatedValues() << '\n';
        }

    return 0;
}
