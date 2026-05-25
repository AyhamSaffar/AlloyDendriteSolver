// script used to test this library's implemention against existing published results for different alloy systems

#include <string>
#include <fstream>
#include <array>
#include <cmath>
#include "solver.h"
#include "alloys.h"
#include "models.h"


int main()
{
    std::string dataPath{DATA_PATH};


    // https://doi.org/10.1016/0025-5416(84)90199-X Fig. 3 & 4
    std::ofstream outfSucAce{dataPath + "/SucAce_LGK.csv"};
    outfSucAce << solver::Result::commaSeparatedColumns << '\n';

    double SucMr{80.09}, AceMr{58.08};
    for (double dT{0.5}; dT<=0.9; dT+=0.4)
        for(double C0MolPercent{0.01}; C0MolPercent<=1; C0MolPercent+=0.01)
        {
            double C0{ (C0MolPercent*AceMr) / (C0MolPercent*AceMr + (100-C0MolPercent)*SucMr) };
            outfSucAce << solver::solve<models::LGK>(dT, C0, alloys::SucAce).commaSeparatedValues() << '\n';
        }

    
    // https://doi.org/10.1007/BF02643853 Fig. 14
    std::ofstream outfAlFe{dataPath + "/AlFe_LGK.csv"};
    outfAlFe << solver::Result::commaSeparatedColumns << '\n';
    
    for (double C0: std::array{0.1, 4.0, 8.0}) // solver cannot handle 0 C0 value
        for(double dTPower{0}; dTPower<=2.7; dTPower+=0.01)
        {
            double dT{std::pow(10, dTPower)};
            outfAlFe << solver::solve<models::LGK>(dT, C0, alloys::AlFe).commaSeparatedValues() << '\n';
        }
        

    return 0;
}

