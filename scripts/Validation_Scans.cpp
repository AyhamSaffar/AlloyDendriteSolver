// script used to test this library's implemention against existing published results for different alloy systems

#include <string>
#include <fstream>
#include <array>
#include <cmath>
#include <tuple>
#include "solvers.h"
#include "alloys.h"
#include "models.h"
#include "approximators.h"


int main()
{
    std::string dataPath{DATA_PATH};


    // https://doi.org/10.1016/0025-5416(84)90199-X Fig. 3 & 4
    std::ofstream outfSucAce{dataPath + "/SucAce_LGK.csv"};
    outfSucAce << solvers::Result::commaSeparatedColumns << '\n';

    double SucMr{80.09}, AceMr{58.08};
    for (double dT{0.5}; dT<=0.9; dT+=0.4)
        for(double C0MolPercent{0.005}; C0MolPercent<=1; C0MolPercent+=0.005)
        {
            double C0{ 100 * (C0MolPercent*AceMr) / (C0MolPercent*AceMr + (100-C0MolPercent)*SucMr) }; // wt.% required
            outfSucAce << solvers::solve<models::LGK>(dT, C0, alloys::SucAce).commaSeparatedValues() << '\n';
        }


    // https://doi.org/10.1007/BF02643853 Fig. 14
    std::ofstream outfAlFe{dataPath + "/AlFe_LGK.csv"};
    outfAlFe << solvers::Result::commaSeparatedColumns << '\n';
    
    for (double C0: std::array{0.1, 0.5, 4.0, 8.0}) // approx module cannot handle 0 C0 value for initial V guess
        for(double dTPower{0}; dTPower<=2.7; dTPower+=0.01)
        {
            double dT{std::pow(10, dTPower)};
            outfAlFe << solvers::solve<models::LGK>(dT, C0, alloys::AlFe).commaSeparatedValues() << '\n';
        }


    // https://doi.org/10.1007/BF02646933 Fig. 12 & 13 (early LKT model skipped as this library doesn't support it)
    std::ofstream outfNiSn{dataPath + "NiSn_LGK.csv"};
    outfNiSn << solvers::Result::commaSeparatedColumns << '\n';

    for (double dT{1}, C0{25}; dT<=1000; ++dT)
        outfNiSn << solvers::solve<models::LGK>(dT, C0, alloys::NiSn).commaSeparatedValues() << '\n';


    // https://doi.org/10.1016/j.actamat.2016.09.047 Fig. 3, 4, & 5
    std::ofstream outfFeCoGamma{dataPath + "FeCoGamma_LKT_BCT.csv"};
    std::ofstream outfFeCoDelta{dataPath + "FeCoDelta_LKT_BCT.csv"};
    outfFeCoGamma << solvers::Result::commaSeparatedColumns << '\n';
    outfFeCoDelta << solvers::Result::commaSeparatedColumns << '\n';

    for (double C0{30}; C0<=50; C0+=10)
    {
        double V0Gamma{approx::getTipVelocity(1.0, C0, alloys::FeCoGamma)};
        double R0Gamma{approx::getTipRadius(1.0, C0, alloys::FeCoGamma)};
        for (double dT{1}; dT<=350; ++dT)
        {
            // model diverges for Gamma at higher dT if approx funcs always used as initial guess for V and R
            solvers::Result result{solvers::solve<models::LKT_BCT>(dT, C0, alloys::FeCoGamma, V0Gamma, R0Gamma)};
            outfFeCoGamma << result.commaSeparatedValues() << '\n';
            std::tie(V0Gamma, R0Gamma) = std::tie(result.V, result.R);
            outfFeCoDelta << solvers::solve<models::LKT_BCT>(dT, C0, alloys::FeCoDelta).commaSeparatedValues() << '\n';
        }
    }


    // https://doi.org/10.1007/s10854-025-14979-6 Fig. 11d
    std::ofstream outfSnAgLGK{dataPath + "SnAg_LGK.csv"};
    std::ofstream outfSnAgLKTBCT{dataPath + "SnAg_LKT_BCT.csv"};
    outfSnAgLGK << solvers::Result::commaSeparatedColumns << '\n';
    outfSnAgLKTBCT << solvers::Result::commaSeparatedColumns << '\n';

    for (double C0{3.5}; C0<=5.0; C0+=1.5)
        for (double dT{1}; dT<=50; ++dT)
        {
            outfSnAgLGK << solvers::solve<models::LGK>(dT, C0, alloys::SnAg).commaSeparatedValues() << '\n';
            outfSnAgLKTBCT << solvers::solve<models::LKT_BCT>(dT, C0, alloys::SnAg).commaSeparatedValues() << '\n';
        }


    return 0;
}

