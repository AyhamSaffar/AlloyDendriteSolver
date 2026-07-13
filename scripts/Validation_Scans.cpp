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
        for(double C0{0.005}; C0<=1; C0+=0.005)
            outfSucAce << solvers::newton<models::LGK>(dT, C0, alloys::SucAce_atp).commaSeparatedValues() << '\n';


    // https://doi.org/10.1007/BF02643853 Fig. 14
    std::ofstream outfAlFe{dataPath + "/AlFe_LGK.csv"};
    outfAlFe << solvers::Result::commaSeparatedColumns << '\n';
    
    // approx module used for initial V, R guess assumes some solute present, so custom guesses needed here
    {
        double C0{0}, dT0{1.0}; 
        double V0{approx::getV(dT0, C0+0.1, alloys::AlFe_wtp)}, R0(approx::getR(dT0, C0+0.1, alloys::AlFe_wtp));
        for(double dTPower{0}; dTPower<=2.7; dTPower+=0.01)
        {
            double dT{std::pow(10, dTPower)};
            solvers::Result result{solvers::newton<models::LGK>(dT, C0, alloys::AlFe_wtp, V0, R0)};
            outfAlFe << result.commaSeparatedValues() << '\n';
            std::tie(V0, R0) = std::tie(result.V, result.R);
        }
    }

    for (double C0{4}; C0<=8; C0+=4)
        for(double dTPower{0}; dTPower<=2.7; dTPower+=0.01)
        {
            double dT{std::pow(10, dTPower)};
            outfAlFe << solvers::newton<models::LGK>(dT, C0, alloys::AlFe_wtp).commaSeparatedValues() << '\n';
        }

    // https://doi.org/10.1007/BF02646933 Fig. 12 & 13 (early LKT model skipped as this library doesn't support it)
    std::ofstream outfNiSn{dataPath + "NiSn_LGK.csv"};
    outfNiSn << solvers::Result::commaSeparatedColumns << '\n';

    for (double dT{1}, C0{25}; dT<=1000; ++dT)
        outfNiSn << solvers::newton<models::LGK>(dT, C0, alloys::NiSn_wtp).commaSeparatedValues() << '\n';


    // https://doi.org/10.1016/j.actamat.2016.09.047 Fig. 3, 4, & 5
    std::ofstream outfFeCoGamma{dataPath + "FeCoGamma_LKT_BCT.csv"};
    std::ofstream outfFeCoDelta{dataPath + "FeCoDelta_LKT_BCT.csv"};
    outfFeCoGamma << solvers::Result::commaSeparatedColumns << '\n';
    outfFeCoDelta << solvers::Result::commaSeparatedColumns << '\n';

    for (double C0{30}; C0<=50; C0+=10)
    {
        double V0Gamma{approx::getV(1.0, C0, alloys::FeCoGamma)};
        double R0Gamma{approx::getR(1.0, C0, alloys::FeCoGamma)};
        for (double dT{1}; dT<=350; ++dT)
        {
            // model diverges for Gamma at higher dT if approx funcs always used as initial guess for V and R
            solvers::Result result{solvers::newton<models::LKT_BCT>(dT, C0, alloys::FeCoGamma, V0Gamma, R0Gamma)};
            outfFeCoGamma << result.commaSeparatedValues() << '\n';
            std::tie(V0Gamma, R0Gamma) = std::tie(result.V, result.R);
            outfFeCoDelta << solvers::newton<models::LKT_BCT>(dT, C0, alloys::FeCoDelta).commaSeparatedValues() << '\n';
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
            constexpr bool legacy{false}; // template arguement for LGK model that makes it consistent with LKT_BCT
            outfSnAgLGK << solvers::newton<models::LGK<legacy>>(dT, C0, alloys::SnAg_wtp).commaSeparatedValues() << '\n';
            outfSnAgLKTBCT << solvers::newton<models::LKT_BCT>(dT, C0, alloys::SnAg_wtp).commaSeparatedValues() << '\n';
        }


    // https://doi.org/10.1103/PhysRevB.45.5019 Fig. 1 & 2b
    //! currently does not quite match published results at higher undercoolings
    std::ofstream outfNiB{dataPath + "NiB_LKT_BCT.csv"};
    outfNiB << solvers::Result::commaSeparatedColumns << ",Cl,Cs\n";

    const alloys::Alloy A{alloys::NiB1997_atp};
    for (double C0: std::array{0.0, 0.7, 1.0})
    {
        double dT0{1}, C00{(C0==0) ? 0.1 : C0}; // approx module cannot handle 0 C0
        double V0{approx::getV(dT0, C00, A)}, R0{approx::getR(dT0, C00, A)};
        for (double dT{dT0}; dT<=400; ++dT)
        {
            solvers::Result result{solvers::newton<models::LKT_BCT>(dT, C0, A, V0, R0)};
            double Pc{result.V*result.R/(2*A.D)}; // solutal Péclet number
            double Ivc{models::ivantsov(Pc)}; // solutal Ivantsov function
            double k{(A.k0+(A.a0*result.V/A.D)) / (1+(A.a0*result.V/A.D)-(1-A.k0)*(C0/100))}; // velocity dependant k
            double Cl{C0/(1-Ivc*(1-k))}; // interface liquid solute conentration
            double Cs{k*Cl}; // interface solid solute concentration

            outfNiB << result.commaSeparatedValues() << ',' << Cl << ',' << Cs << '\n';
            if (result.hasConverged)
                std::tie(V0, R0) = std::tie(result.V, result.R);
        }
    }

    
    // https://doi.org/10.1007/s11433-010-4167-y, Fig.5
    std::ofstream outfCoCu{dataPath + "CoCu_Dynamic.csv"};
    outfCoCu << solvers::Result::commaSeparatedColumns << '\n';

    {
        const alloys::Alloy A{alloys::CoCu_wtp};
        double C0{60}, dT0{1};
        double V0{approx::getV(dT0, C0, A)}, R0{approx::getR(dT0, C0, A)};
        for (double dT{dT0}; dT<=120; ++dT)
        {
            solvers::Result result{solvers::newton<models::dynamic>(dT, C0, A, V0, R0)};
            outfCoCu << result.commaSeparatedValues() << '\n';
            if (result.hasConverged)
                std::tie(V0, R0) = std::tie(result.V, result.R);
        }
    }

    return 0;
}

