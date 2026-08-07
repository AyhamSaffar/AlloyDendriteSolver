// script used to test this library's implemention against existing published Rs for different alloy systems

#include <string>
#include <fstream>
#include <array>
#include <cmath>
#include <tuple>
#include "solvers.h"
#include "alloys.h"
#include "models.h"
#include "approximators.h"
#include "optimiser.h"


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
            solvers::Result R{solvers::newton<models::LGK>(dT, C0, alloys::AlFe_wtp, V0, R0)};
            outfAlFe << R.commaSeparatedValues() << '\n';
            if (R.hasConverged)
                std::tie(V0, R0) = std::tie(R.V, R.R);
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
        // model in paper did not include any kind of kinetic undercooling, which makes a difference at high V
        outfNiSn << solvers::newton<models::LGK>(dT, C0, alloys::NiSn_wtp).commaSeparatedValues() << '\n';

    // Fourth Conference on Rapid Solidification Processing: Principles and Technologies, Application of dendritic
    // growth theory to the interpretation of rapid solidification microstructures, pages 13-25, W.J. Boettinger, S.R.
    // Coriell and R. Trivedi. Fig. 2
    std::ofstream outfAgCu{dataPath + "AgCu_LKT_BCT.csv"};
    outfAgCu << solvers::Result::commaSeparatedColumns << ",k,Cs\n";
    {
        const alloys::Alloy A{alloys::AgCu_wtp};
        double C0{15}, dT0{1};
        
        double V0{approx::getR(dT0, C0, A)}, R0{approx::getR(dT0, C0, A)};
        for (double dT{dT0}; dT<=350; ++dT)
        {
            solvers::Result R{solvers::newton<models::LKT_BCT>(dT, C0, A, V0, R0)};
            if (R.hasConverged) // model diverges if approx funcs always used as initial guess for V and R
                std::tie(V0, R0) = std::tie(R.V, R.R); 

            double Pc{R.V*R.R/(2*A.D)}; // solutal Péclet number
            double Ivc{models::ivantsov(Pc)}; // solutal Ivantsov function
            double k{(A.k0+(A.a0*R.V/A.D)) / (1+(A.a0*R.V/A.D))}; // velocity dependant partition coefficient
            double Cl{C0/(1-Ivc*(1-k))}; // interface liquid solute conentration
            double Cs{k*Cl}; // interface solid solute concentration
                
            outfAgCu << R.commaSeparatedValues() << ',' << k << ',' << Cs << '\n';
        }
    }


    // https://doi.org/10.1016/j.actamat.2016.09.047 Fig. 3, 4, & 5. While paper does not assume the dilute limit for
    // velocity dependent k in LKT-BCT, k0 is so large that it is unlikely to affect the results in a noticable way.
    std::ofstream outfFeCoGamma{dataPath + "FeCoGamma_LKT_BCT.csv"};
    std::ofstream outfFeCoDelta{dataPath + "FeCoDelta_LKT_BCT.csv"};
    outfFeCoGamma << solvers::Result::commaSeparatedColumns << '\n';
    outfFeCoDelta << solvers::Result::commaSeparatedColumns << '\n';

    std::array C0s{30.0, 40.0, 50.0};
    std::array AGammas{alloys::FeCoGamma_30atp, alloys::FeCoGamma_40atp, alloys::FeCoGamma_50atp};
    std::array ADeltas{alloys::FeCoDelta_30atp, alloys::FeCoDelta_40atp, alloys::FeCoDelta_50atp};

    for (std::size_t i{0}; i<=2; ++i)
    {
        double C0{C0s[i]};
        alloys::Alloy AGamma{AGammas[i]}, ADelta{ADeltas[i]};
        
        double V0Gamma{approx::getV(1.0, C0, AGamma)}, R0Gamma{approx::getR(1.0, C0, AGamma)};
        for (double dT{1}; dT<=350; ++dT)
        {
            constexpr bool legacy{false}; // paper uses a slighly more modern form of LKT_BCT
            constexpr models::ModelFunc bct{models::LKT_BCT<legacy>};

            solvers::Result R{solvers::newton<bct>(dT, C0, AGamma, V0Gamma, R0Gamma)};
            outfFeCoGamma << R.commaSeparatedValues() << '\n';
            std::tie(V0Gamma, R0Gamma) = std::tie(R.V, R.R); // model diverges for Gamma if not given a good first guess
            
            outfFeCoDelta << solvers::newton<bct>(dT, C0, ADelta).commaSeparatedValues() << '\n';
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
            constexpr bool legacy{false}; // both LGK and LKT-BCT have a slighly more modern form in this paper
            constexpr models::ModelFunc lgk{models::LGK<legacy>}; // form consistent with LKT-BCT
            constexpr models::ModelFunc bct{models::LKT_BCT<legacy>};

            outfSnAgLGK << solvers::newton<lgk>(dT, C0, alloys::SnAg_wtp).commaSeparatedValues() << '\n';
            outfSnAgLKTBCT << solvers::newton<bct>(dT, C0, alloys::SnAg_wtp).commaSeparatedValues() << '\n';
        }


    // https://doi.org/10.1103/PhysRevB.45.5019 Fig. 1 & 2b
    //! currently does not quite match published Rs at higher undercoolings. Assumed wrong paramters in paper.
    std::ofstream outfNiB{dataPath + "NiB_LKT_BCT.csv"};
    outfNiB << solvers::Result::commaSeparatedColumns << ",Cl,Cs\n";

    const alloys::Alloy A{alloys::NiB1991_atp};
    for (double C0: std::array{0.0, 0.7, 1.0})
    {
        double dT0{1}, C00{(C0==0) ? 0.1 : C0}; // approx module cannot handle 0 C0
        double V0{approx::getV(dT0, C00, A)}, R0{approx::getR(dT0, C00, A)};
        for (double dT{dT0}; dT<=400; ++dT)
        {
            solvers::Result R{solvers::newton<models::LKT_BCT>(dT, C0, A, V0, R0)};
            double Pc{R.V*R.R/(2*A.D)}; // solutal Péclet number
            double Ivc{models::ivantsov(Pc)}; // solutal Ivantsov function
            double k{(A.k0+(A.a0*R.V/A.D)) / (1+(A.a0*R.V/A.D)-(1-A.k0)*(C0/100))}; // velocity dependant k
            double Cl{C0/(1-Ivc*(1-k))}; // interface liquid solute conentration
            double Cs{k*Cl}; // interface solid solute concentration

            outfNiB << R.commaSeparatedValues() << ',' << Cl << ',' << Cs << '\n';
            if (R.hasConverged)
                std::tie(V0, R0) = std::tie(R.V, R.R);
        }
    }

    
    // https://doi.org/10.1007/s11433-010-4167-y Fig. 2, 3, 5 & 6. 20wt.% R values assume linear liquidus and solidus, 
    // which is not true beyond 100K dT. Therefore the CLW model is used instead of LKT-BCT.
    std::ofstream outfCoCuCLW{dataPath + "CoCu_CLW.csv"};
    outfCoCuCLW << solvers::Result::commaSeparatedColumns << ",k0\n";

    for (double C0: {20.0, 60.0})
    {
        const alloys::Alloy A{(C0==20) ? alloys::CoCu_20wtp : alloys::CoCu_60wtp};
        double dT0{1};
        
        solvers::Result R{solvers::bruteForceNewton<models::CLW>(dT0, C0, A)};
        double V0{R.V}, R0{R.R}; // cannot use approx module for non linear phase diagrams

        double Tl{A.TlAtC(C0)};
        for (double dT{dT0}; dT<=400; ++dT) // cannot go below 316K due to phase diagram fit
        {
            R = solvers::newton<models::CLW>(dT, C0, A, V0, R0);
            if (R.hasConverged) // no solution exists = model predicts such high undercooling is impossible
                std::tie(V0, R0) = std::tie(R.V, R.R);

            double k0{A.CsAtT(Tl-dT) / A.ClAtT(Tl-dT)};
            outfCoCuCLW << R.commaSeparatedValues() << ',' << k0 << '\n';
        }
    }

    // https://doi.org/10.1080/09500830903002356 Fig. 2 & 3.
    std::ofstream outfFeSb{dataPath + "FeSb_CLW.csv"};
    outfFeSb << solvers::Result::commaSeparatedColumns << ",k0,kv\n";

    {
        const alloys::Alloy A{alloys::FeSb_wtp};
        const double C0{8.5}, dT0{1};
        solvers::Result R{solvers::bruteForceNewton<models::CLW>(dT0, C0, A)};
        double V0{R.V}, R0{R.R};

        for (double dT{dT0}; dT<=450; ++dT)
        {
            R = solvers::newton<models::CLW>(dT, C0, A, V0, R0);
            if (R.hasConverged)
                std::tie(V0, R0) = std::tie(R.V, R.R);
                            
            double Tl(A.TlAtC(C0)); // liquidus T at bulk C0
            double D{A.DAtT(Tl)}; // diffusivity constant
            double k0{A.CsAtT(Tl)/A.ClAtT(Tl)}; // equilibrium partition coefficient
            double kv{(k0+(A.a0*R.V/D)) / (1+(A.a0*R.V/D))}; // non equilibrium partition coefficient

            outfFeSb << R.commaSeparatedValues() << ',' << k0 << ',' << kv << '\n';
        }      
    }


    return 0;
}

