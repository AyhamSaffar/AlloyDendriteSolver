// sciprt used to sample the optimisation landscape for the space of possible V & R pairs

#include <fstream>
#include <string>
#include <cmath>
#include "alloys.h"
#include "models.h"


int main()
{
    std::string dataPath{DATA_PATH};
    std::ofstream outfError{dataPath + "/scan_data.csv"};
    outfError << "dT,C0,V,R,f1,f2" << '\n';
    
    const alloys::Alloy A{alloys::NiB1997_atp};
    double C0{1}, f1{}, f2{};
    models::DTs _{};
    for (double dT{100}; dT<=300; dT+=100)
        for (double VPower{-3}; VPower<=3.0; VPower+=0.01)
        {
            double V{std::pow(10.0, VPower)};
            for (double RPower{-8}; RPower<=-5; RPower+=0.01)
            {
                double R{std::pow(10.0, RPower)};
                std::tie(f1, f2, _) = models::LKT_BCT(V, R, dT, C0, A);
                outfError << dT << ',' << C0 << ',' << V << ',' << R << ',' << f1 << ',' << f2 << '\n';
            } 
        }


    return 0;
}
