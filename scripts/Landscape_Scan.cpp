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
    
    const alloys::Alloy A{alloys::AgCu};
    double C0{15}, f1{}, f2{};
    for (double dT{50}; dT<=150; dT+=50)
        for (double VPower{-7.0}; VPower<=4.0; VPower+=0.01)
        {
            double V{std::pow(10.0, VPower)};
            for (double RPower{-8}; RPower<=-5.2; RPower+=0.01)
            {
                double R{std::pow(10.0, RPower)};
                std::tie(f1, f2) = models::LKT_BCT(V, R, dT, C0, A);
                outfError << dT << ',' << C0 << ',' << V << ',' << R << ',' << f1 << ',' << f2 << '\n';
            } 
        }


    return 0;
}
