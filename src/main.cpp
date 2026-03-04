#include <iostream>
#include <tuple>
#include "alloy.h"
#include "approximators.h"
// #include "differentials.h"
#include "models.h"
//TODO #include "logger.h"


#include "enzyme.h"
using MyFunc = double (*)(double, double, double, double, const alloy::Alloy&);
double testFunc(double V, double R, double dT, double C0, const alloy::Alloy& A)
{
    return V*R*dT*C0*A.a;
}

template <MyFunc func>
double wrapper(double V, double R, double dT, double C0, const alloy::Alloy& A)
{
    return func(V, R, dT, C0, A);
}


int main()
{
    double V{1e-5};
    double R{1e-6};
    double dT{0.5};
    double C0{5};
    alloy::Alloy alloy{alloy::SnAg};

    // std::tuple f{models::LGK(V, R, dT, C0, alloy)};
    // std::cout << "f1: " << std::get<0>(f) << " f2: " << std::get<1>(f) << '\n';

    double dTestFunc{__enzyme_autodiff<double>(
        (void*)wrapper<testFunc>,
        enzyme_out, V,
        enzyme_const, R, dT, C0, alloy, alloy //? no idea why alloy must be passed twice
    )};

    std::cout << "dTestFunc= "<< dTestFunc << '\n';

    return 0;
}
