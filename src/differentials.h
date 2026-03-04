#ifndef DIFFERENTIALS_H
#define DIFFERENTIALS_H

#include <tuple>
#include "enzyme.h" 
#include "alloy.h"


namespace diff
{
    struct Jacobian{double df1dV{}; double df1dR{}; double df2dV{}; double df2dR{};};
    using ModelFunc =  std::tuple<double, double> (*)(double, double, double, double, const alloy::Alloy&);

    // Enzyme autodiff can only handle functions that return a single value
    template<int fToReturn>
    double ModelWrapper(ModelFunc modelFunc, double V, double R, double dT, double C0, alloy::Alloy A)
    {
        std::tuple<double, double> f{modelFunc(V, R, dT, C0, A)};
        return std::get<fToReturn-1>(f);
    }

    inline Jacobian calculateGrads(ModelFunc modelFunc, double V, double R, double dT, double C0, alloy::Alloy A)
    {
        double dx{1.0}; // enzyme expects to scale calculated gradients by a given value
        Jacobian J{};
        J.df1dV = __enzyme_autodiff<double>(
            (void*)ModelWrapper<1>,
            enzyme_const, modelFunc,
            enzyme_out, V, dx,
            enzyme_const, R, dT, C0, A
        );

        return J;
    }
}

#endif
