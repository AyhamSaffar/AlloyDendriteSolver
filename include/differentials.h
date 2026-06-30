#ifndef DIFFERENTIALS_H
#define DIFFERENTIALS_H

#include <tuple>
#include <iostream>
#include "enzyme.h" 
#include "alloys.h"
#include "models.h"


namespace diff
{    
    // Enzyme autodiff can only handle non out-parameter functions when they return a single value
    template <models::ModelFunc MODEL, int FTORETURN, typename AlloyLike>
    inline double wrapper(double V, double R, double dT, double C0, const AlloyLike& A)
    {
        std::tuple<double, double> f{MODEL(V, R, dT, C0, A)};
        return std::get<FTORETURN-1>(f);
    }
    struct Jacobian{double df1dV{}; double df1dR{}; double df2dV{}; double df2dR{};};
    
    template <models::ModelFunc MODEL, typename AlloyLike>
    inline Jacobian calculateGrads(double V, double R, double dT, double C0, const AlloyLike& A)
    {
        Jacobian J{};
        J.df1dV = __enzyme_autodiff<double>(
            (void*)wrapper<MODEL, 1>,
            enzyme_out, V, enzyme_const, R, enzyme_const, dT, enzyme_const, C0, enzyme_const, &A
        );
        J.df1dR = __enzyme_autodiff<double>(
            (void*)wrapper<MODEL, 1>,
            enzyme_const, V, enzyme_out, R, enzyme_const, dT, enzyme_const, C0, enzyme_const, &A
        );
        J.df2dV = __enzyme_autodiff<double>(
            (void*)wrapper<MODEL, 2>,
            enzyme_out, V, enzyme_const, R, enzyme_const, dT, enzyme_const, C0, enzyme_const, &A
        );
        J.df2dR = __enzyme_autodiff<double>(
            (void*)wrapper<MODEL, 2>,
            enzyme_const, V, enzyme_out, R, enzyme_const, dT, enzyme_const, C0, enzyme_const, &A
        );
        return J;
    }
}

#endif
