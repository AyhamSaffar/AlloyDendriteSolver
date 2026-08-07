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
    template <models::ModelFunc func, int fToReturn>
    inline double wrapper(double V, double R, double dT, double C0, const alloys::Alloy& A)
    {
        std::tuple<double, double, models::DTs> f{func(V, R, dT, C0, A)};
        return std::get<fToReturn-1>(f);
    }
    
    struct Jacobian{double df1dV{}, df1dR{}, df2dV{}, df2dR{};};
    
    struct Diffs{double grad1{}, grad2{};};

    template <models::ModelFunc modelFunc>
    inline Jacobian calculateGrads(double V, double R, double dT, double C0, const alloys::Alloy& A)
    {
        // An Enzyme bug (https://github.com/EnzymeAD/Enzyme/issues/3073) leads to the __enzyme_autodiff always
        // modifying A. As such, temporary copies must be created and passed to this function. 
        static alloys::Alloy ATemp{};
        ATemp = A;
        auto [df1dV, df1dR] = __enzyme_autodiff<Diffs>(
            (void*)wrapper<modelFunc, 1>,
            enzyme_out, V, enzyme_out, R, enzyme_const, dT, enzyme_const, C0, enzyme_const, &ATemp
        );
        ATemp = A;
        auto [df2dV, df2dR] = __enzyme_autodiff<Diffs>(
            (void*)wrapper<modelFunc, 2>,
            enzyme_out, V, enzyme_out, R, enzyme_const, dT, enzyme_const, C0, enzyme_const, &ATemp
        );
        return Jacobian{df1dV, df1dR, df2dV, df2dR};
    }
}

inline std::ostream& operator<<(std::ostream& out, const diff::Jacobian& J)
{
    return out << "Jacobian(δf1/δV=" << J.df1dV << ", δf1/δR=" << J.df1dR << ", δf2/δV=" << J.df2dV << ", δf2/δR=" <<
        J.df2dR << ')';
}

#endif
