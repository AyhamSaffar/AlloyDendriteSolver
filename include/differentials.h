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

    //! As of 25/07/26, Enzyme has a bug where it can sometimes modify arguements passed to it. This means subsequent
    //! calls to __enyme_autodiff evalaute the gradient at the wrong point, giving bad results without any errors being
    //! raised. At the point of writing this, this only occurs when the CLW model is used, where A's m_TlAtCFit private
    //! member variable gets modified. I assume this happens because the CLW model is currently the only model that
    //! calls a member function of A. That is why A is passed by copy to the below function. If passed by const
    //! reference, the modified A will persist and silently cause issues.
    template <models::ModelFunc modelFunc>
    inline Jacobian calculateGrads(double V, double R, double dT, double C0, alloys::Alloy A)
    {
        auto [df1dV, df1dR] = __enzyme_autodiff<Diffs>(
            (void*)wrapper<modelFunc, 1>,
            enzyme_out, V, enzyme_out, R, enzyme_const, dT, enzyme_const, C0, enzyme_const, &A
        );
        //! point where A gets modified by CLW model. Luckily the member of A that gets modified is not used to 
        //! calculate any of the f2 grads, 
        auto [df2dV, df2dR] = __enzyme_autodiff<Diffs>(
            (void*)wrapper<modelFunc, 2>,
            enzyme_out, V, enzyme_out, R, enzyme_const, dT, enzyme_const, C0, enzyme_const, &A
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
