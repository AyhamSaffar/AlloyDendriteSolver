// functionality for using calculated gradients to update predicted parameters to reduce model residuals

#ifndef OPTIMISER_H
#define OPTIMISER_H

#include <tuple>
#include "differentials.h"


namespace optimisers
{
    /// @brief uses Newton Raphson optimisation to find how to update parameters in order to reduce model residuals.
    /// given J∆=−F, ∆=-inv(J)F. This function calculates the inverse of the Jacobian and manually calculates it's dot
    /// product with the model residuals.
    /// @param f1 residual from the f1 model function. 
    /// @param f1 residual from the f2 model function. 
    /// @param J Jacobian of F with respect to V and R.
    /// @return tuple containing ∆V and ∆R required to finimise V and R.
    std::tuple<double, double> newtonRaphson(double f1, double f2, const diff::Jacobian& J)
    {
        double JDet{J.df1dV*J.df2dR - J.df1dR*J.df2dV};
        double dV{-1/JDet * (J.df2dR*f1 - J.df1dR*f2)};
        double dR{-1/JDet * (-J.df2dV*f1 + J.df1dV*f2)};
        return std::tuple(dV, dR);
    }

}

#endif
