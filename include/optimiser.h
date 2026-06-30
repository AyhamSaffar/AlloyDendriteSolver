// functionality for using calculated gradients to update predicted parameters to reduce model residuals

#ifndef OPTIMISER_H
#define OPTIMISER_H

#include <tuple>
#include "differentials.h"
#include "models.h"
#include "alloys.h"


namespace optimisers
{
    inline std::tuple<double, double> newtonRaphson(double f1, double f2, const diff::Jacobian& J);

    // inline double halfLineSearch(
    //     models::ModelFunc model, double V, double R, double dT, double C0, const alloys::Alloy& A,
    //     double dV, double dR, double prevFNorm
    // );

    /// @brief uses Newton Raphson optimisation to find how to update parameters in order to reduce model residuals.
    /// given J∆=−F, ∆=-inv(J)F. This function calculates the inverse of the Jacobian and manually calculates it's dot
    /// product with the model residuals.
    /// @param f1 residual from the f1 model function. 
    /// @param f1 residual from the f2 model function. 
    /// @param J Jacobian of F with respect to V and R.
    /// @return tuple containing ∆V and ∆R required to finimise V and R.
    inline std::tuple<double, double> newtonRaphson(double f1, double f2, const diff::Jacobian& J)
    {
        double JDet{J.df1dV*J.df2dR - J.df1dR*J.df2dV};
        double dV{-1/JDet * (J.df2dR*f1 - J.df1dR*f2)};
        double dR{-1/JDet * (-J.df2dV*f1 + J.df1dV*f2)};
        return std::tuple{dV, dR};
    }

    // using LineSearch = double (*)(
    //     models::ModelFunc, double, double, double, double, const alloys::Alloy&, double, double, double
    // );

    // /// @brief searches along update direction for a step size that gives lower model residuals than the previous
    // /// optimiser step. Starting with a step size of 1, if the step does not give lower model residuals, the step size
    // /// is halved. Fails if step size of 2^{-8} is reached.
    // /// @param model coupled equations that score how consistent a V R pair are with the given Alloy, C0, and dT.
    // /// @param V velocity - m/s
    // /// @param R dendrite tip radius - m
    // /// @param dT undercooling - K
    // /// @param C0 bulk alloy solute concentration - C.%
    // /// @param dV velocity direction to search - m/s
    // /// @param dR dendrite tip radius direction to search - m
    // /// @param prevFNorm euclidean norm of model f1 and f2 from previous optimisation step
    // /// @return optimal step size found. If search fails, returns -1.
    // inline double halfLineSearch(
    //     models::ModelFunc model, double V, double R, double dT, double C0, const alloys::Alloy& A,
    //     double dV, double dR, double prevFNorm
    // )
    // {
    //     double a{1}, f1{}, f2{}; // step size
    //     for (int nAttemps{0}; nAttemps<=8; ++nAttemps)
    //     {
    //         std::tie(f1, f2) = model(V+a*dV, R+a*dR, dT, C0, A);
    //         double fNorm{std::sqrt(f1*f1 + f2*f2)};
    //         if (fNorm < prevFNorm)
    //             return a;
    //         a /= 2;
    //     }
    //     return -1;
    // }
}

#endif
