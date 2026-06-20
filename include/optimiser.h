// functionality for using calculated gradients to update predicted parameters to reduce model residuals

#ifndef OPTIMISER_H
#define OPTIMISER_H

#include <tuple>
#include <exception>
#include "differentials.h"
#include "models.h"
#include "alloys.h"


namespace optimisers
{
    inline std::tuple<double, double> newtonRaphson(double f1, double f2, const diff::Jacobian& J);

    template <models::ModelFunc MODEL>
    using LineSearchFunc = double (*)(double, double, double, double, const alloys::Alloy&, double, double, double);

    template <models::ModelFunc MODEL>
    inline double halfLineSearch(
        double V, double R, double dT, double C0, const alloys::Alloy& A, double dV, double dR, double prevFNorm
    );

    // template <models::ModelFunc MODEL>
    // inline double wolfeLineSearch(
    //     double V, double R, double dT, double C0, const alloys::Alloy& A, double dV, double dR, double prevFNorm
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

    /// @brief searches along update direction for a step size that gives lower model residuals than the previous
    /// optimiser step. Starting with a step size of 1, if the step does not give lower model residuals, the step size
    /// is halved. Fails if step size of 2^{-8} is reached.
    /// @tparam MODEL coupled equations that score how consistent a V R pair are with the given Alloy, C0, and dT.
    /// @param V velocity - m/s
    /// @param R dendrite tip radius - m
    /// @param dT undercooling - K
    /// @param C0 bulk alloy solute concentration - wt.%
    /// @param A struct containing key physical alloy parameters
    /// @param dV velocity direction to search - m/s
    /// @param dR dendrite tip radius direction to search - m
    /// @param prevFNorm euclidean norm of model f1 and f2 from previous optimisation step
    /// @return optimal step size found. If search fails, returns -1.
    template <models::ModelFunc MODEL>
    inline double halfLineSearch(
        double V, double R, double dT, double C0, const alloys::Alloy& A, double dV, double dR, double prevFNorm
    )
    {
        double a{1}, f1{}, f2{}; // step size
        for (int nAttemps{0}; nAttemps<=8; ++nAttemps)
        {
            std::tie(f1, f2) = MODEL(V+a*dV, R+a*dR, dT, C0, A);
            double fNorm{std::sqrt(f1*f1 + f2*f2)};
            if (fNorm < prevFNorm)
                return a;
            a /= 2;
        }
        return -1;
    }

    /* matrix library required for readable implementation
    
    /// @brief binary searches the update direction for a step size that gives a better solution by using the 2 wolfe
    /// conditions. Well explain by Michel Bierlaire (EPFL) in videos 5-9 of the below playlist:
    /// https://www.youtube.com/watch?v=3wh_TLzuiRI&list=PL10NOnsbP5Q7wNrYItE2GhKq05cVov97e. Application to root finding
    /// explained in https://mathsfromnothing.au/line-search-for-root-finding-methods/. Returns -1 if search fails.
    ///
    /// The first wolfe condition sets an upper bound for the step size by making sure that the decrease in F is at
    /// least as big as β1 * the step size * the gradient of F resolved in the update direction, where β1 is typically
    /// = 1e-4. This condition prevents the update from being so large that the linear tangent calculated at the
    /// starting V & R is no longer valid and you have essentially moved to another parabola instead of a lower point on
    /// the parabola you started on.
    ///
    /// The second wolfe condition sets a lower bound for the step size by making sure the ratio of the gradient of F
    /// resolved in the update direction between the end of the step and the start of the step is lower than β2, where
    /// β2 typically = 0.9. Note the gradient decreases as you descend down a parabola, reaching 0 at it's minimum. This
    /// condition ensures you are making some progress in descending the parabola instead of taking many an infinite
    /// number of tiny steps.
    ///
    /// @tparam MODEL coupled equations that score how consistent a V R pair are with the given Alloy, C0, and dT.
    /// @param V velocity - m/s
    /// @param R dendrite tip radius - m
    /// @param dT undercooling - K
    /// @param C0 bulk alloy solute concentration - wt.%
    /// @param A struct containing key physical alloy parameters
    /// @param dV velocity direction to search - m/s
    /// @param dR dendrite tip radius direction to search - m
    /// @param prevFNorm euclidean norm of model f1 and f2 from previous optimisation step
    /// @return optimal step size found. If search fails, returns -1.
    template <models::ModelFunc MODEL>
    inline double wolfeLineSearch(
        double V, double R, double dT, double C0, const alloys::Alloy& A, double dV, double dR, double prevFNorm
    )
    {
        const double b1{1e-4}, b2{0.9}; // wolfe constants
        double aLower{0}, aUpper{1.0}, f1{}, f2{}, dVNew{}, dRNew{};

        for (double step{0}; step<=10; ++step)
        {
            double a{(aUpper-aLower)/2.0}; // step size
            std::tie(f1, f2) = MODEL(V+a*dV, R+a*dR, dT, C0, A);
            double fNorm{std::sqrt(f1*f1 + f2*f2)};
            diff::Jacobian J{diff::calculateGrads<MODEL>(V+a*dV, R+a*dR, dT, C0, A)};
            std::tie(dVNew, dRNew) = newtonRaphson(f1, f2, J);

            // checking if a too high using first wolfe condition
            bool w1F1{};
        }

    };
    */
}

#endif
