#ifndef SOLVER_H
#define SOLVER_H

#include <tuple>
#include <string>
#include <sstream>
#include "alloys.h"
#include "approximators.h"
#include "differentials.h"
#include "models.h"
#include "optimiser.h"


// high level interface for the most common use case - solving V and R for a given Alloy, C0, and dT
namespace solver{

    // struct to hold and log data from a solve attempt
    struct Result{
        bool hasDiverged{};
        bool hasConverged{};
        int steps{};
        double dT{};
        double C0{};
        double V{};
        double R{};
        double f1{};
        double f2{};
        
        inline std::string commaSeparatedValues();
        static inline std::string commaSeparatedColumns{"diverged,converged,steps,dT,C0,V,R,f1,f2"};
    };
    
    /// @brief high level interface for iteratively solving V and R for a given Alloy, C0, and dT.
    /// @tparam MODEL coupled equations that score how consistent a V R pair are with the given Alloy, C0, and dT.
    /// @param dT undercooling - K
    /// @param C0 bulk alloy solute concentration - wt.%
    /// @param A struct containing key physical alloy parameters
    /// @param V0 initial guess for velocity - m/s
    /// @param R0 initial guess for dendrite tip radius - m
    /// @return struct containing V, R, dT, and C0 as well as optimisation flags and parameters
    template <models::ModelFunc MODEL>
    inline Result solve(double dT, double C0, const alloys::Alloy& A, double V0, double R0)
    {
        double f1{}, f2{}, V{V0}, R{R0}, dV{}, dR{};
        diff::Jacobian J{};
        int maxSteps{1000}, step{0};
        bool converged{false}, diverged{false};
    
        for (; step<maxSteps; ++step)
        {
            std::tie(f1, f2) = MODEL(V, R, dT, C0, A);
            J = diff::calculateGrads<MODEL>(V, R, dT, C0, A);
            std::tie(dV, dR) = optimisers::newtonRaphson(f1, f2, J);
            if (std::isnan(dV) || std::isnan(dR))
            {
                diverged = true;
                break;
            }
            V += 0.1*dV; // smaller steps increase the range of starting V and R that don't diverge
            R += 0.1*dR;
            if ((std::abs(f1)<1e-12) && (std::abs(f2)<1e-12))
            {
                converged = true;
                break;
            }
        }
        
        return Result{diverged, converged, step, dT, C0, V, R, f1, f2};
    }

    /// @brief high level interface for iteratively solving V and R for a given Alloy, C0, and dT. Uses approx module
    /// to predict a fair starting guess for V and R.
    /// @tparam MODEL coupled equations that score how consistent a V R pair are with the given Alloy, C0, and dT.
    /// @param dT undercooling - K
    /// @param C0 bulk alloy solute concentration - wt.%
    /// @param A struct containing key physical alloy parameters
    /// @return struct containing V, R, dT, and C0 as well as optimisation flags and parameters
    template <models::ModelFunc MODEL>
    inline Result solve(double dT, double C0, const alloys::Alloy& A)
    {
        return solve<MODEL>(dT, C0, A, approx::getTipVelocity(dT, C0, A), approx::getTipRadius(dT, C0, A));
    }
    
    /// @brief parses struct member variables in a csv compliant form.
    /// @return string of member variables seperated by commas. Does not include trailing newline character.
    inline std::string Result::commaSeparatedValues()
    {
        std::stringstream values{};
        values << std::boolalpha;
        values << hasDiverged << ',' << hasConverged << ',' << steps << ',' <<  dT << ',' << C0 << ',' << V << ',' << R
               << ',' << f1 << ',' << f2 ;
        return values.str();
    }
}

#endif
