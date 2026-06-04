#ifndef SOLVERS_H
#define SOLVERS_H

#include <tuple>
#include <string>
#include <sstream>
#include "alloys.h"
#include "approximators.h"
#include "differentials.h"
#include "models.h"
#include "optimiser.h"


// high level interface for iterative techniques that solve V and R for a given Alloy, C0, and dT.
namespace solvers{

    // struct to hold and log data from a solver attempt
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

    /// @brief Scaled newton method for iteratively solving for V and R. Each iteration approximates the system of
    /// equations with linear tangents at the current V, R pair and updates this pair to the point where those tangents
    /// equal reach zero F. Requires a reasonable initial guess for this pair to converge to the correct solution. 
    /// @tparam MODEL coupled equations that score how consistent a V R pair are with the given Alloy, C0, and dT.
    /// @param dT undercooling - K
    /// @param C0 bulk alloy solute concentration - wt.%
    /// @param A struct containing key physical alloy parameters
    /// @param V0 initial guess for velocity - m/s. Defaults to -1, which uses approx module to get initial guess.
    /// @param R0 initial guess for dendrite tip radius - m. Defaults to -1, which uses approx module to get initial
    /// guess.
    /// @return struct containing V, R, dT, and C0 as well as optimisation flags and parameters
    template <models::ModelFunc MODEL>
    inline Result newton(double dT, double C0, const alloys::Alloy& A, double V0=-1, double R0=-1)
    {
        double V{(V0==-1) ? approx::getV(dT, C0, A): V0};
        double R{(R0==-1) ? approx::getR(dT, C0, A): R0};
        double f1{}, f2{}, dV{}, dR{};
        diff::Jacobian J{};
        int maxSteps{1000}, step{0};
        bool converged{false}, diverged{false};
    
        for (; step<maxSteps; ++step)
        {
            std::tie(f1, f2) = MODEL(V, R, dT, C0, A);
            if ((std::abs(f1)<1e-12) && (std::abs(f2)<1e-12))
            {
                converged = true;
                break;
            }
            J = diff::calculateGrads<MODEL>(V, R, dT, C0, A);
            std::tie(dV, dR) = optimisers::newtonRaphson(f1, f2, J);
            if (std::isnan(dV) || std::isnan(dR))
            {
                diverged = true;
                break;
            }
            V += 0.1*dV; // smaller steps increase the range of starting V and R that don't diverge
            R += 0.1*dR;
        }
        
        return Result{diverged, converged, step, dT, C0, V, R, f1, f2};
    }


    /// @brief Global Newton method for iteratively solving for V and R. Extends Newton method with backtracking at each
    /// iteration. If the Newton step does not give a better V, R pair (measured by the Euclidean norm of
    /// F), the step size is halved until a better solution is reached.
    /// @tparam MODEL coupled equations that score how consistent a V R pair are with the given Alloy, C0, and dT.
    /// @param dT undercooling - K
    /// @param C0 bulk alloy solute concentration - wt.%
    /// @param A struct containing key physical alloy parameters
    /// @param V0 initial guess for velocity - m/s. Defaults to -1, which uses approx module to get initial guess.
    /// @param R0 initial guess for dendrite tip radius - m. Defaults to -1, which uses approx module to get initial
    /// guess.
    /// @return struct containing V, R, dT, and C0 as well as optimisation flags and parameters
    template <models::ModelFunc MODEL>
    inline Result globalNewton(double dT, double C0, const alloys::Alloy& A, double V0=-1, double R0=-1)
    {
        double V{(V0==-1) ? approx::getV(dT, C0, A): V0};
        double R{(R0==-1) ? approx::getR(dT, C0, A): R0};
        double f1{}, f2{}, dV{}, dR{};
        diff::Jacobian J{};
        int maxSteps{100}, step{0};
        bool converged{false}, diverged{false};
    
        for (; step<maxSteps; ++step)
        {
            std::tie(f1, f2) = MODEL(V, R, dT, C0, A);
            if ((std::abs(f1)<1e-12) && (std::abs(f2)<1e-12))
            {
                converged = true;
                break;
            }
            double fNorm{std::sqrt(f1*f1 + f2*f2)};
            J = diff::calculateGrads<MODEL>(V, R, dT, C0, A);
            std::tie(dV, dR) = optimisers::newtonRaphson(f1, f2, J);
            if (std::isnan(dV) || std::isnan(dR))
            {
                diverged = true;
                break;
            }
            double a{1};
            bool searchSucceeded{false};
            for (int nAttemps{0}, nAttemps<10; ++nAttemps)
            {
                std::tie(f1, f2) = MODEL(V+a*dV, R+a*dR, dT, C0, A);
                if (std::sqrt(f1*f1 + f2*f2) < fNorm)
                {
                    searchSucceeded = true;
                    break;
                }
                a /= 2;
            }
            if (!searchSucceeded)
                break;
            V += a*dV;
            R += a*dR;
        }
        
        return Result{diverged, converged, step, dT, C0, V, R, f1, f2};
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
