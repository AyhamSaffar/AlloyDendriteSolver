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
        bool hasDiverged{};     // calculated dV or dR is NaN
        bool hasConverged{};    // undercooling and radius errors within tolerance
        int steps{};            // number of optimisation steps (excluding line search iterations)
        double dT{};            // undercooling - K
        double C0{};            // bulk alloy solute concentration - C.%
        double V{};             // dendrite velocity - m/s
        double R{};             // dendrite tip radius - m
        double f1{};            // undercooling error
        double f2{};            // radius error
        double dTt{};           // calculated thermal undercooling - K
        double dTc{};           // calculated solutal undercooling - K
        double dTr{};           // calculated curvature undercooling - K
        double dTk{};           // calculated kinetic undercooling - K (not calculated in all models)            
        
        inline std::string commaSeparatedValues();
        static inline std::string commaSeparatedColumns{"diverged,converged,steps,dT,C0,V,R,f1,f2,dTt,dTc,dTr,dTk"};
    };

    using SolverFunc = Result (*)(double, double, const alloys::Alloy&);

    /// @brief Scaled newton method for iteratively solving for V and R. Each iteration approximates the system of
    /// equations with linear tangents at the current V, R pair and updates this pair to the point where those tangents
    /// equal reach zero F. Requires a reasonable initial guess for this pair to converge to the correct solution. 
    /// @tparam MODEL coupled equations that score how consistent a V R pair are with the given Alloy, C0, and dT.
    /// @tparam LINESEARCH line search technique to be used. If none given, solver just uses a step size of 0.01.
    /// @param dT undercooling - K
    /// @param C0 bulk alloy solute concentration - C.%
    /// @param A struct containing key physical alloy parameters
    /// @param V0 initial guess for velocity - m/s. Defaults to -1, which uses approx module to get initial guess.
    /// @param R0 initial guess for dendrite tip radius - m. Defaults to -1, which uses approx module to get initial
    /// guess.
    /// @return struct containing V, R, dT, and C0 as well as optimisation flags and parameters
    template <models::ModelFunc MODEL, optimisers::LineSearch LINESEARCH = nullptr>
    inline Result newton(double dT, double C0, const alloys::Alloy& A, double V0=-1, double R0=-1)
    {
        double V{(V0==-1) ? approx::getV(dT, C0, A): V0};
        double R{(R0==-1) ? approx::getR(dT, C0, A): R0};
        double f1{}, f2{}, dV{}, dR{};
        diff::Jacobian J{};
        int maxSteps{1000}, step{0};
        bool converged{false}, diverged{false};
        models::DTs dTs{};

        for (; step<maxSteps; ++step)
        {
            std::tie(f1, f2, dTs) = MODEL(V, R, dT, C0, A);
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
            double a{0.1}; // smaller step size increase the range of starting V and R that don't diverge
            if constexpr (LINESEARCH != nullptr)
            {
                double prevFNorm{std::sqrt(f1*f1 + f2*f2)};
                a = LINESEARCH(MODEL, V, R, dT, C0, A, dV, dR, prevFNorm);
                if (a==-1) // line search failed
                break;
            }
            V += a*dV; 
            R += a*dR;
            if ((V<0) || (R<0)) // required as this gives rise to negative / very high T, breaking phase diagram fits
            {
                diverged = true;
                break;
            }
        }
        
        return Result{diverged, converged, step, dT, C0, V, R, f1, f2, dTs.t, dTs.c, dTs.r, dTs.k};
    }

    /// @brief Brute force newton solver. Searches every possible starting V (1e-6 to 1e3 m/s) and R (1e-9 to 1e-3 m)
    /// pair and attempts to solve pairs that give reasonable initial guesse (f1<1e-3 or f2<1e-3). Returns a solution
    /// if it converges and is physical (R>0 and V>0). Is not especially slow, as the solver has a good chance at
    /// reaching the optimal solution when starting at a point with low f1 or f2.
    /// @tparam MODEL coupled equations that score how consistent a V R pair are with the given Alloy, C0, and dT.
    /// @tparam LINESEARCH line search technique to be used. If none given, solver just uses a step size of 0.01.
    /// @param dT undercooling - K
    /// @param C0 bulk alloy solute concentration - C.%
    /// @param A struct containing key physical alloy parameters
    /// @return struct containing V, R, dT, and C0 as well as optimisation flags and parameters
    template <models::ModelFunc MODEL, optimisers::LineSearch LINESEARCH = nullptr>
    inline Result bruteForceNewton(double dT, double C0, const alloys::Alloy& A)
    {
        double f1{}, f2{}, searchStep{0.2}; // size of each search step in log10
        for (double V0Power{-6}; V0Power<=3; V0Power+=searchStep)
        {
            double V0{std::pow(10.0, V0Power)};
            for (double R0Power{-9}; R0Power<=-3; R0Power+=searchStep)
            {
                double R0{std::pow(10.0, R0Power)};
                models::DTs _{}; // unused struct
                std::tie(f1, f2, _) = MODEL(V0, R0, dT, C0, A);
                double threshold{1e-3};
                if ((f1>threshold) && (f2>threshold))
                    continue;
                Result result{newton<MODEL, LINESEARCH>(dT, C0, A, V0, R0)};
                if(result.hasConverged && (result.V>0) && (result.R>0))
                    return result;
            }
        }
        bool hasDiverged{false}, hasConverged{false};
        return Result{hasDiverged, hasConverged, 0, dT, C0};
    }


    /// @brief parses struct member variables in a csv compliant form.
    /// @return string of member variables seperated by commas. Does not include trailing newline character.
    inline std::string Result::commaSeparatedValues()
    {
        std::stringstream values{};
        values << std::boolalpha;
        values << hasDiverged << ',' << hasConverged << ',' << steps << ',' <<  dT << ',' << C0 << ',' << V << ',' << R
               << ',' << f1 << ',' << f2 << ',' <<  dTt << ',' << dTc << ',' << dTr << ',' << dTk;
        return values.str();
    }
}

#endif
