// high level interface for the most common use case - solving V and R for a given Alloy, C0, and dT

#ifndef SOLVER_H
#define SOLVER_H

#include <tuple>
#include <string>
#include <sstream>
#include "alloy.h"
#include "approximators.h"
#include "differentials.h"
#include "models.h"
#include "optimiser.h"


namespace solver{

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
        
        std::string commaSeparatedValues();
        static inline std::string commaSeparatedColumns{"diverged,converged,steps,C0,dT,V,R,f1,f2"};
    };
    
    template <models::ModelFunc MODEL>
    Result solve(double dT, double C0, const alloy::Alloy& A, double V0, double R0)
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
            V += 0.1*dV; // smaller steps increase the range of starting V and R that dont diverge
            R += 0.1*dR;
            if ((std::abs(f1)<1e-12) && (std::abs(f2)<1e-12))
            {
                converged = true;
                break;
            }
        }
        
        return Result{diverged, converged, step, dT, C0, V, R, f1, f2};
    }
    
    template <models::ModelFunc MODEL>
    Result solve(double dT, double C0, const alloy::Alloy& A)
    {
        return solve<MODEL>(dT, C0, A, approx::getTipVelocity(dT, C0, A), approx::getTipRadius(dT, C0, A));
    }
    
    
    std::string Result::commaSeparatedValues()
    {
        std::stringstream values{};
        values << std::boolalpha;
        values << hasDiverged << ',' << hasConverged << ',' << steps << ',' <<  dT << ',' << C0 << ',' << V << ',' << R
               << ',' << f1 << ',' << f2 ;
        return values.str();
    }
}

#endif
