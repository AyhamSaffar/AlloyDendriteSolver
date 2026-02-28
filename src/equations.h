#ifndef EQUATIONS_H
#define EQUATIONS_H

#include "alloy.h"
#include <cmath> // for std::exp and std::expint

/// @brief Lipton Glicksman Kurz (LGK) model equations that analytically predict how solidification dendrites grow
/// into a molten liquid when the interface is in equilibrium. https://doi.org/10.1016/0025-5416(84)90199-X
namespace LGK
{
    /// @brief calculate the LGK undercooling equation rearranged to equal zero. This equation takes into account the
    /// thermal, constitutional, and curvature undercooling. It uses dimensional analysis to solve for solute and heat
    /// transport across an equilibrium solidification parabaloid dendrite and uses phase diagram constants to calculate
    /// the drop in liquidus temperature ahead of the solidification front due to solute enrichment.
    /// @param V velocity - m/s
    /// @param R dendrite tip radius - m
    /// @param dT undercooling - K
    /// @param C0 bulk alloy composition - wt.%
    /// @param A struct containing key physical alloy parameters
    /// @return the error in the f1 equation. If V, R, dt, and C0 are correct, this equation should equal 0.
    double f1(double V, double R, double dT, double C0, const alloy::Alloy& A)
    {
        double Pt{V*R/(2*A.a)}; // thermal Péclet number
        double Pc{V*R/(2*A.D)}; // solutal Péclet number
        double Ivt{Pt*std::exp(Pt)*std::expint(Pt)}; // thermal Ivantsov function
        double Ivc{Pc*std::exp(Pc)*std::expint(Pc)}; // solutal Ivantsov function
        return A.L*Ivt/A.Cp + A.m*C0*(1 - 1/(1-(1-A.k0)*Ivc)) + 2*A.r/R - dT;
    }

    /// @brief calculate the LGK stability criterion equation rearranged to equal zero. The stability criterion is an
    /// optimal ratio between dendrite tip radius and velocity. If this ratio is too small, the dendrite grows slowly
    /// enough that any small deviation to the solidification front would cause that point grow much faster than the
    /// rest of the solidification front. If this ratio is too large, the dendrite grows fast enough that secondary
    /// dendrites would start growing out from the sides of the initial dendrite, causing the dendrite to grow more
    /// slowly overall.
    /// @param V velocity - m/s
    /// @param R dendrite tip radius - m
    /// @param dT undercooling - K
    /// @param C0 bulk alloy composition - wt.%
    /// @param A struct containing key physical alloy parameters
    /// @return the error in the f2 equation. If V, R, dt, and C0 are correct, this equation should equal 0.
    double f2(double V, double R, double dT, double C0, const alloy::Alloy& A)
    {
        double Pt{V*R/(2*A.a)}; // thermal Péclet number
        double Pc{V*R/(2*A.D)}; // solutal Péclet number
        double Ivc{Pc*std::exp(Pc)*std::expint(Pc)}; // solutal Ivantsov function
        return A.r/A.o / ( Pt*A.L/A.Cp - (Pc*A.m*C0*(1-A.k0))/(1-(1-A.k0)*Ivc) ) - R;
    }
}

#endif
