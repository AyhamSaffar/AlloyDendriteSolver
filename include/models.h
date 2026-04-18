#ifndef MODELS_H
#define MODELS_H

#include <cmath> // for std::exp and std::expint
#include <tuple>
#include "alloy.h"

/// @brief standard models that help calculate solidification parameters from alloy physical parameters.
namespace models
{
    // template for all function headers in this module
    using ModelFunc = std::tuple<double, double> (*)(double, double, double, double, const alloy::Alloy&);

    /// @brief Lipton Glicksman Kurz (LGK) model equations that analytically predict how solidification dendrites grow
    /// into a molten liquid when the interface is in equilibrium. https://doi.org/10.1016/0025-5416(84)90199-X
    ///
    /// The first equation calculates the LGK dendrite undercooling error. This equation takes into
    /// account the thermal, constitutional, and curvature undercooling. It uses dimensional analysis to solve for
    /// solute and heat transport across an equilibrium solidification parabaloid dendrite and uses phase diagram
    /// constants to calculate the drop in liquidus temperature ahead of the solidification front due to solute
    /// enrichment
    ///
    /// The second equation calculates the LGK stability criterion dendrite radius error. The
    /// stability criterion gives an accurate value for dendrite velocity times it's radius squared. Too wide and slow
    /// dendrites split in smaller parallel dendrites. Too narrow and fast dendrites form secondary dendrites that grow
    /// out perpendicularly.
    ///
    /// @param V velocity - m/s
    /// @param R dendrite tip radius - m
    /// @param dT undercooling - K
    /// @param C0 bulk alloy composition - wt.%
    /// @param A struct containing key physical alloy parameters
    /// @return dT and R errors. If V, R, dt, and C0 are perfectly correct, both errors should be zero.
    inline std::tuple<double, double> LGK(double V, double R, double dT, double C0, const alloy::Alloy& A)
    {
        double Pt{V*R/(2*A.a)}; // thermal Péclet number
        double Pc{V*R/(2*A.D)}; // solutal Péclet number
        // <cmath> uses a slightly different form of exponential integral compared to what is needed here
        double Ivt{Pt*std::exp(Pt)*-std::expint(-Pt)}; // thermal Ivantsov function
        double Ivc{Pc*std::exp(Pc)*-std::expint(-Pc)}; // solutal Ivantsov function
        double f1{A.L*Ivt/A.Cp + A.m*C0*(1 - 1/(1-(1-A.k0)*Ivc)) + 2*A.r/R - dT};
        double f2{(A.r/A.o) / ( Pt*A.L/A.Cp - (Pc*A.m*C0*(1-A.k0))/(1-(1-A.k0)*Ivc) ) - R};
        return std::make_tuple(f1, f2);
    }
}

#endif
