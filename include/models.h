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

    // <cmath> uses a slightly different form of exponential integral compared to what is needed here
    double expint(double x) {return -std::expint(-x);}

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
    /// @param C0 bulk alloy solute concentration - wt.%
    /// @param A struct containing key physical alloy parameters
    /// @return dT and R errors. If V, R, dt, and C0 are perfectly correct, both errors should be zero.
    inline std::tuple<double, double> LGK(double V, double R, double dT, double C0, const alloy::Alloy& A)
    {
        double Pt{V*R/(2*A.a)}; // thermal Péclet number
        double Pc{V*R/(2*A.D)}; // solutal Péclet number
        double Ivt{Pt*std::exp(Pt)*expint(Pt)}; // thermal Ivantsov function
        double Ivc{Pc*std::exp(Pc)*expint(Pc)}; // solutal Ivantsov function

        double f1{A.L*Ivt/A.Cp + A.m*C0*(1 - 1/(1-(1-A.k0)*Ivc)) + 2*A.r/R - dT};
        double f2{(A.r/A.o) / ( Pt*A.L/A.Cp - (Pc*A.m*C0*(1-A.k0))/(1-(1-A.k0)*Ivc) ) - R};
        return std::make_tuple(f1, f2);
    }

    /// @brief LKT_BCT explanation
    ///
    /// @param V velocity - m/s
    /// @param R dendrite tip radius - m
    /// @param dT undercooling - K
    /// @param C0 bulk alloy solute concentration - wt.%
    /// @param A struct containing key physical alloy parameters
    /// @return dT and R errors. If V, R, dt, and C0 are perfectly correct, both errors should be zero.
    inline std::tuple<double, double> LKT_BCT(double V, double R, double dT, double C0, const alloy::Alloy& A)
    {
        double Pt{V*R/(2*A.a)}; // thermal Péclet number
        double Pc{V*R/(2*A.D)}; // solutal Péclet number
        double Ivt{Pt*std::exp(Pt)*expint(Pt)}; // thermal Ivantsov function
        double Ivc{Pc*std::exp(Pc)*expint(Pc)}; // solutal Ivantsov function

        double k{(A.k0+(A.a0*V/A.D))/(1+(A.a0*V/A.D))}; // velocity dependant partition coefficient
        double mP{A.m*(1+(A.k0-k*(1-std::log(k/A.k0)))/(1-A.k0))}; // velocity dependant liquidus slope (m prime)
        double Tl{A.Tm + A.m*C0}; // temperature of liquid
        double R0{8.314}; // gas constant
        double mu{A.L*A.V0/(R0*Tl*Tl)}; // interfacial kinetic coefficient
        double xit{1 - 1/std::sqrt(1 + 1/(A.o*Pt*Pt))}; // thermal stability function
        double xic{1 + 2*k/(1-2*k-std::sqrt(1 + 1/(A.o*Pc*Pc)))}; // - solutal stability function

        double f1{A.L*Ivt/A.Cp + A.m*C0*(1 - (mP/A.m)/(1-(1-A.k0)*Ivc)) + 2*A.r/R + V/mu - dT};
        double f2{(A.r/A.o) / ( xit*Pt*A.L/A.Cp - (2*xic*Pc*mP*C0*(1-k))/(1-(1-k)*Ivc) ) - R};
        return std::make_tuple(f1, f2);
    }
}

#endif
