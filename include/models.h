#ifndef MODELS_H
#define MODELS_H

#include <cmath> // for std::exp and std::expint
#include <tuple>
#include "alloy.h"

/// @brief standard models that help calculate solidification parameters from alloy physical parameters. These assume
/// a single nucleation event, for example in small liquid solder balls that don't have any available nucleants.
namespace models
{
    // template for all function headers in this module
    using ModelFunc = std::tuple<double, double> (*)(double, double, double, double, const alloy::Alloy&);

    // cmath module uses a slightly different form of exponential integral compared to what is needed here
    inline double expint(double x) {return -std::expint(-x);}

    /// @brief Lipton, Glicksman, and Kurz model. Useful at moderate undercoolings and velocities (VR/2D << 2π).
    ///
    /// The first equation calculates the LGK dendrite undercooling. It quantifies how the liquid must be cooled below
    /// the temperature of the solid to 1. drive thermal diffusion away from the solid that gives out heat as it
    /// solidifies, 2. reach the lower melting temperature caused by a build up of solute just ahead of the
    /// solidification front, and 3. overcome the energy barrier created by the surface energy of a high curvature
    /// dendrite tip. It uses dimensional analysis to solve for solute and heat transport across an equilibrium
    /// solidification parabaloid dendrite. Phase diagram constants are used to calculate the drop in liquidus
    /// temperature ahead of the solidification front due to solute enrichment.
    ///
    /// The second equation calculates the LGK marginal stability criterion dendrite radius. A planar solidification
    /// front is modified by adding a periodic pertubation. Too small and the curvature will drive the pertubation to
    /// shrink. Too large and purtubation will grow by escaping the cold and solute rich solidification front. The
    /// dendrite radius is approximated as the smallest pertubation that won't shrink. This gives an expression that is
    /// a function of the solute and temperature field gradient, which can be calculated for a parabaloid dendrite using
    /// the same dimensional analysis as the first equation.
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

    /// @brief Lipton, Glicksman, and Kurz - Boettinger Coriell and Trivedi model. Generalises better to higher
    /// undercoolings and velocities.
    ///
    ///The first equation calculates the LGK dendrite undercooling. It quantifies how the liquid must be cooled below
    /// the temperature of the solid to 1. drive thermal diffusion away from the solid that gives out heat as it
    /// solidifies, 2. reach the lower melting temperature caused by a build up of solute just ahead of the
    /// solidification front, and 3. overcome the energy barrier created by the surface energy of a high curvature
    /// dendrite tip. It uses dimensional analysis to solve for solute and heat transport across an equilibrium
    /// solidification parabaloid dendrite. Phase diagram constants are used to calculate the drop in liquidus
    /// temperature ahead of the solidification front due to solute enrichment.
    ///
    /// The second equation calculates the LGK marginal stability criterion dendrite radius. A planar solidification
    /// front is modified by adding a periodic pertubation. Too small and the curvature will drive the pertubation to
    /// shrink. Too large and purtubation will grow by escaping the cold and solute rich solidification front. The
    /// dendrite radius is approximated as the smallest pertubation that won't shrink. This gives an expression that is
    /// a function of the solute and temperature field gradient, which can be calculated for a parabaloid dendrite using
    /// the same dimensional analysis as the first equation.
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

        double k{ (A.k0+(A.a0*V/A.D)) / (1+(A.a0*V/A.D)) }; // velocity dependant partition coefficient
        double mP{A.m*(1+ (A.k0-k*(1-std::log(k/A.k0))) / (1-A.k0) )}; // velocity dependant liquidus slope (m prime)
        double R0{8.314}; // gas constant
        double mu{A.L*A.V0/(R0*A.Tm*A.Tm)}; // interfacial kinetic coefficient
        double xit{1 - 1/std::sqrt(1 + 1/(A.o*Pt*Pt))}; // thermal stability function
        double xic{1 + 2*k/(1-2*k-std::sqrt(1 + 1/(A.o*Pc*Pc)))}; // - solutal stability function

        double f1{A.L*Ivt/A.Cp + A.m*C0*(1 - (mP/A.m)/(1-(1-k)*Ivc)) + 2*A.r/R + V/mu - dT};
        double f2{(A.r/A.o) / ( xit*Pt*A.L/A.Cp + (2*A.m*C0*(k-1)*xic)/(1-(1-k)*Ivc) ) - R};
        return std::make_tuple(f1, f2);
    }
}

#endif
