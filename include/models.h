#ifndef MODELS_H
#define MODELS_H

#include <cmath> // for std::exp and std::expint
#include <tuple>
#include <stdexcept>
#include "alloys.h"

/// @brief standard models that help calculate solidification parameters from alloy physical parameters. These assume
/// a single nucleation event, for example in small liquid solder balls that don't have any available nucleants.
namespace models
{
    // template for all function headers in this module
    using ModelFunc = std::tuple<double, double> (*)(double, double, double, double, const alloys::Alloy&);

    // cmath module uses a slightly different form of exponential integral compared to what is needed here
    inline double expint(double x) {return -std::expint(-x);}

    // The Ivantsov function is numerically unstable for high Peclet numbers (as P grows, std::exp(P) -> ∞ and
    // expint(P) -> 1/∞), so an upper bound check is needed to prevent floating point overflow errors.
    inline double ivantsov(double p)
    {
        if (p<200)
            return p*std::exp(p)*expint(p);
        else
            return 1;
    }

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
    /// Note the extra factor of 2 in the second term of the second equation's denominator. Lipton, Glicksman, & Kurz
    /// remove this factor in their paper in order to coerce this equation into agreeing with a prior published result
    /// for the case where there is zero thermal field gradient and the second equation only depends on the solutal
    /// field gradient. This change is not otherwise justified and is ignored in future iterations of this model such as
    /// LKT-BCT.
    ///
    /// @tparam LEGACY whether to remove the factor of 2 in the f2 solutal field gradient term. If true, the model is
    /// consistent with the original paper by Lipton, Glicksman, & Kurz. If false, the model better matches future
    /// iterations of the model. Defaults to true.
    /// @param V velocity - m/s
    /// @param R dendrite tip radius - m
    /// @param dT undercooling - K
    /// @param C0 bulk alloy solute concentration - C.%
    /// @param A struct containing key physical alloy parameters
    /// @return dT and R errors. If V, R, dt, and C0 are perfectly correct, both errors should be zero.
    template <bool LEGACY=true>
    inline std::tuple<double, double> LGK(double V, double R, double dT, double C0, const alloys::Alloy& A)
    {
        double Pt{V*R/(2*A.a)}; // thermal Péclet number
        double Pc{V*R/(2*A.D)}; // solutal Péclet number
        double Ivt{ivantsov(Pt)}; // thermal Ivantsov function
        double Ivc{ivantsov(Pc)}; // solutal Ivantsov function

        double factor{LEGACY ? 1 : 2}; // solutal field gradient factor
        double f1{A.L*Ivt/A.Cp + A.m*C0*(1 - 1/(1-(1-A.k0)*Ivc)) + 2*A.r/R - dT};
        double f2{(A.r/A.o) / ( Pt*A.L/A.Cp - (factor*Pc*A.m*C0*(1-A.k0))/(1-(1-A.k0)*Ivc) ) - R};
        return std::make_tuple(f1, f2);
    }

    /// @brief Lipton, Glicksman, and Kurz - Boettinger Coriell and Trivedi model. Generalises better to higher
    /// undercoolings and velocities.
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
    /// @tparam NO_PARTITIONING whether to disable the solute from ever crossing the solidification front. This is
    /// usually true at high V, but setting this value to true ensures this at low V aswell. Defaults to false.
    /// @param V velocity - m/s
    /// @param R dendrite tip radius - m
    /// @param dT undercooling - K
    /// @param C0 bulk alloy solute concentration - C.%
    /// @param A struct containing key physical alloy parameters
    /// @return dT and R errors. If V, R, dt, and C0 are perfectly correct, both errors should be zero.
    template <bool NO_PARTITIONING=false>
    inline std::tuple<double, double> LKT_BCT(double V, double R, double dT, double C0, const alloys::Alloy& A)
    {
        if (!A.LKT_BCTCapable)
            throw std::runtime_error("Attempted to pass non LKT-BCT capable Alloy to LKT-BCT model");

        double Pt{V*R/(2*A.a)}; // thermal Péclet number
        double Pc{V*R/(2*A.D)}; // solutal Péclet number
        double Ivt{ivantsov(Pt)}; // thermal Ivantsov function
        double Ivc{ivantsov(Pc)}; // solutal Ivantsov function

        double k{NAN};
        if constexpr (NO_PARTITIONING)
            k = 1;
        else
            k = (A.k0+(A.a0*V/A.D)) / (1+(A.a0*V/A.D)-(1-A.k0)*(C0/100)); // velocity dependant partition coefficient

        double mP{A.m*(1+ (A.k0-k*(1-std::log(k/A.k0))) / (1-A.k0) )}; // velocity dependant liquidus slope (m prime)
        double R0{8.314}; // gas constant
        double mu{A.L*A.V0/(R0*A.Tm*A.Tm)}; // interfacial kinetic coefficient
        double xit{1 - 1/std::sqrt(1 + 1/(A.o*Pt*Pt))}; // thermal stability function
        double xic{1 + 2*k/( 1-2*k-std::sqrt(1 + 1/(A.o*Pc*Pc)) )}; // - solutal stability function
        double Ci{C0/(1-(1-k)*Ivc)}; // interface solute concentration

        double f1{A.L*Ivt/A.Cp + (A.m*C0 - mP*Ci) + 2*A.r/R + V/mu - dT};
        double f2{(A.r/A.o) / ( xit*Pt*A.L/A.Cp - (2*A.m*C0*(1-k)*Pc*xic)/(1-(1-k)*Ivc) ) - R};
        return std::make_tuple(f1, f2);
    }

    inline std::tuple<double, double> dynamic(double V, double R, double dT, double C0, const alloys::Alloy& A)
    {
        if (!A.dynamicCapable)
            throw std::runtime_error("Attempted to pass non dynamic capable Alloy to dynamic model");

        double Tl{A.TlAtC(C0)}; // liquidus temperature
        double D{A.DAtT(Tl-dT)}; // diffusivity constant
        double m{A.mAtC(C0)}; // liquidus gradient
        double k0{A.k0AtT(Tl-dT)}; // equilibrium partition coefficient
    
        double Pt{V*R/(2*A.a)}; // thermal Péclet number
        double Pc{V*R/(2*D)}; // solutal Péclet number
        double Ivt{ivantsov(Pt)}; // thermal Ivantsov function
        double Ivc{ivantsov(Pc)}; // solutal Ivantsov function

        double k{(k0+(A.a0*V/D)) / (1+(A.a0*V/D)-(1-k0)*(C0/100))}; // velocity dependant partition coefficient
        double R0{8.314}; // gas constant
        double mu{A.L*A.V0/(R0*Tl*Tl)}; // interfacial kinetic coefficient //? BCT says to use Tm instead of Tl
        double xit{1 - 1/std::sqrt(1 + 1/(A.o*Pt*Pt))}; // thermal stability function
        double xic{1 + 2*k/( 1-2*k-std::sqrt(1 + 1/(A.o*Pc*Pc)) )}; // solutal stability function
        double Ci{C0/(1-(1-k)*Ivc)}; // interface solute concentration

        double f1{A.L*Ivt/A.Cp + (Tl-A.TlAtC(Ci)) + 2*A.r/R + V/mu - dT};
        double f2{(A.r/A.o) / ( xit*Pt*A.L/A.Cp - (2*m*C0*(1-k)*Pc*xic)/(1-(1-k)*Ivc) ) - R};
        return std::make_tuple(f1, f2);
    }
}

#endif
