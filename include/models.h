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
    // stores different undercooling components
    struct DTs
    {
        double t{}; // thermal undercooling
        double c{}; // solutal undercooling
        double r{}; // curvature undercooling
        double k{}; // kinetic undercooling
    };

    // template for all function headers in this module
    using ModelFunc = std::tuple<double, double, DTs> (*)(double, double, double, double, const alloys::Alloy&);

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

    /// @brief Lipton, Glicksman, and Kurz model. Useful at moderate undercoolings and velocities (VR/2D << 2π) and for
    /// fully linear phase diagrams.
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
    inline std::tuple<double, double, DTs> LGK(double V, double R, double dT, double C0, const alloys::Alloy& A)
    {
        double Pt{V*R/(2*A.a)}; // thermal Péclet number
        double Pc{V*R/(2*A.D)}; // solutal Péclet number
        double Ivt{ivantsov(Pt)}; // thermal Ivantsov function
        double Ivc{ivantsov(Pc)}; // solutal Ivantsov function
        double Ci{C0/(1-(1-A.k0)*Ivc)}; // solute concentration of liquid at interface

        double factor{LEGACY ? 1 : 2}; // solutal field gradient factor
        double dTt{A.L*Ivt/A.Cp}, dTc{A.m*(C0-Ci)}, dTr{2*A.r/R}; // undercooling components
        double f1{dTt+dTc+dTr-dT}; // undercooling error
        double f2{(A.r/A.o) / (Pt*A.L/A.Cp - factor*Pc*A.m*(1-A.k0)*Ci) - R}; // radius error
        return std::make_tuple(f1, f2, DTs{dTt, dTc, dTr});
    }

    /// @brief Lipton, Kurz, and Trivedi - Boettinger Coriell and Trivedi model. Generalises better to higher
    /// undercoolings and velocities for fully linear phase diagrams.
    /// @tparam NO_PARTITIONING whether to disable the solute from ever crossing the solidification front. This is
    /// usually true at high V, but setting this value to true ensures this at low V aswell. Defaults to false.
    /// @param V velocity - m/s
    /// @param R dendrite tip radius - m
    /// @param dT undercooling - K
    /// @param C0 bulk alloy solute concentration - C.%
    /// @param A struct containing key physical alloy parameters
    /// @return dT and R errors. If V, R, dt, and C0 are perfectly correct, both errors should be zero.
    template <bool NO_PARTITIONING=false>
    inline std::tuple<double, double, DTs> LKT_BCT(double V, double R, double dT, double C0, const alloys::Alloy& A)
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
        double Ci{C0/(1-(1-k)*Ivc)}; // solute concentration of liquid at interface

        double dTt{A.L*Ivt/A.Cp}, dTc{A.m*C0 - mP*Ci}, dTr{2*A.r/R}, dTk{V/mu}; // undercooling components
        double f1{dTt+dTc+dTr+dTk-dT}; // undercooling error
        double f2{(A.r/A.o) / (xit*Pt*A.L/A.Cp - 2*A.m*Pc*(1-k)*xic*Ci) - R}; // radius error
        return std::make_tuple(f1, f2, DTs{dTt, dTc, dTr, dTk});
    }

    /// @brief Cao, Wang, Duan, and Bai model. Generalises better to higher undercoolings and velocities for non-linear
    /// phase diagrams.
    /// @param V velocity - m/s
    /// @param R dendrite tip radius - m
    /// @param dT undercooling - K
    /// @param C0 bulk alloy solute concentration - C.%
    /// @param A struct containing key physical alloy parameters
    /// @return dT and R errors. If V, R, dt, and C0 are perfectly correct, both errors should be zero.
    inline std::tuple<double, double, DTs> CLW(double V, double R, double dT, double C0, const alloys::Alloy& A)
    {
        if (!A.CLWCapable)
            throw std::runtime_error("Attempted to pass non CLW capable Alloy to CLW model");

        double Tl{A.TlAtC(C0)}; // liquidus temperature
        double D{A.DAtT(Tl-dT)}; // diffusivity constant
        double m{A.mAtC(C0)}; // liquidus gradient
        double k0{A.k0AtT(Tl-dT)}; // equilibrium partition coefficient
    
        double Pt{V*R/(2*A.a)}; // thermal Péclet number
        double Pc{V*R/(2*D)}; // solutal Péclet number
        double Ivt{ivantsov(Pt)}; // thermal Ivantsov function
        double Ivc{ivantsov(Pc)}; // solutal Ivantsov function
        
        // assumes dilute limit for solute trapping
        double k{(k0+(A.a0*V/D)) / (1+(A.a0*V/D))}; // velocity dependant partition coefficient
        double R0{8.314}; // gas constant
        // BCT paper uses Tm while this model uses Tl
        double mu{A.L*A.V0/(R0*Tl*Tl)}; // interfacial kinetic coefficient
        double xit{1 - 1/std::sqrt(1 + 1/(A.o*Pt*Pt))}; // thermal stability function
        double xic{1 + 2*k/( 1-2*k-std::sqrt(1 + 1/(A.o*Pc*Pc)) )}; // solutal stability function
        double Ci{C0/(1-(1-k)*Ivc)}; // interface solute concentration

        double dTt{A.L*Ivt/A.Cp}, dTc{Tl-A.TlAtC(Ci)}, dTr{2*A.r/R}, dTk{V/mu}; // undercooling components
        double f1{dTt+dTc+dTr+dTk-dT}; // undercooling error
        // Paper divides by xic instead of times by xic, but this must be a missprint.
        double f2{(A.r/A.o) / (xit*Pt*A.L/A.Cp - 2*m*(1-k)*Pc*xic*Ci) - R}; // radius error
        return std::make_tuple(f1, f2, DTs{dTt, dTc, dTr, dTk});
    }
}

#endif
