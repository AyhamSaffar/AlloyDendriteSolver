#ifndef MODELS_H
#define MODELS_H

#include <cmath> // for std::exp and std::expint
#include <tuple>
#include <stdexcept>
#include "alloys.h"
#include "enzyme.h"

/// @brief standard models that help calculate solidification parameters from alloy physical parameters. These assume
/// a single nucleation event, for example in small liquid solder balls that don't have any available nucleants.
namespace models
{
    // stores different undercooling components
    struct DTs
    {
        double t{}; // thermal undercooling
        double c{}; // constitutional (solutal) undercooling
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
    /// @tparam LEGACY whether to use m for marginal stability cretieria. This is what BCT used in their origional
    /// paper. If false, uses mP(V). This form tends to be used in more recent papers. Defaults to true.
    /// @param V velocity - m/s
    /// @param R dendrite tip radius - m
    /// @param dT undercooling - K
    /// @param C0 bulk alloy solute concentration - C.%
    /// @param A struct containing key physical alloy parameters
    /// @return dT and R errors. If V, R, dt, and C0 are perfectly correct, both errors should be zero.
    template <bool LEGACY=true>
    inline std::tuple<double, double, DTs> LKT_BCT(double V, double R, double dT, double C0, const alloys::Alloy& A)
    {
        if (!A.LKT_BCTCapable)
            throw std::runtime_error("Attempted to pass non LKT-BCT capable Alloy to LKT-BCT model");

        double Pt{V*R/(2*A.a)}; // thermal Péclet number
        double Pc{V*R/(2*A.D)}; // solutal Péclet number
        double Ivt{ivantsov(Pt)}; // thermal Ivantsov function
        double Ivc{ivantsov(Pc)}; // solutal Ivantsov function

        double k{(A.k0+(A.a0*V/A.D)) / (1+(A.a0*V/A.D))}; // velocity dependent partition coefficient
        double mP{A.m*(1+ (A.k0-k*(1-std::log(k/A.k0))) / (1-A.k0) )}; // velocity dependent liquidus slope (m prime)

        double R0{8.314}; // molar gas constant
        double mu{A.L*A.V0/(R0*A.Tm*A.Tm)}; // interfacial kinetic coefficient
        double xit{1 - 1/std::sqrt(1 + 1/(A.o*Pt*Pt))}; // thermal stability function
        double xic{1 + 2*k/( 1-2*k-std::sqrt(1 + 1/(A.o*Pc*Pc)) )}; // - solutal stability function
        double Ci{C0/(1-(1-k)*Ivc)}; // solute concentration of liquid at interface

        double dTt{A.L*Ivt/A.Cp}, dTc{A.m*C0 - mP*Ci}, dTr{2*A.r/R}, dTk{V/mu}; // undercooling components
        double f1{dTt+dTc+dTr+dTk-dT}; // undercooling error
        double f2m{LEGACY ? A.m : mP}; // liquidus slope used in f2 expression
        double f2{(A.r/A.o) / (xit*Pt*A.L/A.Cp - 2*f2m*Pc*(1-k)*xic*Ci) - R}; // radius error
        return std::make_tuple(f1, f2, DTs{dTt, dTc, dTr, dTk});
    }

    /// @brief Cao, Wang, Duan, and Bai model. Designed to better generalise to higher undercoolings and velocities for
    /// non-linear phase diagrams, but makes strong assumptions and precise implementation details were never published.
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

        double Pt{V*R/(2*A.a)}; // thermal Péclet number
        double Ivt{ivantsov(Pt)}; // thermal Ivantsov function
        double dTt{A.L*Ivt/A.Cp}; // thermal undercooling
        double T0(A.TlAtC(C0)); // liquidus T at bulk C0

        //! CLW papers never detail which T/C to use for D(T), m(C), or k0(T). This only gets close to published results
        double D{A.DAtT(T0)}; // diffusivity constant
        double Pc{V*R/(2*D)}; // solutal Péclet number
        double Ivc{ivantsov(Pc)}; // solutal Ivantsov function

        double k0{A.CsAtT(T0)/A.ClAtT(T0)}; // equilibrium partition coefficient
        double k{(k0+(A.a0*V/D)) / (1+(A.a0*V/D))}; // non equilibrium partition coefficient
        double Ci{C0/(1-(1-k)*Ivc)}; // interface solute concentration
        double m{A.mlAtC(Ci)}; // liquidus gradient at interface
        
        //! derivation assumes linear liquidus and solidus which is far from true here
        double mP{m * ( 1 + (k0-k*(1-std::log(k/k0))) / (1-k0) )}; // velocity dependent liquidus slope (m prime)
        double R0{8.314}; // gas constant
        double mu{A.L*A.V0/(R0*T0*T0)}; // interfacial kinetic coefficient
        double xit{1 - 1/std::sqrt(1 + 1/(A.o*Pt*Pt))}; // thermal stability function
        double xic{1 + 2*k/( 1-2*k-std::sqrt(1 + 1/(A.o*Pc*Pc)) )}; // solutal stability function

        double dTc{T0-A.TlAtC(Ci)}, dTr{2*A.r/R}, dTk{V/mu}; // undercooling components
        double f1{dTt+dTc+dTr+dTk-dT}; // undercooling error
        //! derivation assumes T(C) = Tm + m*C, which is far from true here. 1998 paper also misplaces xic term 
        double f2{(A.r/A.o) / (xit*Pt*A.L/A.Cp - 2*mP*(1-k)*Pc*xic*Ci) - R}; // radius error
        return std::make_tuple(f1, f2, DTs{dTt, dTc, dTr, dTk});
    }

    // the functions below must be seperate from the WLCYZ function so they can be differentiated within that function.

    // calculates the velocity dependent partition coefficient used in the WLCYZ model
    inline double getkv(double T, double V, const alloys::Alloy& A)
    {
        double ke{A.CsAtT(T)/A.ClAtT(T)}; // equilibrium partition coefficient
        double psi{1 - (V*V)/(A.Vd*A.Vd)}; // diffusion coefficient ψ
        double Vdi{A.D/A.a0}; // maximum speed at interface for diffusion
        return (V<A.Vd) ? ((V/Vdi)+ke*psi) / ((V/Vdi)+psi) : 1;
    }

    // calculates the relaxation effect term N used in the WLCYZ model
    inline double getN(double T, double V, const alloys::Alloy& A)
    {
        double ke{A.CsAtT(T)/A.ClAtT(T)}; // equilibrium partition coefficient
        double psi{1 - (V*V)/(A.Vd*A.Vd)}; // diffusion coefficient ψ
        double Vdi{A.D/A.a0}; // maximum speed at interface for diffusion
        double kv{(V<A.Vd) ? ((V/Vdi)+ke*psi) / ((V/Vdi)+psi) : 1};
        return 1 - kv + std::log(kv/ke) + (1-kv)*(1-kv)*V/A.Vd;
    }

    /// @brief Wang, Liu, Chen, Yang and Zhou  model. Generalises better to higher undercoolings and velocities for non
    /// linear phase diagrams.
    /// @param V velocity - m/s
    /// @param R dendrite tip radius - m
    /// @param dT undercooling - K
    /// @param C0 bulk alloy solute concentration - C.%
    /// @param A struct containing key physical alloy parameters
    /// @return dT and R errors. If V, R, dt, and C0 are perfectly correct, both errors should be zero.
    inline std::tuple<double, double, DTs> WLCYZ(double V, double R, double dT, double C0, const alloys::Alloy& A)
    {
        if (!A.WLCYZCapable)
            throw std::runtime_error("Attempted to pass non WLCYZ capable Alloy to WLCYZ model");

        double Pt{V*R/(2*A.a)}; // thermal Péclet number
        double Ivt{ivantsov(Pt)}; // thermal Ivantsov function
        double dTt{A.L*Ivt/A.Cp}; // thermal undercooling
        //! limit dTt here so it is never greater than dT?
        double Ti{A.TlAtC(C0) - dT + ((dTt<dT) ? dTt : dT)}; // interface temperature. Ti must <= Tl(C0)
        
        double Cle{A.ClAtT(Ti)}, Cse{A.CsAtT(Ti)}; // equilibrium interface solute concentration of liquid & solid
        double ke{Cse/Cle}; // equilibrium partition coefficient
        double dTr{2*A.r/R}; // curvature undercooling
        // P suffix (prime) used to denote a value is curvature adjusted (calculated at T=Ti+dTr)
        double CleP{A.ClAtT(Ti+dTr)}, CseP{A.CsAtT(Ti+dTr)}; // curvature adjusted Cle & Cse
        double keP{CseP/CleP}; // curvature adjusted equilibrium partition coefficient
        
        double psi{1 - (V*V)/(A.Vd*A.Vd)}; // diffusion coefficient ψ
        double Vdi{A.D/A.a0}; // maximum speed at interface for diffusion
        //! is this ever needed?
        double kvP{(V<A.Vd) ? ((V/Vdi)+keP*psi) / ((V/Vdi)+psi) : 1}; // curvature corrected velocity dependent k
        double NP{1 - kvP + std::log(kvP/keP) + (1-kvP)*(1-kvP)*V/A.Vd}; // curvature corrected relaxation term N
        //* can this not just be calculated with dimensional analysis like in all the other models?
        double Cl{(CleP-CseP-(V/A.V0))/NP}; // true interface solute concentration of liquid

        double dTc{A.TlAtC(C0) - A.TlAtC(Cl)}; // constitutional (solutal) undercooling
        double dTk{A.TlAtC(Cl) - A.TlAtC(CleP)}; // kinetic undercooling

        alloys::Alloy ACopy1{A}; // required as __enzye_autodiff sometimes modifies objects passed to it
        //* could not copy assign A to a static object as enyzme would fail to deduce the static object's type?
        double dNdT{  // dN(Ti)/dT
            __enzyme_autodiff<double>((void*)getN, enzyme_out, Ti, enzyme_const, V, enzyme_const, &ACopy1)
        };
        double kv{(V<A.Vd) ? ((V/Vdi)+ke*psi) / ((V/Vdi)+psi) : 1}; // velocity dependent partition coefficient
        double N{1 - kv + std::log(kv/ke) + (1-kv)*(1-kv)*V/A.Vd}; // relaxation term N
        double ml{A.mlAtT(Ti)}, ms{A.msAtT(Ti)}; // solidus and liquidus gradients
        double M{-ml*ms*N/(ml-ms+ml*ms*Cl*dNdT)}; // solutal field gradient coefficient

        double mlP{A.mlAtT(Ti+dTr)}, msP{A.msAtT(Ti+dTr)}; // curvature adjusted solidus and liquidus gradients
        alloys::Alloy ACopy2{A};
        double dNPdT{  // dN(Ti+dTr)/dT
            __enzyme_autodiff<double>((void*)getN, enzyme_out, Ti+dTr, enzyme_const, V, enzyme_const, &ACopy2)
        };
        double MP{-mlP*msP*NP/(mlP-msP+mlP*msP*Cl*dNPdT)}; // curvature adjusted solutal field gradient coefficient
        alloys::Alloy ACopy3{A};
        double dkvPdT{ // dKv(Ti+dTr)/dT
            __enzyme_autodiff<double>((void*)getkv, enzyme_out, Ti+dTr, enzyme_const, V, enzyme_const, &ACopy3)
        };
        double Pc{V*R/(2*A.D)}; // solutal peclet number

        double xic{ // solutal stability function
            (V<A.Vd) ? 1 - (2*kvP+2*MP*Cl*dkvPdT) / ( std::sqrt(1+(psi/(A.o*Pc*Pc))) + 2*kvP - 1 + 2*MP*Cl*dkvPdT) : 0
        };
        double xiL{1 - 1/std::sqrt(1 + 1/(A.o*Pt*Pt))}; // thermal stability function
        double RPred{(A.r/A.o) / (Pt*A.L*xiL/A.Cp + 2*MP*Pc*Cl*(kvP-1)*xic/psi)}; // calculated dendrite radius
     
        return std::make_tuple(RPred-R, dTt+dTc+dTr+dTk-dT, DTs{dTt, dTc, dTr, dTk});
    }

}
#endif
