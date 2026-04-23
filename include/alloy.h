#ifndef ALLOY_H
#define ALLOY_H

#include <iostream>
#include <numbers>
#include <array>
#include <cmath>


/// @brief datastructures needed to track alloy physical constants in SI units
namespace alloy
{
    /// @brief contains key physical constants for a given alloy system in SI units
    struct Alloy
    {
        // LGK parameters
        double L{};     // latent heat of fusion - J/kg
        double Cp{};    // specific heat capacity - J/(Kg K)
        double m{};     // equilibrium liquidus slope - K/wt%
        double k0{};    // partition coefficient - unitless
        double r{};     // Gibbs-Thomson coefficient - K m
        double D{};     // solute diffusion coefficient - m2/s
        double a{};     // thermal diffusivity in liquid - m2/s
        double o{};     // stability constant - unitless

        // LKT-BCT parameters
        double a0{};    // atomic spacing in pure metal - m
        double V0{};    // speed of sound in liquid - m/s
        double Tm{};    // melting point of pure metal - K
        //TODO figure out whether to make second alloy for LKT-BCT or use boolean flag member
    };
    
    /// @brief extends Alloy by adjusting diffusivity parameters with C0 and dT
    class AlloyTDependant: public Alloy
    {
        public:
            AlloyTDependant(const Alloy& A, double TmPure, double DA0, double DEa);
            void updateDiffusivity(double dT, double C0);

        private:
            double m_TmPure{}; // melting point of alloy with zero solute concentration - K
            double m_DA0{}; // Arrhenius constant of diffusivity - m2/s
            double m_DEa{}; // activation energy for diffusion - J/mol
    };
    
    using std::numbers::pi;
    // Taken from ThermoCalc TCSLD 4.1 database as in https://doi.org/10.1007/s10854-025-14979-6
    const Alloy SnAg{61'810.62, 249.0, -3.14, 0.0191, 8.54e-8, 1.82e-9, 1.5e-5, 1/(4*pi*pi), 3.07e-10, 2.47e3, 505.1};

    // Tm taking from https://periodic-table.rsc.org/element/50/tin and diffusion constants taken from
    // https://doi.org/10.1063/1.1708821 using slower c axis. //! note these seem to high compared to TermoCalc numbers
    inline AlloyTDependant SnAgTDependant{SnAg, 505.1, 7.1e-7, 12300.0};
    
    // Succinonitrile Acetone mixture taken from https://doi.org/10.1016/0025-5416(84)90199-X. This polymer system is
    // often used in place of a molten alloy to test solidification models more easily in the lab. Note the equilibrium
    // liquidus scope coversion from its K/mol% value in the paper to K/wt% only holds for small wt% values.
    const Alloy SucAce{46'260, 1937.5, -297.855, 0.103, 6.62e-8, 1.27e-9, 1.14e-7, 1/(4*pi*pi)};
}

/// @brief Configure parameters needed to calculate diffusivity at any T using D = A0*exp(-Ea/RT)
/// @param A alloy base object
/// @param TmPure melting point of alloy with zero solute concentration - K
/// @param DA0 Arrhenius constant of diffusivity - m2/s
/// @param DEa activation energy for diffusion - J/mol
inline alloy::AlloyTDependant::AlloyTDependant(const alloy::Alloy& A, double TmPure, double DA0, double DEa)
    :Alloy{A}, m_TmPure{TmPure}, m_DA0{DA0}, m_DEa{DEa} {};


/// @brief Update alloy diffusivity by calculating T at liquid interface and using Arrhenius relationship
/// @param dT undercooling - K
/// @param C0 bulk alloy solute concentration - wt.%
inline void alloy::AlloyTDependant::updateDiffusivity(double dT, double C0)
{
    double T{m_TmPure + m*C0 - dT}; // note undercooling is how much colder liquid at interface is
    const double R{8.31446262}; // molar gas constant
    D = m_DA0 * std::exp(-m_DEa / (R*T));
}

inline std::ostream& operator<<(std::ostream& out, const alloy::Alloy& alloy)
{
    return out << "Alloy(" <<
    "latent heat of fusion=" << alloy.L <<
    ", specific heat capacity=" << alloy.Cp <<
    ", equilibrium liquidus slope=" << alloy.m <<
    ", partition coefficient=" << alloy.k0 <<
    ", Gibbs-Thomson coefficient=" << alloy.r <<
    ", solute diffusion coefficient=" << alloy.D <<
    ", thermal diffusivity in liquid=" << alloy.a <<
    ", stability constant=" << alloy.o <<
    ')';
}


#endif
