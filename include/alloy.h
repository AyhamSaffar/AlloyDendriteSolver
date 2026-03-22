#ifndef ALLOY_H
#define ALLOY_H

#include <iostream>
#include <numbers>


/// @brief datastructures needed to track alloy physical constants in SI units
namespace alloy
{
    /// @brief contains key physical constants for a given alloy system in SI units
    struct Alloy
    {
        double L{};     // latent heat of fusion - J/kg
        double Cp{};    // specific heat capacity - J/(Kg K)
        double m{};     // equilibrium liquidus slope - K/wt%
        double k0{};    // partition coefficient - unitless
        double r{};     // Gibbs-Thomson coefficient - K m
        double D{};     // solute diffusion coefficient - m2/s
        double a{};     // thermal diffusivity in liquid - m2/s
        double o{};     // stability constant - unitless
    };
    
    using std::numbers::pi;
    
    // Taken from ThermoCalc TCSLD 4.1 database
    constexpr Alloy SnAg{61'810.62, 249.0, -3.14, 0.0191, 8.54e-8, 1.82e-9, 1.5e-5, 1/(4*pi*pi)};

    // Succinonitrile Acetone mixture taken from https://doi.org/10.1016/0025-5416(84)90199-X. This polymer system is
    // often used in place of a molten alloy to test solidification models more easily in the lab. Note the equilibrium
    // liquidus scope coversion from its K/mol% value in the paper to K/wt% only holds for small wt% values.
    constexpr Alloy SucAce{46'260, 1937.5, -297.855, 0.103, 6.62e-8, 1.27e-9, 1.14e-7, 1/(4*pi*pi)};
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
