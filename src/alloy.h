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
    
    //TODO should reference where these come from!
    constexpr Alloy SnAg{61'810.62, 249.0, -3.14, 0.0191, 8.54e-8, 1.82e-9, 1.5e-5, 1/(4*pi*pi)};
}

std::ostream& operator<<(std::ostream& out, const alloy::Alloy& alloy)
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
