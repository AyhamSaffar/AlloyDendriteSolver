#ifndef ALLOY_H
#define ALLOY_H

#include <numbers>


namespace alloy
{
    using namespace std::numbers; // for value of pi

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

    // should reference where these come from!
    Alloy SnAg{61'810.62, 249.0, -3.14, 0.0191, 8.54e-8, 1.82e-9, 1.5e-5, 1/(4*pi*pi)};
}

#endif
