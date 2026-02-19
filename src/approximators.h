#ifndef APPROXIMATORS_H
#define APPROXIMATORS_H
#include <numbers>
#include <cmath>

// approximate analytical solutions to solidification parameters following "Solidification’ by Dantzig & Rappaz
// (1st Ed)" textbook. assumes small undercooling and solutal dendrites.
namespace approx
{
    /// @brief given as equation 8.91 in cited texbook.
    /// @param r Gibbs-Thomson coefficient - K m
    /// @param m quilibrium liquidus slope - K/wt%
    /// @param k0 partition coefficient - unitless
    /// @param C0 bulk alloy composition - wt%
    /// @param dT undercooling - K
    /// @return approximate tip radius - m
    double getTipRaius(double r, double m, double k0, double C0, double dT)
    {
        using std::numbers::pi;
        using std::pow;
        return 6.64*pi*pi*r * pow(-m*(1-k0), 0.25) * (pow(C0, 0.25)/pow(dT, 1.25));
    }

    /// @brief given as equation 8.92 in cited textbook.
    /// @param D solute diffusion coefficient - m2/s
    /// @param m equilibrium liquidus slope - K/wt%
    /// @param k0 partition coefficient - unitless
    /// @param r Gibbs-Thomson coefficient - K m
    /// @param dT undercooling - K
    /// @param C0 bulk alloy composition - wt%
    /// @return approximate tip velocity - m/s
    double getTipVelocity(double D, double m, double k0, double r, double dT, double C0)
    {
        using std::numbers::pi;
        using std::pow;
        return (pow(dT, 2.5) / pow(C0, 1.5)) * D / (5.51*pi*pi*pow(-m*(1-k0), 1.5)*r);
    }
}

#endif
