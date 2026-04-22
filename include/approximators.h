#ifndef APPROXIMATORS_H
#define APPROXIMATORS_H
#include <numbers>
#include <cmath>
#include "alloy.h"

// approximate analytical solutions to solidification parameters following "Solidification" by Dantzig & Rappaz
// (1st Ed)" textbook. assumes small undercooling and solutal dendrites.
namespace approx
{
    /// @brief given as equation 8.91 in cited texbook.
    /// @param dT undercooling - K
    /// @param C0 bulk alloy composition - wt%
    /// @param A struct containing key physical alloy parameters
    /// @return approximate tip radius - m
    inline double getTipRadius(double dT, double C0, const alloy::Alloy& A)
    {
        using std::numbers::pi;
        using std::pow;
        return 6.64*pi*pi*A.r * pow(-A.m*(1-A.k0), 0.25) * (pow(C0, 0.25)/pow(dT, 1.25));
    }

    /// @brief given as equation 8.92 in cited textbook.
    /// @param dT undercooling - K
    /// @param C0 bulk alloy solute concentration - wt%
    /// @param A struct containing key physical alloy parameters
    /// @return approximate tip velocity - m/s
    inline double getTipVelocity(double dT, double C0, const alloy::Alloy& A)
    {
        using std::numbers::pi;
        using std::pow;
        return (pow(dT, 2.5) / pow(C0, 1.5)) * A.D / (5.51*pi*pi*pow(-A.m*(1-A.k0), 1.5)*A.r);
    }
}

#endif
