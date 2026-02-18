#ifndef APPROXIMATORS_H
#define APPROXIMATORS_H


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
    double getTipRaius(double r, double m, double k0, double C0, double dT);

    /// @brief given as equation 8.92 in cited textbook.
    /// @param D solute diffusion coefficient - m2/s
    /// @param m equilibrium liquidus slope - K/wt%
    /// @param k0 partition coefficient - unitless
    /// @param r Gibbs-Thomson coefficient - K m
    /// @param dT undercooling - K
    /// @param C0 bulk alloy composition - wt%
    /// @return approximate tip velocity - m/s
    double getTipVelocity(double D, double m, double k0, double r, double dT, double C0);
}

#endif
