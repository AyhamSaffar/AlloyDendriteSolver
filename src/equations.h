#ifndef EQUATIONS_H
#define EQUATIONS_H


/// @brief Lipton Glicksman Kurz (LGK) model equations that analytically predict how solidification dendrites grow
/// into a molten liquid when the interface is in equilibrium. https://doi.org/10.1016/0025-5416(84)90199-X
namespace equation
{
    /// @brief calculate thermal Peclet number. This is the ratio between thermal transport from fluid movement
    /// (advective) and thermal transport from the temperature gradient (diffusive).
    /// @param V dendrite tip velocity - m/s
    /// @param R dendrite tip radius - m
    /// @param a thermal diffusivity in liquid - m2/s
    /// @return thermal Peclet number - unitless
    inline double getPt(double V, double R, double a)
    {
        return V*R/(2*a);
    }

    /// @brief calculate solutal Peclet number. This is the ratio between mass transport from fluid movement (advective)
    /// and mass transport from the concentration gradient (diffusive).
    /// @param V dendrite tip velocity - m/s
    /// @param R dendrite tip radius - m
    /// @param D solute diffusion coefficient - m2/s
    /// @return solutal peclet number - unitless
    inline double getPc(double V, double R, double D)
    {
        return V*R/(2*D);
    }

    /// @brief calculate thermal Ivantsov function. This equation comes from solving the rate of thermal transport
    /// across the surface of a 3D parabola shaped dendrite.
    /// @param Pt thermal Peclet number - unitless
    /// @param E1Pt 
    /// @return 
    inline double getIvPt(double Pt, double E1Pt);
    inline double getIvPc(double Pc, double E1Pc);
    inline double getF1(double L, double Cp, double IvPt, double m, double C0, double k0, double IvPc, double r,
                        double R, double dT);
    inline double getQ(double Pt, double L, double Cp, double Pc, double m, double C0, double k0,double IvPc);
    inline double getF2(double r, double o, double Q, double R);
}

#endif
