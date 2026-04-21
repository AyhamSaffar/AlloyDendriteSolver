#ifndef ALLOY_H
#define ALLOY_H

#include <iostream>
#include <numbers>
#include <array>


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
    
    /// @brief extends Alloy by adjusting thermodynamic parameters with dT
    class AlloyTDependant: public Alloy
    {
        public:
            AlloyTDependant(const Alloy& A, bool mVaries, bool k0Varies, bool Dvaries)
            :Alloy{A}, m_mVaries{mVaries}, m_k0Varies{k0Varies}, m_DVaries{Dvaries} {};
        
            void updateT(double dT, double C0);
            void addLiquidusFit(double m0, double m1, double m2=0, double m3=0, double m4=0, double m5=0);
            void addSolidusFit(double m0, double m1, double m2=0, double m3=0, double m4=0, double m5=0);
            void addDiffusivityFit(double A0, double Ea);
        
        private:
            bool m_mVaries{false}, m_k0Varies{false}, m_DVaries{false};
            bool m_liquidusFitSet{false}, m_solidusFitSet{false}, m_diffusivityFitSet{false};
            std::array<double, 6> m_solidusFitParams{}, m_liquidusFitParams{};
            double m_diffusivityA0{}, m_diffusivityEa{};
    };

    using std::numbers::pi;
    // Taken from ThermoCalc TCSLD 4.1 database
    constexpr Alloy SnAg{61'810.62, 249.0, -3.14, 0.0191, 8.54e-8, 1.82e-9, 1.5e-5, 1/(4*pi*pi)};
    
    // Succinonitrile Acetone mixture taken from https://doi.org/10.1016/0025-5416(84)90199-X. This polymer system is
    // often used in place of a molten alloy to test solidification models more easily in the lab. Note the equilibrium
    // liquidus scope coversion from its K/mol% value in the paper to K/wt% only holds for small wt% values.
    constexpr Alloy SucAce{46'260, 1937.5, -297.855, 0.103, 6.62e-8, 1.27e-9, 1.14e-7, 1/(4*pi*pi)};
}

/// @brief set a polynomial fit that maps solute concentration (C) to liquidus temperature (Tl), where Tl = m0 + m1 * C
/// + m2 * C^2 + m3 * C^3 + m4 * C^4 + m5 * C^5.
inline void alloy::AlloyTDependant::addLiquidusFit(double m0, double m1, double m2, double m3, double m4, double m5)
{
    m_liquidusFitSet = true;
    m_liquidusFitParams = {m0, m1, m2, m3, m4, m5};
}

/// @brief set a polynomial fit that maps solute concentration (C) to solidus temperature (Ts), where Ts = m0 + m1 * C
/// + m2 * C^2 + m3 * C^3 + m4 * C^4 + m5 * C^5.
inline void alloy::AlloyTDependant::addSolidusFit(double m0, double m1, double m2, double m3, double m4, double m5)
{
    m_solidusFitSet = true;
    m_solidusFitParams = {m0, m1, m2, m3, m4, m5};
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
