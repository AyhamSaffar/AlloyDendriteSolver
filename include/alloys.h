#ifndef ALLOYS_H
#define ALLOYS_H

#include <iostream>
#include <numbers>
#include <array>
#include <cmath>


/// @brief datastructures needed to track alloy physical constants in SI units
namespace alloys
{
    /// @brief contains key physical constants for a given alloy system in SI units
    struct Alloy
    {
        // LGK parameters
        double L{};     // Fusion enthalpy - J/mol
        double Cp{};    // Melt heat capacity - J/(mol K)
        double m{};     // Equilibrium liquidus slope - K/wt.%
        double k0{};    // Partition coefficient - wt.% / wt.%
        double r{};     // Gibbs-Thomson coefficient - K m
        double D{};     // Diffusion coefficient of solute in liquid - m2/s
        double a{};     // Thermal diffusivity of liquid - m2/s
        double o{};     // Stability constant - m/m

        // LKT-BCT parameters 
        // set to nan to ensure the object cannot be used with the LKT-BCT model if values are not given for these
        double a0{NAN};    // Solid atomic spacing - m
        double V0{NAN};    // speed of sound in liquid - m/s
        double Tm{NAN};    // Solid melting point - K
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
    
    // standard solution to marginal stability criterion for a planar interace. //! Could vary with crystal structure.
    static constexpr double o{1/(4*std::numbers::pi*std::numbers::pi)};

    static constexpr double FeMeltDensity{7352.53}, FeAr{55.845e-3}; // Ar in Kg/mol
    // Taken from https://doi.org/10.1016/j.actamat.2016.09.047. Gamma and Delta refer to different crystal phases that
    // form during solidification. The paper lists slighly different parameters for 30, 40, and 50 atom.% Co. The 40
    // atom.% Co parameters are used here as an average of the similar values. Also as Fe and Co have such similar
    // atomic masses (55.845 & 58.993), atom.% is assumed to be equal to wt.%.
    const Alloy FeCoGamma{
        14083, 5796451*FeAr/FeMeltDensity, -0.45, 0.989, 0.319/1032396, 4.7e-9, 5.36e-06, o, 2.354e-10, 550, 1757
    }; //! should use C0 dependant m as is -0.69, 0.45, and 0.13 at 30, 40, and 50 atom.% Co repectively.
    const Alloy FeCoDelta{
        10767, 5704510*FeAr/FeMeltDensity, -1.98, 0.96, 0.206/801030, 4.7e-9, 5.36e-06, o, 2.354e-10, 350, 1733
    };

    static constexpr double NiAr{58.693e-3}; // Ar in Kg/mol
    // Taken from https://doi.org/10.1007/BF02646933 (m and k0 taken from Appendix A liquidus and solidus fits. To
    // remove the T^2 term, a least squares linear trend was fit to the liquidus. k_0 was set to the average calculated
    // values over the 1510K (0K undercooling) to 1210K (300K undercooling) range, as this is where most results are.
    const Alloy NiSn{1.5e5*NiAr, 500*NiAr, -16.2, 0.61, 0.25/8.4e5, 5e-9, 5e-6, o};

    static constexpr double AlMeltDensity{2375}, AlAr{26.982e-3}; // Ar in Kg/mol
    // Taken from https://doi.org/10.1007/BF02643853
    const Alloy AlFe{971e6*AlAr/AlMeltDensity, 2.67e6*AlAr/AlMeltDensity, -3.7, 0.038, 1e-7, 2e-9, 0.34e-4, o};

    // Taken from Fourth Conference on Rapid Solidification Processing: Principles and Technologies, Application of
    // dendritic growth theory to the interpretation of rapid solidification microstructures, pages 13-25, W.J.
    // Boettinger, S.R. Coriell and R. Trivedi.
    const Alloy AgCu{11'900, 31.8, -8.0, 0.37, 1.53e-7, 2.1e-9, 6.6e-5, o, 1.05e-9, 2e3, 1234};

    static constexpr double SnAr{0.11871}; // Ar in Kg/mol 
    // Taken from ThermoCalc TCSLD 4.1 database as in https://doi.org/10.1007/s10854-025-14979-6
    const Alloy SnAg{61'810.62*SnAr, 249.0*SnAr, -3.14, 0.0191, 8.54e-8, 1.82e-9, 1.5e-5, o, 3.07e-10, 2.47e3, 505.1};

    // Tm taking from https://periodic-table.rsc.org/element/50/tin and diffusion constants taken from
    // https://doi.org/10.1063/1.1708821 using slower c axis. //! note these seem to high compared to TermoCalc numbers
    inline AlloyTDependant SnAgTDependant{SnAg, 505.1, 7.1e-7, 12300.0};
    
    static constexpr double SucMr{0.080090}; // relative molecular mass of succinonitrile in Kg/mol
    // Succinonitrile Acetone mixture taken from https://doi.org/10.1016/0025-5416(84)90199-X. This polymer system is
    // often used in place of a molten alloy to test solidification models more easily in the lab. Note the equilibrium
    // liquidus scope coversion from its K/mol% value in the paper to K/wt% only holds for small wt% values.
    const Alloy SucAce{46'260*SucMr, 1937.5*SucMr, -297.855, 0.103, 6.62e-8, 1.27e-9, 1.14e-7, o};
}

/// @brief Configure parameters needed to calculate diffusivity at any T using D = A0*exp(-Ea/RT)
/// @param A alloy base object
/// @param TmPure melting point of alloy with zero solute concentration - K
/// @param DA0 Arrhenius constant of diffusivity - m2/s
/// @param DEa activation energy for diffusion - J/mol
inline alloys::AlloyTDependant::AlloyTDependant(const alloys::Alloy& A, double TmPure, double DA0, double DEa)
    :Alloy{A}, m_TmPure{TmPure}, m_DA0{DA0}, m_DEa{DEa} {};


/// @brief Update alloy diffusivity by calculating T at liquid interface and using Arrhenius relationship
/// @param dT undercooling - K
/// @param C0 bulk alloy solute concentration - wt.%
inline void alloys::AlloyTDependant::updateDiffusivity(double dT, double C0)
{
    double T{m_TmPure + m*C0 - dT}; // note undercooling is how much colder liquid at interface is
    const double R{8.31446262}; // molar gas constant
    D = m_DA0 * std::exp(-m_DEa / (R*T));
}

inline std::ostream& operator<<(std::ostream& out, const alloys::Alloy& alloy)
{
    return out << "Alloy(" <<
    "Fusion enthalpy=" << alloy.L <<
    ", Melt heat capacity=" << alloy.Cp <<
    ", Equilibrium liquidus slope=" << alloy.m <<
    ", Partition coefficient=" << alloy.k0 <<
    ", Gibbs-Thomson coefficient=" << alloy.r <<
    ", Diffusion coefficient of solute in liquid=" << alloy.D <<
    ", Thermal conductivity liquid=" << alloy.a <<
    ", Stability constant=" << alloy.o <<
    ')';
}


#endif
