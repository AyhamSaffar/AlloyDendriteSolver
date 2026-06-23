#ifndef ALLOYS_H
#define ALLOYS_H

#include <iostream>
#include <numbers>
#include <array>
#include <cmath>


/// @brief datastructures needed to track alloy physical constants in SI units
namespace alloys
{
    /// @brief contains key physical constants for a given alloy system in SI units. Concentration units (C%) can either
    /// be atom percent or weight percent.
    struct Alloy
    {
        // LGK parameters
        double L{};     // Fusion enthalpy - J/mol
        double Cp{};    // Melt heat capacity - J/(mol K)
        double m{};     // Equilibrium liquidus slope - K/C%
        double k0{};    // Partition coefficient - C% / C%
        double r{};     // Gibbs-Thomson coefficient - K m
        double D{};     // Diffusion coefficient of solute in liquid - m2/s
        double a{};     // Thermal diffusivity of liquid - m2/s
        double o{};     // Stability constant - m/m

        // LKT-BCT parameters 
        double a0{};    // Solid atomic spacing - m
        double V0{};    // speed of sound in liquid - m/s
        double Tm{};    // Solid melting point - K

        Alloy(double L, double Cp, double m, double k0, double r, double D, double a, double o,
            double a0=NAN, double V0=NAN, double Tm=NAN
        );
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
    
}

/// @brief create a new Alloy object. The final 3 parameters (a0, V0, and Tm) are only required if the alloy is to be
/// used with the LKT-BCT model and so are optional.
/// @param L Fusion enthalpy - J/mol
/// @param Cp Melt heat capacity - J/(mol K)
/// @param m Equilibrium liquidus slope - K/C%
/// @param k0 Partition coefficient - C% / C%
/// @param r Gibbs-Thomson coefficient - K m
/// @param D Diffusion coefficient of solute in liquid - m2/s
/// @param a Thermal diffusivity of liquid - m2/s
/// @param o Stability constant - m/m
/// @param a0 Solid atomic spacing - m. Only needed for LKT-BCT, so default to NAN.
/// @param V0 speed of sound in liquid - m/s. Only needed for LKT-BCT, so default to NAN.
/// @param Tm Solid melting point - K. Only needed for LKT-BCT, so default to NAN.
inline alloys::Alloy::Alloy(double L, double Cp, double m, double k0, double r, double D, double a, double o,
    double a0, double V0, double Tm
): L{L}, Cp{Cp}, m{m}, k0{k0}, r{r}, D{D}, a{a}, o{o}, a0{a0}, V0{V0}, Tm{Tm} {}


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

// bank of known alloy systems

namespace alloys
{
    // standard solution to marginal stability criterion for a planar interace. Could vary with crystal structure.
    static constexpr double o{1.0/(4*std::numbers::pi*std::numbers::pi)};

    static constexpr double NiAr{58.693e-3}, NiDensity{8.907e3}; // In Kg/mol and Kg/m3 respectively
    static constexpr double NiS{(1.72e4*NiDensity/NiAr)/1726}; // S = L/Tm converted to J/m3K
    // Nickel Boron system in at.%. Taken from https://doi.org/10.1103/PhysRevB.45.5019.
    const Alloy NiB_atp{1.72e4, 36.39, -14.3, 8e-6, 0.464/NiS, 2.42e-9, 1e-5, o, (2.42e-9)/7.6, 2e3, 1726};

    static constexpr double FeMeltDensity{7352.53}, FeAr{55.845e-3}; // Ar in Kg/mol
    // Iron Cobalt system in both wt.% and at.%, as Fe and Co have such similar atomic masses (55.845 & 58.993). Taken
    // from https://doi.org/10.1016/j.actamat.2016.09.047. Gamma and Delta refer to different crystal phases that
    // form during solidification. The paper lists slighly different parameters for 30, 40, and 50 atom.% Co. The 40
    // atom.% Co parameters are used here as an average of the similar values.
    const Alloy FeCoGamma{
        14083, 5796451*FeAr/FeMeltDensity, -0.45, 0.989, 0.319/1032396, 4.7e-9, 5.36e-06, o, 2.354e-10, 550, 1757
    }; //! should use C0 dependant m as is -0.69, -0.45, and -0.13 at 30, 40, and 50 atom.% Co repectively.
    const Alloy FeCoDelta{
        10767, 5704510*FeAr/FeMeltDensity, -1.98, 0.96, 0.206/801030, 4.7e-9, 5.36e-06, o, 2.354e-10, 350, 1733
    };

    // Nickel Tin system in wt.%. Taken from https://doi.org/10.1007/BF02646933 (m and k0 taken from Appendix A liquidus
    // and solidus fits. To remove the T^2 term, a least squares linear trend was fit to the liquidus. k_0 was set to
    // the average calculated values over the 1510K (0K undercooling) to 1210K (300K undercooling) range, as this is
    // where most results are.
    const Alloy NiSn_wtp{1.5e5*NiAr, 500*NiAr, -16.2, 0.61, 0.25/8.4e5, 5e-9, 5e-6, o};

    static constexpr double AlMeltDensity{2375}, AlAr{26.982e-3}; // Ar in Kg/mol
    // Aluminum Iron system in wt.%. Taken from https://doi.org/10.1007/BF02643853
    const Alloy AlFe_wtp{971e6*AlAr/AlMeltDensity, 2.67e6*AlAr/AlMeltDensity, -3.7, 0.038, 1e-7, 2e-9, 0.34e-4, o};

    // Silver Copper system in wt.%. Taken from Fourth Conference on Rapid Solidification Processing: Principles and
    // Technologies, Application of dendritic growth theory to the interpretation of rapid solidification
    // microstructures, pages 13-25, W.J. Boettinger, S.R. Coriell and R. Trivedi. 
    const Alloy AgCu_wtp{11'900, 31.8, -8.0, 0.37, 1.53e-7, 2.1e-9, 6.6e-5, o, 1.05e-9, 2e3, 1234};

    static constexpr double SnAr{0.11871}; // Ar in Kg/mol 
    // Tin Silver system in wt.%. Taken from ThermoCalc TCSLD 4.1 database as in
    // https://doi.org/10.1007/s10854-025-14979-6. //! LKT-BCT alloys must be in at.%
    const Alloy SnAg_wtp{61'810.62*SnAr, 249.0*SnAr, -3.14, 0.0191, 8.54e-8, 1.82e-9, 1.5e-5, o, 3.07e-10, 2.47e3, 505.1};

    // Tm taking from https://periodic-table.rsc.org/element/50/tin and diffusion constants taken from
    // https://doi.org/10.1063/1.1708821 using slower c axis. //! note these seem to high compared to TermoCalc numbers
    inline AlloyTDependant SnAgTDependant{SnAg_wtp, 505.1, 7.1e-7, 12300.0};
    
    static constexpr double SucMr{80.090e-3}, AceMr{58.08e-3}; // relative molecular mass of succinonitrile in Kg/mol
    // Succinonitrile Acetone system in wt.%. Taken from https://doi.org/10.1016/0025-5416(84)90199-X. This polymer
    // system is often used in place of a molten alloy to test solidification models more easily in the lab. The
    // equilibrium liquidus scope coversion from its K/mol% value in the paper to K/wt% only holds for small wt% values.
    // Also the equivalence of k0 from mol.%/mol.% value in the paper to its wt.%/wt.% value also only holds for small
    // wt% values.
    const Alloy SucAce_wtp{46'260*SucMr, 1937.5*SucMr, -2.16*SucMr/AceMr, 0.103, 6.62e-8, 1.27e-9, 1.14e-7, o};
}


#endif
