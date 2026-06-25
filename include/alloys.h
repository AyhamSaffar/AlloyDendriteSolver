#ifndef ALLOYS_H
#define ALLOYS_H

#include <iostream>
#include <numbers>
#include <vector>
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
        double Tm{};    // Pure solid melting point - K

        Alloy(double L, double Cp, double m, double k0, double r, double D, double a, double o,
            double a0=NAN, double V0=NAN, double Tm=NAN
        );
    };
    
    /// @brief extends Alloy by adjusting D and k0 with T as well as m with C0
    class DynamicAlloy: public Alloy
    {
        public:
            DynamicAlloy(const Alloy& alloy,
                double DA0, double DEa, double Tm, std::vector<double> mFit, std::vector<double> k0Fit
            );
            void updateParams(double dT, double C0);

        private:
            double m_DA0{}; // Arrhenius constant of diffusivity - m2/s
            double m_DEa{}; // activation energy for diffusion - J/mol
            double m_Tm{}; // Pure solid melting point - K
            std::vector<double> m_mFit{}; // polynomial fit of m for a given C0 (0th to nth order coefficient)
            std::vector<double> m_k0Fit{}; // polynomial fit of k0 for a given T (0th to nth order coefficient)
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
/// @param Tm Pure solid melting point - K. Only needed for LKT-BCT, so default to NAN.
inline alloys::Alloy::Alloy(double L, double Cp, double m, double k0, double r, double D, double a, double o,
    double a0, double V0, double Tm
): L{L}, Cp{Cp}, m{m}, k0{k0}, r{r}, D{D}, a{a}, o{o}, a0{a0}, V0{V0}, Tm{Tm} {}


/// @brief create an T dependant Alloy object whose m, k0, and D parameters can be updated at each C0 and dT
/// @param alloy container for key physical constants of a given alloy system
/// @param DA0 Arrhenius constant of diffusivity - m2/s
/// @param DEa activation energy for diffusion - J/mol
/// @param Tm Pure solid melting point - K
/// @param liquidusFit polynomial liquidus line fit (0th order coefficient to 5th)
/// @param solidusFit polynomial solidus line fit (0th order coefficient to 5th)
inline alloys::DynamicAlloy::DynamicAlloy(const Alloy& alloy,
    double DA0, double DEa, double Tm, std::vector<double> mFit, std::vector<double> k0Fit
): Alloy{alloy}, m_DA0{DA0}, m_DEa{DEa}, m_Tm{Tm}, m_mFit{mFit}, m_k0Fit{k0Fit} {}

/// @brief updates D and k0 with T as well as m with C0
/// @param dT undercooling - K
/// @param C0 bulk alloy solute concentration - C.%
inline void alloys::DynamicAlloy::updateParams(double dT, double C0)
{
    double T{m_Tm-dT}; // dendrite tip temperature
    constexpr double R{8.3145}; // gas constant in J/molK
    D = m_DA0*std::exp(m_DEa/(R*T));
    m = 0;
    for (std::size_t i{0}; i<=std::size(m_mFit); ++i)
        m += m_mFit[i]*std::pow(C0, i);
    k0 = 0;
    for (std::size_t i{0}; i<=std::size(m_k0Fit); ++i)
        k0 += m_k0Fit[i]*std::pow(T, i);
}


// bank of known alloy systems

namespace alloys
{
    // standard solution to marginal stability criterion for a planar interace. Could vary with crystal structure.
    static constexpr double o{1.0/(4*std::numbers::pi*std::numbers::pi)};
    static constexpr double R{8.3145}; // gas constant in J/molK

    // Cobalt Copper system in wt.%. Taken from https://doi.org/10.1007/s11433-010-4167-y. Values used are averages
    // between the 20wt.% Cu and the 60wt.% Cu where values differ. The default D used is the result of the Arrhenius
    // fit at 100K dT. //! The a & a0 varies significantly between 60wt.% the 20wt.%. LKT-BCT alloys must be in at.%
    static constexpr double DA0{(1.58e-7+2.04e-7)/2}, DEa{(55060.0+54096.0)/2}, Tl{(1701.0+1663.0)/2};
    static inline double D{DA0*std::exp( DEa/(R*(Tl-100)) )}; 
    const Alloy CoCu_wtp{
        (15033.0+14057.0)/2, (39.06+36.05)/2, -3.3, 0.67, (3.4e-7+3.33e-7)/2, D, (1.424e-5+2.9e-5)/2, o,
        (1.697e-10+4.294e-10)/2, 4000.0, Tl
    };
    // Cobalt Copper system in wt.%. Taken from https://doi.org/10.1007/s11433-010-4167-y. Values used are averages
    // between the 20wt.% Cu and the 60wt.% Cu where values differ.
    DynamicAlloy CoCuD_wtp{CoCu_wtp, DA0, DEa, 1768};
    // Nickel Borom system in at.%. Taken from https://doi.org/10.1016/j.actamat.2006.08.042. m and k0 were fit using
    // ThermoCalc database TCNI8. Default m and k0 used were the average over the first 100K dT.
    const Alloy NiB2007_atp{1.72e4, 36.39, -16.3, 3.45e-2, 3.42e-7, 3e-9, 8.5e-6, o, (3e-9)/18.9, 425, 1726};

    static constexpr double NiAr{58.693e-3}, NiDensity{8.907e3}; // In Kg/mol and Kg/m3 respectively
    static constexpr double NiS{(1.72e4*NiDensity/NiAr)/1726}; // S = L/Tm and converted to J/m3K. Not supplied in paper
    // Nickel Boron system in at.%. Taken from https://doi.org/10.1103/PhysRevB.45.5019.
    const Alloy NiB1997_atp{1.72e4, 36.39, -14.3, 8e-6, 0.464/NiS, 2.42e-9, 1e-5, o, (2.42e-9)/7.6, 2e3, 1726};

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
    
    static constexpr double SucMr{80.090e-3}, AceMr{58.08e-3}; // relative molecular mass of succinonitrile in Kg/mol
    // Succinonitrile Acetone system in wt.%. Taken from https://doi.org/10.1016/0025-5416(84)90199-X. This polymer
    // system is often used in place of a molten alloy to test solidification models more easily in the lab. The
    // equilibrium liquidus scope coversion from its K/mol% value in the paper to K/wt% only holds for small wt% values.
    // Also the equivalence of k0 from mol.%/mol.% value in the paper to its wt.%/wt.% value also only holds for small
    // wt% values.
    const Alloy SucAce_wtp{46'260*SucMr, 1937.5*SucMr, -2.16*SucMr/AceMr, 0.103, 6.62e-8, 1.27e-9, 1.14e-7, o};
}


#endif
