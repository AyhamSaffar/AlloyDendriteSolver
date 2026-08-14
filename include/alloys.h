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
    /// be atom percent or weight percent if used with LGK or must be in atom percent with any other model.
    class Alloy
    {
        public:
            // LGK parameters
            double L{};     // Fusion enthalpy - J/mol
            double Cp{};    // Melt heat capacity - J/(mol K)
            double m{};     // Equilibrium liquidus slope - K/C%
            double k0{};    // Equilibrium partition coefficient - C% / C%
            double r{};     // Gibbs-Thomson coefficient - K m
            double D{};     // Diffusion coefficient of solute in liquid - m2/s
            double a{};     // Thermal diffusivity of liquid - m2/s
            double o{};     // Stability constant - m/m

            // LKT-BCT parameters
            bool LKT_BCTCapable{false};
            double a0{};    // Solid atomic spacing - m
            double V0{};    // speed of sound in liquid - m/s
            double Tm{};    // Pure solid melting point - K

            // CLW parameters
            bool CLWCapable{false};
            inline double DAtT(double T) const;
            inline double TlAtC(double C) const;
            inline double mlAtC(double C) const;
            inline double ClAtT(double T) const;
            inline double CsAtT(double T) const;

            // WLCYZ parameters
            bool WLCYZCapable{false};
            double Vd{};    // maximum speed of diffusion in the bulk liquid - m/s
            inline double mlAtT(double T) const;
            inline double msAtT(double T) const;

            inline Alloy(
                double L, double Cp, double m, double k0, double r, double D, double a, double o, // LGK
                double a0=-1, double V0=-1, double Tm=-1, // LKT-BCT
                double DA0=-1, double DEa=-1, std::vector<double> TlAtCFit={}, std::vector<double> ClAtTFit={},
                    std::vector<double> CsAtTFit={}, // CLW
                double Vd=-1 // WLCYZ
            );
            Alloy() = default;
            bool operator<=>(const Alloy&) const = default;

        private:
            double m_DA0{}; // Arrhenius constant of diffusivity - m2/s
            double m_DEa{}; // activation energy for diffusion - J/mol
            std::vector<double> m_TlAtCFit{}; // polynomial fit of Tl for a given C (0th to nth order coefficient)
            std::vector<double> m_ClAtTFit{}; // polynomial fit of Cl for a given T (0th to nth order coefficient)
            std::vector<double> m_CsAtTFit{}; // polynomial fit of Cs for a given T (0th to nth order coefficient)
    }; 
}

/// @brief create a new Alloy object.
/// @param L Fusion enthalpy - J/mol
/// @param Cp Melt heat capacity - J/(mol K)
/// @param m Equilibrium liquidus slope - K/C%
/// @param k0 Equilibrium partition coefficient - C% / C%
/// @param r Gibbs-Thomson coefficient - K m
/// @param D Diffusion coefficient of solute in liquid - m2/s
/// @param a Thermal diffusivity of liquid - m2/s
/// @param o Stability constant - m/m
/// @param a0 Solid atomic spacing - m. Only needed for LKT-BCT capable alloys.
/// @param V0 speed of sound in liquid - m/s. Only needed for LKT-BCT capable alloys.
/// @param Tm Pure solid melting point - K. Only needed for LKT-BCT capable alloys.
/// @param DA0 Arrhenius constant of diffusivity - m2/s. Only needed for CLW alloys.
/// @param DEa activation energy for diffusion - J/mol. Only needed for CLW alloys.
/// @param TlAtCFit polynomial fit of Tl for a given C (0th to nth order coefficient). Only needed for CLW alloys.
/// @param mlAtCFit polynomial fit of m for a given C (0th to nth order coefficient). Only needed for CLW alloys.
/// @param ClAtTFit polynomial fit of Cl for a given T (0th to nth order coefficient). Only needed for CLW alloys.
/// @param CsAtTFit polynomial fit of Cs for a given T (0th to nth order coefficient). Only needed for CLW alloys.
/// @param Vd maximum speed of diffusion in the bulk liquid - m/s. Only needed for WLCYZ alloys.
inline alloys::Alloy::Alloy(
    double L, double Cp, double m, double k0, double r, double D, double a, double o, double a0, double V0, double Tm,
    double DA0, double DEa, std::vector<double> TlAtCFit, std::vector<double> ClAtTFit, std::vector<double> CsAtTFit,
    double Vd  
): L{L}, Cp{Cp}, m{m}, k0{k0}, r{r}, D{D}, a{a}, o{o}, a0{a0}, V0{V0}, Tm{Tm}, m_DA0{DA0}, m_DEa{DEa},
    m_TlAtCFit{TlAtCFit}, m_ClAtTFit{ClAtTFit}, m_CsAtTFit{CsAtTFit}, Vd{Vd}
{
    if ((a0!=-1) && (V0!=-1) && (Tm!=-1))
        LKT_BCTCapable = true;
    if (LKT_BCTCapable && (DA0!=-1) && (DEa!=-1) && (!TlAtCFit.empty()) && (!ClAtTFit.empty()) && (!CsAtTFit.empty()))
        CLWCapable = true;
    if (LKT_BCTCapable && (!TlAtCFit.empty()) && (!ClAtTFit.empty()) && (!CsAtTFit.empty()) && (Vd!=-1))
        WLCYZCapable = true;
}                  

/// @brief Calculates the solute liquid diffusivity at a given T
/// @param T Liquid temperature - K
/// @return Diffusion coefficient of solute in liquid - m2/s
inline double alloys::Alloy::DAtT(double T) const
{
    constexpr double R{8.3145}; // gas constant in J/molK
    return m_DA0*std::exp(-m_DEa/(R*T));
}

/// @brief Calculates the equilibrium liquidus temperature at a given C
/// @param C Liquid solute concentration - C%
/// @return Equilibrium liquidus temperature - K
inline double alloys::Alloy::TlAtC(double C) const
{
    double Tl{0};
    for (std::size_t i{0}; i<std::size(m_TlAtCFit); ++i)
        Tl += m_TlAtCFit[i] * std::pow(C, i);
    return Tl;
}

/// @brief Calculates the equilibrium liquidus slope at a given C
/// @param C Liquid solute concentration - C%
/// @return Equilibrium liquidus slope - K/C%
inline double alloys::Alloy::mlAtC(double C) const
{
    double ml{0};
    for (std::size_t i{1}; i<std::size(m_TlAtCFit); ++i) // i must start at 1 as otherwise uint{0}-1 gives underflow
        ml += i * m_TlAtCFit[i] * std::pow(C, i-1); // ml(C) = dTl(C)/dC
    return ml;
}

/// @brief Calculates the equilibrium liquidus concentration at a given T
/// @param T Liquid temperature - K
/// @return Equilibrium liquidus concentration - C% / C%
inline double alloys::Alloy::ClAtT(double T) const
{
    double Cl{0};
    for (std::size_t i{0}; i<std::size(m_ClAtTFit); ++i)
        Cl += m_ClAtTFit[i] * std::pow(T, i);
    return Cl;
}

/// @brief Calculates the equilibrium solidus concentration at a given T
/// @param T Liquid temperature - K
/// @return Equilibrium solidus concentration - C% / C%
inline double alloys::Alloy::CsAtT(double T) const
{
    double Cs{0};
    for (std::size_t i{0}; i<std::size(m_CsAtTFit); ++i)
        Cs += m_CsAtTFit[i] * std::pow(T, i);
    return Cs;
}

/// @brief Calculates the equilibrium liquidus gradient at a given T
/// @param T Liquid temperature - K
/// @return Equilibrium liquidus gradient - K/C%
inline double alloys::Alloy::mlAtT(double T) const
{
    double dCldT{0};
    for (std::size_t i{1}; i<std::size(m_ClAtTFit); ++i) // i must start at 1 as otherwise uint{0}-1 gives underflow
        dCldT += i * m_ClAtTFit[i] * std::pow(T, i-1);
    return 1/dCldT; // ml = dTl/dC = 1 / (dCl/dT)
}

/// @brief Calculates the equilibrium solidus gradient at a given T
/// @param T Liquid temperature - K
/// @return Equilibrium solidus gradient - K/C%
inline double alloys::Alloy::msAtT(double T) const
{
    double dCsdT{0};
    for (std::size_t i{1}; i<std::size(m_CsAtTFit); ++i) // i must start at 1 as otherwise uint{0}-1 gives underflow
        dCsdT += i * m_CsAtTFit[i] * std::pow(T, i-1);
    return 1/dCsdT; // ml = dTs/dC = 1 / (dCs/dT)
}


// bank of known alloy systems

namespace alloys
{
    // standard solution to marginal stability criterion for a planar interace. Could vary with crystal structure.
    static constexpr double o{1.0/(4*std::numbers::pi*std::numbers::pi)};
    static constexpr double R{8.3145}; // gas constant in J/molK
    static constexpr double NA{0}; //* used E.G. for m and k0 in alloys with non linear phase diagrams.

    // calculated using least squares fitting of the phase diagrams from the TCNI8 database. Valid down to 1343K.
    static std::vector<double> NiBTlAtCFit{1728.310506759,-13.74052637939,-0.2444339920279,-0.01972205626106};
    static std::vector<double> NiBClAtTFit{280.856238051,-0.537495568016,0.000381791884713,-9.53611701497e-08};
    static std::vector<double> NiBCsAtTFit{
        -23.86104359553,0.0647957147612,-6.61489819253e-05,3.083066146706e-08,-5.57017810813e-12
    };

    // Nickel Borom system in at.%. Taken from https://doi.org/10.1016/j.actamat.2006.08.042. The constant m & k0 values
    // are taken from https://www.sciencedirect.com/science/article/pii/S1359646207003302.
    const Alloy NiB2007_atp{
        1.72e4, 36.39, -14.3, 0.0155, 3.42e-7, 3e-9, 8.5e-6, o, (3e-9)/18.9, 425, 1728, -1, -1,
        NiBTlAtCFit, NiBClAtTFit, NiBCsAtTFit, 18.9
    };

    // Derived from Binary Alloy Phase Diagrams, Vol. 1104, American Society for Metals, Metals Park, OH, 1986.
    static std::vector<double> FeSbTlAtCFit{1811.15, -5.7341, -0.02588, -0.00103};
    static std::vector<double> FeSbClAtTFit{393.39, -0.6802, 4.8758e-4, -1.2803e-7};
    static std::vector<double> FeSbCsAtTFit{33.559, -0.01853};

    // Iron Antimony system in wt.%. Take from https://doi.org/10.1080/09500830903002356.
    const Alloy FeSb_wtp{15'027, 43.77, NA, NA, 3.56e-7, NA, 7.3e-6, o, 2.5e-10, 3000, NA, 4.11e-7, 5e4, 
        FeSbTlAtCFit, FeSbClAtTFit, FeSbCsAtTFit};

    /// ThermoCalc 2026b with the TCBIN v1.1 database was used to get phase diagram data and this was fit to 7th order
    // polynomials using the least squares method in Python's Numpy library. This fit is only valid above 1385K (316K
    // dT), as below this T, the phase diagram transitions from an FCC Cu & Liquid mix to an FCC Cu & FCC Co mix. 
    static std::vector<double> CoCuTlAtCFit{
        1768.4309105217317, -2.4830585090005366, -0.06569466381811043, -0.002782299761566121, 0.00024348143492583415,
        -5.255801655884799e-06, 4.880250303891831e-08, -1.7583823392354547e-10
    };
    static std::vector<double> CoCuClAtTFit{
        779838079.9241728, -3507066.1641644835, 6748.775967182457, -7.203433488153699, 0.00460573309595906,
        -1.7639643327132657e-06, 3.746931408209084e-10, -3.4051708146153894e-14
    };
    static std::vector<double> CoCuCsAtTFit{
        -35635062.73621782, 163520.05663544996, -321.2683140293004, 0.35032250139780624, -0.0002289764005790414,
        8.970805214138605e-08, -1.9505707064548234e-11, 1.8158181787312304e-15
    };
    // below paper uses noticably different a0 and a for C0=20wt.% and C0=60wt.%, so these must be 2 seperate Alloys.

    // Cobalt Copper system for 20wt.% Cu. Taken from https://doi.org/10.1007/s11433-010-4167-y. Phase diagram fits used
    // in paper was not used as it lacked decimal places in coefficients as well as a Tl(C) fit.
    const Alloy CoCu_20wtp{15033, 39.06, NA, NA, 3.4e-7, NA, 1.424e-5, o, 1.697e-10, 4000, NA, 1.58e-7, 55060,
        CoCuTlAtCFit, CoCuClAtTFit, CoCuCsAtTFit};

    // Cobalt Copper system for 60wt.% Cu. Taken from https://doi.org/10.1007/s11433-010-4167-y. Phase diagram fits used
    // in paper was not used as it lacked decimal places in coefficients as well as a Tl(C) fit.
    const Alloy CoCu_60wtp{14057, 36.05, NA, NA, 3.33e-7, NA, 2.9e-5, o, 4.294e-10, 4000, NA, 2.04e-7, 54069,
        CoCuTlAtCFit, CoCuClAtTFit, CoCuCsAtTFit};

    static constexpr double NiAr{58.693e-3}, NiDensity{8.907e3}; // In Kg/mol and Kg/m3 respectively
    static constexpr double NiS{(1.72e4*NiDensity/NiAr)/1726}; // S = L/Tm and converted to J/m3K. Not supplied in paper
    // Nickel Boron system in at.%. Taken from https://doi.org/10.1103/PhysRevB.45.5019.
    const Alloy NiB1991_atp{1.72e4, 36.39, -14.3, 8e-6, 0.464/NiS, 2.42e-9, 1e-5, o, (2.42e-9)/7.6, 2e3, 1726};

    static constexpr double FeMeltDensity{7352.53}, FeAr{55.845e-3}; // Ar in Kg/mol
    // Iron Cobalt system in both wt.% and at.%, as Fe and Co have such similar atomic masses (55.845 & 58.993). Taken
    // from https://doi.org/10.1016/j.actamat.2016.09.047. Gamma and Delta refer to different crystal phases that
    // form during solidification. In order to model a non-linear phase diagram with a model that assumes a linear phase
    // diagram (LKT-BCT), the alloy has different values for m and k0 within each concentration range. Also uses C
    // dependent Tm.
    const Alloy FeCoGamma_30atp{
        14098, 5749190*FeAr/7612, -0.69, 0.977, 0.319/1020485, 4.7e-9, 5.46e-06, o, 2.358e-10, 550, 1763
    };
    const Alloy FeCoDelta_30atp{
        10999, 5712909*FeAr/7242, -1.99, 0.949, 0.206/801219, 4.7e-9, 5.46e-06, o, 2.358e-10, 350, 1753
    };
    const Alloy FeCoGamma_40atp{ // m value in paper missing minus sign. See reference [1] in paper.
        14083, 5796451*FeAr/7718, -0.45, 0.989, 0.319/1032396, 4.7e-9, 5.36e-06, o, 2.354e-10, 550, 1757
    };
    const Alloy FeCoDelta_40atp{
        10767, 5704510*FeAr/7352.53, -1.98, 0.96, 0.206/801030, 4.7e-9, 5.36e-06, o, 2.354e-10, 350, 1733
    };
    const Alloy FeCoGamma_50atp{ // m value in paper also missing minus sign. See reference [1] in paper.
        14154, 5822432*FeAr/7824, -0.13, 0.997, 0.319/1043547, 4.7e-9, 5.29e-06, o, 2.35e-10, 550, 1754
    };
    const Alloy FeCoDelta_50atp{
        10795, 5666976*FeAr/7423.47, -1.85, 0.969, 0.206/815419, 4.7e-9, 5.29e-06, o, 2.35e-10, 350, 1714
    };

    // Nickel Tin system in wt.%. Taken from https://doi.org/10.1007/BF02646933 (m and k0 taken from Appendix A liquidus
    // and solidus fits. To remove the T^2 term, a least squares linear trend was fit to the liquidus. k_0 was set to
    // the average calculated values over the 1510K (0K undercooling) to 1210K (300K undercooling) range, as this is
    // where most results are.
    const Alloy NiSn_wtp{1.5e5*NiAr, 500*NiAr, -16.2, 0.61, 0.25/8.4e5, 5e-9, 5e-6, o};

    static constexpr double AlMeltDensity{2375}, AlAr{26.982e-3}; // Ar in Kg/mol
    // Aluminum Iron system in wt.%. Taken from https://doi.org/10.1007/BF02643853
    const Alloy AlFe_wtp{971e6*AlAr/AlMeltDensity, 2.67e6*AlAr/AlMeltDensity, -3.7, 0.038, 1e-7, 2e-9, 0.34e-4, o};

    // Silver Copper system in at.%. Taken from Fourth Conference on Rapid Solidification Processing: Principles and
    // Technologies, Application of dendritic growth theory to the interpretation of rapid solidification
    // microstructures, pages 13-25, W.J. Boettinger, S.R. Coriell and R. Trivedi. Converted to at.% by converting
    // wt.% line to at.% at each point and fitting the result.
    const Alloy AgCu_atp{11'900, 31.8, -11.349, 0.35, 1.53e-7, 2.1e-9, 6.6e-5, o, 1.05e-9, 2e3, 1234};

    // Silver Copper system in wt.%. Taken from Fourth Conference on Rapid Solidification Processing: Principles and
    // Technologies, Application of dendritic growth theory to the interpretation of rapid solidification
    // microstructures, pages 13-25, W.J. Boettinger, S.R. Coriell and R. Trivedi. 
    const Alloy AgCu_wtp{11'900, 31.8, -8.0, 0.37, 1.53e-7, 2.1e-9, 6.6e-5, o, 1.05e-9, 2e3, 1234};

    static constexpr double SnAr{0.11871}; // Ar in Kg/mol 
    // Tin Silver system in wt.%. Taken from ThermoCalc TCSLD 4.1 database as in
    // https://doi.org/10.1007/s10854-025-14979-6.
    const Alloy SnAg_wtp{61'810.62*SnAr, 249.0*SnAr, -3.14, 0.0191, 8.54e-8, 1.82e-9, 1.5e-5, o, 3.07e-10, 2470, 505.1};
    
    static constexpr double SucMr{80.090e-3}; // relative molecular mass of succinonitrile in Kg/mol
    // Succinonitrile Acetone system in at.%. Taken from https://doi.org/10.1016/0025-5416(84)90199-X. This polymer
    // system is often used in place of a molten alloy to test solidification models more easily in the lab. Conversion
    // of L & Cp from /kg to /mol assume a small fraction of solute in succinonitrile. 
    const Alloy SucAce_atp{46'260*SucMr, 1937.5*SucMr, -2.16, 0.103, 6.62e-8, 1.27e-9, 1.14e-7, o};
}

#endif
