#include <numbers>
#include <cmath>
#include <tuple>
#include <string>
#include <catch2/catch_test_macros.hpp>
#include "alloys.h"
#include "models.h"
#include "solvers.h"


// numericallly derived LGK fit from https://doi.org/10.1007/s10854-025-14979-6
double getPublishedLGKSnAgVFit(double dT, double C0)
{
    return std::pow(dT*std::pow(C0, -0.59)/28.6, 1/0.35);
}

TEST_CASE("LGK model V prediction agrees with published LGK SnAg numerical fit and gives positive R", "[solvers]")
{
    for (double dT{2.5}; dT<=60; dT+=2.5)
        for (double C0{3.0}; C0<=6.0; C0+=1.0)
        {
            INFO("dT = " + std::to_string(dT) + ", and C0 = " + std::to_string(C0));
            solvers::Result R{solvers::newton<models::LGK>(dT, C0, alloys::SnAg_wtp)};
            REQUIRE(R.hasConverged);
            REQUIRE(R.R > 0);

            double VFit{getPublishedLGKSnAgVFit(dT, C0)};
            REQUIRE(std::abs(R.V-VFit)/VFit < 0.20); // maximum of 20% error as numerical fit
        }
}


TEST_CASE("LKT-BCT model V prediction agrees with LGK at low undercooling", "[solvers]")
{
    // low undercoolings -> lower V and R -> k does not vary from k0, stability functions equal 1, and negligible
    // kinetic undercooling. Therefore LKT_BCT reduces to LGK 
    for (double C0{1}; C0<=6; C0+=1)
        for (double dT{1}; dT<=200; dT+=1)
        {
            INFO("dT = " + std::to_string(dT) + ", and C0 = " + std::to_string(C0));

            const alloys::Alloy A{alloys::SnAg_wtp};
            constexpr bool legacy{false}; // ensures LGK form is consistent with LKT_BCT
            solvers::Result LGKR{solvers::newton<models::LGK<legacy>>(dT, C0, A)};
            solvers::Result BCTR{solvers::newton<models::LKT_BCT>(dT, C0, A)};

            const solvers::Result& R{BCTR};
            double Pt{R.V*R.R/(2*A.a)}; // thermal Péclet number
            double Pc{R.V*R.R/(2*A.D)}; // solutal Péclet number
            double k{(A.k0+(A.a0*R.V/A.D)) / (1+(A.a0*R.V/A.D))}; // non equilibrium partition coefficient
            double xit{1 - 1/std::sqrt(1 + 1/(A.o*Pt*Pt))}; // thermal stability function
            double xic{1 + 2*k/( 1-2*k-std::sqrt(1 + 1/(A.o*Pc*Pc)) )}; // - solutal stability function
            double dTk{R.V/(A.L*A.V0/(8.314*A.Tm*A.Tm))}; // thermal undercooling

            if((xit<0.9) || (xic<0.9) || ((std::abs(k-A.k0)/A.k0)>0.1) || ((dTk/dT)>0.1))
                break; // LKT-BCT no longer equivalent to LGK

            REQUIRE(BCTR.hasConverged);
            REQUIRE(std::abs(BCTR.R-LGKR.R)/LGKR.R < 0.01); // maximum of 1% error
            REQUIRE(std::abs(BCTR.V-LGKR.V)/LGKR.V < 0.01); // maximum of 1% error
        }
}

TEST_CASE("Linearised CLW model agrees with LKT-BCT at pre solute trapping undercoolings", "[Solvers]")
{
    // CLW starts to diverge from LKT-BCT at higher V as dTk expression slighly different. As such, test only goes up
    // to V = Vd/20.
    const alloys::Alloy A{alloys::AgCu_wtp}; // LKT-BCT capable Alloy
    const alloys::Alloy ALin{ // CLW model with fixed D, k0 and linear Tl
        A.L, A.Cp, A.m, A.k0, A.r, A.D, A.a, A.o, A.a0, A.V0, A.Tm, A.D, 0, 
        {alloys::Fit{{A.Tm, A.m}}}, {alloys::Fit{{1}}}, {alloys::Fit{{A.k0}}}, 
    };
    const double C0{3}, dT0{1}, Vd{A.D/A.a0}; // C0 low as CLW assumes dilute limit in solute trapping expression
    double V0{approx::getV(dT0, C0, A)}, R0{approx::getR(dT0, C0, A)};
    
    for (double dT{1}; dT<500; ++dT)
    {
        INFO("dT = " + std::to_string(dT));
        solvers::Result BCTR{solvers::newton<models::LKT_BCT>(dT, C0, A, V0, R0)};
        solvers::Result CLWR{solvers::newton<models::CLW>(dT, C0, ALin, V0, R0)};
        INFO("V LKT-BCT = " + std::to_string(BCTR.V) + ", V CLW = " + std::to_string(CLWR.V) + '\n');
        INFO("R LKT-BCT = " + std::to_string(BCTR.R) + ", R CLW = " + std::to_string(CLWR.R) + '\n');

        REQUIRE(BCTR.hasConverged);
        REQUIRE(CLWR.hasConverged);
        REQUIRE((std::abs(BCTR.R - CLWR.R)/BCTR.R) < 0.05); // maximum of 5% error
        REQUIRE((std::abs(BCTR.V - CLWR.V)/BCTR.V) < 0.05);

        if (CLWR.V > (Vd/20))
            break;
        std::tie(V0, R0) = std::tie(CLWR.V, CLWR.R);
    }

}

// currently fails at intermediate dTs (~33K) meaning there is likely an issue with WLCYZ
TEST_CASE("Linearised WLCYZ model agrees with GD for Alloy with linear phase diagram", "[Solvers]")
{
    // WLCYZ paper states this model reduces to non-equilibrium bulk diffusion adjusted LKT-BCT for linear liquidus and
    // solidus.
    const alloys::Alloy A{alloys::NiB2007_atp}; // WLCYZ capable Alloy

    // linearised phase diagram identities
    std::vector<alloys::Fit> TlAtC{alloys::Fit{{A.Tm, A.m}}}; // Tl = Tm + ml*C
    std::vector<alloys::Fit> ClAtT{alloys::Fit{{-A.Tm/A.m, 1/A.m}}}; // Cl = (T-Tm)/ml
    double ms{A.m/A.k0}; // derived from combining Cl = (T-Tm)/ml and Cs = (T-Tm)/ms given k0 = Cs/Cl
    std::vector<alloys::Fit> CsAtT{alloys::Fit{{-A.Tm/ms, 1/ms}}}; // Cs = (T-Tm)/ms
    const alloys::Alloy ALin{ // linearised version
        A.L, A.Cp, A.m, A.k0, A.r, A.D, A.a, A.o, A.a0, A.V0, A.Tm, -1, -1, TlAtC, ClAtT, CsAtT, A.Vd
    };

    const double C0{0.7}, dT0{1};
    double V0{approx::getV(dT0, C0, A)}, R0{approx::getR(dT0, C0, A)};
    
    for (double dT{dT0}; dT<300; ++dT)
    {
        INFO("dT = " + std::to_string(dT));
        solvers::Result GDR{solvers::newton<models::GD>(dT, C0, ALin, V0, R0)};
        solvers::Result WLCYZR{solvers::newton<models::WLCYZ>(dT, C0, ALin, V0, R0)};
        INFO("V LKT-BCT = " + std::to_string(GDR.V) + ", V WLCYZ = " + std::to_string(WLCYZR.V) + '\n');
        INFO("R LKT-BCT = " + std::to_string(GDR.R) + ", R WLCYZ = " + std::to_string(WLCYZR.R) + '\n');

        REQUIRE(GDR.hasConverged);
        REQUIRE(WLCYZR.hasConverged);

        double f1{}, f2{};
        models::DTs DTs{};
        if ((std::abs(GDR.V - WLCYZR.V)/GDR.V) > 0.05)
        {
            std::tie(f1, f2, DTs) = models::GD(GDR.V, GDR.R, dT, C0, ALin);
            std::tie(f1, f2, DTs) = models::WLCYZ(GDR.V, GDR.R, dT, C0, ALin);
        }

        REQUIRE((std::abs(GDR.R - WLCYZR.R)/GDR.R) < 0.05); // maximum of 5% error
        REQUIRE((std::abs(GDR.V - WLCYZR.V)/GDR.V) < 0.05);

        std::tie(V0, R0) = std::tie(WLCYZR.V, WLCYZR.R);
    }
}
