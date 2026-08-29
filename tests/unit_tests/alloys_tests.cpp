#include <catch2/catch_test_macros.hpp>
#include <vector>
#include "alloys.h"

TEST_CASE("Alloy phase diagram methods give correct results for linear phase diagram", "[alloys]")
{
    double Tm{500}, ml{-10}, ms{-20}; // melting temperature, liquidus gradient, and solidus gradient
    // fits given as coefficients of polynomials. 1st num is y intercept and second num is gradient
    std::vector<alloys::Fit> TlAtC{alloys::Fit{{Tm, ml}}};
    std::vector<alloys::Fit> ClAtT{alloys::Fit{{-Tm/ml, 1/ml}}};
    std::vector<alloys::Fit> CsAtT{alloys::Fit{{-Tm/ms, 1/ms}}};
       
    alloys::Alloy A{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, TlAtC, ClAtT, CsAtT, 1};

    double C{10}, T{50};
    REQUIRE(A.TlAtC(C) == (Tm+ml*C));
    REQUIRE(A.mlAtC(C) == ml);
    REQUIRE(A.ClAtT(T) == (T-Tm)/ml);
    REQUIRE(A.CsAtT(T) == (T-Tm)/ms);
    REQUIRE(A.mlAtT(T) == ml);
    REQUIRE(A.msAtT(T) == ms);
}
