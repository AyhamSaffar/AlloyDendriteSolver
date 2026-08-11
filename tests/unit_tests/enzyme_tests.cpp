#include <catch2/catch_test_macros.hpp>
#include "enzyme.h"

double cube(double x) {return x*x*x;}

double cubeDiffed(double x) {return __enzyme_autodiff<double>((void*)cube, enzyme_out, x);}

TEST_CASE("functions that call differentiated functions can be differentiated correctly", "[enzyme]")
{
    double x{10};
    double cubeDoubleDiffed{__enzyme_autodiff<double>((void*)cubeDiffed, enzyme_out, x)};
    REQUIRE(cubeDoubleDiffed == 6*x); // d^2(x^3)/dx^2  = d(3x^2)/dx = 6x 
}
