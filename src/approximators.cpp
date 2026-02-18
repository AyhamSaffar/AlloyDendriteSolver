#include <numbers>
#include <cmath>
#include "approximators.h"

namespace approx
{
    double getTipRaius(double r, double m, double k0, double C0, double dT)
    {
        using std::numbers::pi;
        using std::pow;
        return 6.64*pi*pi*r * pow(-m*(1-k0), 0.25) * (pow(C0, 0.25)/pow(dT, 1.25));
    }

    double getTipVelocity(double D, double m, double k0, double r, double dT, double C0)
    {
        using std::numbers::pi;
        using std::pow;
        return (pow(dT, 2.5) / pow(C0, 1.5)) * D / (5.51*pi*pi*pow(-m*(1-k0), 1.5)*r);
    }
}
