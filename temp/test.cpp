#include <iostream>
#include "enzyme.h"

double f(double x) {return x * x;}

int main()
{
    double x = 5.0;
    double dx = 1.0;
    double df_dx = __enzyme_fwddiff<double>((void*)f, enzyme_dup, x, dx); 
    std::cout << "f(x) = " << f(x) << " f'(x) = " << df_dx << '\n';
}
