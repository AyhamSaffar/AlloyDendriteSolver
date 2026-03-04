#include <iostream>

int enzyme_out;
int enzyme_const;

template < typename return_type, typename ... T >
return_type __enzyme_autodiff(void*, T ... );

struct MyStruct{double A{}; double B{}; double C{};};

double myFunc(double x, const MyStruct& s) {return x * s.A;}

int main()
{
    MyStruct s{2.0, 1.0, 4.0};
    double x{3.0};
    double DmyFuncDx {__enzyme_autodiff<double>((void*)myFunc, enzyme_out, x, enzyme_const, s)};

    std::cout << "DmyFuncDx: " << DmyFuncDx << '\n';

    return 0;
}




