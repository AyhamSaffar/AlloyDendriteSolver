// #ifndef DIFFERENTIALS_H
// #define DIFFERENTIALS_H

// #include <tuple>
// #include <iostream>
// #include "enzyme.h" 
// #include "alloy.h"
// #include "models.h"


// namespace diff
// {    
//     // Enzyme autodiff can only handle non out-parameter functions when they return a single value
//     template <models::ModelFunc func, int fToReturn>
//     double wrapper(double V, double R, double dT, double C0, const alloy::Alloy& A)
//     {
//         std::tuple<double, double> f{func(V, R, dT, C0, A)};
//         return std::get<fToReturn-1>(f);
//     }
//     struct Jacobian{double df1dV{}; double df1dR{}; double df2dV{}; double df2dR{};};
    
//     template <models::ModelFunc modelFunc>
//     inline Jacobian calculateGrads(double V, double R, double dT, double C0, alloy::Alloy A)
//     {
//         Jacobian J{};
//         // alloy object must be passed twice as enzye thinks all object pointers (e.g. struct reference) is a pointer
//         // to a decayed pointer parameter whose gradients should be stored in a second pointer to an identical type?
//         // here however, A doesn't seem to get modified, maybe because it is labelled with "enzyme_const" meaning its
//         // derivates are never calculated.
//         J.df1dV = __enzyme_autodiff<double>(
//             (void*)wrapper<modelFunc, 1>, enzyme_out, V, enzyme_const, R, dT, C0, A, A 
//         );
//         J.df1dR = __enzyme_autodiff<double>(
//             (void*)wrapper<modelFunc, 1>, enzyme_const, V, enzyme_out, R, enzyme_const, dT, C0, A, A
//         );
//         J.df2dV = __enzyme_autodiff<double>(
//             (void*)wrapper<modelFunc, 2>, enzyme_out, V, enzyme_const, R, dT, C0, A, A
//         );
//         J.df2dR = __enzyme_autodiff<double>(
//             (void*)wrapper<modelFunc, 2>, enzyme_const, V, enzyme_out, R, enzyme_const, dT, C0, A, A
//         );
//         return J;
//     }
// }

// std::ostream& operator<<(std::ostream& out, const diff::Jacobian& J)
// {
//     return out << "Jacobian(δf1/δV=" << J.df1dV << ", δf1/δR=" << J.df1dR << ", δf2/δV=" << J.df2dV << ", δf2/δR=" <<
//         J.df2dR << ')';
// }

// #endif
