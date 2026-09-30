// forward declarations for common Enzyme automatic differentation functionality

#ifndef ENZYME_H
#define ENZYME_H

// label for object pointer arguement whose derivative should be stored to a second object pointer of the same type.
extern int enzyme_dup;
// same as above but for when you dont need the value of the object, just its derivative. Should be faster.
extern int enzyme_dupnoneed;
extern int enzyme_out; // label for arguements whose derivate should be calculated and returned.
extern int enzyme_const; // label for arguements whose derivate should not be calculated. May still modify arugements.

// should be used when autodiffing functions that take const object arguements. This lets Enzyme check to see if the
// object arguements need to be differentiated at runtime, and if so allocate memory for a "shadow pointer" to store
// these gradients. See https://enzymead.github.io/Enzyme.jl/stable/faq/#faq-runtime-activity for more detail.
extern int enzyme_runtime_activity;

template < typename return_type, typename ... T >
return_type __enzyme_fwddiff(void*, T ... );

template < typename return_type, typename ... T >
return_type __enzyme_autodiff(void*, T ... );

#endif
