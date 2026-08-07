// forward declarations for common Enzyme automatic differentation functionality

#ifndef ENZYME_H
#define ENZYME_H

// label for object pointer arguement whose derivative should be stored to a second object of the same type
extern int enzyme_dup;

// same as above but for when you dont need the return value of an out parameter, just the gradient. Should be faster.
extern int enzyme_dupnoneed;

extern int enzyme_out; // label for arguements whose derivate should be calculated
extern int enzyme_const; // label for arguements whose derivate should not be calculated

template < typename return_type, typename ... T >
return_type __enzyme_fwddiff(void*, T ... );

template < typename return_type, typename ... T >
return_type __enzyme_autodiff(void*, T ... );

#endif
