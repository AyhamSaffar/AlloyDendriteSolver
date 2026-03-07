// forward declarations for common Enzyme automatic differentation functionality

#ifndef ENZYME_H
#define ENZYME_H

// label for (decayed) object pointer arguement whose derivative should be returned to a second pointer of the same type
int enzyme_dup;

// same as above but for when you dont need the return value of an out parameter, just the gradient. Should be faster.
int enzyme_dupnoneed;

int enzyme_out; // label for arguements whose derivate should be calculated
int enzyme_const; // label for arguements whose derivate should not be calculated

template < typename return_type, typename ... T >
return_type __enzyme_fwddiff(void*, T ... );

template < typename return_type, typename ... T >
return_type __enzyme_autodiff(void*, T ... );

#endif
