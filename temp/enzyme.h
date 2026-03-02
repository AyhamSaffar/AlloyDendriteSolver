// forward declarations for common Enzyme automatic differentation functionality

#ifndef ENZYME_H
#define ENZYME_H

int enzyme_dup;
int enzyme_dupnoneed;
int enzyme_out;
int enzyme_const;

template < typename return_type, typename ... T >
return_type __enzyme_fwddiff(void*, T ... );

template < typename return_type, typename ... T >
return_type __enzyme_autodiff(void*, T ... );

#endif
