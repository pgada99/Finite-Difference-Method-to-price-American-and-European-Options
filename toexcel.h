// Pietro Gadaleta / Francesco Stampini
// Finite Differences Option Pricing
// C++ Project

/*! \file toexcel.h contains the definition of the function "get_results" which we pass to Excel VBA. */

#ifndef toexcel_h
#define toexcel_h

#include "Option.h"

#pragma GCC visibility push(default)

//!By declaring extern C, the compiler won't mangle the function names, and it's exported with the understanding of C linkage.
extern "C"
{
//!Function that returns all the results in the class Option, in order to export them in VBA.
double get_results(double strike, double S_0, double r, double sigma,
                    double t_fin, double Nd, double Md, double calld,
                   double americand, double intervald, double extra_plotsd, double * results, double * price, double * greek_vect );
}

#pragma GCC visibility pop
#endif /* toexcel_h */
