#ifndef PARAMS_H_
#define PARAMS_H_

#include <cmath>
#include <string>
#include <iostream>
#include <boost/multiprecision/cpp_dec_float.hpp>
using HighPrecision = boost::multiprecision::cpp_dec_float_100;

struct ModelParameters
{
    HighPrecision t0 = HighPrecision("1.51926764e-8"); // Initial time in GeV^-1
    HighPrecision m;                                   // Mass of the particle in GeV
    HighPrecision lambda;                              // Dimensionless coupling constant << 1.
    HighPrecision b;                                   // Dimensionless expansion parameter.
    HighPrecision xi;                                  // Gravitational coupling.
    HighPrecision G_N = HighPrecision("1.683e-37");    // Gravitational constant in GeV^-2
    // Universe matter content; n = 0 Minkowskian, n = 1 stiff, n = 2 radiation, n = 4 matter.
    HighPrecision alpha(HighPrecision n);
};


inline HighPrecision ModelParameters::alpha(HighPrecision n)
{
    return sqrt(HighPrecision(1.0) - n*(n - HighPrecision(2.0))
           * (HighPrecision(6.0)*xi - HighPrecision(1.0))) / (HighPrecision(2.0) + n);
}

#endif
