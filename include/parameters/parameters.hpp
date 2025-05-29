#ifndef PARAMS_H_
#define PARAMS_H_

#include <cmath>
#include <boost/multiprecision/cpp_dec_float.hpp>
using HighPrecision = boost::multiprecision::cpp_dec_float_100;

struct ModelParameters
{
    HighPrecision t0; // = HighPrecision(1e-32);  // Initial time in seconds.
    HighPrecision m;            // Mass of the particle in GeV
    HighPrecision lambda;       // Dimensionless coupling constant << 1.
    HighPrecision b;            // Dimensionless expansion parameter.
    HighPrecision xi;           // Gravitational coupling.
    HighPrecision G_N = HighPrecision("6.708e-39"); // Gravitational constant in GeV^-2
    // Universe matter content; n = 0 Minkowskian, n = 1 stiff, n = 2 radiation, n = 4 matter.
    HighPrecision alpha(HighPrecision n);

};


inline HighPrecision ModelParameters::alpha(HighPrecision n)
{
    return sqrt(1 - n*(n - 2) * (6*xi - 1)) / (2 + n);
}

#endif
