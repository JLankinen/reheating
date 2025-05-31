#ifndef PARAMS_H_
#define PARAMS_H_

#include <cmath>
#include <string>
#include <iostream>
#include <boost/multiprecision/cpp_dec_float.hpp>
using HighPrecision = boost::multiprecision::cpp_dec_float_100;

struct ModelParameters
{
    HighPrecision t0 = HighPrecision("1e-32"); // = HighPrecision(1e-32);  // Initial time in seconds or 1.86e11 in Planck units.
    HighPrecision m = HighPrecision("1e30");            // Mass of the particle in GeV
    HighPrecision lambda = 0.01;       // Dimensionless coupling constant << 1.
    HighPrecision b = 1;            // Dimensionless expansion parameter.
    HighPrecision xi = HighPrecision(1.0 / 6.0);          // Gravitational coupling.
    HighPrecision G_N = HighPrecision("1.683e-37"); // Gravitational constant in GeV^-2
    // Universe matter content; n = 0 Minkowskian, n = 1 stiff, n = 2 radiation, n = 4 matter.
    HighPrecision alpha(HighPrecision n);
    
    void logParameters() const {
        std::string separator(94, '-');
        std::cout << "[ModelParameters] "
                  << "t0: " << t0 << " | "
                  << "m: " << m << " | "
                  << "lambda: " << lambda << " | "
                  << "b: " << b << " | "
                  << "xi: " << xi << " | "
                  << "G_N: " << G_N << " | " << "\n"
                  << separator << std::endl;
    }
};


inline HighPrecision ModelParameters::alpha(HighPrecision n)
{
    return sqrt(HighPrecision(1.0) - n*(n - HighPrecision(2.0))
           * (HighPrecision(6.0)*xi - HighPrecision(1.0))) / (HighPrecision(2.0) + n);
}

#endif
