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
    HighPrecision m = HighPrecision("1e32");           // Mass of the particle in GeV
    HighPrecision lambda = HighPrecision(0.01);        // Dimensionless coupling constant << 1.
    HighPrecision b = HighPrecision(1.0);              // Dimensionless expansion parameter.
    HighPrecision xi = HighPrecision(1.0 / 6.0);       // Gravitational coupling.
    HighPrecision G_N = HighPrecision("1.683e-37");    // Gravitational constant in GeV^-2
    // Universe matter content; n = 0 Minkowskian, n = 1 stiff, n = 2 radiation, n = 4 matter.
    HighPrecision alpha(HighPrecision n);
    
    void logParameters() const {
        std::string separator(54, '-');
        std::cout << "\n🌌 Cosmological Simulation Run Begins 🌌\n";
        std::cout << "---------------------------------------\n";
        std::cout << "🚀 Initializing model parameters...\n";
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
