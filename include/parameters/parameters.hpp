#ifndef PARAMS_H_
#define PARAMS_H_

#include <cmath>
#include <string>
#include <iostream>
#include <utils/types.hpp>

struct ModelParameters
{
    double t0 = double(1.51926764e-8); // Initial time in GeV^-1
    double m;                                   // Mass of the particle in GeV
    double lambda;                              // Coupling constant << 1.
    double b;                                   // Dimensionless expansion parameter.
    double xi;                                  // Gravitational coupling.
    double G_N = double(1.683e-37);    // Gravitational constant in GeV^-2
    // Universe matter content; n = 0 Minkowskian, n = 1 stiff, n = 2 radiation, n = 4 matter.
    double alpha(double n);
};


inline double ModelParameters::alpha(double n)
{
    return sqrt(double(1.0) - n*(n - double(2.0))
           * (double(6.0)*xi - double(1.0))) / (double(2.0) + n);
}


struct State
{
    double minkowski = 0.0;
    double stiff = 1.0;
    double radiation = 2.0;
    double matter = 4.0;
};


#endif
