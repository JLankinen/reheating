// Contains type definitions used throughout the program.

#ifndef TYPES_H_
#define TYPES_H_

#include <functional>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "parameters/parameters.hpp"

using EnergyDensity = std::function<double(double t)>;

#endif