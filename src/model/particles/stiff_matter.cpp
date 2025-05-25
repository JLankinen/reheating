#include <cmath>
#include <numbers>

#include "stiff_matter.hpp"


EnergyDensity StiffMatter::energyDensity()
{
    return [this](HighPrecision t) -> HighPrecision
    {
        return HighPrecision(1) / (HighPrecision(24.0) * std::numbers::pi * p.G_N * pow(t, HighPrecision(2.0)));
    };
};