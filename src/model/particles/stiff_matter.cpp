#include <cmath>
#include <numbers>

#include "model/particles/stiff_matter.hpp"


EnergyDensity StiffMatter::energyDensity()
{
    return [this](HighPrecision t) -> HighPrecision
    {
        return HighPrecision(1.0) / (HighPrecision(24.0) * std::numbers::pi * this->p.G_N * pow(t, HighPrecision(2.0)));
    };
};