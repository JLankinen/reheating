#include <cmath>
#include <numbers>

#include "model/particles/stiff_matter.hpp"


EnergyDensity StiffMatter::energyDensity()
{
    return [this](double t) -> double
    {
        return 1.0 / (24.0 * std::numbers::pi * this->p.G_N * pow(t, 2.0));
    };
};