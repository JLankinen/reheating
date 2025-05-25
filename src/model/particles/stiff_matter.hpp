#ifndef STIFF_MATTER_H_
#define STIFF_MATTER_H_

#include "parameters/parameters.hpp"
#include "utils/types.hpp"


class StiffMatter final
{
    private:
        ModelParameters p;
    public:
        explicit StiffMatter(const ModelParameters& _p) : p{_p} {};

        EnergyDensity energyDensity();
};

#endif