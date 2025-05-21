#ifndef CREATION_DECAY_H_
#define CREATION_DECAY_H_

#include <boost/multiprecision/cpp_dec_float.hpp>
#include "parameters/parameters.hpp"

/**
 * Decay rate for massive phi particle to decay into massless chi particles.
 */
HighPrecision ChiDecayRate(ModelParameters p, HighPrecision n, HighPrecision t0, HighPrecision t);

/**
 * Creation rate for massive phi particles in a stiff matter dominated universe.
 */
HighPrecision PhiCreationRate(ModelParameters p, HighPrecision t);

/*struct ChiDecayRate
{
    ModelParameters p;
    HighPrecision n;
    HighPrecision initialBessel;

    ChiDecayRate(ModelParameters p_, HighPrecision n_);
    HighPrecision operator()(HighPrecision t);
};
*/

#endif