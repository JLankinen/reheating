#ifndef CREATION_DECAY_H_
#define CREATION_DECAY_H_

#include <boost/multiprecision/cpp_dec_float.hpp>
#include "parameters/parameters.hpp"

/**
 * Creation rate for massive phi particles in a stiff matter dominated universe.
 */
HighPrecision PhiCreationRate(ModelParameters& p, HighPrecision t);

/**
 * @brief Computes the decay rate of the massless chi particle as a function of time.
 *
 * Precomputes and stores the time-independent part of the decay rate 
 * integral (at the lower bound t0) to avoid redundant computation during simulation runs.
 * Speeds up the computation effectively 715x.
 * Provides a callable interface via operator() to evaluate the full decay rate 
 * at arbitrary time t.
 */
class ChiDecayRate
{
    private:
        const ModelParameters& p;
        HighPrecision n;
        HighPrecision t0;
        HighPrecision alpha;
        HighPrecision initialBessel;
    public:
        ChiDecayRate(ModelParameters& p_, HighPrecision n_, HighPrecision t0_);

        HighPrecision operator()(HighPrecision t) const;
};
#endif