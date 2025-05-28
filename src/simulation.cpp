#include "model/energy/creation_decay.hpp"
#include "model/energy/energy_densities.hpp"
#include "parameters/parameters.hpp"


/**
 * Runs the simulation with given parameters. Results are stored in a parquet file.
 */
void RunSimulation(ModelParameters p)
{
    // Calculation in stiff matter phase:
    auto t_eq = EqualTime(RhoStiff, RhoPhiStiff, RhoChiStiff, p, p.t0, 1e-10);

};


class PhiParticle {
    private:
        ModelParameters params;
    
    public:
        PhiParticle(const ModelParameters& p) : params(p) {}
    
        HighPrecision RhoStiff(HighPrecision t) const {
            // Use params directly without passing it as argument
            return /* computation using params and t */;
        }
    
        // Other member functions using params...
    };