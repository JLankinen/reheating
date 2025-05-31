/* General integration function used throughout the simulation for consistency.*/

#ifndef INTEGRATION_H_
#define INTEGRATION_H_

#include <boost/math/quadrature/gauss_kronrod.hpp>

using boost::math::quadrature::gauss_kronrod;

namespace IntegrationUtils{

template<typename F, typename T>
T integrate(F&& f, T lower, T upper, double tol = 1e-15)
{
    static gauss_kronrod<T, 31> integrator;
    return integrator.integrate(std::forward<F>(f), lower, upper, tol);
}

};


#endif