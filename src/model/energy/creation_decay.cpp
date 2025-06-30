#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/math/special_functions/hankel.hpp>
#include <boost/math/special_functions/airy.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <cmath>
#include "model/energy/creation_decay.hpp"

using boost::math::cyl_hankel_1;
using boost::math::cyl_hankel_2;
using boost::math::airy_ai;
using boost::math::airy_bi;
using boost::math::cyl_bessel_j;
using boost::math::cyl_neumann;
using HighPrecision = boost::multiprecision::cpp_dec_float_100;

HighPrecision PhiCreationRate(ModelParameters& p, HighPrecision t)
{
    HighPrecision arg = -pow((HighPrecision(3.0) * p.m * t / HighPrecision(2.0)), HighPrecision(2.0) / HighPrecision(3.0));

    HighPrecision airyAi = boost::math::airy_ai(arg);
    HighPrecision airyBi = boost::math::airy_bi(arg);

    HighPrecision airySum = pow(airyAi, 2.0) + pow(airyBi, 2.0);

    HighPrecision factor = HighPrecision(3.0) * pow(p.m * p.b, HighPrecision(13.0) / HighPrecision(3.0)) / (HighPrecision(32.0) * p.b);

    return factor * t * airySum;
}


ChiDecayRate::ChiDecayRate(ModelParameters& p_, HighPrecision n_, HighPrecision t0_):
        p{p_}, n{n_}, t0{t0_}, alpha{p_.alpha(n_)}
        {
            HighPrecision argt0 = HighPrecision(p.m * t0);
            HighPrecision factor2 = pow(p.lambda * t0, HighPrecision(2.0)) / HighPrecision(64.0);


            HighPrecision J_alpha_t0  = boost::math::cyl_bessel_j(alpha, argt0);
            HighPrecision J_alpha1_t0 = boost::math::cyl_bessel_j(alpha - HighPrecision(1), argt0);
            HighPrecision J_alphaP1_t0 = boost::math::cyl_bessel_j(alpha + HighPrecision(1), argt0);
            HighPrecision N_alpha_t0  = boost::math::cyl_neumann(alpha, argt0);
            HighPrecision N_alpha1_t0 = boost::math::cyl_neumann(alpha - HighPrecision(1), argt0);
            HighPrecision N_alphaP1_t0 = boost::math::cyl_neumann(alpha + HighPrecision(1), argt0);
        
            HighPrecision bessel2 = pow(J_alpha_t0, 2) 
                                  - J_alpha1_t0 * J_alphaP1_t0 
                                  - N_alphaP1_t0 * N_alpha1_t0 
                                  + pow(N_alpha_t0, 2);
            initialBessel = factor2 * bessel2;        
        };

// Use as function.
HighPrecision ChiDecayRate::operator()(HighPrecision t) const
        {
            HighPrecision arg = HighPrecision(p.m * t);
            HighPrecision factor1 = pow(p.lambda * t, HighPrecision(2.0)) / HighPrecision(64.0);
            // Bessel functions for t
            HighPrecision J_alpha_t  = boost::math::cyl_bessel_j(alpha, arg);
            HighPrecision J_alpha1_t = boost::math::cyl_bessel_j(alpha - HighPrecision(1), arg);
            HighPrecision J_alphaP1_t = boost::math::cyl_bessel_j(alpha + HighPrecision(1), arg);
            HighPrecision N_alpha_t  = boost::math::cyl_neumann(alpha, arg);
            HighPrecision N_alpha1_t = boost::math::cyl_neumann(alpha - HighPrecision(1), arg);
            HighPrecision N_alphaP1_t = boost::math::cyl_neumann(alpha + HighPrecision(1), arg);

            HighPrecision bessel1 = pow(J_alpha_t, 2) 
                                - J_alpha1_t * J_alphaP1_t 
                                - N_alphaP1_t * N_alpha1_t 
                                + pow(N_alpha_t, 2);
            return factor1 * bessel1 - initialBessel;
        }