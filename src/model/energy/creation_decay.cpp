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

double PhiCreationRate(ModelParameters& p, double t)
{
    double arg = -pow((3.0 * p.m * t / 2.0), 2.0 / 3.0);

    double airyAi = boost::math::airy_ai(arg);
    double airyBi = boost::math::airy_bi(arg);

    double airySum = pow(airyAi, 2.0) + pow(airyBi, 2.0);

    double factor = 3.0 * pow(p.m * p.b, 13.0 / 3.0) / (32.0 * p.b);

    return factor * t * airySum;
}


ChiDecayRate::ChiDecayRate(ModelParameters& p_, double n_, double t0_):
        p{p_}, n{n_}, t0{t0_}, alpha{p_.alpha(n_)}
        {
            double argt0 = p.m * t0;
            double factor2 = pow(p.lambda * t0, 2.0) / 64.0;


            double J_alpha_t0  = boost::math::cyl_bessel_j(alpha, argt0);
            double J_alpha1_t0 = boost::math::cyl_bessel_j(alpha - 1.0, argt0);
            double J_alphaP1_t0 = boost::math::cyl_bessel_j(alpha + 1.0, argt0);
            double N_alpha_t0  = boost::math::cyl_neumann(alpha, argt0);
            double N_alpha1_t0 = boost::math::cyl_neumann(alpha - 1.0, argt0);
            double N_alphaP1_t0 = boost::math::cyl_neumann(alpha + 1.0, argt0);
        
            double bessel2 = pow(J_alpha_t0, 2) 
                                  - J_alpha1_t0 * J_alphaP1_t0 
                                  - N_alphaP1_t0 * N_alpha1_t0 
                                  + pow(N_alpha_t0, 2);
            initialBessel = factor2 * bessel2;        
        };

// Use as function.
double ChiDecayRate::operator()(double t) const
        {
            double arg = p.m * t;
            double factor1 = pow(p.lambda * t, 2.0) / 64.0;
            // Bessel functions for t
            double J_alpha_t  = boost::math::cyl_bessel_j(alpha, arg);
            double J_alpha1_t = boost::math::cyl_bessel_j(alpha - 1.0, arg);
            double J_alphaP1_t = boost::math::cyl_bessel_j(alpha + 1.0, arg);
            double N_alpha_t  = boost::math::cyl_neumann(alpha, arg);
            double N_alpha1_t = boost::math::cyl_neumann(alpha - 1.0, arg);
            double N_alphaP1_t = boost::math::cyl_neumann(alpha + 1.0, arg);

            double bessel1 = pow(J_alpha_t, 2) 
                                - J_alpha1_t * J_alphaP1_t 
                                - N_alphaP1_t * N_alpha1_t 
                                + pow(N_alpha_t, 2);
            return factor1 * bessel1 - initialBessel;
        }