#include <gtest/gtest.h>
#include <model/particles/phi_particle.hpp>

TEST(PhiParticleTest, EnergyDensityMatterMonotonicDecrease) {
    ModelParameters p;
    p.t0 = HighPrecision("1e-32");
    p.m = HighPrecision("1e36");
    p.lambda = HighPrecision("0.001");
    p.b = HighPrecision("1.0");
    PhiParticle phi{p};
    phi.setInitialRhoMatter(HighPrecision("1e-10"));

    auto rho = phi.energyDensityMatter(p.t0);

    HighPrecision t1("1e-28");
    HighPrecision t2("2e-28");

    HighPrecision val1 = rho(t1);
    HighPrecision val2 = rho(t2);

    EXPECT_GT(val1, val2);  // Should decrease with time
}

TEST(PhiParticleTest, EnergyDensityRadiationMonotonicDecrease) {
    ModelParameters p;
    p.t0 = HighPrecision("1e-32");
    p.m = HighPrecision("1e36");
    p.lambda = HighPrecision("0.001");
    p.b = HighPrecision("1.0");
    PhiParticle phi{p};
    phi.setInitialRhoRadiation(HighPrecision("1e-10"));

    auto rho = phi.energyDensityRadiation(p.t0);

    HighPrecision t1("1e-28");
    HighPrecision t2("2e-28");

    HighPrecision val1 = rho(t1);
    HighPrecision val2 = rho(t2);

    EXPECT_GT(val1, val2);  // Should decrease with time
}

