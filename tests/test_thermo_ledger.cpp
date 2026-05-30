#include <gtest/gtest.h>

#include "../LIB/statmech.h"

#include <cmath>

using statmech::StatMechEngine;
using statmech::kB_kcal;

namespace {
constexpr double kT = 300.0;
constexpr double kTol = 1e-9;
}

TEST(ThermoLedger, SingleStateIdentities) {
    StatMechEngine engine(kT);
    engine.add_sample(-12.0);

    const auto b = engine.compute_breakdown();

    EXPECT_NEAR(b.G_config_kcal_mol, -kB_kcal * kT * b.logZ_config, kTol);
    EXPECT_NEAR(b.H_eff_kcal_mol, -12.0, kTol);
    EXPECT_NEAR(b.S_config_kcal_mol_K, 0.0, kTol);
    EXPECT_NEAR(b.minus_T_S_config_kcal_mol, b.G_config_kcal_mol - b.H_eff_kcal_mol, kTol);
    EXPECT_NEAR(b.Cv_kcal_mol_K, 0.0, kTol);
    EXPECT_NEAR(b.sigma_E_kcal_mol, 0.0, kTol);
    EXPECT_NEAR(b.G_total_kcal_mol, b.G_config_kcal_mol, kTol);
}

TEST(ThermoLedger, TwoEqualStatesHaveConfigEntropy) {
    StatMechEngine engine(kT);
    engine.add_sample(-10.0);
    engine.add_sample(-10.0);

    const auto b = engine.compute_breakdown();

    const double expected_G = -10.0 - kB_kcal * kT * std::log(2.0);
    EXPECT_NEAR(b.logZ_config, std::log(2.0) - (-10.0 / (kB_kcal * kT)), 1e-8);
    EXPECT_NEAR(b.G_config_kcal_mol, expected_G, 1e-8);
    EXPECT_NEAR(b.H_eff_kcal_mol, -10.0, kTol);
    EXPECT_NEAR(b.S_config_kcal_mol_K, kB_kcal * std::log(2.0), 1e-10);
    EXPECT_NEAR(b.Cv_kcal_mol_K, 0.0, kTol);
}

TEST(ThermoLedger, CorrectionsSumIntoTotalOnly) {
    StatMechEngine engine(kT);
    engine.add_sample(-10.0);
    engine.add_sample(-8.0);

    const auto legacy = engine.compute();
    const auto b = engine.compute_breakdown(-0.3, 0.2, -0.1, true, true, true);

    EXPECT_NEAR(b.G_config_kcal_mol, legacy.free_energy, kTol);
    EXPECT_NEAR(b.H_eff_kcal_mol, legacy.mean_energy, kTol);
    EXPECT_NEAR(b.S_config_kcal_mol_K, legacy.entropy, kTol);
    EXPECT_NEAR(b.G_total_kcal_mol, legacy.free_energy - 0.3 + 0.2 - 0.1, kTol);
    EXPECT_TRUE(b.has_vib);
    EXPECT_TRUE(b.has_natural);
    EXPECT_TRUE(b.has_other);
}

TEST(ThermoLedger, SigmaEnergyMatchesVariance) {
    StatMechEngine engine(kT);
    engine.add_sample(-11.0);
    engine.add_sample(-8.0);
    engine.add_sample(-7.0, 2.0);

    const auto legacy = engine.compute();
    const auto b = engine.compute_breakdown();
    const double variance = legacy.mean_energy_sq - legacy.mean_energy * legacy.mean_energy;

    EXPECT_NEAR(b.sigma_E_kcal_mol, std::sqrt(std::max(0.0, variance)), kTol);
}
