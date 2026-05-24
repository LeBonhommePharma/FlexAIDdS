// tests/test_dift.cpp — unit tests for the DiFT torsional parametrization engine
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
//
// Covers: forward transform / round-trip exactness, Shannon-collapse
// truncation, QM–MM iterative refinement, per-bond torsional thermodynamics,
// Boltzmann inversion, circular mean, and the GA scoring adapter.

#include <gtest/gtest.h>

#include "../LIB/DiFT/DiFT.h"
#include "../LIB/DiFT/DiFTGAAdapter.h"

#include <cmath>
#include <numeric>
#include <vector>

using namespace dift;

namespace {

constexpr double kTwoPi = 6.283185307179586476925286766559;

// Build an M-point sample of  V(φ) = mean + Σ Aₙ cos(nφ − ωₙ).
std::vector<double> make_profile(int M, double mean,
                                 const std::vector<int>& mult,
                                 const std::vector<double>& amp,
                                 const std::vector<double>& phase) {
    std::vector<double> p(static_cast<std::size_t>(M));
    for (int k = 0; k < M; ++k) {
        const double phi = kTwoPi * k / M;
        double v = mean;
        for (std::size_t t = 0; t < mult.size(); ++t)
            v += amp[t] * std::cos(mult[t] * phi - phase[t]);
        p[static_cast<std::size_t>(k)] = v;
    }
    return p;
}

// Look up the amplitude of a given multiplicity in a spectrum.
double amp_of(const std::vector<FourierTerm>& spec, int n) {
    for (const auto& t : spec) if (t.multiplicity == n) return t.amplitude;
    return 0.0;
}

} // namespace

// ─── forward transform: power-of-two grid (radix-2 FFT path) ─────────────────
TEST(DiFTTransform, RecoversAmplitudesPow2) {
    DiFTEngine engine(300.0);
    const int M = 64;
    auto profile = make_profile(M, 1.5, {1, 3, 5},
                                {2.0, 1.0, 0.5}, {0.0, 0.7, -1.2});
    double mean = 0.0;
    auto spec = engine.transform(profile, mean);

    EXPECT_NEAR(mean, 1.5, 1e-9);
    EXPECT_NEAR(amp_of(spec, 1), 2.0, 1e-9);
    EXPECT_NEAR(amp_of(spec, 3), 1.0, 1e-9);
    EXPECT_NEAR(amp_of(spec, 5), 0.5, 1e-9);
    EXPECT_NEAR(amp_of(spec, 2), 0.0, 1e-9);   // absent frequency
    EXPECT_NEAR(amp_of(spec, 4), 0.0, 1e-9);
}

// ─── forward transform: non-power-of-two grid (direct DFT path) ──────────────
TEST(DiFTTransform, RecoversAmplitudesNonPow2) {
    DiFTEngine engine(300.0);
    const int M = 36;                          // typical 10°-resolution scan
    auto profile = make_profile(M, -0.4, {2, 4},
                                {1.3, 0.8}, {0.3, 2.1});
    double mean = 0.0;
    auto spec = engine.transform(profile, mean);

    EXPECT_NEAR(mean, -0.4, 1e-9);
    EXPECT_NEAR(amp_of(spec, 2), 1.3, 1e-8);
    EXPECT_NEAR(amp_of(spec, 4), 0.8, 1e-8);
}

// ─── round-trip: full spectrum reconstructs the input exactly ────────────────
TEST(DiFTTransform, RoundTripReconstruction) {
    DiFTEngine engine(300.0);
    const int M = 72;
    auto profile = make_profile(M, 0.9, {1, 2, 6},
                                {1.1, 0.6, 0.25}, {-0.5, 1.0, 2.8});
    double mean = 0.0;
    auto spec = engine.transform(profile, mean);

    TorsionalPotential pot;
    pot.terms = spec;
    pot.mean  = mean;
    auto model = pot.sample(M);

    for (int k = 0; k < M; ++k)
        EXPECT_NEAR(model[static_cast<std::size_t>(k)],
                    profile[static_cast<std::size_t>(k)], 1e-7);
    EXPECT_NEAR(DiFTEngine::r_squared(profile, model), 1.0, 1e-9);
}

TEST(DiFTTransform, RejectsShortProfile) {
    DiFTEngine engine;
    std::vector<double> tiny{1.0};
    double mean = 0.0;
    EXPECT_THROW(engine.transform(tiny, mean), std::invalid_argument);
}

// ─── spectral Shannon entropy ────────────────────────────────────────────────
TEST(DiFTSpectralEntropy, SingleModeCollapsesToOne) {
    DiFTEngine engine;
    auto profile = make_profile(64, 0.0, {2}, {3.0}, {0.4});
    double mean = 0.0;
    auto spec = engine.transform(profile, mean);

    const double H = spectral_entropy(spec);
    EXPECT_NEAR(H, 0.0, 1e-9);                 // one mode → zero entropy
    EXPECT_NEAR(std::exp(H), 1.0, 1e-9);       // N_eff ≈ 1
}

TEST(DiFTSpectralEntropy, EqualPowerSpreadsToLogN) {
    // Three equal-amplitude modes → H = ln 3, N_eff = 3.
    DiFTEngine engine;
    auto profile = make_profile(64, 0.0, {1, 2, 3},
                                {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0});
    double mean = 0.0;
    auto spec = engine.transform(profile, mean);
    EXPECT_NEAR(spectral_entropy(spec), std::log(3.0), 1e-9);
}

// ─── Shannon-collapse truncation in parametrize() ────────────────────────────
TEST(DiFTParametrize, KeepsDominantModesOnly) {
    DiFTEngine engine;
    // Two strong modes + a tiny ripple the collapse should discard.
    auto profile = make_profile(128, 0.2, {1, 4, 11},
                                {2.5, 1.8, 0.02}, {0.1, -0.6, 1.4});
    auto pot = engine.parametrize(profile);

    EXPECT_GE(pot.n_terms(), 2u);
    EXPECT_LE(pot.n_terms(), 3u);
    EXPECT_GT(pot.effective_modes, 1.0);
    EXPECT_LT(pot.effective_modes, 3.0);
    EXPECT_GT(pot.r_squared, 0.99);            // dominant modes capture the shape
}

TEST(DiFTParametrize, MaxMultiplicityGuard) {
    DiFTEngine engine;
    auto profile = make_profile(64, 0.0, {2, 9}, {1.0, 1.0}, {0.0, 0.0});
    auto pot = engine.parametrize(profile, /*max_multiplicity=*/6);
    for (const auto& t : pot.terms)
        EXPECT_LE(t.multiplicity, 6);
}

// ─── iterative QM–MM refinement (paper eq. 18) ───────────────────────────────
TEST(DiFTRefine, ConvergesToQMTarget) {
    DiFTEngine engine(300.0);
    const int M = 96;
    auto qm = make_profile(M, 0.0, {1, 2, 3},
                           {3.0, 1.5, 0.8}, {0.2, -0.9, 1.7});
    auto mm = make_profile(M, 0.0, {1, 2},
                           {2.0, 1.0}, {0.2, -0.9});   // missing the n=3 term

    auto corrected = engine.refine(qm, mm, /*lambda=*/0.5, /*r2_target=*/0.98);

    EXPECT_GE(corrected.r_squared, 0.98);
    EXPECT_GT(corrected.refinement_iters, 0);

    // mm + correction should now track qm.
    auto corr = corrected.sample(M);
    std::vector<double> total(static_cast<std::size_t>(M));
    for (int k = 0; k < M; ++k)
        total[static_cast<std::size_t>(k)] =
            mm[static_cast<std::size_t>(k)] + corr[static_cast<std::size_t>(k)];
    EXPECT_GE(DiFTEngine::r_squared(qm, total), 0.98);
}

TEST(DiFTRefine, RejectsMismatchedGrids) {
    DiFTEngine engine;
    std::vector<double> qm(64, 0.0), mm(32, 0.0);
    EXPECT_THROW(engine.refine(qm, mm), std::invalid_argument);
}

// ─── per-bond torsional thermodynamics ───────────────────────────────────────
TEST(DiFTThermodynamics, FreeRotorHasZeroExcessEntropy) {
    DiFTEngine engine(300.0);
    std::vector<double> flat(64, 0.0);          // V ≡ 0 → free rotor
    auto pot = engine.parametrize(flat);
    auto th  = engine.thermodynamics(pot);

    EXPECT_NEAR(th.partition_function, 1.0, 1e-6);
    EXPECT_NEAR(th.entropy,  0.0, 1e-6);
    EXPECT_NEAR(th.minus_TS, 0.0, 1e-6);
}

TEST(DiFTThermodynamics, ConfinedBondPaysEntropyPenalty) {
    DiFTEngine engine(300.0);
    // A deep single-well torsional barrier (5 kcal/mol amplitude).
    auto profile = make_profile(128, 0.0, {1}, {5.0}, {0.0});
    auto pot = engine.parametrize(profile);
    auto th  = engine.thermodynamics(pot);

    EXPECT_LT(th.partition_function, 1.0);      // confinement shrinks z
    EXPECT_LT(th.entropy,  0.0);                // entropy LOSS vs free rotor
    EXPECT_GT(th.minus_TS, 0.0);                // → positive ΔG penalty
    EXPECT_GT(th.mean_energy, 0.0);
}

TEST(DiFTThermodynamics, DeeperWellPaysMorePenalty) {
    DiFTEngine engine(300.0);
    auto shallow = engine.parametrize(
        make_profile(128, 0.0, {1}, {1.0}, {0.0}));
    auto deep = engine.parametrize(
        make_profile(128, 0.0, {1}, {6.0}, {0.0}));
    EXPECT_GT(engine.thermodynamics(deep).minus_TS,
              engine.thermodynamics(shallow).minus_TS);
}

// ─── Boltzmann inversion of a CG dihedral histogram (paper eq. 20) ───────────
TEST(DiFTBoltzmannInvert, UniformHistogramIsFlat) {
    DiFTEngine engine(300.0);
    std::vector<double> uniform(36, 100.0);
    auto energy = engine.boltzmann_invert(uniform);
    for (double e : energy) EXPECT_NEAR(e, 0.0, 1e-9);
}

TEST(DiFTBoltzmannInvert, PeakedHistogramBecomesWell) {
    DiFTEngine engine(300.0);
    std::vector<double> hist(36, 1.0);
    hist[18] = 1000.0;                          // sharp population peak
    auto energy = engine.boltzmann_invert(hist);

    EXPECT_NEAR(energy[18], 0.0, 1e-9);         // peak → energy minimum
    for (std::size_t k = 0; k < energy.size(); ++k) {
        EXPECT_GE(energy[k], -1e-12);           // all energies ≥ 0 after shift
        if (k != 18) { EXPECT_GT(energy[k], 0.0); }
    }
}

TEST(DiFTBoltzmannInvert, RejectsEmptyHistogram) {
    DiFTEngine engine;
    std::vector<double> empty(20, 0.0);
    EXPECT_THROW(engine.boltzmann_invert(empty), std::invalid_argument);
}

// ─── circular mean of phase offsets ──────────────────────────────────────────
TEST(DiFTCircularMean, HandlesWraparound) {
    // +170° and −170° straddle ±π; the circular mean is 180°, not 0°.
    std::vector<double> angles{ 170.0 * M_PI / 180.0,
                               -170.0 * M_PI / 180.0 };
    const double m = DiFTEngine::circular_mean(angles);
    EXPECT_NEAR(std::abs(m), M_PI, 1e-6);
}

TEST(DiFTCircularMean, SimpleAverage) {
    std::vector<double> angles{0.1, 0.2, 0.3};
    EXPECT_NEAR(DiFTEngine::circular_mean(angles), 0.2, 1e-6);
}

// ─── GA scoring adapter ──────────────────────────────────────────────────────
TEST(DiFTGAAdapter, ScoresTorsionalEnergyAndEntropy) {
    // Two rotatable bonds, each with a single-well DiFT potential. The π phase
    // makes V(φ) = −A·cos(nφ), so the well minimum sits at φ = 0.
    auto rbt0 = make_bond_torsion(
        make_profile(128, 0.0, {1}, {3.0}, {M_PI}), /*gene_index=*/0);
    auto rbt1 = make_bond_torsion(
        make_profile(128, 0.0, {2}, {2.0}, {M_PI}), /*gene_index=*/1);
    std::vector<RotatableBondTorsion> bonds{rbt0, rbt1};

    // Place both bonds at their well minima (φ = 0).
    std::vector<double> angles{0.0, 0.0};
    auto score = score_torsional(bonds, angles, 300.0);

    EXPECT_EQ(score.n_bonds, 2);
    EXPECT_NEAR(score.energy, 0.0, 1e-3);       // at the minima → ~0 energy
    EXPECT_GT(score.minus_TS, 0.0);             // confinement penalty > 0
    EXPECT_GT(score.total(), 0.0);
}

TEST(DiFTGAAdapter, EnergyRisesAwayFromMinimum) {
    // V(φ) = −4·cos(φ): well minimum at φ = 0, barrier crest at φ = π.
    auto rbt = make_bond_torsion(
        make_profile(128, 0.0, {1}, {4.0}, {M_PI}), 0);
    std::vector<RotatableBondTorsion> bonds{rbt};

    std::vector<double> at_min{0.0};
    std::vector<double> at_top{M_PI};
    const double e_min = score_torsional(bonds, at_min).energy;
    const double e_top = score_torsional(bonds, at_top).energy;
    EXPECT_GT(e_top, e_min + 1.0);              // barrier crest costs energy
}

TEST(DiFTGAAdapter, MismatchedInputReturnsZero) {
    auto rbt = make_bond_torsion(std::vector<double>(64, 0.0), 0);
    std::vector<RotatableBondTorsion> bonds{rbt};
    std::vector<double> angles{0.0, 0.0};       // wrong length
    auto score = score_torsional(bonds, angles);
    EXPECT_EQ(score.n_bonds, 0);
    EXPECT_EQ(score.total(), 0.0);
}
