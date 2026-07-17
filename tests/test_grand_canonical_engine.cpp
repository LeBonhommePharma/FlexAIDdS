// test_grand_canonical_engine.cpp — GrandCanonicalEngine + PartitionFunctionBase
//
// Covers multi-N Ξ(μ,V,T), competitive channel parity with GPF, NVT adapter,
// MOR/naloxone and 5-HT2A toy benchmarks (synthetic log_Z only).
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include "GrandCanonicalEngine.h"
#include "PartitionFunctionBase.h"
#include "GrandPartitionFunction.h"
#include "statmech.h"

#include <cmath>
#include <memory>
#include <vector>

using flexaids::CanonicalPartitionAdapter;
using flexaids::EnsembleKind;
using flexaids::GrandCanonicalEngine;
using flexaids::PartitionFunctionBase;

static constexpr double kT_300 = statmech::kB_kcal * 300.0;

// ════════════════════════════════════════════════════════════════════════
// PartitionFunctionBase / CanonicalPartitionAdapter
// ════════════════════════════════════════════════════════════════════════

TEST(PartitionFunctionBase, CanonicalAdapterFromLogZ) {
    CanonicalPartitionAdapter nvt(300.0, 10.0);
    EXPECT_EQ(nvt.ensemble(), EnsembleKind::NVT);
    EXPECT_NEAR(nvt.log_partition(), 10.0, 1e-15);
    EXPECT_NEAR(nvt.free_energy(), -kT_300 * 10.0, 1e-10);
    EXPECT_NEAR(nvt.temperature(), 300.0, 1e-15);
}

TEST(PartitionFunctionBase, CanonicalAdapterFromEngine) {
    statmech::StatMechEngine eng(300.0);
    eng.add_sample(-5.0, 1.0);
    eng.add_sample(-4.0, 2.0);
    auto th = eng.compute();

    CanonicalPartitionAdapter nvt(eng);
    EXPECT_NEAR(nvt.log_partition(), th.log_Z, 1e-12);
    EXPECT_NEAR(nvt.free_energy(), th.free_energy, 1e-12);
    EXPECT_TRUE(nvt.has_full_thermo());
}

TEST(PartitionFunctionBase, PolymorphicSurface) {
    CanonicalPartitionAdapter nvt(310.0, 5.0);
    GrandCanonicalEngine gce(310.0);
    gce.set_active_channel(GrandCanonicalEngine::ActiveChannel::Competitive);
    gce.add_competitive_ligand("A", 5.0, 1.0);

    PartitionFunctionBase* bases[2] = {&nvt, &gce};
    EXPECT_EQ(bases[0]->ensemble(), EnsembleKind::NVT);
    EXPECT_EQ(bases[1]->ensemble(), EnsembleKind::muVT);
    EXPECT_NEAR(bases[0]->log_partition(), 5.0, 1e-12);
    // Competitive: Ξ = 1 + e^5 → ln Ξ = log(1+e^5)
    EXPECT_NEAR(bases[1]->log_partition(), std::log(1.0 + std::exp(5.0)), 1e-10);
}

// ════════════════════════════════════════════════════════════════════════
// Multi-N channel: Ξ = Σ_N λ^N Z_N
// ════════════════════════════════════════════════════════════════════════

TEST(GrandCanonicalEngine, MultiN_EmptyIsOne) {
    GrandCanonicalEngine gce(300.0);
    gce.set_active_channel(GrandCanonicalEngine::ActiveChannel::MultiN);
    EXPECT_NEAR(gce.log_Xi_multiN(), 0.0, 1e-15);  // Ξ=1
    EXPECT_NEAR(gce.mean_N_multiN(), 0.0, 1e-15);
}

TEST(GrandCanonicalEngine, MultiN_TwoLevelLangmuir) {
    // Classic Langmuir: Z_0 = 1, Z_1 = K (association constant in 1/M units
    // when λ = c/c°). At λ=1, K=1 → half occupancy.
    GrandCanonicalEngine gce(300.0);
    gce.set_active_channel(GrandCanonicalEngine::ActiveChannel::MultiN);
    gce.set_canonical_log_Z(0, 0.0);  // ln Z_0 = 0 → Z_0 = 1
    gce.set_canonical_log_Z(1, 0.0);  // ln Z_1 = 0 → Z_1 = 1
    gce.set_fugacity(1.0);            // λ = 1

    // Ξ = 1 + 1 = 2 → ln Ξ = ln 2
    EXPECT_NEAR(gce.log_Xi_multiN(), std::log(2.0), 1e-12);
    EXPECT_NEAR(gce.mean_N_multiN(), 0.5, 1e-12);
    EXPECT_NEAR(gce.occupancy_probability(0), 0.5, 1e-12);
    EXPECT_NEAR(gce.occupancy_probability(1), 0.5, 1e-12);
    EXPECT_NEAR(gce.var_N_multiN(), 0.25, 1e-12);
}

TEST(GrandCanonicalEngine, MultiN_MeanN_IncreasesWithFugacity) {
    GrandCanonicalEngine gce(300.0);
    gce.set_active_channel(GrandCanonicalEngine::ActiveChannel::MultiN);
    gce.set_canonical_log_Z(0, 0.0);
    gce.set_canonical_log_Z(1, 2.0);  // Z_1 = e^2
    gce.set_fugacity(0.01);
    const double n_lo = gce.mean_N_multiN();
    gce.set_fugacity(10.0);
    const double n_hi = gce.mean_N_multiN();
    EXPECT_GT(n_hi, n_lo);
    EXPECT_GE(n_hi, 0.0);
    EXPECT_LE(n_hi, 1.0);
}

TEST(GrandCanonicalEngine, MultiN_ChemicalPotentialConsistency) {
    GrandCanonicalEngine gce(300.0);
    gce.set_chemical_potential(0.0);  // μ=0 → λ=1
    EXPECT_NEAR(gce.fugacity(), 1.0, 1e-12);
    gce.set_concentration_M(0.01);
    EXPECT_NEAR(gce.fugacity(), 0.01, 1e-12);
    // μ = kT ln(0.01)
    EXPECT_NEAR(gce.chemical_potential(), kT_300 * std::log(0.01), 1e-10);
}

TEST(GrandCanonicalEngine, MultiN_FromStatMechEngine) {
    statmech::StatMechEngine n0(300.0);
    n0.add_sample(0.0, 1.0);  // empty reference-ish
    statmech::StatMechEngine n1(300.0);
    n1.add_sample(-8.0, 5.0);
    n1.add_sample(-6.0, 2.0);

    GrandCanonicalEngine gce(300.0);
    gce.set_active_channel(GrandCanonicalEngine::ActiveChannel::MultiN);
    gce.set_canonical_from_engine(0, n0);
    gce.set_canonical_from_engine(1, n1);
    gce.set_fugacity(1.0);

    EXPECT_GT(gce.log_Xi_multiN(), 0.0);
    EXPECT_GT(gce.mean_N_multiN(), 0.0);
    EXPECT_LE(gce.mean_N_multiN(), 1.0);
    // free_energy via PartitionFunctionBase = −kT ln Ξ
    EXPECT_NEAR(gce.free_energy(), -kT_300 * gce.log_partition(), 1e-12);
}

// ════════════════════════════════════════════════════════════════════════
// Competitive channel parity with GrandPartitionFunction
// ════════════════════════════════════════════════════════════════════════

TEST(GrandCanonicalEngine, CompetitiveParityWithGPF) {
    GrandCanonicalEngine gce(300.0);
    gce.add_competitive_ligand("A", 10.0, 1.0);
    gce.add_competitive_ligand("B", 8.0, 0.1);

    target::GrandPartitionFunction gpf(300.0);
    gpf.add_ligand("A", 10.0, 1.0);
    gpf.add_ligand("B", 8.0, 0.1);

    EXPECT_NEAR(gce.log_Xi_competitive(), gpf.log_Xi(), 1e-12);
    EXPECT_NEAR(gce.mean_N_competitive(), gpf.mean_N(), 1e-12);
    EXPECT_NEAR(gce.selectivity("A", "B"), gpf.selectivity("A", "B"), 1e-12);
    EXPECT_NEAR(gce.mixing_entropy(), gpf.mixing_entropy(), 1e-12);
}

TEST(GrandCanonicalEngine, LigandVectorApply) {
    GrandCanonicalEngine gce(300.0);
    gce.add_competitive_ligand("fentanyl", 12.0, 1e-9);
    gce.add_competitive_ligand("naloxone", 10.0, 1e-9);

    GrandCanonicalEngine::LigandVector lv;
    lv.names = {"fentanyl", "naloxone"};
    lv.concentrations_M = {1e-9, 1e-4};
    gce.set_ligand_vector(lv);
    gce.apply_ligand_vector_to_competitive();

    EXPECT_NEAR(gce.competitive().concentration("naloxone"), 1e-4, 1e-15);
    EXPECT_GT(gce.competitive().binding_probability("naloxone"),
              gce.competitive().binding_probability("fentanyl") * 0.0);
}

// ════════════════════════════════════════════════════════════════════════
// Toy benchmarks: MOR (fentanyl+naloxone) and 5-HT2A (5-MeO-DMT + 5-HT)
// Synthetic log_Z — not experimental affinity claims.
// ════════════════════════════════════════════════════════════════════════

TEST(GrandCanonicalEngine, Benchmark_MOR_FentanylNaloxone) {
    constexpr double T = 310.0;
    // Synthetic NVT ensembles (CF-proxy stand-ins)
    statmech::StatMechEngine fen(T);
    fen.add_sample(-12.0, 20);
    fen.add_sample(-10.5, 8);
    auto fen_th = fen.compute();

    statmech::StatMechEngine nal(T);
    nal.add_sample(-9.0, 15);
    nal.add_sample(-7.0, 10);
    auto nal_th = nal.compute();

    // NVT: fentanyl stronger alone
    EXPECT_LT(fen_th.free_energy, nal_th.free_energy);

    GrandCanonicalEngine gce(T);
    gce.set_active_channel(GrandCanonicalEngine::ActiveChannel::Competitive);
    gce.add_competitive_ligand("fentanyl", fen_th.log_Z, 1e-9);   // 1 nM
    gce.add_competitive_ligand("naloxone", nal_th.log_Z, 1e-6);  // 1 µM

    const double p_f = gce.competitive().binding_probability("fentanyl");
    const double p_n = gce.competitive().binding_probability("naloxone");
    const double p_e = gce.competitive().empty_probability();
    EXPECT_NEAR(p_f + p_n + p_e, 1.0, 1e-10);
    EXPECT_NEAR(gce.mean_N_competitive(), 1.0 - p_e, 1e-12);

    // Apparent Ki diagnostic positive
    EXPECT_GT(gce.competitive().apparent_Ki_M("fentanyl"), 0.0);
    EXPECT_GT(gce.competitive().apparent_Ki_M("naloxone"), 0.0);

    // Raise naloxone → occupancy shifts
    const double p_n0 = p_n;
    gce.set_competitive_concentration("naloxone", 1e-4);
    EXPECT_GT(gce.competitive().binding_probability("naloxone"), p_n0);

    // Collapse in [0,1]
    EXPECT_GE(gce.ligand_entropy_collapse(), 0.0);
    EXPECT_LE(gce.ligand_entropy_collapse(), 1.0);

    // Isotherm monotonic
    auto curve = gce.occupancy_vs_concentration(
        "naloxone", {1e-9, 1e-7, 1e-5, 1e-3});
    ASSERT_GE(curve.size(), 2u);
    for (std::size_t i = 1; i < curve.size(); ++i)
        EXPECT_GE(curve[i].p_species, curve[i - 1].p_species);
}

TEST(GrandCanonicalEngine, Benchmark_5HT2A_5MeODMT_Endogenous) {
    constexpr double T = 310.0;
    // Toy: 5-MeO-DMT high affinity; 5-HT (endogenous) moderate; same site
    const double log_Z_dmt = 14.0;
    const double log_Z_ht = 11.0;

    GrandCanonicalEngine gce(T);
    gce.add_competitive_ligand("5-MeO-DMT", log_Z_dmt, 1e-8);  // 10 nM
    gce.add_competitive_ligand("5-HT", log_Z_ht, 1e-6);        // 1 µM local

    EXPECT_NEAR(
        gce.competitive().binding_probability("5-MeO-DMT")
            + gce.competitive().binding_probability("5-HT")
            + gce.competitive().empty_probability(),
        1.0, 1e-10);

    // Intrinsic selectivity favors DMT (higher Z)
    EXPECT_GT(gce.log_intrinsic_selectivity("5-MeO-DMT", "5-HT"), 0.0);

    // At equal c, apparent = intrinsic
    gce.set_competitive_concentrations({"5-MeO-DMT", "5-HT"}, {1e-6, 1e-6});
    EXPECT_NEAR(gce.log_intrinsic_selectivity("5-MeO-DMT", "5-HT"),
                gce.competitive().log_selectivity("5-MeO-DMT", "5-HT"),
                1e-12);

    // High endogenous 5-HT can compete at equal intrinsic if c is huge
    gce.set_competitive_concentration("5-HT", 1e-3);
    // Not claiming experimental Ki — just probability conservation
    const double sum = gce.competitive().binding_probability("5-MeO-DMT")
                     + gce.competitive().binding_probability("5-HT")
                     + gce.competitive().empty_probability();
    EXPECT_NEAR(sum, 1.0, 1e-10);
    EXPECT_GT(gce.competitive().apparent_Ki_M("5-MeO-DMT"), 0.0);
}

// ════════════════════════════════════════════════════════════════════════
// Multi-N OpenMP path smoke (larger N_max)
// ════════════════════════════════════════════════════════════════════════

TEST(GrandCanonicalEngine, MultiN_LargeNmax_Smoke) {
    GrandCanonicalEngine gce(300.0);
    gce.set_active_channel(GrandCanonicalEngine::ActiveChannel::MultiN);
    // Soft multi-occupancy toy: Z_N = 1 / N!  → ln Z_N = −ln(N!)
    double ln_fact = 0.0;
    gce.set_canonical_log_Z(0, 0.0);
    for (int N = 1; N <= 64; ++N) {
        ln_fact += std::log(static_cast<double>(N));
        gce.set_canonical_log_Z(N, -ln_fact);
    }
    gce.set_fugacity(2.0);  // λ = 2 → Poisson-like mean ≈ 2 for ideal gas
    const double mean = gce.mean_N_multiN();
    EXPECT_TRUE(std::isfinite(mean));
    EXPECT_GT(mean, 0.0);
    EXPECT_LT(mean, 64.0);
    EXPECT_TRUE(std::isfinite(gce.log_Xi_multiN()));
    EXPECT_NEAR(gce.free_energy(), -kT_300 * gce.log_partition(), 1e-10);
}
