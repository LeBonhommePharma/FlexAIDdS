// test_grand_partition.cpp — Unit tests for GrandPartitionFunction
//
// Verifies log-space arithmetic, binding probabilities, selectivity
// ratios, and ranking against analytical hand calculations.
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>
#include "GrandPartitionFunction.h"
#include "MultiSiteGPF.h"
#include "statmech.h"

#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

using namespace target;

static constexpr double kT_300 = statmech::kB_kcal * 300.0;  // ~0.596 kcal/mol

// ════════════════════════════════════════════════════════════════════════
// Basic construction
// ════════════════════════════════════════════════════════════════════════

TEST(GrandPartition, EmptyXi) {
    GrandPartitionFunction gpf(300.0);
    EXPECT_EQ(gpf.num_ligands(), 0);
    // Ξ = 1 (empty site only) → ln Ξ = 0
    EXPECT_NEAR(gpf.log_Xi(), 0.0, 1e-12);
    // p(empty) = 1/1 = 1
    EXPECT_NEAR(gpf.empty_probability(), 1.0, 1e-12);
}

TEST(GrandPartition, InvalidTemperature) {
    EXPECT_THROW(GrandPartitionFunction(0.0), std::invalid_argument);
    EXPECT_THROW(GrandPartitionFunction(-100.0), std::invalid_argument);
}

// ════════════════════════════════════════════════════════════════════════
// Single ligand
// ════════════════════════════════════════════════════════════════════════

TEST(GrandPartition, SingleLigand) {
    GrandPartitionFunction gpf(300.0);

    // Z_A = e^10 (arbitrary large partition function)
    double log_Z_A = 10.0;
    gpf.add_ligand("ligandA", log_Z_A);

    EXPECT_EQ(gpf.num_ligands(), 1);
    EXPECT_TRUE(gpf.has_ligand("ligandA"));

    // Ξ = 1 + Z_A = 1 + e^10
    double expected_log_xi = std::log(1.0 + std::exp(10.0));
    EXPECT_NEAR(gpf.log_Xi(), expected_log_xi, 1e-10);

    // p(A) = Z_A / Ξ = e^10 / (1 + e^10) ≈ 0.999955
    double expected_pA = std::exp(10.0) / (1.0 + std::exp(10.0));
    EXPECT_NEAR(gpf.binding_probability("ligandA"), expected_pA, 1e-8);

    // p(empty) = 1/Ξ ≈ 0.000045
    EXPECT_NEAR(gpf.empty_probability(), 1.0 - expected_pA, 1e-8);

    // F_bound = -kT * ln(Z_A) = -kT * 10
    EXPECT_NEAR(gpf.F_bound("ligandA"), -kT_300 * 10.0, 1e-10);
}

// ════════════════════════════════════════════════════════════════════════
// Two ligands — competitive binding
// ════════════════════════════════════════════════════════════════════════

TEST(GrandPartition, TwoLigandsCompetitive) {
    GrandPartitionFunction gpf(300.0);

    // Z_A = e^10, Z_B = e^8
    gpf.add_ligand("A", 10.0);
    gpf.add_ligand("B",  8.0);

    // Ξ = 1 + e^10 + e^8
    double xi = 1.0 + std::exp(10.0) + std::exp(8.0);
    EXPECT_NEAR(std::exp(gpf.log_Xi()), xi, xi * 1e-10);

    // selectivity A/B = Z_A/Z_B = e^(10-8) = e^2 ≈ 7.389
    EXPECT_NEAR(gpf.selectivity("A", "B"), std::exp(2.0), 1e-10);
    EXPECT_NEAR(gpf.selectivity("B", "A"), std::exp(-2.0), 1e-10);

    // Probabilities must sum to 1
    double sum = gpf.binding_probability("A")
               + gpf.binding_probability("B")
               + gpf.empty_probability();
    EXPECT_NEAR(sum, 1.0, 1e-10);
}

// ════════════════════════════════════════════════════════════════════════
// Equal ligands
// ════════════════════════════════════════════════════════════════════════

TEST(GrandPartition, EqualLigands) {
    GrandPartitionFunction gpf(300.0);

    // Three identical ligands with Z = e^5
    gpf.add_ligand("X", 5.0);
    gpf.add_ligand("Y", 5.0);
    gpf.add_ligand("Z", 5.0);

    // All should have equal binding probability
    double pX = gpf.binding_probability("X");
    double pY = gpf.binding_probability("Y");
    double pZ = gpf.binding_probability("Z");
    EXPECT_NEAR(pX, pY, 1e-12);
    EXPECT_NEAR(pY, pZ, 1e-12);

    // Selectivity between equal ligands = 1.0
    EXPECT_NEAR(gpf.selectivity("X", "Y"), 1.0, 1e-12);

    // All F_bound should be equal
    EXPECT_NEAR(gpf.F_bound("X"), gpf.F_bound("Y"), 1e-12);
}

// ════════════════════════════════════════════════════════════════════════
// Ranking
// ════════════════════════════════════════════════════════════════════════

TEST(GrandPartition, Ranking) {
    GrandPartitionFunction gpf(300.0);

    gpf.add_ligand("weak",    2.0);
    gpf.add_ligand("strong", 20.0);
    gpf.add_ligand("medium", 10.0);

    auto ranks = gpf.rank();
    ASSERT_EQ(ranks.size(), 3u);

    // Sorted by ΔG ascending (most negative = strongest binder first)
    EXPECT_EQ(ranks[0].name, "strong");
    EXPECT_EQ(ranks[1].name, "medium");
    EXPECT_EQ(ranks[2].name, "weak");

    // Check ΔG values
    EXPECT_NEAR(ranks[0].dG, -kT_300 * 20.0, 1e-10);
    EXPECT_NEAR(ranks[1].dG, -kT_300 * 10.0, 1e-10);
    EXPECT_NEAR(ranks[2].dG, -kT_300 *  2.0, 1e-10);

    // Probabilities must sum to < 1 (leaving room for empty)
    double psum = 0;
    for (const auto& r : ranks) psum += r.p_bound;
    EXPECT_LT(psum, 1.0);
    EXPECT_NEAR(psum + gpf.empty_probability(), 1.0, 1e-10);
}

// ════════════════════════════════════════════════════════════════════════
// Update (re-docking merge)
// ════════════════════════════════════════════════════════════════════════

TEST(GrandPartition, OverwriteLigand) {
    GrandPartitionFunction gpf(300.0);

    // Initial: Z_A = e^5
    gpf.add_ligand("A", 5.0);
    double dG_before = gpf.F_bound("A");

    // Re-dock with a better estimate: overwrite replaces the value
    gpf.overwrite_ligand("A", 8.0);

    // Z should now be e^8 (not e^5 + e^8)
    EXPECT_NEAR(gpf.F_bound("A"), -kT_300 * 8.0, 1e-10);

    // F_bound should be more negative after overwrite with larger Z
    EXPECT_LT(gpf.F_bound("A"), dG_before);
}

TEST(GrandPartition, MergeLigand) {
    GrandPartitionFunction gpf(300.0);

    // Initial: Z_A = e^5
    gpf.add_ligand("A", 5.0);
    double dG_before = gpf.F_bound("A");

    // Merge independent ensemble: Z_merged = 2·e^5
    gpf.merge_ligand("A", 5.0);

    // ln(Z_merged) = ln(2·e^5) = 5 + ln(2)
    double expected_dG = -kT_300 * (5.0 + std::log(2.0));
    EXPECT_NEAR(gpf.F_bound("A"), expected_dG, 1e-10);

    // F_bound should be more negative (more favorable) after merging
    EXPECT_LT(gpf.F_bound("A"), dG_before);
}

// ════════════════════════════════════════════════════════════════════════
// Remove ligand
// ════════════════════════════════════════════════════════════════════════

TEST(GrandPartition, RemoveLigand) {
    GrandPartitionFunction gpf(300.0);
    gpf.add_ligand("A", 10.0);
    gpf.add_ligand("B",  5.0);
    EXPECT_EQ(gpf.num_ligands(), 2);

    double log_xi_before = gpf.log_Xi();
    gpf.remove_ligand("A");
    double log_xi_after = gpf.log_Xi();

    EXPECT_EQ(gpf.num_ligands(), 1);
    EXPECT_FALSE(gpf.has_ligand("A"));
    EXPECT_TRUE(gpf.has_ligand("B"));
    EXPECT_LT(log_xi_after, log_xi_before);  // Ξ decreases after removal

    // Remove remaining ligand — Ξ must return to 1 (empty site only).
    gpf.remove_ligand("B");
    // log_Xi must equal log(1) = 0 when all ligands removed at zero concentration
    ASSERT_NEAR(gpf.log_Xi(), 0.0, 1e-12);
}

// ════════════════════════════════════════════════════════════════════════
// Error handling
// ════════════════════════════════════════════════════════════════════════

TEST(GrandPartition, DuplicateAddThrows) {
    GrandPartitionFunction gpf(300.0);
    gpf.add_ligand("A", 5.0);
    EXPECT_THROW(gpf.add_ligand("A", 10.0), std::invalid_argument);
}

TEST(GrandPartition, QueryMissingThrows) {
    GrandPartitionFunction gpf(300.0);
    EXPECT_THROW((void)gpf.binding_probability("X"), std::invalid_argument);
    EXPECT_THROW((void)gpf.F_bound("X"), std::invalid_argument);
    EXPECT_THROW(gpf.overwrite_ligand("X", 5.0), std::invalid_argument);
    EXPECT_THROW(gpf.merge_ligand("X", 5.0), std::invalid_argument);
    EXPECT_THROW(gpf.remove_ligand("X"), std::invalid_argument);
}

// ════════════════════════════════════════════════════════════════════════
// Numerical stability — extreme values
// ════════════════════════════════════════════════════════════════════════

TEST(GrandPartition, ExtremeLogZ) {
    GrandPartitionFunction gpf(300.0);

    // Very large Z (ln Z = 500)
    gpf.add_ligand("huge", 500.0);
    // Very small Z (ln Z = -500)
    gpf.add_ligand("tiny", -500.0);

    // Ξ ≈ e^500 (dominated by "huge")
    EXPECT_NEAR(gpf.log_Xi(), 500.0, 1.0);  // within 1.0 of 500

    // p(huge) ≈ 1.0
    EXPECT_NEAR(gpf.binding_probability("huge"), 1.0, 1e-6);
    // p(tiny) ≈ 0.0
    EXPECT_NEAR(gpf.binding_probability("tiny"), 0.0, 1e-6);

    // No NaN or Inf
    EXPECT_TRUE(std::isfinite(gpf.log_Xi()));
    EXPECT_TRUE(std::isfinite(gpf.binding_probability("huge")));
    EXPECT_TRUE(std::isfinite(gpf.binding_probability("tiny")));
    EXPECT_TRUE(std::isfinite(gpf.empty_probability()));
}

// ════════════════════════════════════════════════════════════════════════
// StatMechEngine integration
// ════════════════════════════════════════════════════════════════════════

TEST(GrandPartition, FromStatMechEngine) {
    GrandPartitionFunction gpf(300.0);

    statmech::StatMechEngine engine(300.0);
    engine.add_sample(-5.0, 1);   // E = -5 kcal/mol
    engine.add_sample(-3.0, 2);   // E = -3 kcal/mol (2 counts)
    engine.add_sample(-1.0, 1);

    gpf.add_ligand("from_engine", engine);
    EXPECT_TRUE(gpf.has_ligand("from_engine"));

    auto thermo = engine.compute();
    EXPECT_NEAR(gpf.F_bound("from_engine"), thermo.free_energy, 1e-10);
}

// ════════════════════════════════════════════════════════════════════════
// Concentration-dependent binding
// ════════════════════════════════════════════════════════════════════════

TEST(GrandPartition, ConcentrationAffectsProbability) {
    GrandPartitionFunction gpf(300.0);

    // Same partition function, different concentrations
    // A at 1 M, B at 0.01 M (10 mM)
    gpf.add_ligand("A", 10.0, 1.0);      // ln(z·Z) = ln(1) + 10 = 10
    gpf.add_ligand("B", 10.0, 0.01);     // ln(z·Z) = ln(0.01) + 10 = 10 - ln(100)

    // A should have higher binding probability (higher concentration)
    EXPECT_GT(gpf.binding_probability("A"), gpf.binding_probability("B"));

    // Selectivity A/B should reflect concentration ratio
    // z_A·Z_A / (z_B·Z_B) = (1.0 / 0.01) = 100
    EXPECT_NEAR(gpf.selectivity("A", "B"), 100.0, 0.1);
}

TEST(GrandPartition, DefaultConcentrationIsMolar) {
    GrandPartitionFunction gpf(300.0);

    // Both at default 1 M → should have equal probability
    gpf.add_ligand("X", 5.0);   // default c = 1 M
    gpf.add_ligand("Y", 5.0, 1.0);  // explicit 1 M

    EXPECT_NEAR(gpf.binding_probability("X"), gpf.binding_probability("Y"), 1e-12);
    EXPECT_NEAR(gpf.selectivity("X", "Y"), 1.0, 1e-12);
}

TEST(GrandPartition, InvalidConcentrationThrows) {
    GrandPartitionFunction gpf(300.0);
    EXPECT_THROW(gpf.add_ligand("bad", 5.0, 0.0), std::domain_error);
    EXPECT_THROW(gpf.add_ligand("bad2", 5.0, -1.0), std::domain_error);
}

TEST(GrandPartition, ImpossibleConcentrationThrows) {
    GrandPartitionFunction gpf(300.0);
    EXPECT_THROW(gpf.add_ligand("bad", 5.0, 1e6), std::invalid_argument);
}

// ════════════════════════════════════════════════════════════════════════
// mean_occupancy() and occupancy_variance()
// ════════════════════════════════════════════════════════════════════════

TEST(GrandPartition, MeanOccupancyEmpty) {
    // No ligands → receptor is always apo → ⟨n⟩ = 0
    GrandPartitionFunction gpf(300.0);
    EXPECT_NEAR(gpf.mean_occupancy(), 0.0, 1e-12);
    EXPECT_NEAR(gpf.occupancy_variance(), 0.0, 1e-12);
}

TEST(GrandPartition, MeanOccupancyStrongBinder) {
    // Z ≫ 1 → receptor almost always occupied → ⟨n⟩ ≈ 1
    GrandPartitionFunction gpf(300.0);
    gpf.add_ligand("A", 500.0);  // overwhelming partition function
    EXPECT_NEAR(gpf.mean_occupancy(), 1.0, 1e-6);
    // Variance ≈ 0 at near-saturation
    EXPECT_NEAR(gpf.occupancy_variance(), 0.0, 1e-6);
}

TEST(GrandPartition, MeanOccupancyHalfSaturation) {
    // At half-saturation: log_Z = 0, concentration = 1 M → Z = 1
    // Ξ = 1 + 1 = 2 → p(empty) = 0.5 → ⟨n⟩ = 0.5
    GrandPartitionFunction gpf(300.0);
    gpf.add_ligand("A", 0.0, 1.0);
    EXPECT_NEAR(gpf.mean_occupancy(), 0.5, 1e-12);
    // Variance is maximised at 0.25 when ⟨n⟩ = 0.5
    EXPECT_NEAR(gpf.occupancy_variance(), 0.25, 1e-12);
}

TEST(GrandPartition, MeanOccupancyPlusPEmptyEqualsOne) {
    GrandPartitionFunction gpf(300.0);
    gpf.add_ligand("A", 5.0);
    gpf.add_ligand("B", 3.0);
    EXPECT_NEAR(gpf.mean_occupancy() + gpf.empty_probability(), 1.0, 1e-12);
}

// ════════════════════════════════════════════════════════════════════════
// Log-selectivity (overflow-safe)
// ════════════════════════════════════════════════════════════════════════

TEST(GrandPartition, LogSelectivity) {
    GrandPartitionFunction gpf(300.0);
    gpf.add_ligand("A", 500.0);
    gpf.add_ligand("B", 3.0);

    // ln(Z_A/Z_B) = 500 - 3 = 497
    EXPECT_NEAR(gpf.log_selectivity("A", "B"), 497.0, 1e-10);
    EXPECT_NEAR(gpf.log_selectivity("B", "A"), -497.0, 1e-10);

    // selectivity() returns a finite value (exp(497) is representable)
    EXPECT_GT(gpf.selectivity("A", "B"), 0.0);
    EXPECT_TRUE(std::isfinite(gpf.selectivity("A", "B")));
    EXPECT_NEAR(gpf.selectivity("B", "A"), 0.0, 1e-200);

    // Extreme values: selectivity returns DBL_MAX sentinel (not infinity —
    // avoids UB under -ffast-math/-ffinite-math-only)
    GrandPartitionFunction gpf2(300.0);
    gpf2.add_ligand("X", 800.0);
    gpf2.add_ligand("Y", 0.0);
    EXPECT_EQ(gpf2.selectivity("X", "Y"), std::numeric_limits<double>::max());
    EXPECT_EQ(gpf2.selectivity("Y", "X"), 0.0);
}

TEST(GrandPartition, IntrinsicVsApparentSelectivity) {
    GrandPartitionFunction gpf(300.0);
    gpf.add_ligand("A", 10.0, 1.0);       // 1 M
    gpf.add_ligand("B", 10.0, 0.01);      // 10 mM

    // Intrinsic selectivity (concentration-independent): Z_A/Z_B = 1
    EXPECT_NEAR(gpf.log_intrinsic_selectivity("A", "B"), 0.0, 1e-12);

    // Apparent selectivity (concentration-weighted): z_A·Z_A / (z_B·Z_B) = 100
    EXPECT_NEAR(gpf.log_selectivity("A", "B"), std::log(100.0), 1e-10);
}

// ─────────────────────────────────────────────────────────────────────
// Regression guards for log_intrinsic_selectivity()
//
// These tests exist to prevent a very specific class of "helpful
// refactor" — namely, replacing the implementation with
//   return log_c(a) − log_c(b);
// which would silently compute ln(c_A/c_B) (the concentration ratio)
// instead of ln(Z_A/Z_B) (the intrinsic affinity ratio).
//
// Each test below would fail under that broken implementation.
// ─────────────────────────────────────────────────────────────────────

TEST(GrandPartition, IntrinsicSelectivityIndependentOfConcentration) {
    // Same partition functions, wildly different concentrations.
    // Intrinsic selectivity must be zero for ALL concentration pairs.
    // A log_c-based implementation would return ln(c_A/c_B) = ±N, not 0.
    const double concentrations[] = {1e-9, 1e-6, 1e-3, 1.0, 100.0};
    for (double c_a : concentrations) {
        for (double c_b : concentrations) {
            GrandPartitionFunction gpf(300.0);
            gpf.add_ligand("A", 7.0, c_a);
            gpf.add_ligand("B", 7.0, c_b);
            EXPECT_NEAR(gpf.log_intrinsic_selectivity("A", "B"), 0.0, 1e-12)
                << "c_a=" << c_a << " c_b=" << c_b;
        }
    }
}

TEST(GrandPartition, IntrinsicSelectivityTracksLogZAtAnyConcentration) {
    // Different Z values, equal concentrations: intrinsic = log_Z_A − log_Z_B.
    // Repeat at several concentrations — the value must not shift.
    for (double c : {1e-9, 1e-3, 1.0, 100.0}) {
        GrandPartitionFunction gpf(300.0);
        gpf.add_ligand("A", 12.5, c);
        gpf.add_ligand("B",  4.5, c);
        EXPECT_NEAR(gpf.log_intrinsic_selectivity("A", "B"), 8.0, 1e-12)
            << "c=" << c;
    }
}

TEST(GrandPartition, IntrinsicSelectivityAntisymmetryAndReflexivity) {
    GrandPartitionFunction gpf(300.0);
    gpf.add_ligand("A", 9.0,  0.5);
    gpf.add_ligand("B", 2.0,  0.02);

    // Reflexive: ln(Z_A/Z_A) = 0
    EXPECT_NEAR(gpf.log_intrinsic_selectivity("A", "A"), 0.0, 1e-12);
    EXPECT_NEAR(gpf.log_intrinsic_selectivity("B", "B"), 0.0, 1e-12);

    // Antisymmetric: ln(Z_A/Z_B) = −ln(Z_B/Z_A)
    double s_ab = gpf.log_intrinsic_selectivity("A", "B");
    double s_ba = gpf.log_intrinsic_selectivity("B", "A");
    EXPECT_NEAR(s_ab, -s_ba, 1e-12);
    EXPECT_NEAR(s_ab,  7.0, 1e-12);  // 9 − 2
}

TEST(GrandPartition, IntrinsicAndApparentSelectivityDifferByLogConcentrationRatio) {
    // Algebraic identity: log_selectivity − log_intrinsic_selectivity = ln(c_A/c_B)
    //   apparent = (log_c_A + log_Z_A) − (log_c_B + log_Z_B)
    //   intrinsic = log_Z_A − log_Z_B
    //   apparent − intrinsic = log_c_A − log_c_B = ln(c_A/c_B)
    // This would collapse to 0 under a log_c-based implementation.
    GrandPartitionFunction gpf(300.0);
    gpf.add_ligand("A", 6.5,  2.5);   // 2.5 M
    gpf.add_ligand("B", 1.25, 1e-4);  // 0.1 mM

    double apparent  = gpf.log_selectivity("A", "B");
    double intrinsic = gpf.log_intrinsic_selectivity("A", "B");
    double expected_conc_ratio = std::log(2.5 / 1e-4);

    EXPECT_NEAR(apparent - intrinsic, expected_conc_ratio, 1e-10);
    // And independently: intrinsic must equal log_Z_A − log_Z_B exactly.
    EXPECT_NEAR(intrinsic, 6.5 - 1.25, 1e-12);
}

// ════════════════════════════════════════════════════════════════════════
// F_bound vs delta_G_bind
// ════════════════════════════════════════════════════════════════════════

TEST(GrandPartition, FBoundVsDeltaGBind) {
    GrandPartitionFunction gpf(300.0);
    gpf.add_ligand("A", 10.0);

    // F_bound = −kT ln Z
    double F = gpf.F_bound("A");
    EXPECT_NEAR(F, -kT_300 * 10.0, 1e-10);

    // delta_G_bind with F_ref = 0 → same as F_bound
    EXPECT_NEAR(gpf.delta_G_bind("A", 0.0), F, 1e-10);

    // delta_G_bind with a realistic reference state
    double F_ref = 2.0;  // unbound ligand free energy
    EXPECT_NEAR(gpf.delta_G_bind("A", F_ref), F - F_ref, 1e-10);
}

// ════════════════════════════════════════════════════════════════════════
// add_or_overwrite (atomic insert-or-replace)
// ════════════════════════════════════════════════════════════════════════

TEST(GrandPartition, AddOrOverwriteInsert) {
    GrandPartitionFunction gpf(300.0);
    gpf.add_or_overwrite("A", 5.0);
    EXPECT_TRUE(gpf.has_ligand("A"));
    EXPECT_NEAR(gpf.F_bound("A"), -kT_300 * 5.0, 1e-10);
}

TEST(GrandPartition, AddOrOverwriteReplace) {
    GrandPartitionFunction gpf(300.0);
    gpf.add_or_overwrite("A", 5.0);
    EXPECT_NEAR(gpf.F_bound("A"), -kT_300 * 5.0, 1e-10);

    // Overwrite with new value
    gpf.add_or_overwrite("A", 10.0);
    EXPECT_NEAR(gpf.F_bound("A"), -kT_300 * 10.0, 1e-10);
    EXPECT_EQ(gpf.num_ligands(), 1);  // still just one ligand
}

TEST(GrandPartition, AddOrOverwriteConcentration) {
    GrandPartitionFunction gpf(300.0);
    gpf.add_or_overwrite("A", 10.0, 1.0);   // 1 M
    gpf.add_or_overwrite("B", 10.0, 0.01);  // 10 mM

    // A should have higher binding probability
    EXPECT_GT(gpf.binding_probability("A"), gpf.binding_probability("B"));

    // Now overwrite B with higher concentration
    gpf.add_or_overwrite("B", 10.0, 1.0);  // upgrade to 1 M
    EXPECT_NEAR(gpf.binding_probability("A"), gpf.binding_probability("B"), 1e-12);
}

// ════════════════════════════════════════════════════════════════════════
// Log-Xi caching (dirty flag)
// ════════════════════════════════════════════════════════════════════════

TEST(GrandPartition, CacheConsistencyAfterMutation) {
    GrandPartitionFunction gpf(300.0);
    gpf.add_ligand("A", 5.0);

    double log_xi_1 = gpf.log_Xi();
    double pA_1 = gpf.binding_probability("A");

    // Mutate: add another ligand
    gpf.add_ligand("B", 8.0);

    // Cache should have been invalidated — new values must differ
    double log_xi_2 = gpf.log_Xi();
    EXPECT_GT(log_xi_2, log_xi_1);  // Ξ grows with more ligands

    double pA_2 = gpf.binding_probability("A");
    EXPECT_LT(pA_2, pA_1);  // A's share decreases with competition
}

TEST(GrandPartition, CacheConsistencyAfterOverwrite) {
    GrandPartitionFunction gpf(300.0);
    gpf.add_ligand("A", 5.0);
    double pA_before = gpf.binding_probability("A");

    gpf.overwrite_ligand("A", 20.0);  // much stronger
    double pA_after = gpf.binding_probability("A");
    EXPECT_GT(pA_after, pA_before);
}

TEST(GrandPartition, CacheConsistencyAfterMerge) {
    GrandPartitionFunction gpf(300.0);
    gpf.add_ligand("A", 5.0);
    double dG_before = gpf.F_bound("A");

    gpf.merge_ligand("A", 5.0);  // double the ensemble
    double dG_after = gpf.F_bound("A");
    EXPECT_LT(dG_after, dG_before);  // more favorable after merge
}

TEST(GrandPartition, CacheConsistencyAfterRemove) {
    GrandPartitionFunction gpf(300.0);
    gpf.add_ligand("A", 10.0);
    gpf.add_ligand("B", 5.0);
    double pA_with_B = gpf.binding_probability("A");

    gpf.remove_ligand("B");
    double pA_alone = gpf.binding_probability("A");
    EXPECT_GT(pA_alone, pA_with_B);  // A's probability increases without B
}

TEST(GrandPartition, RepeatedQueriesUseCache) {
    // Verify that repeated queries return identical results (cache doesn't drift)
    GrandPartitionFunction gpf(300.0);
    gpf.add_ligand("A", 10.0);
    gpf.add_ligand("B", 5.0);

    double log_xi_1 = gpf.log_Xi();
    double log_xi_2 = gpf.log_Xi();
    EXPECT_EQ(log_xi_1, log_xi_2);  // bitwise identical from cache

    double p1 = gpf.binding_probability("A");
    double p2 = gpf.binding_probability("A");
    EXPECT_EQ(p1, p2);
}

// ════════════════════════════════════════════════════════════════════════
// all_log_zZ accessor
// ════════════════════════════════════════════════════════════════════════

TEST(GrandPartition, AllLogZzAccessor) {
    GrandPartitionFunction gpf(300.0);
    gpf.add_ligand("A", 10.0, 1.0);    // log_zZ = 0 + 10 = 10
    gpf.add_ligand("B", 10.0, 0.01);   // log_zZ = ln(0.01) + 10

    auto zz = gpf.all_log_zZ();
    ASSERT_EQ(zz.size(), 2u);

    // Find entries (order not guaranteed)
    double log_zz_A = 0, log_zz_B = 0;
    for (const auto& [name, val] : zz) {
        if (name == "A") log_zz_A = val;
        if (name == "B") log_zz_B = val;
    }
    EXPECT_NEAR(log_zz_A, 10.0, 1e-10);
    EXPECT_NEAR(log_zz_B, std::log(0.01) + 10.0, 1e-10);
}

TEST(GrandPartition, AllLogZAccessor) {
    GrandPartitionFunction gpf(300.0);
    gpf.add_ligand("A", 10.0, 1.0);
    gpf.add_ligand("B", 20.0, 0.5);

    auto zz = gpf.all_log_Z();  // intrinsic log(Z), not concentration-weighted
    ASSERT_EQ(zz.size(), 2u);

    double logZ_A = 0, logZ_B = 0;
    for (const auto& [name, val] : zz) {
        if (name == "A") logZ_A = val;
        if (name == "B") logZ_B = val;
    }
    // all_log_Z returns intrinsic partition function (no concentration weighting)
    EXPECT_NEAR(logZ_A, 10.0, 1e-10);
    EXPECT_NEAR(logZ_B, 20.0, 1e-10);
}

// ════════════════════════════════════════════════════════════════════════
// Performance / stress benchmark for Grand Canonical implementation
// Run with: build/test_grand_partition --gtest_filter="*BenchmarkStress*"
// Reports wall-clock timings for core operations under load.
// These are not strict assertions (timings are environment-dependent) but
// exercise scaling, cache behavior, and numerical paths at volume.
// ════════════════════════════════════════════════════════════════════════

#include <chrono>
#include <iostream>
#include <random>
#include <vector>

TEST(GrandPartition, BenchmarkStress) {
    using clock = std::chrono::high_resolution_clock;
    auto t0 = clock::now();

    const int N = 200;  // scale: 200 competing ligands
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> logz_dist(-30.0, 120.0);
    std::uniform_real_distribution<double> conc_dist(1e-6, 10.0);

    GrandPartitionFunction gpf(310.0);  // physiological-ish temp

    auto t_add_start = clock::now();
    for (int i = 0; i < N; ++i) {
        std::string name = "L" + std::to_string(i);
        double lz = logz_dist(rng);
        double c = conc_dist(rng);
        gpf.add_ligand(name, lz, c);
    }
    auto t_add = std::chrono::duration<double, std::milli>(clock::now() - t_add_start).count();

    // Force cache population
    auto t_xi_start = clock::now();
    double lx = gpf.log_Xi();
    auto t_xi = std::chrono::duration<double, std::milli>(clock::now() - t_xi_start).count();

    auto t_prob_start = clock::now();
    double sum_p = gpf.empty_probability();
    std::vector<double> probs;
    probs.reserve(N);
    for (int i = 0; i < N; ++i) {
        std::string name = "L" + std::to_string(i);
        double p = gpf.binding_probability(name);
        probs.push_back(p);
        sum_p += p;
    }
    auto t_prob = std::chrono::duration<double, std::milli>(clock::now() - t_prob_start).count();

    // Verify thermodynamic closure under load
    EXPECT_NEAR(sum_p, 1.0, 1e-8);

    auto t_rank_start = clock::now();
    auto ranks = gpf.rank();
    auto t_rank = std::chrono::duration<double, std::milli>(clock::now() - t_rank_start).count();
    ASSERT_EQ(ranks.size(), static_cast<size_t>(N));

    // Repeated queries should be near-instant (cache hit)
    auto t_cached_start = clock::now();
    for (int r = 0; r < 100; ++r) {
        (void)gpf.log_Xi();
        (void)gpf.binding_probability("L0");
        (void)gpf.mean_occupancy();
    }
    auto t_cached = std::chrono::duration<double, std::milli>(clock::now() - t_cached_start).count();

    // MultiSite stress
    auto t_ms_start = clock::now();
    MultiSiteGPF msgpf(300.0);
    for (int s = 0; s < 4; ++s) {
        msgpf.add_site("site_" + std::to_string(s));
    }
    for (int s = 0; s < 4; ++s) {
        for (int i = 0; i < 30; ++i) {  // 120 total ligand registrations
            double lz = logz_dist(rng) * 0.6;
            msgpf.add_ligand(s, "M" + std::to_string(i), lz, 1.0);
        }
    }
    msgpf.set_cooperativity(0, 1, 2.5);
    msgpf.set_cooperativity(2, 3, 0.4);
    double mxi = msgpf.log_Xi();
    double p0 = msgpf.binding_probability(0, "M5");
    auto t_ms = std::chrono::duration<double, std::milli>(clock::now() - t_ms_start).count();
    EXPECT_TRUE(std::isfinite(mxi));
    EXPECT_GE(p0, 0.0);

    auto total_ms = std::chrono::duration<double, std::milli>(clock::now() - t0).count();

    std::cout << "\n[GrandPartition BenchmarkStress]\n"
              << "  Ligands: " << N << "\n"
              << "  add_ligand (200): " << t_add << " ms\n"
              << "  first log_Xi: " << t_xi << " ms\n"
              << "  binding probs + empty (200): " << t_prob << " ms\n"
              << "  rank(): " << t_rank << " ms\n"
              << "  100x cached queries: " << t_cached << " ms\n"
              << "  MultiSite (4 sites, ~120 regs + 2 coop + queries): " << t_ms << " ms\n"
              << "  TOTAL: " << total_ms << " ms\n";

    // Soft sanity: operations should complete in reasonable time even for N=200
    // (these are generous to avoid flaky CI on slow CI runners)
    EXPECT_LT(t_add, 500.0);
    EXPECT_LT(t_rank, 200.0);
    EXPECT_LT(t_cached, 50.0);
}
