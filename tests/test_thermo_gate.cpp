// test_thermo_gate.cpp — thermodynamic impossibility gate (LP whiteboard IMG_3696)
//
// Physics: ΔG = ΔH − TΔS. When ΔH > 0 and ΔS < 0, −TΔS > 0 for every T > 0,
// so ΔG is strictly positive at all temperatures — the pose can never bind
// spontaneously. Such poses are penalised to a large positive sentinel.
//
// Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

#include <gtest/gtest.h>
#include "../LIB/ThermodynamicEngine.h"

using thermo_gate::is_impossible;
using thermo_gate::apply_gate;
using thermo_gate::kImpossibleSentinel;

// ── The two cases named in the spec ─────────────────────────────────────────

TEST(ThermoGate, PositiveEnthalpyNegativeEntropyIsImpossible) {
    // dH = +5, dS = −0.01 → both conditions hold → impossible.
    bool flagged = false;
    const float dG = apply_gate(/*dG_eff=*/-42.0f, /*dH=*/5.0f, /*dS=*/-0.01f,
                                &flagged);
    EXPECT_TRUE(flagged);
    EXPECT_GE(dG, 999.0f);
    EXPECT_FLOAT_EQ(dG, kImpossibleSentinel);
    EXPECT_TRUE(is_impossible(5.0f, -0.01f));
}

TEST(ThermoGate, NegativeEnthalpyAloneIsNotFlagged) {
    // dH = −5, dS = −0.01 → only dS < 0; enthalpy is favourable → allowed.
    bool flagged = true;
    const float dG = apply_gate(/*dG_eff=*/-42.0f, /*dH=*/-5.0f, /*dS=*/-0.01f,
                                &flagged);
    EXPECT_FALSE(flagged);
    EXPECT_FLOAT_EQ(dG, -42.0f);   // passed through unchanged
    EXPECT_FALSE(is_impossible(-5.0f, -0.01f));
}

// ── The other two quadrants: the gate is strictly two-sided ─────────────────

TEST(ThermoGate, PositiveEntropyIsNeverFlagged) {
    // dH = +5 but dS = +0.01 → entropy-driven binding is possible at high T.
    EXPECT_FALSE(is_impossible(5.0f, 0.01f));
    bool flagged = true;
    EXPECT_FLOAT_EQ(apply_gate(-42.0f, 5.0f, 0.01f, &flagged), -42.0f);
    EXPECT_FALSE(flagged);
}

TEST(ThermoGate, BothFavourableIsNeverFlagged) {
    EXPECT_FALSE(is_impossible(-5.0f, 0.01f));
}

// ── Boundary behaviour: the comparisons are strict (> 0, < 0) ───────────────

TEST(ThermoGate, ZeroBoundariesAreNotImpossible) {
    EXPECT_FALSE(is_impossible(0.0f, -0.01f));  // dH == 0 is not > 0
    EXPECT_FALSE(is_impossible(5.0f, 0.0f));    // dS == 0 is not < 0
    EXPECT_FALSE(is_impossible(0.0f, 0.0f));
}

// ── Guard against the dead-code trap documented in ThermodynamicEngine.h ────
// Shannon entropy H = −Σ P_i·ln P_i is ≥ 0 by construction, so wiring H in as
// the ΔS source would make the gate unreachable. This pins that reasoning.

TEST(ThermoGate, ShannonEntropySignCanNeverTripTheGate) {
    // Any H produced by a real probability distribution is >= 0.
    for (float H : {0.0f, 0.5f, 2.009579f, 3.508342f, 3.912023f}) {
        EXPECT_FALSE(is_impossible(5.0f, H))
            << "H=" << H << " must not trip the gate; Shannon entropy is >= 0";
    }
}

// ── The sentinel must dominate any realistic ΔG_eff so clustering deranks it ─

TEST(ThermoGate, SentinelDominatesRealisticScores) {
    // Observed ΔG_eff on 1SG0/2GBP/1OF1 was roughly -63 .. -154.
    for (float dG : {-63.296303f, -118.405441f, -154.344559f}) {
        bool flagged = false;
        EXPECT_GT(apply_gate(dG, 5.0f, -0.01f, &flagged), dG);
        EXPECT_TRUE(flagged);
    }
}
