// tests/test_metal_coordination.cpp
// Unit tests for the metal ion coordination scoring potential
// Apache-2.0 (c) 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include "../LIB/metal_coordination.h"
#include <cmath>
#include <cstring>
#include <vector>

// Undefine the single-letter macro from flexaid.h to avoid name collisions
#ifdef E
#undef E
#endif

using namespace metal_coord;

static constexpr double EPS = 1e-6;

// ─── Helper: create a minimal atom with the given SYBYL type ────────────────
static atom_struct make_atom(int sybyl_type, float x = 0.f, float y = 0.f,
                              float z = 0.f) {
    atom_struct a{};
    a.type = sybyl_type;
    a.coor[0] = x;
    a.coor[1] = y;
    a.coor[2] = z;
    a.charge = 0.0f;
    a.type256 = 0;
    a.bond[0] = 0;
    std::memset(a.name, 0, sizeof(a.name));
    std::memset(a.element, 0, sizeof(a.element));
    return a;
}

// ===========================================================================
// GAUSSIAN WELL POTENTIAL — UNIT TESTS
// ===========================================================================

TEST(GaussianWell, MinimumAtIdealDistance) {
    // At r = r0, Gaussian well = well_depth
    double e = gaussian_well(2.36, 2.36, -5.0, 0.45);
    EXPECT_NEAR(e, -5.0, EPS);
}

TEST(GaussianWell, AsymptoteAtLargeDistance) {
    // At large r, Gaussian well -> 0
    double e = gaussian_well(10.0, 2.36, -5.0, 0.45);
    EXPECT_NEAR(e, 0.0, 0.001);
}

TEST(GaussianWell, NoRepulsionAtShortDistance) {
    // At r << r0, Gaussian well -> 0 (NOT repulsive, unlike Morse)
    double e = gaussian_well(0.5, 2.36, -5.0, 0.45);
    EXPECT_NEAR(e, 0.0, 0.01);
    EXPECT_LE(e, 0.0);  // Always <= 0 (or zero) for negative well_depth
}

TEST(GaussianWell, SymmetricWellShape) {
    // Verify that the well is deeper at r0 than at r0 +/- delta
    double e_center = gaussian_well(2.36, 2.36, -5.0, 0.45);
    double e_plus   = gaussian_well(2.36 + 0.3, 2.36, -5.0, 0.45);
    double e_minus  = gaussian_well(2.36 - 0.3, 2.36, -5.0, 0.45);
    EXPECT_LT(e_center, e_plus);
    EXPECT_LT(e_center, e_minus);
}

TEST(GaussianWell, SigmaControlsWidth) {
    // Larger sigma -> broader well -> more energy at same displacement
    double e_narrow = gaussian_well(2.36 + 0.3, 2.36, -5.0, 0.30);
    double e_wide   = gaussian_well(2.36 + 0.3, 2.36, -5.0, 0.60);
    EXPECT_LT(e_wide, e_narrow);  // wider well = more negative (deeper) at +0.3
}

// ===========================================================================
// METAL PARAMETER LOOKUP
// ===========================================================================

TEST(MetalParams, Ca2PlusRecognized) {
    const MetalParams* p = get_metal_params(36);
    ASSERT_NE(p, nullptr);
    EXPECT_EQ(p->ideal_cn, 7);
    EXPECT_EQ(p->geometry, Geometry::PENTAGONAL_BIPYRAMIDAL);
    EXPECT_NEAR(p->angle_primary, 72.0, EPS);
}

TEST(MetalParams, Zn2PlusRecognized) {
    const MetalParams* p = get_metal_params(35);
    ASSERT_NE(p, nullptr);
    EXPECT_EQ(p->ideal_cn, 4);
    EXPECT_EQ(p->geometry, Geometry::TETRAHEDRAL);
}

TEST(MetalParams, Mg2PlusRecognized) {
    const MetalParams* p = get_metal_params(28);
    ASSERT_NE(p, nullptr);
    EXPECT_EQ(p->ideal_cn, 6);
    EXPECT_EQ(p->geometry, Geometry::OCTAHEDRAL);
}

TEST(MetalParams, NonMetalReturnsNull) {
    EXPECT_EQ(get_metal_params(1), nullptr);   // C.1
    EXPECT_EQ(get_metal_params(14), nullptr);  // O.3
    EXPECT_EQ(get_metal_params(0), nullptr);   // out of range
    EXPECT_EQ(get_metal_params(99), nullptr);  // out of range
}

TEST(MetalParams, PlaceholderMetalReturnsNull) {
    EXPECT_EQ(get_metal_params(29), nullptr);  // Sr — no params
}

// ===========================================================================
// DONOR AFFINITY LOOKUP — with scientific corrections
// ===========================================================================

TEST(DonorAffinity, Ca2PlusCarboxylateO) {
    auto aff = get_donor_affinity(36, 15); // Ca2+ - O.CO2
    ASSERT_TRUE(aff.has_value());
    EXPECT_NEAR(aff->ideal_dist, 2.36, EPS);
    EXPECT_NEAR(aff->well_depth, -5.0, EPS);
}

TEST(DonorAffinity, Ca2PlusCarbonylO) {
    auto aff = get_donor_affinity(36, 13); // Ca2+ - O.2
    ASSERT_TRUE(aff.has_value());
    EXPECT_NEAR(aff->ideal_dist, 2.38, EPS);  // corrected from 2.36
    EXPECT_NEAR(aff->well_depth, -4.0, EPS);
}

TEST(DonorAffinity, Ca2PlusHydroxylO_CorrectedDistance) {
    auto aff = get_donor_affinity(36, 14); // Ca2+ - O.3
    ASSERT_TRUE(aff.has_value());
    EXPECT_NEAR(aff->ideal_dist, 2.43, EPS);  // corrected from 2.38 (Harding CSD)
    EXPECT_NEAR(aff->well_depth, -3.5, EPS);
}

TEST(DonorAffinity, N4QuaternaryExcluded) {
    // N.4 (quaternary ammonium, SYBYL 9) has NO lone pair and CANNOT
    // coordinate metals — must return nullopt
    auto aff = get_donor_affinity(36, 9);  // Ca2+ - N.4
    EXPECT_FALSE(aff.has_value()) << "N.4 must not coordinate metals (no lone pair)";
    auto aff2 = get_donor_affinity(35, 9); // Zn2+ - N.4
    EXPECT_FALSE(aff2.has_value());
}

TEST(DonorAffinity, NAR_FullStrength) {
    auto aff = get_donor_affinity(35, 10); // Zn2+ - N.AR (His imidazole)
    ASSERT_TRUE(aff.has_value());
    // N.AR gets full strength (n_strength = 1.0)
    EXPECT_NEAR(aff->well_depth, -5.0, EPS);  // -5.0 * 1.0
}

TEST(DonorAffinity, NAM_WeakStrength) {
    auto aff = get_donor_affinity(35, 11); // Zn2+ - N.AM (amide)
    ASSERT_TRUE(aff.has_value());
    // N.AM gets 20% strength (n_strength = 0.2)
    EXPECT_NEAR(aff->well_depth, -1.0, EPS);  // -5.0 * 0.2
}

TEST(DonorAffinity, Zn2PlusSulfur) {
    auto aff = get_donor_affinity(35, 18); // Zn2+ - S.3
    ASSERT_TRUE(aff.has_value());
    EXPECT_NEAR(aff->ideal_dist, 2.30, EPS);
    EXPECT_NEAR(aff->well_depth, -6.0, EPS);  // reduced from -16.0
}

TEST(DonorAffinity, SulfoxideExcluded) {
    // S.O (19), S.O2 (20), S.AR (21) should return nullopt
    EXPECT_FALSE(get_donor_affinity(36, 19).has_value()) << "S.O should not coordinate";
    EXPECT_FALSE(get_donor_affinity(35, 20).has_value()) << "S.O2 should not coordinate";
    EXPECT_FALSE(get_donor_affinity(35, 21).has_value()) << "S.AR should not coordinate";
}

TEST(DonorAffinity, UnknownMetalDonorPair) {
    // Ca2+ with a halogen (no affinity defined)
    EXPECT_FALSE(get_donor_affinity(36, 23).has_value()); // Ca2+ - F
    EXPECT_FALSE(get_donor_affinity(36, 24).has_value()); // Ca2+ - Cl
}

TEST(DonorAffinity, NonMetalReturnsNullopt) {
    EXPECT_FALSE(get_donor_affinity(1, 14).has_value());  // C.1 is not a metal
    EXPECT_FALSE(get_donor_affinity(14, 36).has_value()); // O.3 is not a metal
}

// ===========================================================================
// TYPE CLASSIFICATION
// ===========================================================================

TEST(TypeClassification, MetalTypes) {
    EXPECT_TRUE(is_metal_type(36));  // Ca
    EXPECT_TRUE(is_metal_type(35));  // Zn
    EXPECT_TRUE(is_metal_type(28));  // Mg
    EXPECT_TRUE(is_metal_type(37));  // Fe
    EXPECT_FALSE(is_metal_type(1));  // C.1
    EXPECT_FALSE(is_metal_type(14)); // O.3
    EXPECT_FALSE(is_metal_type(29)); // Sr (placeholder)
}

TEST(TypeClassification, DonorTypes) {
    EXPECT_TRUE(is_coord_donor_type(13));   // O.2
    EXPECT_TRUE(is_coord_donor_type(14));   // O.3
    EXPECT_TRUE(is_coord_donor_type(15));   // O.CO2
    EXPECT_TRUE(is_coord_donor_type(18));   // S.3
    EXPECT_TRUE(is_coord_donor_type(6));    // N.1
    EXPECT_TRUE(is_coord_donor_type(22));   // P.3
    EXPECT_FALSE(is_coord_donor_type(9));   // N.4 — EXCLUDED (no lone pair)
    EXPECT_FALSE(is_coord_donor_type(1));   // C.1
    EXPECT_FALSE(is_coord_donor_type(3));   // C.3
    EXPECT_FALSE(is_coord_donor_type(36));  // Ca (metal, not donor)
    EXPECT_FALSE(is_coord_donor_type(19));  // S.O — excluded
    EXPECT_FALSE(is_coord_donor_type(20));  // S.O2 — excluded
    EXPECT_FALSE(is_coord_donor_type(21));  // S.AR — excluded
}

// ===========================================================================
// FULL SCORING FUNCTION (compute_metal_coord_energy)
// ===========================================================================

TEST(MetalCoordEnergy, Ca2PlusO_CO2_IdealDistance) {
    atom_struct atoms[2];
    atoms[0] = make_atom(36, 0.f, 0.f, 0.f);       // Ca2+
    atoms[1] = make_atom(15, 2.36f, 0.f, 0.f);      // O.CO2 at ideal distance
    double e = compute_metal_coord_energy(atoms, 0, 1, 2.36, 1.0, 0.45);
    // At ideal distance: E = weight * well_depth = 1.0 * -5.0
    EXPECT_NEAR(e, -5.0, EPS);
}

TEST(MetalCoordEnergy, Ca2PlusO_CO2_FarAway) {
    atom_struct atoms[2];
    atoms[0] = make_atom(36, 0.f, 0.f, 0.f);
    atoms[1] = make_atom(15, 8.0f, 0.f, 0.f);
    double e = compute_metal_coord_energy(atoms, 0, 1, 8.0, 1.0, 0.45);
    EXPECT_NEAR(e, 0.0, 0.001);
}

TEST(MetalCoordEnergy, NoRepulsionAtShortDistance) {
    // Gaussian well produces NO repulsion (unlike old Morse)
    atom_struct atoms[2];
    atoms[0] = make_atom(36, 0.f, 0.f, 0.f);
    atoms[1] = make_atom(15, 1.0f, 0.f, 0.f);
    double e = compute_metal_coord_energy(atoms, 0, 1, 1.0, 1.0, 0.45);
    EXPECT_LE(e, 0.0) << "Gaussian well should never be repulsive";
}

TEST(MetalCoordEnergy, ReverseAtomOrder) {
    // Should work regardless of which atom is metal vs donor
    atom_struct atoms[2];
    atoms[0] = make_atom(15, 2.36f, 0.f, 0.f);      // O.CO2
    atoms[1] = make_atom(36, 0.f, 0.f, 0.f);         // Ca2+
    double e = compute_metal_coord_energy(atoms, 0, 1, 2.36, 1.0, 0.45);
    EXPECT_NEAR(e, -5.0, EPS);
}

TEST(MetalCoordEnergy, NonMetalPairReturnsZero) {
    atom_struct atoms[2];
    atoms[0] = make_atom(3, 0.f, 0.f, 0.f);   // C.3
    atoms[1] = make_atom(14, 2.0f, 0.f, 0.f);  // O.3
    double e = compute_metal_coord_energy(atoms, 0, 1, 2.0, 1.0, 0.45);
    EXPECT_EQ(e, 0.0);
}

TEST(MetalCoordEnergy, MetalMetalReturnsZero) {
    atom_struct atoms[2];
    atoms[0] = make_atom(36, 0.f, 0.f, 0.f);  // Ca2+
    atoms[1] = make_atom(35, 3.0f, 0.f, 0.f);  // Zn2+
    double e = compute_metal_coord_energy(atoms, 0, 1, 3.0, 1.0, 0.45);
    EXPECT_EQ(e, 0.0);
}

TEST(MetalCoordEnergy, WeightScaling) {
    atom_struct atoms[2];
    atoms[0] = make_atom(36, 0.f, 0.f, 0.f);
    atoms[1] = make_atom(15, 2.36f, 0.f, 0.f);
    double e1 = compute_metal_coord_energy(atoms, 0, 1, 2.36, 1.0, 0.45);
    double e2 = compute_metal_coord_energy(atoms, 0, 1, 2.36, 2.0, 0.45);
    EXPECT_NEAR(e2, 2.0 * e1, EPS);
}

TEST(MetalCoordEnergy, Zn2PlusS3) {
    atom_struct atoms[2];
    atoms[0] = make_atom(35, 0.f, 0.f, 0.f);       // Zn2+
    atoms[1] = make_atom(18, 2.30f, 0.f, 0.f);      // S.3
    double e = compute_metal_coord_energy(atoms, 0, 1, 2.30, 1.0, 0.45);
    EXPECT_NEAR(e, -6.0, EPS);
}

TEST(MetalCoordEnergy, Mg2PlusO_CO2) {
    atom_struct atoms[2];
    atoms[0] = make_atom(28, 0.f, 0.f, 0.f);       // Mg2+
    atoms[1] = make_atom(15, 2.06f, 0.f, 0.f);      // O.CO2
    double e = compute_metal_coord_energy(atoms, 0, 1, 2.06, 1.0, 0.45);
    EXPECT_NEAR(e, -4.5, EPS);
}

TEST(MetalCoordEnergy, N4QuaternaryReturnsZero) {
    atom_struct atoms[2];
    atoms[0] = make_atom(35, 0.f, 0.f, 0.f);       // Zn2+
    atoms[1] = make_atom(9, 2.05f, 0.f, 0.f);       // N.4
    double e = compute_metal_coord_energy(atoms, 0, 1, 2.05, 1.0, 0.45);
    EXPECT_EQ(e, 0.0) << "N.4 should not produce metal coordination energy";
}

TEST(MetalCoordEnergy, SulfoxideReturnsZero) {
    atom_struct atoms[2];
    atoms[0] = make_atom(35, 0.f, 0.f, 0.f);       // Zn2+
    atoms[1] = make_atom(19, 2.30f, 0.f, 0.f);      // S.O
    double e = compute_metal_coord_energy(atoms, 0, 1, 2.30, 1.0, 0.45);
    EXPECT_EQ(e, 0.0) << "S.O should not coordinate metals";
}

TEST(MetalCoordEnergy, DistanceCutoff) {
    atom_struct atoms[2];
    atoms[0] = make_atom(36, 0.f, 0.f, 0.f);       // Ca2+
    atoms[1] = make_atom(15, 5.0f, 0.f, 0.f);       // O.CO2 far away
    // ideal_dist=2.36, cutoff = 2.36 + 2.5*0.45 = 3.485
    double e = compute_metal_coord_energy(atoms, 0, 1, 5.0, 1.0, 0.45);
    EXPECT_EQ(e, 0.0) << "Should return 0 beyond cutoff";
}

// ===========================================================================
// COORDINATION NUMBER PENALTY — CORRECTED SIGN
// ===========================================================================

TEST(CNPenalty, IdealCNZeroPenalty) {
    EXPECT_NEAR(cn_penalty(7, 7, 0.5), 0.0, EPS);
    EXPECT_NEAR(cn_penalty(4, 4, 0.5), 0.0, EPS);
    EXPECT_NEAR(cn_penalty(6, 6, 0.5), 0.0, EPS);
}

TEST(CNPenalty, OneOffPenalty_Positive) {
    // delta=1 -> 0.5 * 1 = +0.5 (POSITIVE = unfavorable)
    double pen = cn_penalty(6, 7, 0.5);
    EXPECT_NEAR(pen, +0.5, EPS);
    EXPECT_GT(pen, 0.0) << "CN penalty must be POSITIVE (unfavorable)";
}

TEST(CNPenalty, TwoOffPenalty_Positive) {
    // delta=2 -> 0.5 * 4 = +2.0
    double pen = cn_penalty(5, 7, 0.5);
    EXPECT_NEAR(pen, +2.0, EPS);
}

TEST(CNPenalty, SymmetricPenalty) {
    // Under-coordination and over-coordination penalized equally
    EXPECT_EQ(cn_penalty(5, 7, 0.5), cn_penalty(9, 7, 0.5));
}

TEST(CNPenalty, WeightScaling) {
    double pen1 = cn_penalty(3, 7, 0.5);   // delta=4 -> 0.5*16 = 8.0
    double pen2 = cn_penalty(3, 7, 1.0);   // delta=4 -> 1.0*16 = 16.0
    EXPECT_NEAR(pen2, 2.0 * pen1, EPS);
}

// ===========================================================================
// COORDINATION CUTOFF
// ===========================================================================

TEST(CoordCutoff, WithinCutoff) {
    EXPECT_TRUE(is_coordinating(2.36, 2.36));
    EXPECT_TRUE(is_coordinating(3.0, 2.36));   // 2.36 + 0.8 = 3.16
}

TEST(CoordCutoff, OutsideCutoff) {
    EXPECT_FALSE(is_coordinating(3.5, 2.36));  // > 3.16
}

// ===========================================================================
// MULTIPLE METALS — SIMULTANEOUS SCORING
// ===========================================================================

TEST(MetalCoordEnergy, FeNitrogen_HisImidazole) {
    atom_struct atoms[2];
    atoms[0] = make_atom(37, 0.f, 0.f, 0.f);       // Fe
    atoms[1] = make_atom(10, 2.15f, 0.f, 0.f);      // N.AR (His imidazole)
    double e = compute_metal_coord_energy(atoms, 0, 1, 2.15, 1.0, 0.45);
    // Fe N.AR: -5.0 * n_strength(10) = -5.0 * 1.0 = -5.0
    EXPECT_NEAR(e, -5.0, EPS);
}

TEST(MetalCoordEnergy, Cu2PlusSulfur) {
    atom_struct atoms[2];
    atoms[0] = make_atom(30, 0.f, 0.f, 0.f);       // Cu2+
    atoms[1] = make_atom(18, 2.15f, 0.f, 0.f);      // S.3
    double e = compute_metal_coord_energy(atoms, 0, 1, 2.15, 1.0, 0.45);
    EXPECT_NEAR(e, -5.5, EPS);
}

TEST(MetalCoordEnergy, UnknownDonorForMetal) {
    atom_struct atoms[2];
    atoms[0] = make_atom(36, 0.f, 0.f, 0.f);       // Ca2+
    atoms[1] = make_atom(25, 3.0f, 0.f, 0.f);       // BR — no affinity
    double e = compute_metal_coord_energy(atoms, 0, 1, 3.0, 1.0, 0.45);
    EXPECT_EQ(e, 0.0);
}
