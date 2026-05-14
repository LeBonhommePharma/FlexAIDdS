// test_tencom_entropy_diff.cpp — Unit tests for TENCoM entropy differential
//
// Tests:
//   1. build_from_ca() API on TorsionalENM
//   2. Eigenvalue differential computation between two structures
//   3. Eigenvector overlap calculation
//   4. Vibrational entropy differential (ΔS_vib) via ENCoM
//   5. Self-differential is zero (identity test)
//   6. Perturbation produces non-zero differential

#include <gtest/gtest.h>
#include "tENCoM/tencm.h"
#include "encom.h"
#include "statmech.h"

#include <array>
#include <cmath>
#include <numbers>
#include <vector>

// ─── Helper: generate an ideal α-helix Cα trace ────────────────────────────

static std::vector<std::array<float,3>> make_helix(
    int n_residues,
    float rise    = 1.5f,
    float radius  = 2.3f,
    float turn_deg = 100.0f)
{
    std::vector<std::array<float,3>> ca;
    ca.reserve(n_residues);
    const float turn_rad = turn_deg * std::numbers::pi_v<float> / 180.0f;
    for (int r = 0; r < n_residues; ++r) {
        float angle = r * turn_rad;
        ca.push_back({
            radius * std::cos(angle),
            radius * std::sin(angle),
            r * rise
        });
    }
    return ca;
}

// Perturb a helix: uniformly shift all Cα by a small random displacement
static std::vector<std::array<float,3>> perturb_helix(
    const std::vector<std::array<float,3>>& ca,
    float amplitude = 0.5f,
    unsigned seed = 123)
{
    std::mt19937 rng(seed);
    std::normal_distribution<float> gauss(0.0f, amplitude);
    auto out = ca;
    for (auto& c : out) {
        c[0] += gauss(rng);
        c[1] += gauss(rng);
        c[2] += gauss(rng);
    }
    return out;
}

// ─── Vibrational entropy from TENCoM modes (same logic as standalone) ───────

static encom::VibrationalEntropy tencom_svib(
    const std::vector<tencm::NormalMode>& modes,
    double T,
    const encom::FrequencyCalibration& calibration =
        encom::FrequencyCalibration::model_scale(),
    int skip = 6)
{
    std::vector<encom::NormalMode> em;
    for (int m = skip; m < static_cast<int>(modes.size()); ++m) {
        if (modes[m].eigenvalue < 1e-8) continue;
        encom::NormalMode e;
        e.index = m + 1;
        e.eigenvalue = modes[m].eigenvalue;
        e.frequency  = std::sqrt(std::abs(e.eigenvalue));
        em.push_back(e);
    }
    if (em.empty()) {
        encom::VibrationalEntropy z{};
        z.temperature = T;
        z.eigenvalue_to_omega = calibration.eigenvalue_to_omega;
        z.calibrated = calibration.calibrated;
        z.calibration_label = calibration.label;
        z.calibration_provenance = calibration.provenance;
        return z;
    }
    return encom::ENCoMEngine::compute_vibrational_entropy(em, T, calibration);
}

// ─── Tests ──────────────────────────────────────────────────────────────────

class TENCoMEntropyDiffTest : public ::testing::Test {
protected:
    void SetUp() override {
        ca_ref_ = make_helix(50);
        ca_pert_ = perturb_helix(ca_ref_, 0.3f);
    }

    std::vector<std::array<float,3>> ca_ref_;
    std::vector<std::array<float,3>> ca_pert_;
};

TEST_F(TENCoMEntropyDiffTest, BuildFromCaWorks) {
    tencm::TorsionalENM enm;
    enm.build_from_ca(ca_ref_);
    ASSERT_TRUE(enm.is_built());
    EXPECT_EQ(enm.n_residues(), 50);
    EXPECT_EQ(enm.n_bonds(), 49);
    EXPECT_FALSE(enm.modes().empty());
}

TEST_F(TENCoMEntropyDiffTest, BuildFromCaTooFew) {
    tencm::TorsionalENM enm;
    std::vector<std::array<float,3>> tiny = {{0,0,0}, {1,0,0}};
    enm.build_from_ca(tiny);
    EXPECT_FALSE(enm.is_built());
}

TEST_F(TENCoMEntropyDiffTest, SelfDifferentialIsZero) {
    tencm::TorsionalENM enm1, enm2;
    enm1.build_from_ca(ca_ref_);
    enm2.build_from_ca(ca_ref_);

    ASSERT_TRUE(enm1.is_built());
    ASSERT_TRUE(enm2.is_built());

    // Same structure → eigenvalues should be identical
    const auto& m1 = enm1.modes();
    const auto& m2 = enm2.modes();
    ASSERT_EQ(m1.size(), m2.size());

    for (std::size_t i = 0; i < m1.size(); ++i) {
        EXPECT_NEAR(m1[i].eigenvalue, m2[i].eigenvalue, 1e-10)
            << "Mode " << i << " eigenvalue mismatch";
    }

    // ΔS_vib should be zero
    double T = 300.0;
    auto s1 = tencom_svib(m1, T);
    auto s2 = tencom_svib(m2, T);
    EXPECT_NEAR(s1.S_vib_kcal_mol_K - s2.S_vib_kcal_mol_K, 0.0, 1e-12);
}

TEST_F(TENCoMEntropyDiffTest, PerturbationProducesNonZeroDelta) {
    tencm::TorsionalENM ref_enm, pert_enm;
    ref_enm.build_from_ca(ca_ref_);
    pert_enm.build_from_ca(ca_pert_);

    ASSERT_TRUE(ref_enm.is_built());
    ASSERT_TRUE(pert_enm.is_built());

    const auto& mr = ref_enm.modes();
    const auto& mp = pert_enm.modes();

    // At least some eigenvalues should differ
    double sum_delta = 0.0;
    int n = std::min(mr.size(), mp.size());
    for (int i = 6; i < n; ++i)
        sum_delta += std::abs(mp[i].eigenvalue - mr[i].eigenvalue);

    EXPECT_GT(sum_delta, 0.0) << "Perturbation should change eigenvalues";

    // ΔS_vib should be non-zero
    double T = 300.0;
    auto sr = tencom_svib(mr, T);
    auto sp = tencom_svib(mp, T);
    double delta_s = sp.S_vib_kcal_mol_K - sr.S_vib_kcal_mol_K;
    EXPECT_NE(delta_s, 0.0);
}

TEST_F(TENCoMEntropyDiffTest, EigenvectorOverlapIdentity) {
    tencm::TorsionalENM enm1, enm2;
    enm1.build_from_ca(ca_ref_);
    enm2.build_from_ca(ca_ref_);

    const auto& m1 = enm1.modes();
    const auto& m2 = enm2.modes();

    // Same structure → eigenvector overlap should be 1.0 for all modes
    for (std::size_t i = 6; i < m1.size() && i < 26; ++i) {
        const auto& v1 = m1[i].eigenvector;
        const auto& v2 = m2[i].eigenvector;
        ASSERT_EQ(v1.size(), v2.size());

        double dot = 0, n1 = 0, n2 = 0;
        for (std::size_t j = 0; j < v1.size(); ++j) {
            dot += v1[j] * v2[j];
            n1  += v1[j] * v1[j];
            n2  += v2[j] * v2[j];
        }
        double overlap = std::abs(dot / (std::sqrt(n1) * std::sqrt(n2)));
        EXPECT_NEAR(overlap, 1.0, 1e-10)
            << "Mode " << i << " self-overlap should be 1.0";
    }
}

TEST_F(TENCoMEntropyDiffTest, BfactorDifferentialSign) {
    tencm::TorsionalENM ref_enm, pert_enm;
    ref_enm.build_from_ca(ca_ref_);
    pert_enm.build_from_ca(ca_pert_);

    auto bf_ref  = ref_enm.bfactors(300.0f);
    auto bf_pert = pert_enm.bfactors(300.0f);

    ASSERT_EQ(bf_ref.size(), bf_pert.size());

    // At least some B-factors should differ
    bool any_diff = false;
    for (std::size_t i = 0; i < bf_ref.size(); ++i) {
        if (std::abs(bf_pert[i] - bf_ref[i]) > 1e-6) {
            any_diff = true;
            break;
        }
    }
    EXPECT_TRUE(any_diff) << "Perturbed structure should have different B-factors";
}

TEST_F(TENCoMEntropyDiffTest, CaPositionsAccessor) {
    tencm::TorsionalENM enm;
    enm.build_from_ca(ca_ref_);
    const auto& pos = enm.ca_positions();
    ASSERT_EQ(pos.size(), ca_ref_.size());
    for (std::size_t i = 0; i < pos.size(); ++i) {
        EXPECT_FLOAT_EQ(pos[i][0], ca_ref_[i][0]);
        EXPECT_FLOAT_EQ(pos[i][1], ca_ref_[i][1]);
        EXPECT_FLOAT_EQ(pos[i][2], ca_ref_[i][2]);
    }
}

TEST_F(TENCoMEntropyDiffTest, DifferentSizeStructures) {
    // Two structures with different residue counts should both build OK
    auto ca_small = make_helix(30);
    auto ca_large = make_helix(60);

    tencm::TorsionalENM enm_s, enm_l;
    enm_s.build_from_ca(ca_small);
    enm_l.build_from_ca(ca_large);

    ASSERT_TRUE(enm_s.is_built());
    ASSERT_TRUE(enm_l.is_built());
    EXPECT_EQ(enm_s.n_residues(), 30);
    EXPECT_EQ(enm_l.n_residues(), 60);

    // Mode counts differ — comparison uses min
    int n = std::min(enm_s.modes().size(), enm_l.modes().size());
    EXPECT_GT(n, 6) << "Should have enough modes for comparison";
}

TEST_F(TENCoMEntropyDiffTest, VibrationalEntropyPositive) {
    tencm::TorsionalENM enm;
    enm.build_from_ca(ca_ref_);

    auto vs = tencom_svib(enm.modes(), 300.0);
    // Vibrational entropy should be positive at finite temperature
    EXPECT_GT(vs.S_vib_kcal_mol_K, 0.0);
    EXPECT_GT(vs.n_modes, 0);
}

TEST_F(TENCoMEntropyDiffTest, DefaultEntropyMetadataIsModelScale) {
    tencm::TorsionalENM enm;
    enm.build_from_ca(ca_ref_);

    auto vs = tencom_svib(enm.modes(), 300.0);
    EXPECT_FALSE(vs.calibrated);
    EXPECT_FALSE(vs.absolute_claim_allowed());
    EXPECT_EQ(vs.calibration_label, "model-scale");
    EXPECT_STREQ(vs.calibration_status(), "model_scale_heuristic");
}

TEST_F(TENCoMEntropyDiffTest, CalibratedScalePropagatesInEntropyDiffHelper) {
    tencm::TorsionalENM enm;
    enm.build_from_ca(ca_ref_);
    auto calibration = encom::FrequencyCalibration::calibrated_scale(
        1.0e12, "synthetic-test", "unit-test provenance");

    auto base = tencom_svib(enm.modes(), 300.0);
    auto scaled = tencom_svib(enm.modes(), 300.0, calibration);

    EXPECT_TRUE(scaled.calibrated);
    EXPECT_TRUE(scaled.absolute_claim_allowed());
    EXPECT_EQ(scaled.calibration_label, "synthetic-test");
    EXPECT_NEAR(scaled.omega_eff,
                base.omega_eff * 1.0e12,
                std::abs(base.omega_eff * 1.0e12) * 1e-12);
}
