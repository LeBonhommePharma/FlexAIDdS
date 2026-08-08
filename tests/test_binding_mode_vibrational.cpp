// tests/test_binding_mode_vibrational.cpp
// Unit tests for BindingMode vibrational correction (Phase 3)
// Validates ENCoM-based -T*S_vib integration into BindingMode free energy
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include "../LIB/BindingMode.h"
#include "../LIB/statmech.h"
#include "../LIB/encom.h"
#include "../LIB/gaboom.h"
#include <cmath>
#include <vector>

// ===========================================================================
// Test-accessible subclass to reach protected methods
// ===========================================================================
class TestableBindingMode : public BindingMode {
public:
    using BindingMode::BindingMode;
    using BindingMode::compute_vibrational_correction;
};

// ===========================================================================
// TEST FIXTURE
// ===========================================================================

class BindingModeVibrationalTest : public ::testing::Test {
protected:
    FA_Global* mock_fa;
    GB_Global* mock_gb;
    VC_Global* mock_vc;
    chromosome* mock_chroms;
    genlim* mock_gene_lim;
    atom* mock_atoms;
    resid* mock_residue;
    gridpoint* mock_cleftgrid;
    BindingPopulation* test_population;

    static constexpr int NUM_ATOMS = 20;
    static constexpr int NUM_RESIDUES = 2;
    static constexpr double TEST_TEMPERATURE = 300.0;
    static constexpr double EPSILON = 1e-6;

    void SetUp() override {
        mock_fa = new FA_Global();
        mock_fa->temperature = static_cast<uint>(TEST_TEMPERATURE);
        mock_fa->normal_modes = 0;  // No vibrational correction by default
        // Vibrational correction lives on the physical ranking/ledger path.
        // Classic soft-β global-Z ranking is covered in test_classic_entropy_ranking.
        mock_fa->force_cf_rank_emission = true;

        mock_gb = new GB_Global();
        mock_gb->num_genes = 6;

        mock_vc = new VC_Global();

        mock_chroms = new chromosome[5];
        mock_gene_lim = new genlim[mock_gb->num_genes];
        mock_atoms = new atom[NUM_ATOMS];
        mock_residue = new resid[NUM_RESIDUES];
        mock_cleftgrid = new gridpoint[100];

        // Initialize atoms with null eigen pointers
        for (int i = 0; i < NUM_ATOMS; ++i) {
            mock_atoms[i].eigen = nullptr;
        }

        test_population = new BindingPopulation(
            mock_fa, mock_gb, mock_vc,
            mock_chroms, mock_gene_lim,
            mock_atoms, mock_residue, mock_cleftgrid,
            5
        );
    }

    void TearDown() override {
        // Clean up eigenvalue arrays if allocated
        if (mock_atoms[0].eigen) {
            for (int m = 0; m < mock_fa->normal_modes; ++m) {
                delete[] mock_atoms[0].eigen[m];
            }
            delete[] mock_atoms[0].eigen;
            mock_atoms[0].eigen = nullptr;
        }

        delete test_population;
        delete[] mock_cleftgrid;
        delete[] mock_residue;
        delete[] mock_atoms;
        delete[] mock_gene_lim;
        delete[] mock_chroms;
        delete mock_vc;
        delete mock_gb;
        delete mock_fa;
    }

    Pose create_mock_pose(double cf_value, int index) {
        std::vector<float> empty_vec;
        Pose p(&mock_chroms[index], index, 0, 0.0, TEST_TEMPERATURE, empty_vec);
        p.CF = cf_value;
        return p;
    }

    // Reproduce the INVALID atom-0 layout these tests used to depend on.
    //
    // This layout does not exist in production and never did:
    //   * `atom::eigen` holds normal-mode eigenVECTORS, not eigenvalues.
    //     assign_eigen.cpp allocates eigen[m] as THREE floats (the x/y/z
    //     displacement components), whereas this helper allocates one.
    //   * assign_eigen.cpp only populates real atoms — its residue loop starts
    //     at index 1 — and read_pdb.cpp:38 explicitly sets atoms[0].eigen to
    //     NULL, so the sentinel atom never carries mode data at all.
    //
    // The helper is retained precisely so the tests below can prove that
    // BindingMode REFUSES to read this fabricated channel. Treating
    // eigen[m][0] as an eigenvalue would manufacture a vibrational entropy the
    // elastic-network model never produced.
    void setup_invalid_atom0_eigen_layout(int n_modes, double base_value = 1.0) {
        mock_fa->normal_modes = n_modes;
        mock_atoms[0].eigen = new float*[n_modes];
        for (int m = 0; m < n_modes; ++m) {
            mock_atoms[0].eigen[m] = new float[1];
            mock_atoms[0].eigen[m][0] = static_cast<float>(base_value * (m + 1));
        }
    }
};

// ===========================================================================
// VIBRATIONAL CORRECTION TESTS
// ===========================================================================

TEST_F(BindingModeVibrationalTest, NoModesZeroCorrection) {
    // When normal_modes == 0, vibrational correction should be 0
    TestableBindingMode mode(test_population);
    Pose p = create_mock_pose(-10.0, 0);
    mode.add_Pose(p);

    double energy_no_vib = mode.compute_energy();
    double correction = mode.compute_vibrational_correction();

    EXPECT_NEAR(correction, 0.0, EPSILON);

    // Free energy should equal StatMech-only free energy
    auto thermo = mode.get_thermodynamics();
    // With zero correction, compute_energy == statmech free energy
    EXPECT_NEAR(energy_no_vib, thermo.free_energy, EPSILON);
}

TEST_F(BindingModeVibrationalTest, InvalidAtomZeroEigenLayoutIsIgnored) {
    // The fabricated atom-0 "eigenvalue" channel must NOT be consumed.
    setup_invalid_atom0_eigen_layout(5, 0.5);

    TestableBindingMode mode(test_population);
    Pose p = create_mock_pose(-10.0, 0);
    mode.add_Pose(p);

    double correction = mode.compute_vibrational_correction();

    // Fail closed: no valid eigenvalue channel exists, so there is no
    // vibrational correction to report.
    EXPECT_TRUE(std::isfinite(correction));
    EXPECT_NEAR(correction, 0.0, EPSILON);
}

TEST_F(BindingModeVibrationalTest, ThermodynamicsLedgerStaysInternallyCoherent) {
    setup_invalid_atom0_eigen_layout(5, 0.5);

    TestableBindingMode mode(test_population);
    for (int i = 0; i < 3; ++i) {
        Pose p = create_mock_pose(-10.0 - i * 2.0, i);
        mode.add_Pose(p);
    }

    double correction = mode.compute_vibrational_correction();
    EXPECT_NEAR(correction, 0.0, EPSILON);

    // get_thermodynamics() must return the configurational ensemble unmodified,
    // so the struct satisfies its own identity F = H - T*S. Shifting F alone
    // (the previous behaviour) silently violated this.
    auto thermo = mode.get_thermodynamics();
    ASSERT_GT(thermo.temperature, 0.0);
    EXPECT_NEAR(thermo.free_energy,
                thermo.mean_energy - thermo.temperature * thermo.entropy,
                1e-9);

    // A CF/contact-function ensemble is never promotable to a physical claim.
    EXPECT_TRUE(thermo.is_proxy_only());
}

TEST_F(BindingModeVibrationalTest, ModeCountCannotResurrectTheDisabledCorrection) {
    // Previously "more modes -> more negative correction". With the invalid
    // channel disabled, the declared mode count is irrelevant: both are zero.
    setup_invalid_atom0_eigen_layout(3, 0.5);

    TestableBindingMode mode_few(test_population);
    Pose p1 = create_mock_pose(-10.0, 0);
    mode_few.add_Pose(p1);
    double correction_few = mode_few.compute_vibrational_correction();

    for (int m = 0; m < 3; ++m) delete[] mock_atoms[0].eigen[m];
    delete[] mock_atoms[0].eigen;
    mock_atoms[0].eigen = nullptr;

    setup_invalid_atom0_eigen_layout(10, 0.5);

    TestableBindingMode mode_many(test_population);
    Pose p2 = create_mock_pose(-10.0, 0);
    mode_many.add_Pose(p2);
    double correction_many = mode_many.compute_vibrational_correction();

    EXPECT_NEAR(correction_few, 0.0, EPSILON);
    EXPECT_NEAR(correction_many, 0.0, EPSILON);
    EXPECT_NEAR(correction_many, correction_few, EPSILON);
}

TEST_F(BindingModeVibrationalTest, NullEigenReturnsZero) {
    // normal_modes > 0 but eigen pointer is null → should return 0 safely
    mock_fa->normal_modes = 5;
    // mock_atoms[0].eigen is already nullptr from SetUp

    TestableBindingMode mode(test_population);
    Pose p = create_mock_pose(-10.0, 0);
    mode.add_Pose(p);

    double correction = mode.compute_vibrational_correction();
    EXPECT_NEAR(correction, 0.0, EPSILON);
}

TEST_F(BindingModeVibrationalTest, TemperatureCannotResurrectTheDisabledCorrection) {
    setup_invalid_atom0_eigen_layout(5, 0.5);

    // Test at 300K
    TestableBindingMode mode_300(test_population);
    Pose p1 = create_mock_pose(-10.0, 0);
    mode_300.add_Pose(p1);
    double correction_300 = mode_300.compute_vibrational_correction();

    // Change temperature to 600K
    mock_fa->temperature = 600;
    test_population->Temperature = 600;

    TestableBindingMode mode_600(test_population);
    Pose p2 = create_mock_pose(-10.0, 0);
    mode_600.add_Pose(p2);
    double correction_600 = mode_600.compute_vibrational_correction();

    // -T*S_vib scaled with T only because S_vib was fabricated. With the
    // channel disabled the correction is identically zero at every T.
    EXPECT_NEAR(correction_300, 0.0, EPSILON);
    EXPECT_NEAR(correction_600, 0.0, EPSILON);

    // Restore
    mock_fa->temperature = static_cast<uint>(TEST_TEMPERATURE);
    test_population->Temperature = static_cast<unsigned int>(TEST_TEMPERATURE);
}

// ===========================================================================
// EDGE CASES — ZERO POSES
// ===========================================================================

TEST_F(BindingModeVibrationalTest, ZeroPosesVibrationalCorrectionIsZero) {
    // A BindingMode with no poses: correction must be 0 (don't access uninitialised memory)
    TestableBindingMode mode(test_population);
    double correction = mode.compute_vibrational_correction();
    EXPECT_NEAR(correction, 0.0, EPSILON);
}

TEST_F(BindingModeVibrationalTest, ZeroPosesFreeEnergyThrows) {
    // get_thermodynamics on empty mode should throw (statmech requires at least 1 sample)
    TestableBindingMode mode(test_population);
    EXPECT_THROW(mode.get_thermodynamics(), std::runtime_error);
}

// ===========================================================================
// EDGE CASES — EXTREME EIGENVALUES
// ===========================================================================

TEST_F(BindingModeVibrationalTest, ExtremeStiffValuesStillProduceNoCorrection) {
    // A huge fabricated value must not be laundered into a large "stabilising"
    // free-energy term. Magnitude is irrelevant when the channel is invalid.
    setup_invalid_atom0_eigen_layout(5, 1e6);

    TestableBindingMode mode(test_population);
    Pose p = create_mock_pose(-10.0, 0);
    mode.add_Pose(p);

    double correction = mode.compute_vibrational_correction();
    EXPECT_TRUE(std::isfinite(correction));
    EXPECT_NEAR(correction, 0.0, EPSILON);
}

TEST_F(BindingModeVibrationalTest, ExtremeFloppyValuesStillProduceNoCorrection) {
    setup_invalid_atom0_eigen_layout(5, 1e-4);

    TestableBindingMode mode(test_population);
    Pose p = create_mock_pose(-10.0, 0);
    mode.add_Pose(p);

    double correction = mode.compute_vibrational_correction();
    EXPECT_TRUE(std::isfinite(correction));
    EXPECT_NEAR(correction, 0.0, EPSILON);
}

TEST_F(BindingModeVibrationalTest, ZeroValuedEntryHandledSafely) {
    // A zero entry must not produce log(0) NaN. It cannot: the channel is
    // never read at all.
    mock_fa->normal_modes = 3;
    mock_atoms[0].eigen = new float*[3];
    mock_atoms[0].eigen[0] = new float[1]; mock_atoms[0].eigen[0][0] = 0.0f; // zero!
    mock_atoms[0].eigen[1] = new float[1]; mock_atoms[0].eigen[1][0] = 1.0f;
    mock_atoms[0].eigen[2] = new float[1]; mock_atoms[0].eigen[2][0] = 2.0f;

    TestableBindingMode mode(test_population);
    Pose p = create_mock_pose(-10.0, 0);
    mode.add_Pose(p);

    double correction = mode.compute_vibrational_correction();
    EXPECT_TRUE(std::isfinite(correction));
    EXPECT_NEAR(correction, 0.0, EPSILON);
}

// ===========================================================================
// EDGE CASES — SINGLE POSE
// ===========================================================================

TEST_F(BindingModeVibrationalTest, SinglePoseGetsNoFabricatedCorrection) {
    setup_invalid_atom0_eigen_layout(5, 0.5);

    TestableBindingMode mode(test_population);
    Pose p = create_mock_pose(-12.0, 0);
    mode.add_Pose(p);

    double correction = mode.compute_vibrational_correction();
    EXPECT_TRUE(std::isfinite(correction));
    EXPECT_NEAR(correction, 0.0, EPSILON);
}

TEST_F(BindingModeVibrationalTest, SinglePoseFreeEnergyEqualsPoseEnergy) {
    setup_invalid_atom0_eigen_layout(5, 0.5);

    TestableBindingMode mode(test_population);
    Pose p = create_mock_pose(-12.0, 0);
    mode.add_Pose(p);

    double correction = mode.compute_vibrational_correction();
    double total      = mode.compute_energy();
    // For a single state statmech F = E, and the disabled correction adds
    // nothing, so the ranking energy is exactly the pose energy.
    EXPECT_NEAR(correction, 0.0, EPSILON);
    EXPECT_NEAR(total, -12.0, EPSILON);
}

// ===========================================================================
// EDGE CASES — VIBRATIONAL CORRECTION PRESERVES RELATIVE MODE RANKING
// ===========================================================================

TEST_F(BindingModeVibrationalTest, RankingPreservedWhenBothModesHaveSameModes) {
    // If both modes see the same eigenvalue set, the mode with lower CF
    // should still have lower total free energy after correction.
    setup_invalid_atom0_eigen_layout(5, 0.5);

    TestableBindingMode stable(test_population), weak(test_population);
    Pose ps = create_mock_pose(-15.0, 0);
    Pose pw = create_mock_pose(-8.0, 0);
    stable.add_Pose(ps);
    weak.add_Pose(pw);

    double F_stable = stable.compute_energy();
    double F_weak   = weak.compute_energy();
    EXPECT_LT(F_stable, F_weak);
}

TEST_F(BindingModeVibrationalTest, DeltaGRelativeToAnotherModeConsistent) {
    setup_invalid_atom0_eigen_layout(4, 0.5);

    TestableBindingMode mode_a(test_population), mode_b(test_population);
    for (double e : {-14.0, -13.0, -12.0}) {
        Pose pa = create_mock_pose(e, 0);
        mode_a.add_Pose(pa);
    }
    for (double e : {-9.0,  -8.0,  -7.0}) {
        Pose pb = create_mock_pose(e, 1);
        mode_b.add_Pose(pb);
    }

    double F_a = mode_a.compute_energy();
    double F_b = mode_b.compute_energy();

    // ΔG(a relative to b) = F_a - F_b (how much higher energy a is vs b)
    double dG = mode_a.delta_G_relative_to(mode_b);
    EXPECT_NEAR(dG, F_a - F_b, EPSILON);
}

// ===========================================================================
// EDGE CASES — MANY POSES
// ===========================================================================

TEST_F(BindingModeVibrationalTest, ManyPosesVibrationalCorrectionFinite) {
    setup_invalid_atom0_eigen_layout(8, 0.3);

    TestableBindingMode mode(test_population);
    for (int i = 0; i < 8; ++i) {
        Pose p = create_mock_pose(-10.0 - i * 1.5, i % 5);  // mock_chroms has 5 slots
        mode.add_Pose(p);
    }

    double correction = mode.compute_vibrational_correction();
    EXPECT_TRUE(std::isfinite(correction));
    EXPECT_NEAR(correction, 0.0, EPSILON);
}

// ===========================================================================
// MAIN
// ===========================================================================

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
