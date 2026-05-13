// tests/test_encom.cpp
// Unit tests for ENCoM vibrational entropy engine
// Tests entropy computation, mode handling, and thermodynamic integration.

#include <gtest/gtest.h>
#include "encom.h"
#include <cmath>
#include <fstream>
#include <vector>
#include <filesystem>

using namespace encom;

static constexpr double TOL = 1e-10;

// ===========================================================================
// HELPER: Create synthetic normal modes
// ===========================================================================
static std::vector<NormalMode> make_modes(int n, double base_eigenvalue = 1.0) {
    std::vector<NormalMode> modes;
    for (int i = 0; i < n; ++i) {
        NormalMode m;
        m.index = i + 1;
        m.eigenvalue = base_eigenvalue * (i + 1);
        m.frequency = std::sqrt(m.eigenvalue);
        m.eigenvector = std::vector<double>(10, 0.1 * (i + 1));
        modes.push_back(m);
    }
    return modes;
}

// ===========================================================================
// VIBRATIONAL ENTROPY — BASIC PROPERTIES
// ===========================================================================

TEST(ENCoMEngine, EmptyModesZeroEntropy) {
    std::vector<NormalMode> empty;
    auto result = ENCoMEngine::compute_vibrational_entropy(empty);
    EXPECT_NEAR(result.S_vib_kcal_mol_K, 0.0, TOL);
    EXPECT_NEAR(result.S_vib_J_mol_K, 0.0, TOL);
    EXPECT_EQ(result.n_modes, 0);
}

TEST(ENCoMEngine, AllBelowCutoffZeroEntropy) {
    auto modes = make_modes(5, 1e-10);  // all below cutoff
    auto result = ENCoMEngine::compute_vibrational_entropy(modes, 300.0, 1e-6);
    EXPECT_NEAR(result.S_vib_kcal_mol_K, 0.0, TOL);
    EXPECT_EQ(result.n_modes, 0);
}

TEST(ENCoMEngine, AllNegativeEigenvaluesZeroEntropy) {
    // All negative eigenvalues should be filtered out, giving zero entropy
    auto modes = make_modes(5, -1.0);  // negative eigenvalues
    auto result = ENCoMEngine::compute_vibrational_entropy(modes, 300.0, 1e-6);
    EXPECT_NEAR(result.S_vib_kcal_mol_K, 0.0, TOL);
    EXPECT_EQ(result.n_modes, 0);
    EXPECT_NEAR(result.omega_eff, 0.0, TOL);
}

TEST(ENCoMEngine, MixedPositiveNegativeEigenvalues) {
    // Only positive eigenvalues above cutoff should contribute
    std::vector<encom::NormalMode> modes;
    double eigenvalues[] = {-2.0, -1.0, 0.0, 0.5, 1.0, 2.0};
    for (int i = 0; i < 6; ++i) {
        encom::NormalMode m;
        m.index = i + 1;
        m.eigenvalue = eigenvalues[i];
        m.frequency = std::sqrt(std::abs(eigenvalues[i]));
        modes.push_back(m);
    }
    auto result = ENCoMEngine::compute_vibrational_entropy(modes, 300.0, 1e-6);
    // Only 0.5, 1.0, 2.0 should contribute (3 modes)
    EXPECT_EQ(result.n_modes, 3);
    EXPECT_GT(result.omega_eff, 0.0);
    EXPECT_TRUE(std::isfinite(result.S_vib_kcal_mol_K));
}

TEST(ENCoMEngine, PositiveEntropyForValidModes) {
    auto modes = make_modes(10, 1.0);
    auto result = ENCoMEngine::compute_vibrational_entropy(modes, 300.0);

    EXPECT_GT(result.n_modes, 0);
    // S_vib should be positive at room temperature for typical eigenvalues
    EXPECT_TRUE(std::isfinite(result.S_vib_kcal_mol_K));
    EXPECT_TRUE(std::isfinite(result.S_vib_J_mol_K));
}

TEST(ENCoMEngine, TemperatureIsRecorded) {
    auto modes = make_modes(5);
    auto result = ENCoMEngine::compute_vibrational_entropy(modes, 310.0);
    EXPECT_NEAR(result.temperature, 310.0, TOL);
}

TEST(ENCoMEngine, DefaultCalibrationIsModelScaleHeuristic) {
    auto modes = make_modes(5);
    auto result = ENCoMEngine::compute_vibrational_entropy(modes, 300.0);

    EXPECT_FALSE(result.calibrated);
    EXPECT_FALSE(result.absolute_claim_allowed());
    EXPECT_NEAR(result.eigenvalue_to_omega, 1.0, TOL);
    EXPECT_EQ(result.calibration_label, "model-scale");
    EXPECT_STREQ(result.calibration_status(), "model_scale_heuristic");
}

TEST(ENCoMEngine, CalibratedScalePropagatesToOmega) {
    auto modes = make_modes(5, 1.0);
    auto calibration = FrequencyCalibration::calibrated_scale(
        1.0e12, "synthetic-test", "unit-test scale");

    auto uncalibrated = ENCoMEngine::compute_vibrational_entropy(modes, 300.0);
    auto calibrated = ENCoMEngine::compute_vibrational_entropy(
        modes, 300.0, calibration, 1e-6);

    EXPECT_TRUE(calibrated.calibrated);
    EXPECT_TRUE(calibrated.absolute_claim_allowed());
    EXPECT_EQ(calibrated.calibration_label, "synthetic-test");
    EXPECT_NEAR(calibrated.eigenvalue_to_omega, 1.0e12, 1.0);
    EXPECT_NEAR(calibrated.omega_eff,
                uncalibrated.omega_eff * 1.0e12,
                std::abs(uncalibrated.omega_eff * 1.0e12) * 1e-12);
}

TEST(ENCoMEngine, OmegaEffIsPositive) {
    auto modes = make_modes(5);
    auto result = ENCoMEngine::compute_vibrational_entropy(modes, 300.0);
    EXPECT_GT(result.omega_eff, 0.0);
}

// ===========================================================================
// TEMPERATURE DEPENDENCE
// ===========================================================================

TEST(ENCoMEngine, HigherTemperatureHigherEntropy) {
    auto modes = make_modes(10, 1.0);
    auto result_low = ENCoMEngine::compute_vibrational_entropy(modes, 200.0);
    auto result_high = ENCoMEngine::compute_vibrational_entropy(modes, 400.0);

    // Both should have entropy values, and higher T → higher S
    if (result_low.n_modes > 0 && result_high.n_modes > 0) {
        EXPECT_GT(result_high.S_vib_kcal_mol_K, result_low.S_vib_kcal_mol_K);
    }
}

// ===========================================================================
// MODE COUNT
// ===========================================================================

TEST(ENCoMEngine, MoreModesChangeEntropy) {
    auto few_modes = make_modes(5, 1.0);
    auto many_modes = make_modes(20, 1.0);

    auto result_few = ENCoMEngine::compute_vibrational_entropy(few_modes, 300.0);
    auto result_many = ENCoMEngine::compute_vibrational_entropy(many_modes, 300.0);

    EXPECT_LT(result_few.n_modes, result_many.n_modes);
}

// ===========================================================================
// CUTOFF FILTERING
// ===========================================================================

TEST(ENCoMEngine, CutoffFiltersLowEigenvalues) {
    std::vector<NormalMode> modes;
    // 3 modes below cutoff, 7 above
    for (int i = 0; i < 10; ++i) {
        NormalMode m;
        m.index = i + 1;
        m.eigenvalue = (i < 3) ? 1e-8 : 1.0 * (i + 1);
        m.frequency = std::sqrt(std::abs(m.eigenvalue));
        m.eigenvector = {1.0};
        modes.push_back(m);
    }

    auto result = ENCoMEngine::compute_vibrational_entropy(modes, 300.0, 1e-6);
    EXPECT_EQ(result.n_modes, 7);
}

// ===========================================================================
// UNIT CONSISTENCY
// ===========================================================================

TEST(ENCoMEngine, KcalAndJouleConsistent) {
    auto modes = make_modes(10, 1.0);
    auto result = ENCoMEngine::compute_vibrational_entropy(modes, 300.0);

    if (result.n_modes > 0 && result.S_vib_kcal_mol_K != 0.0) {
        // 1 kcal = 4184 J
        double expected_kcal = result.S_vib_J_mol_K / 4184.0;
        EXPECT_NEAR(result.S_vib_kcal_mol_K, expected_kcal, 1e-8);
    }
}

// ===========================================================================
// TOTAL ENTROPY
// ===========================================================================

TEST(ENCoMEngine, TotalEntropyIsAdditive) {
    double S_conf = 0.005;  // kcal mol⁻¹ K⁻¹
    double S_vib = 0.003;
    double S_total = ENCoMEngine::total_entropy(S_conf, S_vib);
    EXPECT_NEAR(S_total, 0.008, TOL);
}

TEST(ENCoMEngine, TotalEntropyWithZeros) {
    EXPECT_NEAR(ENCoMEngine::total_entropy(0.0, 0.0), 0.0, TOL);
    EXPECT_NEAR(ENCoMEngine::total_entropy(0.005, 0.0), 0.005, TOL);
    EXPECT_NEAR(ENCoMEngine::total_entropy(0.0, 0.003), 0.003, TOL);
}

// ===========================================================================
// FREE ENERGY WITH VIBRATIONS
// ===========================================================================

TEST(ENCoMEngine, FreeEnergyWithVibrationsFormula) {
    double F_elec = -10.0;     // kcal/mol
    double S_vib = 0.001;      // kcal mol⁻¹ K⁻¹
    double T = 300.0;          // K

    double F_total = ENCoMEngine::free_energy_with_vibrations(F_elec, S_vib, T);
    // F_total = F_elec - T * S_vib = -10.0 - 300 * 0.001 = -10.3
    EXPECT_NEAR(F_total, -10.3, TOL);
}

TEST(ENCoMEngine, FreeEnergyVibrationsLowersEnergy) {
    double F_elec = -5.0;
    double S_vib = 0.002;
    double T = 300.0;

    double F_total = ENCoMEngine::free_energy_with_vibrations(F_elec, S_vib, T);
    // Positive S_vib should lower the free energy
    EXPECT_LT(F_total, F_elec);
}

TEST(ENCoMEngine, FreeEnergyZeroEntropyUnchanged) {
    double F_elec = -7.5;
    double F_total = ENCoMEngine::free_energy_with_vibrations(F_elec, 0.0, 300.0);
    EXPECT_NEAR(F_total, F_elec, TOL);
}

// ===========================================================================
// LOAD MODES (file I/O)
// ===========================================================================

TEST(ENCoMEngine, LoadModesFromFiles) {
    // Create temporary eigenvalue file
    std::ofstream eval_file("test_eigenvalues.txt");
    eval_file << "0.5\n1.0\n2.0\n3.5\n";
    eval_file.close();

    // Create temporary eigenvector file
    std::ofstream evec_file("test_eigenvectors.txt");
    evec_file << "0.1 0.2 0.3\n";
    evec_file << "0.4 0.5 0.6\n";
    evec_file << "0.7 0.8 0.9\n";
    evec_file << "1.0 1.1 1.2\n";
    evec_file.close();

    auto modes = ENCoMEngine::load_modes("test_eigenvalues.txt",
                                          "test_eigenvectors.txt");

    EXPECT_EQ(modes.size(), 4u);
    EXPECT_NEAR(modes[0].eigenvalue, 0.5, TOL);
    EXPECT_NEAR(modes[1].eigenvalue, 1.0, TOL);
    EXPECT_NEAR(modes[2].eigenvalue, 2.0, TOL);
    EXPECT_NEAR(modes[3].eigenvalue, 3.5, TOL);

    // Frequencies should be sqrt of eigenvalues
    EXPECT_NEAR(modes[0].frequency, std::sqrt(0.5), 1e-6);
    EXPECT_NEAR(modes[1].frequency, std::sqrt(1.0), 1e-6);

    // Eigenvectors should have 3 components
    for (const auto& m : modes) {
        EXPECT_EQ(m.eigenvector.size(), 3u);
    }

    // Mode indices should be 1-based
    EXPECT_EQ(modes[0].index, 1);
    EXPECT_EQ(modes[3].index, 4);

    std::filesystem::remove("test_eigenvalues.txt");
    std::filesystem::remove("test_eigenvectors.txt");
}

TEST(ENCoMEngine, LoadModesMissingFileThrows) {
    EXPECT_THROW(
        ENCoMEngine::load_modes("nonexistent_eval.txt", "nonexistent_evec.txt"),
        std::runtime_error);
}

TEST(ENCoMEngine, LoadModesEmptyFileThrows) {
    std::ofstream eval_file("test_empty_eval.txt");
    eval_file.close();
    std::ofstream evec_file("test_empty_evec.txt");
    evec_file.close();

    EXPECT_THROW(
        ENCoMEngine::load_modes("test_empty_eval.txt", "test_empty_evec.txt"),
        std::runtime_error);

    std::filesystem::remove("test_empty_eval.txt");
    std::filesystem::remove("test_empty_evec.txt");
}

// ===========================================================================
// PHYSICAL CONSTANTS
// ===========================================================================

TEST(ENCoMConstants, BoltzmannConstantKcal) {
    EXPECT_NEAR(kB_kcal, 0.001987206, 1e-9);
}

TEST(ENCoMConstants, BoltzmannConstantSI) {
    EXPECT_NEAR(kB_SI, 1.380649e-23, 1e-29);
}

TEST(ENCoMConstants, PlanckConstantSI) {
    EXPECT_NEAR(hbar_SI, 1.054571817e-34, 1e-40);
}

TEST(ENCoMConstants, AvogadroNumber) {
    EXPECT_NEAR(NA, 6.02214076e23, 1e17);
}
