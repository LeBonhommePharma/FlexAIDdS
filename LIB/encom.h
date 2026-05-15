// encom.h — ENCoM (Elastic Network Contact Model) Integration for FlexAID∆S
//
// Computes vibrational entropy contributions from normal mode analysis:
//   – Parse ENCoM eigenvector/eigenvalue files
//   – Calculate quasi-harmonic vibrational entropy S_vib
//   – Combine with configurational entropy S_conf for total entropy
//
// Calibration status:
//   ENCoM eigenvalues are model-scale quantities, not SI frequencies by
//   default. Absolute S_vib and -T*S_vib values are therefore heuristic unless
//   the eigenvalue-to-frequency scale is calibrated for the benchmark system.
//   Differential comparisons under one fixed protocol are the supported use.
//
// Reference:
//   Frappier et al. (2015). *Proteins* 83(11):2073-82.
//   DOI: 10.1002/prot.24922
//
// Mathematical framework:
//   S_vib = (3N - 6) × k_B × [1 + ln(kT/ħω_eff)]
//   ω_eff = calibrated geometric mean frequency from non-zero eigenvalues

#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <numeric>

namespace encom {

// ──────────────────────────────────────────────────────────────────────────────
// Physical constants
// ──────────────────────────────────────────────────────────────────────────────

inline constexpr double kB_kcal    = 0.001987206;      // kcal mol⁻¹ K⁻¹
inline constexpr double kB_SI      = 1.380649e-23;     // J K⁻¹
inline constexpr double hbar_SI    = 1.054571817e-34;  // J·s
inline constexpr double NA         = 6.02214076e23;    // mol⁻¹
inline constexpr double amu_to_kg  = 1.66053906660e-27;// kg

// ──────────────────────────────────────────────────────────────────────────────
// Data structures
// ──────────────────────────────────────────────────────────────────────────────

struct NormalMode {
    int     index;               // Mode number (1-based)
    double  eigenvalue;          // λ_i (arbitrary units from ENCoM)
    double  frequency;           // ω_i = sqrt(λ_i) before any calibration scale
    std::vector<double> eigenvector; // Displacement vector (3N components)
};

struct VibrationalEntropy {
    double S_vib_kcal_mol_K;     // Heuristic vibrational entropy unless calibrated
    double S_vib_J_mol_K;        // Same value in J mol⁻¹ K⁻¹
    double omega_eff;            // Effective model frequency
    int    n_modes;              // Number of non-zero modes (3N - 6)
    double temperature;          // K
};

// ──────────────────────────────────────────────────────────────────────────────
// ENCoM mode reader and vibrational entropy calculator
// ──────────────────────────────────────────────────────────────────────────────

class ENCoMEngine {
public:
    /// Load eigenvalues and eigenvectors from ENCoM output files
    /// Format: plain text, one eigenvalue per line, eigenvectors in separate file
    static std::vector<NormalMode> load_modes(
        const std::string& eigenvalue_file,
        const std::string& eigenvector_file
    );
    
    /// Compute quasi-harmonic vibrational entropy from normal modes.
    /// Absolute magnitudes require a calibrated eigenvalue-to-frequency scale.
    /// Without calibration, use only as a relative heuristic.
    static VibrationalEntropy compute_vibrational_entropy(
        const std::vector<NormalMode>& modes,
        double temperature_K = 300.0,
        double eigenvalue_cutoff = 1e-6  // Skip modes below this threshold
    );
    
    /// Combine configurational entropy (from StatMechEngine) with vibrational
    static double total_entropy(
        double S_conf_kcal_mol_K,
        double S_vib_kcal_mol_K
    ) noexcept {
        return S_conf_kcal_mol_K + S_vib_kcal_mol_K;
    }
    
    /// Compute free energy including vibrational correction:
    /// F_total = F_elec + F_vib = (H_elec - T·S_conf) + (-T·S_vib)
    static double free_energy_with_vibrations(
        double F_electronic,        // from BindingMode::compute_energy()
        double S_vib_kcal_mol_K,
        double temperature_K
    ) noexcept {
        return F_electronic - temperature_K * S_vib_kcal_mol_K;
    }

private:
    /// Helper: geometric mean of eigenvalues
    static double geometric_mean(const std::vector<double>& values);
};

}  // namespace encom
