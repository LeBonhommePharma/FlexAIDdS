// encom.h — ENCoM (Elastic Network Contact Model) Integration for FlexAID∆S
//
// Computes vibrational entropy contributions from normal mode analysis:
//   – Parse ENCoM eigenvector/eigenvalue files
//   – Calculate quasi-harmonic vibrational entropy S_vib
//   – Combine with configurational entropy S_conf for total entropy
//
// ENCoM/tENCoM eigenvalues are stiffness-like model quantities. Absolute
// vibrational entropy/free-energy claims require either a mass/inertia model or
// an empirical eigenvalue-to-angular-frequency calibration. Without that,
// S_vib is a model-scale heuristic suitable for relative comparisons only.
//
// Reference:
//   Frappier et al. (2015). *Proteins* 83(11):2073-82.
//   DOI: 10.1002/prot.24922
//
// Mathematical framework:
//   S_vib = (3N - 6) × k_B × [1 + ln(kBT/ħω_eff)]
//   ω_eff = scale × geometric_mean(sqrt(non-zero eigenvalues))

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
    double  frequency;           // sqrt(λ_i) in model scale unless calibrated
    std::vector<double> eigenvector; // Displacement vector (3N components)
};

struct FrequencyCalibration {
    double eigenvalue_to_omega = 1.0;  // rad s⁻¹ per sqrt(model eigenvalue)
    bool calibrated = false;
    std::string label = "model-scale";
    std::string provenance =
        "No mass/inertia or empirical eigenvalue-to-frequency calibration supplied";

    static FrequencyCalibration model_scale();

    static FrequencyCalibration calibrated_scale(
        double eigenvalue_to_omega,
        const std::string& label,
        const std::string& provenance = "");

    bool valid() const noexcept {
        return std::isfinite(eigenvalue_to_omega) && eigenvalue_to_omega > 0.0;
    }

    const char* status() const noexcept {
        return calibrated ? "calibrated" : "model_scale_heuristic";
    }
};

struct VibrationalEntropy {
    double S_vib_kcal_mol_K = 0.0; // Vibrational entropy (kcal mol⁻¹ K⁻¹)
    double S_vib_J_mol_K = 0.0;    // Vibrational entropy (J mol⁻¹ K⁻¹)
    double omega_eff = 0.0;        // Effective frequency (rad/s)
    int    n_modes = 0;            // Number of non-zero modes (3N - 6)
    double temperature = 300.0;    // K
    double eigenvalue_to_omega = 1.0; // calibration scale used for omega_eff
    bool   calibrated = false;     // true only when scale has physical/empirical provenance
    std::string calibration_label = "model-scale";
    std::string calibration_provenance =
        "No mass/inertia or empirical eigenvalue-to-frequency calibration supplied";

    bool absolute_claim_allowed() const noexcept { return calibrated; }
    const char* calibration_status() const noexcept {
        return calibrated ? "calibrated" : "model_scale_heuristic";
    }
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
    /// Default uses model-scale sqrt(eigenvalue) and is not an absolute entropy
    /// claim. Use the overload with FrequencyCalibration for calibrated output.
    static VibrationalEntropy compute_vibrational_entropy(
        const std::vector<NormalMode>& modes,
        double temperature_K = 300.0,
        double eigenvalue_cutoff = 1e-6  // Skip modes below this threshold
    );

    /// Compute quasi-harmonic vibrational entropy with an explicit
    /// eigenvalue-to-angular-frequency calibration.
    static VibrationalEntropy compute_vibrational_entropy(
        const std::vector<NormalMode>& modes,
        double temperature_K,
        const FrequencyCalibration& calibration,
        double eigenvalue_cutoff = 1e-6
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
