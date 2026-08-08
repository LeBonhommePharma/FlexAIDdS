// ShannonThermoStack.h — Shannon Entropy + Torsional ENCoM thermodynamic stack
//
// Combines:
//   – Shannon configurational entropy over GA ensemble (binned into 256 mega-clusters)
//   – Torsional ENCoM vibrational entropy from NormalMode fluctuations (protein + nucleotide backbones)
//   – Hardware-accelerated histogram computation (Metal on Apple Silicon, OpenMP/Eigen on other platforms)
//
// Reuses StatMechEngine (statmech.h) and TorsionalENM (tencm.h) without modification.
// BindingPopulation is untouched; SugarPuckerGene flexibility is separate.
#pragma once

#include "../statmech.h"
#include "../tENCoM/tencm.h"
#include <vector>
#include <string>
#include <cmath>
#include <mutex>
#include <atomic>

namespace shannon_thermo {

// ─── constants ───────────────────────────────────────────────────────────────
inline constexpr int   SHANNON_BINS      = 256;    // mega-cluster discretisation
inline constexpr double kB_kcal          = 0.001987206; // kcal mol⁻¹ K⁻¹  (= R/4184)
inline constexpr double kB_SI            = 1.380649e-23; // J K⁻¹
inline constexpr double hbar_SI          = 1.054571817e-34; // J·s
inline constexpr double TEMPERATURE_K    = 298.15;
inline constexpr int   DEFAULT_HIST_BINS = 20;
inline constexpr int   GPU_DISPATCH_THRESHOLD = 500000; // only use GPU for N > 500K

// ─── Shannon Energy Collapse thresholds ──────────────────────────────────────
// All internal entropy APIs return nats (natural log). Convert to bits at
// reporting/convergence boundaries only: H_bits = H_nats / ln(2).
//
// These named constants are the single source of truth for H(X) < threshold
// comparisons. Never compare a nats value against the raw bits constant or
// vice versa — always use the matching _nats or _bits form.
//
//   Soft collapse:  H < 2.0 bits  -> effective support < 4 clusters
//   Hard collapse:  H < 1.0 bit   -> one cluster has >50% probability
//
// Derivation: H_nats = H_bits × ln(2);  ln(2) = 0.693147...
inline constexpr double kHSC_soft_bits = 2.0;
inline constexpr double kHSC_hard_bits = 1.0;
inline constexpr double kHSC_soft_nats = kHSC_soft_bits * 0.6931471805599453; // 2 × ln(2)
inline constexpr double kHSC_hard_nats = kHSC_hard_bits * 0.6931471805599453; // 1 × ln(2)

// ⚠ SUPPORT-SIZE CAVEAT — the absolute thresholds above were derived against
// the SHANNON_BINS = 256 "mega-cluster" discretisation described at the top of
// this header (ceiling H_max = ln 256 = 8 bits), but every live caller computes
// H with DEFAULT_HIST_BINS = 20 (ceiling H_max = ln 20 = 4.32 bits). As
// fractions of the reachable maximum the intended 25% / 12.5% lines therefore
// land at 46% / 23%, i.e. both gates are considerably more eager than their
// derivation implies.
//
// The constants are left at their shipped values because they gate GA
// termination and changing them moves docking results — that belongs in a
// benchmarked A/B, not a drive-by edit. New code that wants a gate which means
// the same thing at any bin count should use collapse_threshold_nats() below.
inline constexpr double kHSC_soft_frac_of_max = 0.25;  // "support < 1/4 of bins"
inline constexpr double kHSC_hard_frac_of_max = 0.125; // "one bin dominates"

// Collapse threshold scaled to the support actually used by the estimator:
// returns frac_of_max × ln(num_bins) nats. Bin-count independent by
// construction, unlike the absolute kHSC_*_nats constants above.
inline double collapse_threshold_nats(int num_bins,
                                      double frac_of_max = kHSC_soft_frac_of_max) noexcept {
    if (num_bins <= 1) return 0.0;
    return frac_of_max * std::log(static_cast<double>(num_bins));
}

// ─── result struct ───────────────────────────────────────────────────────────
struct FullThermoResult {
    double deltaG;              // base ΔG plus calibrated entropy terms (kcal/mol)
    double shannonEntropy;      // dimensionless nats (conformational, natural log)
    double torsionalVibEntropy; // kcal/mol·K; heuristic unless calibrated elsewhere
    double entropyContribution; // applied -T*S term (kcal/mol), excludes heuristics
    std::string report;
};

// ─── 256×256 precomputed Shannon energy lookup ───────────────────────────────
// E[i][j] = -kT * p_i * ln(p_j)
// Generated at startup with seed 42 + Gaussian perturbation of uniform priors.
class ShannonEnergyMatrix {
public:
    static ShannonEnergyMatrix& instance();

    // Initialise matrix from uniform priors (seed = 42)
    void initialise();

    // O(1) pairwise entropy contribution lookup
    double lookup(int bin_i, int bin_j) const noexcept {
        return matrix_[bin_i * SHANNON_BINS + bin_j];
    }

    // Load trained 256×256 matrix from binary blob (SHNN header).
    // Returns false on I/O error or magic mismatch.
    bool initialise_from_file(const std::string& path);

    // Load matrix directly from float data (256*256 float32 values).
    void initialise_from_data(const float* data, int count);

    bool is_initialised() const noexcept { return initialised_; }

    // Raw access to underlying data (for pybind11 zero-copy views)
    const double* data() const noexcept { return matrix_.data(); }
    int size() const noexcept { return static_cast<int>(matrix_.size()); }

private:
    ShannonEnergyMatrix() = default;
    std::vector<double> matrix_; // SHANNON_BINS × SHANNON_BINS
    std::atomic<bool> initialised_{false};
    mutable std::once_flag init_once_;   // serialises initialise() across threads
    mutable std::mutex     mtx_;         // guards initialise_from_*() races
};

// ─── Shannon entropy computation ──────────────────────────────────────────────
// Bins a vector of continuous values into numBins and computes Shannon entropy H.
// Uses OpenMP parallelism when available; Metal GPU on Apple Silicon.
double compute_shannon_entropy(const std::vector<double>& values,
                               int num_bins = DEFAULT_HIST_BINS);

// Same for integer state labels (discrete)
double compute_shannon_entropy_discrete(const std::vector<int>& states);

// ─── torsional vibrational entropy from ENCoM modes ─────────────────────────
// Sums classical harmonic oscillator entropy using model-scale torsional
// eigenvalues:
//   S_mode = kB * [1 + ln(kBT/(hbar*omega))]
//   omega = sqrt(lambda) unless an external calibration is introduced.
// Uncalibrated values are useful as relative flexibility heuristics only.
double compute_torsional_vibrational_entropy(
    const std::vector<tencm::NormalMode>& modes,
    double temperature_K = TEMPERATURE_K);

// ─── full stack entry point ───────────────────────────────────────────────────
// Runs the complete ShannonThermoStack on a GA population ensemble.
//
// Parameters:
//   stat_engine    – populated StatMechEngine from the GA run
//   tencm_model    – built TorsionalENM (may be default-constructed if backbone
//                    flexibility is disabled)
//   base_deltaG    – enthalpy-dominated ΔG from scoring function (kcal/mol)
//   temperature_K  – simulation temperature
FullThermoResult run_shannon_thermo_stack(
    const statmech::StatMechEngine& stat_engine,
    const tencm::TorsionalENM&      tencm_model,
    double                          base_deltaG,
    double                          temperature_K = TEMPERATURE_K);

// ─── entropy plateau detection ──────────────────────────────────────────────
// Returns true if the last `window` entries in `history` all have relative
// change < `rel_threshold` from the first entry in the window.
// Used for GA early-termination when Shannon entropy stabilises.
bool detect_entropy_plateau(const std::vector<double>& history,
                            int window, double rel_threshold);

// ─── entropy event classification (unified framework) ──────────────────────
// Mirrors Shannon's EntropyEvent enum and NATURaL's EntropyEvent enum.
// Collapse  = entropy drops (binding lock-in, ordering)
// Expansion = entropy rises (solvation release, disordering)
// Oscillation = rapid alternation (unstable binding site, external perturbation)
enum class EntropyEventType { None = 0, Collapse = 1, Expansion = 2, Oscillation = 3 };

struct EntropyEventResult {
    EntropyEventType event    = EntropyEventType::None;
    double            delta   = 0.0;
    double            z_score = 0.0;
    double            entropy = 0.0;
};

// ─── sliding-window entropy event detector ──────────────────────────────────
// Tracks entropy over a window of recent measurements and classifies events.
// Used for GA diversity monitoring and binding trajectory analysis.
class EntropyEventDetector {
public:
    explicit EntropyEventDetector(
        int window_size = 8,
        double collapse_threshold = -3.2,
        double expansion_threshold = +3.2,
        int oscillation_window = 5);

    EntropyEventResult push(double entropy);
    void reset();

    const std::vector<double>& history() const noexcept { return history_; }

private:
    int window_size_;
    double collapse_threshold_;
    double expansion_threshold_;
    int oscillation_window_;
    std::vector<double> history_;
    std::vector<EntropyEventType> event_history_;

    bool detect_oscillation() const;
};

} // namespace shannon_thermo
