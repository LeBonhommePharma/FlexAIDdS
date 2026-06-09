// statmech.h — Statistical Mechanics Engine for FlexAIDdS
//
// Full thermodynamic analysis of the GA conformational ensemble:
//   – Partition function Z(T) with log-sum-exp numerical stability
//   – Helmholtz free energy  F = −kT ln Z
//   – Average energy ⟨E⟩, variance ⟨E²⟩−⟨E⟩², heat capacity C_v
//   – Conformational entropy  S = (⟨E⟩ − F)/T
//   – Boltzmann-weighted probability of each sampled state
//   – Parallel tempering (replica exchange) swap acceptance
//   – Boltzmann-reweighted PMF (free energy profiles along an arbitrary coordinate)
//   – Thermodynamic integration (TI) via trapezoidal rule
//   – Fast Boltzmann lookup table for inner-loop evaluation
#pragma once

#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <span>
#include <random>
#include <stdexcept>

namespace statmech {

// ─── physical constants ──────────────────────────────────────────────────────
inline constexpr double kB_kcal = 0.001987206;   // kcal mol⁻¹ K⁻¹
inline constexpr double kB_SI   = 1.380649e-23;  // J K⁻¹

// ─── data structures ─────────────────────────────────────────────────────────

struct State {
    double energy;     // CF value (kcal/mol; negative = favourable)
    double count;      // degeneracy / sampling multiplicity (double so Boltzmann weights pass without truncation)
};

struct Thermodynamics {
    double temperature;       // K
    double log_Z;             // ln Z  (for numerical stability we store the log)
    double free_energy;       // F = −kT ln Z  (kcal/mol)
    double mean_energy;       // ⟨E⟩
    double mean_energy_sq;    // ⟨E²⟩
    double heat_capacity;     // C_v = (⟨E²⟩ − ⟨E⟩²) / (kT²)
    double entropy;           // S = (⟨E⟩ − F) / T  (kcal mol⁻¹ K⁻¹)
    double std_energy;        // σ_E = sqrt(C_v kT²)
};

// ─── THERMODYNAMIC LEDGER (Task 1 — auditable breakdown) ─────────────────────
// Single source of truth for all thermodynamic quantities exposed by the engine.
// All fields carry explicit units in their names per architectural principles.
// This struct aggregates the canonical ensemble result (G_config etc.) plus
// optional additive corrections (vibrational, NATURaL, other) WITHOUT changing
// any legacy ranking or public API behaviour.
//
// G_total = G_config + G_vib + G_natural + G_other  (always)
// Legacy Thermodynamics (free_energy, mean_energy, entropy, ...) remain the
// source of truth for the configurational part; this ledger is derived from it.
//
// DO NOT use for ranking, pose selection, or optimization in early phases.
// All new fields are additive and optional. has_* flags indicate presence.

// EnergyComponents must be defined before ThermodynamicBreakdown (which contains it)
enum class ComponentStatus {
    Available,
    IncludedInOther,
    NotComputed,
    Experimental
};

struct EnergyComponents {
    double total = 0.0;
    double cf = 0.0;
    double receptor_strain = 0.0;
    double ligand_internal = 0.0;
    double hbond = 0.0;
    double gist = 0.0;
    double metal = 0.0;
    double water = 0.0;
    double other = 0.0;
    bool complete = false;

    ComponentStatus cf_status = ComponentStatus::Available;
    ComponentStatus receptor_strain_status = ComponentStatus::NotComputed;
    ComponentStatus ligand_internal_status = ComponentStatus::NotComputed;
    ComponentStatus hbond_status = ComponentStatus::NotComputed;
    ComponentStatus gist_status = ComponentStatus::NotComputed;
    ComponentStatus metal_status = ComponentStatus::NotComputed;
    ComponentStatus water_status = ComponentStatus::NotComputed;
    ComponentStatus other_status = ComponentStatus::Available;

    bool has_meaningful_components() const {
        return cf_status == ComponentStatus::Available ||
               receptor_strain_status == ComponentStatus::Available;
    }
};

struct ComponentAverages {
    double mean_CF_kcal_mol = 0.0;
    double mean_receptor_strain_kcal_mol = 0.0;
    double mean_ligand_internal_kcal_mol = 0.0;
    double mean_hbond_kcal_mol = 0.0;
    double mean_gist_kcal_mol = 0.0;
    double mean_metal_kcal_mol = 0.0;
    double mean_water_kcal_mol = 0.0;
    double mean_other_kcal_mol = 0.0;
    double component_sum_kcal_mol = 0.0;
    bool component_completeness_flag = false;
    ComponentStatus component_status = ComponentStatus::NotComputed;
};

struct ThermodynamicBreakdown {
    double temperature_K = 300.0;

    // Configurational ensemble (from StatMechEngine / GA poses)
    double logZ_config = 0.0;                 // ln Z (dimensionless)
    double G_config_kcal_mol = 0.0;           // F_config = -kB T logZ
    double H_eff_kcal_mol = 0.0;              // ⟨E⟩ Boltzmann-weighted mean
    double S_config_kcal_mol_K = 0.0;         // (H_eff - G_config) / T
    double minus_T_S_config_kcal_mol = 0.0;   // G_config - H_eff
    double Cv_kcal_mol_K = 0.0;               // variance(E) / (kB T²)
    double sigma_E_kcal_mol = 0.0;            // sqrt(variance(E))

    // Additive corrections (populated by callers: BindingMode, tENCoM, NATURaL, ...)
    double G_vib_kcal_mol = 0.0;              // ENCoM / tENCoM vibrational free energy correction
    double G_natural_kcal_mol = 0.0;          // NATURaL co-translational / receptor strain correction
    double G_other_kcal_mol = 0.0;            // Future: explicit GIST, custom terms, etc.
    double G_total_kcal_mol = 0.0;            // G_config + G_vib + G_natural + G_other

    // Presence flags (true only when the corresponding correction was intentionally supplied)
    bool has_vib = false;
    bool has_natural = false;
    bool has_other = false;

    // ═══ Task 4: Diagnostic metrics (never for ranking) ═══
    double entropy_fraction() const;
    double enthalpy_fraction() const;
    double compensation_score() const;

    // ═══ COMPONENT-WISE BOLTZMANN AVERAGES (Task 3) ═══
    // These are ensemble averages: <X> = Σ_i p_i * X_i using the same Boltzmann weights
    // as the rest of the ledger. They are populated when component data is available.
    //
    // IMPORTANT: H_eff is the weighted total energy. component_sum may differ from H_eff
    // when not all energy terms are tracked in EnergyComponents (common case).
    // The completeness flag tells consumers whether they can treat component_sum ≈ H_eff.

    EnergyComponents component_means;   // all fields are <X> = Σ p_i X_i
    double component_sum_kcal_mol = 0.0; // sum of the mean components (for diagnostics)
    bool   components_complete = false;  // true only if every significant term was tracked
    ComponentAverages components;         // legacy BindingMode component diagnostic surface
    bool   has_components = false;

    // ─── Enthalpy-Entropy Index (Williams et al. 2017, Drug Discov. Today) ────
    // I_EE = (ΔH + T·ΔS) / ΔG
    //
    // Sign convention (binding context):
    //   ΔH    = H_eff_kcal_mol           (negative = exothermic, favourable)
    //   T·ΔS  = -minus_T_S_config_kcal_mol  (positive = entropic driving force)
    //   ΔG    = G_total_kcal_mol         (negative = spontaneous binding)
    //
    // Interpretation (ΔG < 0 for all real binders):
    //   I_EE > 1  → entropy-assisted   (TΔS adds to enthalpy to drive ΔG)
    //   I_EE = 1  → pure enthalpy      (TΔS ≈ 0)
    //   I_EE < 1  → entropy-opposed    (enthalpy must overcome unfavourable ΔS)
    //   I_EE < 0  → enthalpy-opposed   (entropy-driven binding, rare)
    //
    // has_I_EE = false when |ΔG| < 1e-6 kcal/mol (numerically undefined).
    // DIAGNOSTIC ONLY — never use for ranking.
    double I_EE     = std::numeric_limits<double>::quiet_NaN();
    bool   has_I_EE = false;
};

struct Replica {
    int    id;
    double temperature;
    double beta;              // 1/(kT)
    double current_energy;
};

// ─── Diagnostic Enthalpy–Entropy Metrics (Task 4) ────────────────────────────
// These functions are **diagnostic only**.
// They must never be used for ranking, pose selection, optimization,
// or any affinity claim.
//
// compensation_score high → strong enthalpy-entropy compensation (G small relative to parts)
// compensation_score low  → one term dominates
//
// All functions are safe for near-zero denominators (return well-defined values).

inline constexpr double kDiagnosticEpsilon = 1e-12;

inline double entropy_fraction(double H_eff_kcal_mol, double minus_T_S_config_kcal_mol) {
    const double denom = std::abs(H_eff_kcal_mol) + std::abs(minus_T_S_config_kcal_mol) + kDiagnosticEpsilon;
    return std::abs(minus_T_S_config_kcal_mol) / denom;
}

inline double enthalpy_fraction(double H_eff_kcal_mol, double minus_T_S_config_kcal_mol) {
    const double denom = std::abs(H_eff_kcal_mol) + std::abs(minus_T_S_config_kcal_mol) + kDiagnosticEpsilon;
    return std::abs(H_eff_kcal_mol) / denom;
}

inline double compensation_score(double G_config_kcal_mol,
                                 double H_eff_kcal_mol,
                                 double minus_T_S_config_kcal_mol) {
    const double denom = std::abs(H_eff_kcal_mol) + std::abs(minus_T_S_config_kcal_mol) + kDiagnosticEpsilon;
    double score = 1.0 - (std::abs(G_config_kcal_mol) / denom);
    if (score < 0.0) score = 0.0;
    if (score > 1.0) score = 1.0;
    return score;
}

// ─── Enthalpy-Entropy Index (Williams et al. 2017) ───────────────────────────
// compute_IEE(delta_H, T_delta_S, delta_G)
//
//   I_EE = (ΔH + T·ΔS) / ΔG
//
// All arguments in kcal/mol. T·ΔS is the full T×entropy product, not ΔS alone,
// so no temperature factor is needed here.
//
// Returns quiet_NaN when |delta_G| < 1e-6 kcal/mol.
//
// DIAGNOSTIC ONLY — never use for ranking or pose selection.
inline double compute_IEE(double delta_H,
                          double T_delta_S,
                          double delta_G) noexcept
{
    if (std::abs(delta_G) < 1e-6)
        return std::numeric_limits<double>::quiet_NaN();
    return (delta_H + T_delta_S) / delta_G;
}

// Convenience overload that reads directly from a ThermodynamicBreakdown ledger.
// Mutates the passed ledger: sets bd.I_EE and bd.has_I_EE in place.
inline void fill_IEE(ThermodynamicBreakdown& bd) noexcept
{
    const double T_dS = -bd.minus_T_S_config_kcal_mol; // T·ΔS (kcal/mol, positive = favourable)
    bd.I_EE     = compute_IEE(bd.H_eff_kcal_mol, T_dS, bd.G_total_kcal_mol);
    bd.has_I_EE = std::isfinite(bd.I_EE);
}

// ─── Joint Receptor–Ligand Ensemble (Task 5 — EXPERIMENTAL) ──────────────────
// Formalizes the joint microstate analysis over receptor conformers (r) and
// ligand poses (i):  Z = Σ_r Σ_i exp[-β E(r,i)]
//
// This is marked EXPERIMENTAL until properly benchmarked.
// If receptor_conformer_id is not available, fallback mode sets
// S_receptor = 0 and mutual_information = 0.

struct JointMicrostate {
    int receptor_conformer_id = -1;   // -1 means unknown / single conformer
    int ligand_pose_id = -1;
    int binding_mode_id = -1;
    EnergyComponents energy;          // decomposed energy for this microstate
    double log_multiplicity = 0.0;    // log(n) for degeneracy
};

struct JointEnsembleResult {
    double temperature_K = 300.0;

    double logZ = 0.0;
    double G_kcal_mol = 0.0;
    double H_kcal_mol = 0.0;
    double S_joint_kcal_mol_K = 0.0;
    double S_receptor_kcal_mol_K = 0.0;
    double S_ligand_kcal_mol_K = 0.0;
    double mutual_information_dimensionless = 0.0;

    std::vector<double> receptor_population;  // p(r)
    std::vector<double> ligand_population;    // p(i)

    bool experimental = true;                 // always true for now
    bool fallback_single_receptor = false;    // true if no receptor conformer info was available
};

// ─── Standard-State Affinity Calibration (Task 6 — EXPERIMENTAL / SAFE ONLY) ─
// This provides utilities to convert between standard-state ΔG° and Kd (in molar)
// while strictly enforcing safety rules.
//
// Key invariants:
// - Never output a real Kd unless calibrated == true.
// - Relative free energies (ΔΔG) are allowed but must be clearly labelled "relative".
// - All functions reject invalid inputs (T <= 0, Kd <= 0).
// - This is **not** true experimental affinity unless a calibration benchmark exists.

struct AffinityCalibration {
    double temperature_K = 300.0;

    // Bound and unbound reference free energies (if available)
    double F_bound_kcal_mol = 0.0;
    double F_unbound_receptor_kcal_mol = 0.0;
    double F_unbound_ligand_kcal_mol = 0.0;

    double standard_state_correction_kcal_mol = 0.0;  // RT ln(c° / 1M) etc.
    double deltaG_standard_kcal_mol = 0.0;            // ΔG° at standard state
    double predicted_Kd_M = 0.0;                      // Only valid if calibrated == true

    bool calibrated = false;   // Must be true before using predicted_Kd_M as real affinity
    bool experimental = true;  // Always true until a proper calibration benchmark suite exists
};

// Safe conversion utilities (Task 6)
double deltaG_standard_to_Kd_M(double deltaG_kcal_mol, double T_K, double c0_M = 1.0);
double Kd_M_to_deltaG_standard(double Kd_M, double T_K, double c0_M = 1.0);

struct WHAMBin {
    double coord_center;
    double count;
    double free_energy;       // kcal/mol
};

struct TIPoint {
    double lambda;            // coupling parameter [0,1]
    double dV_dlambda;        // ⟨∂V/∂λ⟩_λ
};

// ─── main engine ─────────────────────────────────────────────────────────────

class StatMechEngine {
public:
    explicit StatMechEngine(double temperature_K = 300.0);

    // Add a sampled configuration.
    // multiplicity is double so Boltzmann weights (0.0–1.0) can be passed directly
    // without silent int-truncation to zero (which caused log(0) = -inf).
    void add_sample(double energy, double multiplicity = 1.0);

    // Compute full thermodynamics over the current ensemble
    Thermodynamics compute() const;

    // Re-evaluate the existing ensemble at a different temperature WITHOUT
    // re-running the GA. Only β changes; the sample set is unchanged.
    // Used for finite-difference ΔCp computation and temperature scanning.
    Thermodynamics compute_at_temperature(double T_K) const;

    // Compute an auditable thermodynamic ledger without changing legacy fields.
    ThermodynamicBreakdown compute_breakdown(
        double G_vib_kcal_mol = 0.0,
        double G_natural_kcal_mol = 0.0,
        double G_other_kcal_mol = 0.0,
        bool has_vib = false,
        bool has_natural = false,
        bool has_other = false) const;

    // Boltzmann-weight arbitrary energy components over the current ensemble.
    ComponentAverages component_averages(
        std::span<const EnergyComponents> components) const;

    // Boltzmann weight vector (same order as insertion).
    // Uses the PHYSICAL β = 1/(kB·T). This is the thermodynamically correct
    // weighting and the public-API contract exercised by the Python tests —
    // do NOT change its β convention.
    std::vector<double> boltzmann_weights() const;

    // Selection-temperature weights for GA/cluster *selection* (NOT physics).
    // Uses β_sel = 1/T — the SAME convention as FA->beta in the clustering
    // code (read_input.cpp:251) — instead of the kB-folded β = 1/(kB·T).
    //
    // Rationale: CF "energies" are in arbitrary CF units, not kcal/mol, so
    // folding the kcal-specific kB (≈0.001987) into the selection softmax
    // inflates β to ≈1.68 at T=300. A CF spread of ~100 units then yields a
    // Boltzmann ratio of e^{168} → effectively a zero-temperature argmax,
    // collapsing the thermal diversity SMFREE is meant to inject. Using
    // β_sel = 1/T matches the clustering convention (β≈0.0033) so selection
    // and clustering weight poses consistently. compute()/boltzmann_weights()
    // physics are left untouched. See P1 fix.
    std::vector<double> selection_weights() const;

    double beta_selection() const noexcept { return beta_selection_; }
    void   set_selection_beta(double b) noexcept { beta_selection_ = b; }
    void   set_selection_temperature(double T_sel) {
        if (T_sel <= 0.0)
            throw std::invalid_argument("StatMechEngine::set_selection_temperature: T must be > 0");
        beta_selection_ = 1.0 / T_sel;
    }

    // ΔG relative to another engine's ensemble
    double delta_G(const StatMechEngine& reference) const;

    // Parallel tempering: set up replicas at given temperatures
    static std::vector<Replica> init_replicas(
        std::span<const double> temperatures);

    // Attempt Metropolis swap between replicas a and b.
    // Returns true if accepted.
    static bool attempt_swap(Replica& a, Replica& b, std::mt19937& rng);

    // Boltzmann-reweighted free energy profile along a 1D collective
    // coordinate. For each bin b:
    //
    //     F_b = -kT · ln( Σ_{i∈b} exp(-β E_i) / N_b )
    //
    // This is NOT multi-window WHAM (Kumar et al. 1992). It is a single-
    // window post-hoc reweighting of an existing ensemble — useful for
    // building a PMF from a converged GA trajectory along an arbitrary
    // reaction coordinate. Multi-window WHAM requires biased simulations
    // and per-window offsets, neither of which are provided here.
    //
    // Use this when you have a single biased/unbiased ensemble and want
    // a 1D free energy curve. For umbrella-sampling unbiasing, use a
    // dedicated multi-window WHAM implementation.
    static std::vector<WHAMBin> boltzmann_pmf(
        std::span<const double> energies,
        std::span<const double> coordinates,
        double temperature,
        int    n_bins,
        int    max_iter  = 1000,
        double tolerance = 1e-6);

    // Backward-compatible alias for the historical (misleading) name.
    [[deprecated("This is single-window Boltzmann reweighting, not multi-window WHAM. "
                 "Use boltzmann_pmf() — same arguments, accurate name.")]]
    static std::vector<WHAMBin> wham(
        std::span<const double> energies,
        std::span<const double> coordinates,
        double temperature,
        int    n_bins,
        int    max_iter  = 1000,
        double tolerance = 1e-6) {
        return boltzmann_pmf(energies, coordinates, temperature, n_bins, max_iter, tolerance);
    }

    // Thermodynamic integration via trapezoidal rule
    static double thermodynamic_integration(std::span<const TIPoint> points);

    // ── Ensemble merging (for parallel grid-decomposed docking) ────────────
    // Merge another engine's ensemble into this one.
    // Thermodynamically correct: Z_merged = Σ_all exp(-βE_i).
    void merge(const StatMechEngine& other);

    // Merge from raw arrays (for MPI deserialization).
    // multiplicities is double to round-trip with serialize_multiplicities()
    // (which returns vector<double> after the C-1 fix) and to allow fractional
    // weights without silent int-truncation.
    void merge_samples(std::span<const double> energies,
                       std::span<const double> multiplicities);

    // Serialize ensemble for transport (MPI, socket, etc.)
    std::vector<double> serialize_energies() const;
    std::vector<double> serialize_multiplicities() const;

    // Accessors
    double temperature() const noexcept { return T_; }
    double beta()        const noexcept { return beta_; }
    size_t size()        const noexcept { return ensemble_.size(); }
    void   clear()               { ensemble_.clear(); }

    // Read-only access to ensemble (for serialization/inspection)
    const std::vector<State>& ensemble() const noexcept { return ensemble_; }

    // Convenience: Helmholtz free energy from a raw energy vector
    static double helmholtz(std::span<const double> energies, double T);

    // ─── Thermodynamic ledger factory (Task 1) ──────────────────────────────
    // Builds a fully-audited breakdown from an existing engine result.
    // Corrections are additive and optional. When a correction is supplied,
    // the corresponding has_* flag must be set by the caller.
    // This function performs no I/O, no ranking, and has no side effects on
    // the engine or any global state.
    static ThermodynamicBreakdown make_breakdown(
        const StatMechEngine& engine,
        double G_vib_kcal_mol = 0.0,     bool has_vib = false,
        double G_natural_kcal_mol = 0.0, bool has_natural = false,
        double G_other_kcal_mol = 0.0,   bool has_other = false);

    // ─── Component-wise ensemble averages (Task 3) ──────────────────────────
    // Given a vector of Boltzmann weights (from boltzmann_weights()) and a
    // parallel vector of EnergyComponents (one per microstate), returns the
    // properly weighted averages:  <X> = Σ (w_i * X_i) / Σ w_i
    //
    // This is the function that implements Σ_i p_i * CF_i etc.
    // It does NOT modify ranking or total_energy.
    static EnergyComponents compute_weighted_components(
        std::span<const double> weights,
        std::span<const EnergyComponents> components);

    // Convenience overload: compute both the ledger and the component averages
    // in one call when you have the raw data.
    static ThermodynamicBreakdown make_breakdown_with_components(
        const StatMechEngine& engine,
        std::span<const EnergyComponents> components,
        double G_vib_kcal_mol = 0.0,     bool has_vib = false,
        double G_natural_kcal_mol = 0.0, bool has_natural = false,
        double G_other_kcal_mol = 0.0,   bool has_other = false);

    // ─── Joint Receptor–Ligand Ensemble (Task 5 — EXPERIMENTAL) ─────────────
    static JointEnsembleResult compute_joint_ensemble(
        std::span<const JointMicrostate> microstates,
        double temperature_K = 300.0);

private:
    double T_;
    double beta_;            // physical β = 1/(kB·T) — thermodynamics
    double beta_selection_;  // selection β = 1/T — GA/cluster softmax (P1)
    std::vector<State> ensemble_;

    // Numerically stable log(Σ exp(x_i))
    static double log_sum_exp(std::span<const double> x);
};

// ─── ΔCp of binding — finite-difference temperature derivative ───────────────
//
// ΔCp = (∂ΔH/∂T)p ≈ [ΔH(T+dT) − ΔH(T−dT)] / (2·dT)     (central diff)
//
// The same GA ensemble is re-evaluated at T±dT via compute_at_temperature().
// No additional simulation is required. dT = 10 K is safe: truncation error
// is O(dT²); below 1 K floating-point cancellation dominates.
//
// Consistency cross-check: ΔCp ≈ T_ref × (ΔS_hi − ΔS_lo)/(2·dT)
// Both paths should agree within ~5%.
//
// IMPORTANT: distinct from Cv = (⟨E²⟩−⟨E⟩²)/(kBT²).
//   Cv   = fluctuation at one temperature — what StatMechEngine::compute() gives.
//   ΔCp  = finite diff of binding enthalpy across two temperatures — this.

struct DeltaCpResult {
    double T_ref_K = 0.0;               // central temperature (K)
    double dT_K    = 0.0;               // step used (K)
    double delta_H_lo = 0.0;            // ΔH at T_ref − dT  (kcal/mol)
    double delta_H_hi = 0.0;            // ΔH at T_ref + dT  (kcal/mol)
    double delta_S_lo = 0.0;            // ΔS at T_ref − dT  (kcal/mol·K)
    double delta_S_hi = 0.0;            // ΔS at T_ref + dT  (kcal/mol·K)
    double delta_Cp = 0.0;              // kcal/(mol·K) — enthalpy path
    double delta_Cp_from_entropy = 0.0; // kcal/(mol·K) — entropy path (cross-check)
    double consistency_check = 0.0;     // |ΔCp_H − ΔCp_S| / (|mean| + 1e-9)
    bool   consistent = false;          // true when consistency_check < 0.05
};

// Compute ΔCp of binding via central finite difference.
// bound:   engine for the ligand-bound GA ensemble.
// unbound: engine for the apo / reference ensemble.
// Both must carry populations sampled at T_ref_K.
DeltaCpResult compute_delta_Cp(
    const StatMechEngine& bound,
    const StatMechEngine& unbound,
    double T_ref_K = 298.15,
    double dT_K    = 10.0);

// ─── fast Boltzmann lookup table ─────────────────────────────────────────────
//  Pre-tabulates exp(−β E) over [E_min, E_max] for O(1) inner-loop access.

class BoltzmannLUT {
public:
    BoltzmannLUT(double beta, double e_min, double e_max, int n_bins = 10000);
    double operator()(double energy) const noexcept;

private:
    [[maybe_unused]] double beta_;
    double e_min_, inv_bin_width_;
    int    n_bins_;
    std::vector<double> table_;
};

}  // namespace statmech
