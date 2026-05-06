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

struct Replica {
    int    id;
    double temperature;
    double beta;              // 1/(kT)
    double current_energy;
};

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

    // Boltzmann weight vector (same order as insertion)
    std::vector<double> boltzmann_weights() const;

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

private:
    double T_;
    double beta_;
    std::vector<State> ensemble_;

    // Numerically stable log(Σ exp(x_i))
    static double log_sum_exp(std::span<const double> x);
};

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
