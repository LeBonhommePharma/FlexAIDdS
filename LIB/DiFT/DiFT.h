// DiFT.h — Discrete Fourier Transform torsional parametrization for FlexAID∆S
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
//
// ─────────────────────────────────────────────────────────────────────────────
// DiFT: automated, unbiased parametrization of dihedral (torsional) potentials.
//
// Reference:
//   Flores-Trujillo, Rodríguez-Segura, Amador-Bedolla & Domínguez (2026).
//   "Fast Fourier Transform Enables Automated Parametrization of Complex
//    Dihedral Potentials in All-Atom and Coarse-Grained Force Fields."
//   J. Chem. Inf. Model.  DOI: 10.1021/acs.jcim.6c00123
//
// ─────────────────────────────────────────────────────────────────────────────
// WHY THIS LIVES IN AN ENTROPY-DRIVEN DOCKING ENGINE
//
// A torsional potential V(φ) and its Fourier spectrum {Aₙ, ωₙ} are conjugate
// representations of the SAME object. FlexAID∆S scores binding by free energy,
// ΔG = ΔH − TΔS — yet ligand rotatable-bond entropy has classically been a
// crude per-bond count penalty. DiFT fixes that by first principles:
//
//   ONE Fast Fourier Transform → TWO payoffs
//     (1) ENERGY  — a truncated cosine series V_tors(φ) = Σ Aₙ cos(nφ − ωₙ)
//                   feeds the GA fitness as a real analytical potential.
//     (2) ENTROPY — the SAME spectrum yields, in closed form, the 1-D torsional
//                   partition function z = ⟨exp(−βV)⟩ and hence the per-bond
//                   torsional entropy S_tors = (⟨V⟩ − F)/T. This is a rigorous
//                   statistical-mechanical ΔS contribution, not a heuristic.
//
// SHANNON-COLLAPSE TRUNCATION (the unbiased stopping criterion)
//   The paper increments the term count F by brute force and tests R². We
//   replace that with a parameter-free criterion: the spectral Shannon entropy
//   of the normalized power spectrum  pₙ = Aₙ² / Σ Aₘ².  The effective number
//   of modes is the spectral participation ratio  N_eff = exp(H_spec).  A
//   profile dominated by one cosine has H_spec → 0 (N_eff → 1); a flat noisy
//   profile has H_spec → ln N (all modes). The spectrum itself tells you how
//   many frequencies it needs — no user-defined amplitude threshold. This is
//   Shannon's energy collapse applied to the dihedral spectrum.
// ─────────────────────────────────────────────────────────────────────────────
#pragma once

#include <vector>
#include <span>
#include <cstddef>

namespace dift {

// ─── physical constants ──────────────────────────────────────────────────────
inline constexpr double kB_kcal = 0.001987206;   // kcal mol⁻¹ K⁻¹  (matches statmech)

// ─── data structures ─────────────────────────────────────────────────────────

// A single cosine term of the torsional potential.
//   contribution(φ) = amplitude · cos(multiplicity·φ − phase)
// The GAFF/CHARMM "1 + cos" constant is an additive offset folded into
// TorsionalPotential::mean; here we keep the pure oscillatory part.
struct FourierTerm {
    int    multiplicity = 0;   // n  — integer frequency
    double amplitude    = 0.0; // Aₙ — force constant (kcal/mol), ≥ 0 by convention
    double phase        = 0.0; // ωₙ — phase shift (radians), in (−π, π]
    double power        = 0.0; // Aₙ² — spectral power (cached for truncation)
};

// A DiFT-parametrized torsional potential: the analytical representation
// V(φ) = mean + Σ Aₙ cos(nφ − ωₙ) over the retained (truncated) terms.
struct TorsionalPotential {
    std::vector<FourierTerm> terms;            // retained terms, n ≥ 1, sorted by n
    double mean             = 0.0;             // A₀ — DC offset (kcal/mol)
    double v_min            = 0.0;             // min of V(φ) over a fine grid
    int    n_samples        = 0;               // M — input grid resolution
    double r_squared        = 0.0;             // fit quality vs the target profile
    double spectral_entropy = 0.0;             // H_spec (nats) of the FULL spectrum
    double effective_modes  = 0.0;             // N_eff = exp(H_spec)
    int    refinement_iters = 0;               // QM–MM refinement iterations used

    // Absolute potential at an arbitrary angle φ (radians).
    double evaluate(double phi) const noexcept;
    // Potential relative to its global minimum — the form used for scoring,
    // since only relative torsional energies are physically meaningful.
    double relative(double phi) const noexcept { return evaluate(phi) - v_min; }
    // Sample the analytical potential on a uniform grid of n points over [0,2π).
    std::vector<double> sample(int n) const;
    // Total number of retained cosine terms.
    std::size_t n_terms() const noexcept { return terms.size(); }
};

// Per-bond torsional thermodynamics, derived from the SAME Fourier spectrum.
// Entropy is reported as an EXCESS quantity relative to a free rotor (V ≡ 0):
//   free rotor      → S_tors = 0     (no confinement, no penalty)
//   confined bond   → S_tors < 0     (entropy loss; −TS > 0 = ΔG penalty)
struct TorsionalThermo {
    double temperature        = 300.0; // K
    double partition_function = 1.0;   // z = ⟨exp(−βV)⟩  (1 for a free rotor)
    double free_energy        = 0.0;   // F = −kT ln z            (kcal/mol)
    double mean_energy        = 0.0;   // ⟨V⟩                     (kcal/mol)
    double entropy            = 0.0;   // S_tors = (⟨V⟩ − F)/T    (kcal mol⁻¹ K⁻¹)
    double minus_TS           = 0.0;   // −T·S_tors               (kcal/mol)
};

// ─── core engine ─────────────────────────────────────────────────────────────

class DiFTEngine {
public:
    explicit DiFTEngine(double temperature_K = 300.0);

    // Forward DiFT: M equally-spaced samples of a torsional profile over [0,2π)
    // → the full Fourier spectrum (terms n = 1 … ⌊M/2⌋). The DC mean is
    // returned via the out-parameter. Uses a radix-2 FFT for power-of-two M
    // and an exact direct DFT otherwise.
    std::vector<FourierTerm> transform(std::span<const double> profile,
                                       double& mean_out) const;

    // Parametrize a target torsional profile: forward transform followed by
    // Shannon-collapse spectral truncation. If max_multiplicity > 0, terms with
    // n above it are discarded first (anti-overfit guard, paper uses n ≤ 6).
    TorsionalPotential parametrize(std::span<const double> profile,
                                   int max_multiplicity = 0) const;

    // Iterative QM–MM refinement (paper eq. 18): V_{i+1} = V_i + λ·D_i, with
    // D_i = V_QM − V_MM,i. Both profiles must share the SAME M-point grid.
    // Each step is FFT band-limited; stops when R² ≥ r2_target.
    TorsionalPotential refine(std::span<const double> qm,
                              std::span<const double> mm_initial,
                              double lambda           = 0.5,
                              double r2_target        = 0.98,
                              int    max_iter         = 50,
                              int    max_multiplicity = 6) const;

    // Boltzmann-invert a coarse-grained dihedral-angle histogram into an energy
    // profile (paper eq. 20):  E(φ) = −kT ln p(φ), shifted so min = 0.
    // Empty bins are assigned the maximum finite energy (capped well wall).
    std::vector<double> boltzmann_invert(std::span<const double> histogram) const;

    // Per-bond torsional thermodynamics from a parametrized potential.
    // Integrates z = ⟨exp(−βV)⟩ on a fine grid using the analytical V(φ).
    TorsionalThermo thermodynamics(const TorsionalPotential& pot) const;

    // Circular (directional) mean of a set of phase offsets — required when
    // averaging dihedral phase shifts across distinct torsional pathways
    // (paper eq. 14). Returns atan2(Σ sin, Σ cos) in (−π, π].
    static double circular_mean(std::span<const double> angles) noexcept;

    // Coefficient of determination R² between an observed profile and a model
    // profile of equal length.
    static double r_squared(std::span<const double> observed,
                            std::span<const double> model) noexcept;

    double temperature() const noexcept { return T_; }
    void   set_temperature(double T_K) noexcept;

private:
    double T_;
    double beta_;   // 1 / (kB·T)
};

// Spectral Shannon entropy H_spec (nats) of a power spectrum: the normalized
// powers pₙ = Aₙ²/ΣAₘ² (n ≥ 1) form a distribution; H = −Σ pₙ ln pₙ.
// exp(H) is the effective number of contributing frequencies (N_eff).
double spectral_entropy(std::span<const FourierTerm> spectrum) noexcept;

} // namespace dift
