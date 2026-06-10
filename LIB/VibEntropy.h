// VibEntropy.h — H(ω) Vibrational-Mode Shannon Entropy (FlexAIDdS Level 3)
//
// A convergence diagnostic for the genetic algorithm: instead of measuring
// Shannon collapse on the *energy* landscape (H(E), the classic SEC metric),
// this module measures it on the *vibrational* landscape — the distribution of
// normal-mode frequencies (ENCoM / tENCoM Hessian eigenvalues) computed across
// the population of cluster representatives in a GA generation.
//
//   • H_pop  (H(ω))   — Shannon entropy of the pooled, log-binned frequency
//                       distribution across ALL cluster reps. Collapse ⇒ the
//                       population has converged onto a single vibrational
//                       fingerprint (conformational collapse in ω-space).
//   • H_rep_mean      — mean per-rep Shannon entropy (mode diversity *within*
//                       an individual binding pose).
//   • D_vib           — inter-rep divergence: mean symmetric-KL between each
//                       rep's frequency distribution and the population mean.
//                       Large D_vib ⇒ reps occupy distinct vibrational basins.
//
// This is LP's Shannon Energy Collapse (SEC) idea applied to the vibrational
// spectrum rather than the binding-energy spectrum. It is a *relative*
// diagnostic over generations; the eigenvalues are ENCoM model-scale stiffness
// quantities (see encom.h) and carry no absolute thermodynamic meaning here.
//
// Design notes:
//   • Eigenvalues span many orders of magnitude (soft global modes ≪ stiff
//     local modes), so binning is done in LOG-frequency space. Uniform bins
//     would crush all the low-frequency structure into a single bin.
//   • All entropies are reported in NATS (natural log), consistent with the
//     internal convention in ShannonThermoStack. Convert to bits at reporting
//     boundaries via H_bits = H_nats / ln(2) if desired.
//   • Pure functions, no global/static mutable state — thread-safe.
//
// Env overrides:
//   FLEXAIDDS_VIB_ENTROPY_BINS — integer bin count (default 32). Values < 2 or
//                                non-numeric are ignored (default retained).

#pragma once

#include <vector>
#include <cstddef>

namespace vibentropy {

// ──────────────────────────────────────────────────────────────────────────────
// Constants
// ──────────────────────────────────────────────────────────────────────────────

inline constexpr int    kDefaultBins = 32;        // log-spaced frequency bins
inline constexpr double kEpsilon     = 1e-10;     // KL / log(0) guard
inline constexpr double kLn2         = 0.6931471805599453; // ln(2), nats→bits

// ──────────────────────────────────────────────────────────────────────────────
// Result struct
// ──────────────────────────────────────────────────────────────────────────────

struct VibEntropyResult {
    double H_pop        = 0.0;  // H(ω): entropy of pooled log-binned distribution (nats)
    double H_rep_mean   = 0.0;  // mean per-rep entropy (nats)
    double D_vib        = 0.0;  // mean symmetric-KL(rep ‖ population) (nats)
    int    n_reps       = 0;    // number of cluster representatives contributing
    int    n_modes_per_rep = 0; // mean modes per contributing rep (rounded)
};

// ──────────────────────────────────────────────────────────────────────────────
// API
// ──────────────────────────────────────────────────────────────────────────────

// Resolve the active bin count, honoring FLEXAIDDS_VIB_ENTROPY_BINS.
int resolve_bin_count();

// Main entry point.
//
// `eigenvalues` is one eigenvalue array per cluster representative:
//   eigenvalues[r][m] = m-th normal-mode eigenvalue (λ ≥ 0) of representative r.
//
// Behavior / edge cases (never throws):
//   • Empty input, or fewer than 1 rep with usable modes  → zeroed result.
//   • Non-positive / non-finite eigenvalues are skipped (log requires λ > 0).
//   • If all usable eigenvalues are equal (degenerate spectrum), every
//     distribution is a single occupied bin → H_pop = H_rep_mean = D_vib = 0,
//     which is the correct "fully collapsed" reading.
//   • With a single rep, D_vib is 0 by definition (rep == population mean).
VibEntropyResult compute_vib_entropy_collapse(
    const std::vector<std::vector<double>>& eigenvalues);

// Same as above but with an explicit bin count (bypasses the env override).
// `n_bins` is clamped to >= 2.
VibEntropyResult compute_vib_entropy_collapse(
    const std::vector<std::vector<double>>& eigenvalues, int n_bins);

} // namespace vibentropy
