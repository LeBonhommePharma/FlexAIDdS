// VibEntropy.cpp — H(ω) Vibrational-Mode Shannon Entropy (FlexAIDdS Level 3)
//
// See VibEntropy.h for the conceptual overview. All routines here are pure and
// operate only on their arguments + one read-only env lookup; no shared mutable
// state, safe to call concurrently across generations / threads.

#include "VibEntropy.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <numeric>

namespace vibentropy {

namespace {

// Map an eigenvalue λ > 0 to a frequency proxy in log space. ENCoM eigenvalues
// are stiffness-like; ω ∝ sqrt(λ), so log(ω) = 0.5·log(λ). Working directly in
// log(λ) is monotonically equivalent for binning purposes and avoids the sqrt.
inline double log_freq(double lambda) {
    return std::log(lambda);
}

inline bool usable(double lambda) {
    return std::isfinite(lambda) && lambda > 0.0;
}

// Shannon entropy (nats) of a probability vector assumed to sum to 1.
double shannon_entropy(const std::vector<double>& p) {
    double H = 0.0;
    for (double pi : p) {
        if (pi > 0.0) H -= pi * std::log(pi);
    }
    // Numerical guard: clamp tiny negatives from round-off.
    return H > 0.0 ? H : 0.0;
}

// Symmetric KL divergence (nats): KL(p‖q) + KL(q‖p), with epsilon smoothing so
// zero-probability bins do not produce log(0) / division-by-zero.
double symmetric_kl(const std::vector<double>& p, const std::vector<double>& q) {
    const std::size_t n = p.size();
    double d = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double pi = p[i] + kEpsilon;
        const double qi = q[i] + kEpsilon;
        d += pi * std::log(pi / qi);  // KL(p‖q)
        d += qi * std::log(qi / pi);  // KL(q‖p)
    }
    return d > 0.0 ? d : 0.0;
}

// Normalize a histogram (counts) into a probability distribution summing to 1.
// If the total mass is zero, returns an all-zero vector (caller treats that rep
// as contributing no entropy).
void normalize(std::vector<double>& h) {
    const double total = std::accumulate(h.begin(), h.end(), 0.0);
    if (total <= 0.0) return;
    const double inv = 1.0 / total;
    for (double& v : h) v *= inv;
}

} // namespace

// ──────────────────────────────────────────────────────────────────────────────

int resolve_bin_count() {
    if (const char* env = std::getenv("FLEXAIDDS_VIB_ENTROPY_BINS")) {
        char* end = nullptr;
        const long v = std::strtol(env, &end, 10);
        // Require full, valid parse and a sane lower bound.
        if (end && *end == '\0' && v >= 2 && v <= 100000) {
            return static_cast<int>(v);
        }
    }
    return kDefaultBins;
}

VibEntropyResult compute_vib_entropy_collapse(
    const std::vector<std::vector<double>>& eigenvalues) {
    return compute_vib_entropy_collapse(eigenvalues, resolve_bin_count());
}

VibEntropyResult compute_vib_entropy_collapse(
    const std::vector<std::vector<double>>& eigenvalues, int n_bins) {

    VibEntropyResult out;
    if (n_bins < 2) n_bins = 2;

    // ── Pass 1: find the usable reps and the global log-frequency range ──────
    double lo =  std::numeric_limits<double>::infinity();
    double hi = -std::numeric_limits<double>::infinity();
    int    n_reps = 0;
    long   total_modes = 0;

    for (const auto& rep : eigenvalues) {
        bool rep_has_mode = false;
        for (double lam : rep) {
            if (!usable(lam)) continue;
            const double lf = log_freq(lam);
            lo = std::min(lo, lf);
            hi = std::max(hi, lf);
            rep_has_mode = true;
            ++total_modes;
        }
        if (rep_has_mode) ++n_reps;
    }

    // Nothing usable → fully collapsed / undefined; return zeroed result.
    if (n_reps == 0 || !std::isfinite(lo) || !std::isfinite(hi)) {
        return out;
    }

    out.n_reps = n_reps;
    out.n_modes_per_rep =
        static_cast<int>(std::llround(static_cast<double>(total_modes) / n_reps));

    // Degenerate spectrum: every usable eigenvalue is identical. All mass lands
    // in one bin for every rep and the population → H = 0, D_vib = 0. That is
    // the correct "vibrationally collapsed" reading, so short-circuit cleanly
    // (also avoids a zero-width bin span).
    const double span = hi - lo;
    if (span <= kEpsilon) {
        out.H_pop = 0.0;
        out.H_rep_mean = 0.0;
        out.D_vib = 0.0;
        return out;
    }

    // Log-spaced bin edges: bin index for log-freq lf is
    //   idx = floor((lf - lo) / span * n_bins), clamped to [0, n_bins-1].
    const double inv_bin = static_cast<double>(n_bins) / span;
    auto bin_of = [&](double lf) -> int {
        int idx = static_cast<int>((lf - lo) * inv_bin);
        if (idx < 0) idx = 0;
        if (idx >= n_bins) idx = n_bins - 1;
        return idx;
    };

    // ── Pass 2: build per-rep histograms + accumulate the population pool ─────
    std::vector<std::vector<double>> rep_dists;
    rep_dists.reserve(n_reps);
    std::vector<double> pool(n_bins, 0.0);

    for (const auto& rep : eigenvalues) {
        std::vector<double> h(n_bins, 0.0);
        bool rep_has_mode = false;
        for (double lam : rep) {
            if (!usable(lam)) continue;
            const int b = bin_of(log_freq(lam));
            h[b] += 1.0;
            pool[b] += 1.0;     // pooled distribution uses raw mode counts
            rep_has_mode = true;
        }
        if (!rep_has_mode) continue;     // skip empty reps (already excluded from n_reps)
        normalize(h);
        rep_dists.push_back(std::move(h));
    }

    // Population distribution = normalized pooled counts (this is the H(ω) of
    // the whole generation, and the reference for D_vib).
    normalize(pool);
    out.H_pop = shannon_entropy(pool);

    // ── Per-rep entropy + divergence from the population mean ────────────────
    double sum_H_rep = 0.0;
    double sum_D_vib = 0.0;
    for (const auto& d : rep_dists) {
        sum_H_rep += shannon_entropy(d);
        sum_D_vib += symmetric_kl(d, pool);
    }
    const double denom = static_cast<double>(rep_dists.size());
    out.H_rep_mean = denom > 0.0 ? sum_H_rep / denom : 0.0;
    // Single rep ⇒ rep == population ⇒ D_vib ≈ 0 (epsilon-floor handles it).
    out.D_vib = denom > 0.0 ? sum_D_vib / denom : 0.0;

    return out;
}

} // namespace vibentropy
