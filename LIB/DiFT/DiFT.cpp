// DiFT.cpp — Discrete Fourier Transform torsional parametrization for FlexAID∆S
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
//
// Implementation of the DiFT engine. Self-contained: depends only on the
// C++ standard library (no GPL, no external FFT dependency — the radix-2
// transform is implemented here for licensing cleanliness).
//
// See DiFT.h for the method, references, and the energy/entropy rationale.

#include "DiFT.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <numeric>
#include <stdexcept>

namespace dift {

namespace {

constexpr double kTwoPi = 6.283185307179586476925286766559;

// True iff n is a power of two and n ≥ 2.
inline bool is_pow2(std::size_t n) noexcept {
    return n >= 2 && (n & (n - 1)) == 0;
}

// In-place iterative radix-2 Cooley–Tukey FFT.
// Sign convention:  X_n = Σ_k x_k · exp(−i·2π·n·k/N)   — O(N log N).
void fft_radix2(std::vector<std::complex<double>>& a) {
    const std::size_t n = a.size();
    if (n < 2) return;

    // Bit-reversal permutation.
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }

    for (std::size_t len = 2; len <= n; len <<= 1) {
        const double ang = -kTwoPi / static_cast<double>(len);
        const std::complex<double> wlen(std::cos(ang), std::sin(ang));
        for (std::size_t i = 0; i < n; i += len) {
            std::complex<double> w(1.0, 0.0);
            for (std::size_t k = 0; k < (len >> 1); ++k) {
                const std::complex<double> u = a[i + k];
                const std::complex<double> v = a[i + k + (len >> 1)] * w;
                a[i + k]               = u + v;
                a[i + k + (len >> 1)]  = u - v;
                w *= wlen;
            }
        }
    }
}

// Exact direct DFT for arbitrary (non-power-of-two) lengths — O(N²).
// Torsional grids are small (typically 24–72 points), so this is negligible
// and keeps the transform mathematically exact rather than zero-padding a
// periodic signal (which would distort its spectrum).
std::vector<std::complex<double>> dft_direct(std::span<const double> x) {
    const std::size_t M = x.size();
    std::vector<std::complex<double>> X(M);
    for (std::size_t n = 0; n < M; ++n) {
        std::complex<double> s(0.0, 0.0);
        for (std::size_t k = 0; k < M; ++k) {
            const double ang = -kTwoPi * static_cast<double>(n) *
                               static_cast<double>(k) / static_cast<double>(M);
            s += x[k] * std::complex<double>(std::cos(ang), std::sin(ang));
        }
        X[n] = s;
    }
    return X;
}

// Wrap an angle into (−π, π].
inline double wrap_pi(double a) noexcept {
    a = std::fmod(a + M_PI, kTwoPi);
    if (a < 0.0) a += kTwoPi;
    return a - M_PI;
}

} // namespace

// ─────────────────────────────────────────────────────────────────────────────
DiFTEngine::DiFTEngine(double temperature_K) {
    set_temperature(temperature_K);
}

void DiFTEngine::set_temperature(double T_K) noexcept {
    T_    = (T_K > 0.0) ? T_K : 300.0;
    beta_ = 1.0 / (kB_kcal * T_);
}

// ─────────────────────────────────────────────────────────────────────────────
std::vector<FourierTerm>
DiFTEngine::transform(std::span<const double> profile, double& mean_out) const {
    const std::size_t M = profile.size();
    if (M < 2)
        throw std::invalid_argument("DiFT::transform: profile needs ≥ 2 samples");
    for (double v : profile) {
        if (!std::isfinite(v))
            throw std::invalid_argument("DiFT::transform: profile contains non-finite values");
    }

    // Forward transform → full complex spectrum X_0 … X_{M−1}.
    std::vector<std::complex<double>> X;
    if (is_pow2(M)) {
        X.resize(M);
        for (std::size_t k = 0; k < M; ++k) X[k] = std::complex<double>(profile[k], 0.0);
        fft_radix2(X);
    } else {
        X = dft_direct(profile);
    }

    const double invM = 1.0 / static_cast<double>(M);

    // DC component → mean offset.
    mean_out = X[0].real() * invM;

    // Real-signal reconstruction:
    //   x_k = mean + Σ_{n=1}^{M/2} Aₙ·cos(n·φ_k − ωₙ),  φ_k = 2πk/M
    // with, for 1 ≤ n < M/2:  Aₙ = (2/M)|Xₙ|,  ωₙ = −arg(Xₙ).
    // The Nyquist term (n = M/2, M even) carries weight 1/M.
    const std::size_t n_max = M / 2;
    std::vector<FourierTerm> spectrum;
    spectrum.reserve(n_max);

    for (std::size_t n = 1; n <= n_max; ++n) {
        const bool   nyquist = (M % 2 == 0) && (n == n_max);
        const double scale   = nyquist ? invM : (2.0 * invM);
        double amp           = scale * std::abs(X[n]);
        double phase         = -std::atan2(X[n].imag(), X[n].real());

        // Normalize to amplitude ≥ 0 (fold a negative sign into a π phase flip).
        if (amp < 0.0) { amp = -amp; phase = wrap_pi(phase + M_PI); }
        else            { phase = wrap_pi(phase); }

        FourierTerm t;
        t.multiplicity = static_cast<int>(n);
        t.amplitude    = amp;
        t.phase        = phase;
        t.power        = amp * amp;
        spectrum.push_back(t);
    }
    return spectrum;
}

// ─────────────────────────────────────────────────────────────────────────────
TorsionalPotential
DiFTEngine::parametrize(std::span<const double> profile, int max_multiplicity) const {
    double mean = 0.0;
    std::vector<FourierTerm> spectrum = transform(profile, mean);

    // Spectral Shannon entropy of the FULL spectrum (before any truncation).
    const double H_spec = spectral_entropy(spectrum);
    const double N_eff  = std::exp(H_spec);

    // Anti-overfit guard: drop frequencies above the requested cap.
    if (max_multiplicity > 0) {
        std::erase_if(spectrum, [max_multiplicity](const FourierTerm& t) {
            return t.multiplicity > max_multiplicity;
        });
    }

    // Shannon-collapse truncation: keep the ⌈N_eff⌉ highest-power terms.
    // The spectrum's own entropy decides the term count — no user threshold.
    std::size_t keep = static_cast<std::size_t>(std::ceil(N_eff));
    keep = std::clamp<std::size_t>(keep, 1, spectrum.size());

    std::sort(spectrum.begin(), spectrum.end(),
              [](const FourierTerm& a, const FourierTerm& b) {
                  return a.power > b.power;
              });
    if (spectrum.size() > keep) spectrum.resize(keep);

    // Retained terms back in ascending-multiplicity order (canonical form).
    std::sort(spectrum.begin(), spectrum.end(),
              [](const FourierTerm& a, const FourierTerm& b) {
                  return a.multiplicity < b.multiplicity;
              });

    TorsionalPotential pot;
    pot.terms            = std::move(spectrum);
    pot.mean             = mean;
    pot.n_samples        = static_cast<int>(profile.size());
    pot.spectral_entropy = H_spec;
    pot.effective_modes  = N_eff;
    pot.refinement_iters = 0;

    // Fit quality of the truncated model vs the input profile.
    std::vector<double> model = pot.sample(static_cast<int>(profile.size()));
    pot.r_squared = r_squared(profile, model);

    // Cache the global minimum (used for relative energies / scoring).
    std::vector<double> fine = pot.sample(720);
    pot.v_min = *std::min_element(fine.begin(), fine.end());

    return pot;
}

// ─────────────────────────────────────────────────────────────────────────────
TorsionalPotential
DiFTEngine::refine(std::span<const double> qm,
                   std::span<const double> mm_initial,
                   double lambda, double r2_target,
                   int max_iter, int max_multiplicity) const {
    const std::size_t M = qm.size();
    if (M < 2 || mm_initial.size() != M)
        throw std::invalid_argument(
            "DiFT::refine: qm and mm_initial must share a grid of ≥ 2 samples");

    // The torsional correction satisfies  mm_initial + V_t ≈ qm.
    // We approach it through the damped, FFT band-limited loop of eq. 18.
    std::vector<double> correction(M, 0.0);
    std::vector<double> mm_current(M, 0.0);
    TorsionalPotential  pot;
    int    iters = 0;
    double r2    = 0.0;

    for (int it = 1; it <= max_iter; ++it) {
        iters = it;
        // D_i = V_QM − V_MM,i
        for (std::size_t k = 0; k < M; ++k)
            mm_current[k] = mm_initial[k] + correction[k];

        // V_{i+1} = V_i + λ·D_i
        for (std::size_t k = 0; k < M; ++k)
            correction[k] += lambda * (qm[k] - mm_current[k]);

        // Band-limit the correction through the FFT (paper: FFT applied to V_i).
        pot        = parametrize(correction, max_multiplicity);
        correction = pot.sample(static_cast<int>(M));

        // Convergence: does mm_initial + correction reproduce qm?
        for (std::size_t k = 0; k < M; ++k)
            mm_current[k] = mm_initial[k] + correction[k];
        r2 = r_squared(qm, mm_current);
        if (r2 >= r2_target) break;
    }

    pot.r_squared        = r2;
    pot.refinement_iters = iters;
    return pot;
}

// ─────────────────────────────────────────────────────────────────────────────
std::vector<double>
DiFTEngine::boltzmann_invert(std::span<const double> histogram) const {
    const std::size_t M = histogram.size();
    if (M < 2)
        throw std::invalid_argument("DiFT::boltzmann_invert: need ≥ 2 bins");

    double total = 0.0;
    for (double c : histogram) {
        if (!std::isfinite(c))
            throw std::invalid_argument("DiFT::boltzmann_invert: non-finite count");
        if (c < 0.0)
            throw std::invalid_argument("DiFT::boltzmann_invert: negative count");
        total += c;
    }
    if (total <= 0.0)
        throw std::invalid_argument("DiFT::boltzmann_invert: empty histogram");

    // E(φ) = −kT ln p(φ). No Jacobian: dihedral angles are uniformly
    // distributed for a free rotor (unlike bond lengths or bond angles).
    const double kT = kB_kcal * T_;
    std::vector<double> energy(M);
    double e_min =  1e300;
    double e_max = -1e300;
    for (std::size_t k = 0; k < M; ++k) {
        if (histogram[k] > 0.0) {
            energy[k] = -kT * std::log(histogram[k] / total);
            e_min = std::min(e_min, energy[k]);
            e_max = std::max(e_max, energy[k]);
        } else {
            energy[k] = std::nan("");   // flag empty bins; filled below
        }
    }
    // Empty bins → capped well wall at the highest observed energy.
    for (double& e : energy)
        if (std::isnan(e)) e = e_max;
    // Shift so the global minimum sits at 0 (relative energies).
    for (double& e : energy) e -= e_min;
    return energy;
}

// ─────────────────────────────────────────────────────────────────────────────
TorsionalThermo
DiFTEngine::thermodynamics(const TorsionalPotential& pot) const {
    // Integrate the 1-D torsional partition function on a fine grid using the
    // analytical potential. Energies are taken relative to the global minimum
    // so the result is the EXCESS thermodynamics vs an (unconfined) free rotor.
    constexpr int Nf = 1440;
    double Z = 0.0, EZ = 0.0;
    for (int j = 0; j < Nf; ++j) {
        const double phi = kTwoPi * static_cast<double>(j) / static_cast<double>(Nf);
        const double V   = pot.evaluate(phi) - pot.v_min;   // ≥ 0
        const double w   = std::exp(-beta_ * V);
        Z  += w;
        EZ += w * V;
    }
    const double invNf = 1.0 / static_cast<double>(Nf);
    const double z     = Z * invNf;                 // ⟨exp(−βV)⟩; 1 for a free rotor
    const double mean  = (Z > 0.0) ? (EZ / Z) : 0.0;
    const double F     = -kB_kcal * T_ * std::log(z);

    TorsionalThermo th;
    th.temperature        = T_;
    th.partition_function = z;
    th.free_energy        = F;
    th.mean_energy        = mean;
    th.entropy            = (mean - F) / T_;        // ≤ 0: confinement = entropy loss
    th.minus_TS           = -T_ * th.entropy;       // ≥ 0: the ΔG penalty
    return th;
}

// ─────────────────────────────────────────────────────────────────────────────
double DiFTEngine::circular_mean(std::span<const double> angles) noexcept {
    double s = 0.0, c = 0.0;
    for (double a : angles) { s += std::sin(a); c += std::cos(a); }
    if (s == 0.0 && c == 0.0) return 0.0;
    return std::atan2(s, c);
}

double DiFTEngine::r_squared(std::span<const double> observed,
                             std::span<const double> model) noexcept {
    const std::size_t n = observed.size();
    if (n == 0 || model.size() != n) return 0.0;

    double mean = 0.0;
    for (double v : observed) mean += v;
    mean /= static_cast<double>(n);

    double ss_res = 0.0, ss_tot = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double d_res = observed[i] - model[i];
        const double d_tot = observed[i] - mean;
        ss_res += d_res * d_res;
        ss_tot += d_tot * d_tot;
    }
    if (ss_tot <= 0.0)                     // flat target
        return (ss_res <= 1e-18) ? 1.0 : 0.0;
    return 1.0 - ss_res / ss_tot;
}

// ─────────────────────────────────────────────────────────────────────────────
double TorsionalPotential::evaluate(double phi) const noexcept {
    double v = mean;
    for (const FourierTerm& t : terms)
        v += t.amplitude *
             std::cos(static_cast<double>(t.multiplicity) * phi - t.phase);
    return v;
}

std::vector<double> TorsionalPotential::sample(int n) const {
    std::vector<double> out;
    if (n < 1) return out;
    out.reserve(static_cast<std::size_t>(n));
    for (int j = 0; j < n; ++j)
        out.push_back(evaluate(kTwoPi * static_cast<double>(j) /
                                static_cast<double>(n)));
    return out;
}

// ─────────────────────────────────────────────────────────────────────────────
double spectral_entropy(std::span<const FourierTerm> spectrum) noexcept {
    double total = 0.0;
    for (const FourierTerm& t : spectrum) total += t.power;
    if (total <= 0.0) return 0.0;

    double H = 0.0;
    for (const FourierTerm& t : spectrum) {
        const double p = t.power / total;
        if (p > 0.0) H -= p * std::log(p);
    }
    return H;
}

} // namespace dift
