// statmech.cpp — Statistical Mechanics Engine implementation
//
// Notation:
//   β  = 1/(kB T)
//   Z  = Σ_i  n_i exp(−β E_i)        (canonical partition function)
//   F  = −kT ln Z                      (Helmholtz free energy)
//   ⟨E⟩ = (1/Z) Σ_i  n_i E_i exp(−β E_i)
//   C_v = (⟨E²⟩ − ⟨E⟩²) / (kT²)      (heat capacity)
//   S  = (⟨E⟩ − F) / T                 (entropy)
//
// All sums use log-sum-exp for numerical stability when energies span
// hundreds of kcal/mol (common in docking).
//
// Hardware dispatch (runtime via UnifiedHardwareDispatch layer):
//   1. AVX-512 16-wide SIMD (+ OpenMP)
//   2. Eigen3 vectorised array ops (auto-vectorises to AVX2/AVX-512)
//   3. OpenMP parallel reductions for large ensembles
//   4. Scalar fallback (always available)

#include "statmech.h"
#include "UnifiedHardwareDispatch.h"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <stdexcept>

#include <Eigen/Dense>

#ifdef _OPENMP
#  include <omp.h>
#endif

#if defined(__AVX512F__) && defined(__AVX512DQ__)
#  include <immintrin.h>
#  define STATMECH_HAS_AVX512 1
#else
#  define STATMECH_HAS_AVX512 0
#endif

namespace statmech {

// Threshold above which OpenMP parallelisation pays off for reductions.
[[maybe_unused]] static constexpr std::size_t OMP_THRESHOLD = 4096;

// ─── construction ────────────────────────────────────────────────────────────

StatMechEngine::StatMechEngine(double temperature_K)
    : T_(temperature_K)
    , beta_(1.0 / (kB_kcal * temperature_K))
{
    if (temperature_K <= 0.0)
        throw std::invalid_argument("StatMechEngine: temperature must be > 0");
}

// ─── add_sample ──────────────────────────────────────────────────────────────

void StatMechEngine::add_sample(double energy, double multiplicity) {
    ensemble_.push_back({energy, multiplicity});
}

// ─── log_sum_exp ─────────────────────────────────────────────────────────────

double StatMechEngine::log_sum_exp(std::span<const double> x) {
    // Delegate to the unified hardware dispatch layer which handles
    // AVX-512, Eigen, OpenMP, and scalar paths with runtime selection.
    return flexaids::log_sum_exp_dispatch(x);
}

// ─── compute ─────────────────────────────────────────────────────────────────

Thermodynamics StatMechEngine::compute() const {
    if (ensemble_.empty())
        throw std::runtime_error("StatMechEngine::compute: empty ensemble");

    const std::size_t N = ensemble_.size();

    // Build array of log-weights:  w_i = ln(n_i) − β E_i
    std::vector<double> log_w(N);

    // Eigen-vectorised log-weight construction
    {
        Eigen::ArrayXd counts(static_cast<Eigen::Index>(N));
        Eigen::ArrayXd energies(static_cast<Eigen::Index>(N));
        for (std::size_t i = 0; i < N; ++i) {
            counts(static_cast<Eigen::Index>(i))   = ensemble_[i].count;
            energies(static_cast<Eigen::Index>(i)) = ensemble_[i].energy;
        }
        Eigen::ArrayXd lw = counts.log() - beta_ * energies;
        Eigen::Map<Eigen::ArrayXd>(log_w.data(), static_cast<Eigen::Index>(N)) = lw;
    }

    double lnZ = log_sum_exp(log_w);

    // ⟨E⟩  = (1/Z) Σ n_i E_i exp(−β E_i)
    //       = exp(−lnZ) Σ E_i exp(log_w_i)
    // To keep stability: ⟨E⟩ = Σ E_i exp(log_w_i − lnZ)
    double E_avg  = 0.0;
    double E2_avg = 0.0;

    // Build contiguous energy array for vectorised paths.
    std::vector<double> energies_vec(N);
    for (std::size_t i = 0; i < N; ++i)
        energies_vec[i] = ensemble_[i].energy;

#if STATMECH_HAS_AVX512
    // AVX-512 path: 8-wide fused probability × energy moment accumulation.
    if (N >= 16) {
        __m512d v_Eavg  = _mm512_setzero_pd();
        __m512d v_E2avg = _mm512_setzero_pd();
        __m512d v_lnZ   = _mm512_set1_pd(lnZ);

        std::size_t i = 0;
        for (; i + 7 < N; i += 8) {
            __m512d v_lw = _mm512_loadu_pd(log_w.data() + i);
            __m512d v_E  = _mm512_loadu_pd(energies_vec.data() + i);

            // p_i = exp(log_w_i - lnZ)
            __m512d v_arg = _mm512_sub_pd(v_lw, v_lnZ);
            alignas(64) double tmp_exp[8];
            _mm512_storeu_pd(tmp_exp, v_arg);
            for (int k = 0; k < 8; ++k) tmp_exp[k] = std::exp(tmp_exp[k]);
            __m512d v_p = _mm512_loadu_pd(tmp_exp);

            // Accumulate p_i * E_i and p_i * E_i^2
            v_Eavg  = _mm512_fmadd_pd(v_p, v_E, v_Eavg);
            v_E2avg = _mm512_fmadd_pd(v_p, _mm512_mul_pd(v_E, v_E), v_E2avg);
        }
        E_avg  = _mm512_reduce_add_pd(v_Eavg);
        E2_avg = _mm512_reduce_add_pd(v_E2avg);

        // Scalar tail
        for (; i < N; ++i) {
            double p_i = std::exp(log_w[i] - lnZ);
            double Ei  = energies_vec[i];
            E_avg  += p_i * Ei;
            E2_avg += p_i * Ei * Ei;
        }
    } else
#endif

    // Eigen vectorised path: auto-vectorises to AVX2/AVX-512 via Eigen's backend.
    if (N >= 16) {
        Eigen::Map<const Eigen::ArrayXd> lw(log_w.data(), static_cast<Eigen::Index>(N));
        Eigen::Map<const Eigen::ArrayXd> E(energies_vec.data(), static_cast<Eigen::Index>(N));

        Eigen::ArrayXd probs = (lw - lnZ).exp();
        E_avg  = (probs * E).sum();
        E2_avg = (probs * E * E).sum();
    } else
    {
        for (std::size_t i = 0; i < N; ++i) {
            double p_i = std::exp(log_w[i] - lnZ);
            double Ei  = energies_vec[i];
            E_avg  += p_i * Ei;
            E2_avg += p_i * Ei * Ei;
        }
    }

    double kT  = kB_kcal * T_;
    double var = E2_avg - E_avg * E_avg;

    Thermodynamics th;
    th.temperature    = T_;
    th.log_Z          = lnZ;
    th.free_energy    = -kT * lnZ;
    th.mean_energy    = E_avg;
    th.mean_energy_sq = E2_avg;
    th.heat_capacity  = var / (kB_kcal * T_ * T_);
    th.entropy        = (E_avg - th.free_energy) / T_;
    th.std_energy     = std::sqrt(std::max(0.0, var));
    return th;
}

ThermodynamicBreakdown StatMechEngine::compute_breakdown(
    double G_vib_kcal_mol,
    double G_natural_kcal_mol,
    double G_other_kcal_mol,
    bool has_vib,
    bool has_natural,
    bool has_other) const
{
    const Thermodynamics th = compute();

    ThermodynamicBreakdown b;
    b.temperature_K = th.temperature;
    b.logZ_config = th.log_Z;
    b.G_config_kcal_mol = th.free_energy;
    b.H_eff_kcal_mol = th.mean_energy;
    b.S_config_kcal_mol_K = th.entropy;
    b.minus_T_S_config_kcal_mol = th.free_energy - th.mean_energy;
    b.Cv_kcal_mol_K = th.heat_capacity;
    b.sigma_E_kcal_mol = th.std_energy;
    b.G_vib_kcal_mol = G_vib_kcal_mol;
    b.G_natural_kcal_mol = G_natural_kcal_mol;
    b.G_other_kcal_mol = G_other_kcal_mol;
    b.G_total_kcal_mol = b.G_config_kcal_mol
                       + b.G_vib_kcal_mol
                       + b.G_natural_kcal_mol
                       + b.G_other_kcal_mol;
    b.has_vib = has_vib;
    b.has_natural = has_natural;
    b.has_other = has_other;
    return b;
}

// ─── boltzmann_weights ───────────────────────────────────────────────────────

std::vector<double> StatMechEngine::boltzmann_weights() const {
    if (ensemble_.empty()) return {};

    const std::size_t N = ensemble_.size();

    // Build raw energy array and use the unified dispatch layer.
    std::vector<double> energies(N);
    for (std::size_t i = 0; i < N; ++i)
        energies[i] = ensemble_[i].energy;

    auto result = flexaids::compute_boltzmann_batch(energies, beta_);

    // Normalise weights accounting for multiplicities.
    std::vector<double> w(N);
    double Z_with_mult = 0.0;
    for (std::size_t i = 0; i < N; ++i)
        Z_with_mult += ensemble_[i].count * result.weights[i];

    if (Z_with_mult > 0.0) {
        for (std::size_t i = 0; i < N; ++i)
            w[i] = ensemble_[i].count * result.weights[i] / Z_with_mult;
    }
    return w;
}

// ─── delta_G ─────────────────────────────────────────────────────────────────
// ΔG = F_this − F_ref = −kT (ln Z_this − ln Z_ref)

double StatMechEngine::delta_G(const StatMechEngine& reference) const {
    auto this_th = this->compute();
    auto ref_th  = reference.compute();
    double kT = kB_kcal * T_;
    return -kT * (this_th.log_Z - ref_th.log_Z);
}

// ─── Helmholtz convenience ───────────────────────────────────────────────────

double StatMechEngine::helmholtz(std::span<const double> energies, double T) {
    if (energies.empty())
        throw std::invalid_argument("helmholtz: empty energy list");
    double beta = 1.0 / (kB_kcal * T);

    // Use unified dispatch for the Boltzmann batch computation.
    auto result = flexaids::compute_boltzmann_batch(energies, beta);
    // F = -kT * ln(Z) where log_Z already accounts for E_min shift.
    return -(kB_kcal * T) * result.log_Z;
}

// ─── replica exchange ────────────────────────────────────────────────────────

std::vector<Replica>
StatMechEngine::init_replicas(std::span<const double> temperatures) {
    std::vector<Replica> reps;
    reps.reserve(temperatures.size());
    int id = 0;
    for (double T : temperatures) {
        Replica r;
        r.id             = id++;
        r.temperature    = T;
        r.beta           = 1.0 / (kB_kcal * T);
        r.current_energy = 0.0;
        reps.push_back(r);
    }
    return reps;
}

bool StatMechEngine::attempt_swap(Replica& a, Replica& b, std::mt19937& rng) {
    // Metropolis criterion:
    //   Δ = (β_a − β_b)(E_a − E_b)
    //   P_accept = min(1, exp(Δ))
    double delta = (a.beta - b.beta) * (a.current_energy - b.current_energy);
    if (delta >= 0.0) {
        std::swap(a.current_energy, b.current_energy);
        return true;
    }
    std::uniform_real_distribution<double> U(0.0, 1.0);
    if (U(rng) < std::exp(delta)) {
        std::swap(a.current_energy, b.current_energy);
        return true;
    }
    return false;
}

// ─── Boltzmann-reweighted PMF ────────────────────────────────────────────────
// Single-window post-hoc reweighting of an ensemble onto a 1D collective
// coordinate. NOT multi-window WHAM (Kumar et al. 1992) — that requires
// biased simulations with per-window offsets, neither of which are
// available here. The historical name `wham()` survives as a deprecated
// alias in the header (see statmech.h).

std::vector<WHAMBin> StatMechEngine::boltzmann_pmf(
    std::span<const double> energies,
    std::span<const double> coordinates,
    double temperature,
    int    n_bins,
    int    max_iter,
    double tolerance)
{
    if (energies.size() != coordinates.size())
        throw std::invalid_argument("boltzmann_pmf: energies and coordinates size mismatch");
    if (energies.empty() || n_bins <= 0)
        throw std::invalid_argument("boltzmann_pmf: invalid input");

    const std::size_t N = energies.size();
    double beta = 1.0 / (kB_kcal * temperature);

    // Find coordinate range (single pass)
    auto [cmin_it, cmax_it] = std::minmax_element(coordinates.begin(), coordinates.end());
    double cmin = *cmin_it;
    double cmax = *cmax_it;
    double bin_w = (cmax - cmin) / n_bins;
    if (bin_w <= 0.0) bin_w = 1.0;

    // Histogram + Boltzmann-weighted histogram.
    // Use the dispatch layer for the Boltzmann weight computation,
    // then bin the pre-computed weights for O(N) histogramming.
    auto boltz_result = flexaids::compute_boltzmann_batch(energies, beta);

    std::vector<double> raw_count(static_cast<std::size_t>(n_bins), 0.0);
    std::vector<double> boltz_sum(static_cast<std::size_t>(n_bins), 0.0);
    [[maybe_unused]] double inv_bw = 1.0 / bin_w;

#ifdef _OPENMP
    // OpenMP parallel histogram with per-thread private bins
    if (N >= OMP_THRESHOLD) {
        int n_threads = omp_get_max_threads();
        std::vector<std::vector<double>> t_raw(n_threads,
            std::vector<double>(static_cast<std::size_t>(n_bins), 0.0));
        std::vector<std::vector<double>> t_boltz(n_threads,
            std::vector<double>(static_cast<std::size_t>(n_bins), 0.0));

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < static_cast<int>(N); ++i) {
            int tid = omp_get_thread_num();
            int b = static_cast<int>((coordinates[i] - cmin) * inv_bw);
            b = std::min(std::max(b, 0), n_bins - 1);
            t_raw[tid][static_cast<std::size_t>(b)]  += 1.0;
            t_boltz[tid][static_cast<std::size_t>(b)] += boltz_result.weights[i];
        }
        // Reduce thread-private histograms
        for (auto& tr : t_raw)
            for (int b = 0; b < n_bins; ++b)
                raw_count[static_cast<std::size_t>(b)] += tr[static_cast<std::size_t>(b)];
        for (auto& tb : t_boltz)
            for (int b = 0; b < n_bins; ++b)
                boltz_sum[static_cast<std::size_t>(b)] += tb[static_cast<std::size_t>(b)];
    } else
#endif
    {
        for (std::size_t i = 0; i < N; ++i) {
            int b = static_cast<int>((coordinates[i] - cmin) / bin_w);
            if (b < 0) b = 0;
            if (b >= n_bins) b = n_bins - 1;
            raw_count[static_cast<std::size_t>(b)] += 1.0;
            boltz_sum[static_cast<std::size_t>(b)] += boltz_result.weights[i];
        }
    }

    // Free energy per bin: F_b = −kT ln( weighted_count_b / raw_count_b )
    // Iterative self-consistency (single-window simplification)
    std::vector<double> f_old(static_cast<std::size_t>(n_bins), 0.0);
    std::vector<double> f_new(static_cast<std::size_t>(n_bins), 0.0);

    for (int iter = 0; iter < max_iter; ++iter) {
        // Eigen vectorised WHAM self-consistency update
        Eigen::Map<const Eigen::ArrayXd> rc(raw_count.data(), n_bins);
        Eigen::Map<const Eigen::ArrayXd> bs(boltz_sum.data(), n_bins);
        Eigen::Map<Eigen::ArrayXd> fn(f_new.data(), n_bins);
        Eigen::Map<Eigen::ArrayXd> fo(f_old.data(), n_bins);

        // F_b = -kT * ln(boltz_sum_b / raw_count_b) where raw_count_b > 0
        auto occupied = (rc > 0.0);
        Eigen::ArrayXd safe_rc = occupied.select(rc, Eigen::ArrayXd::Ones(n_bins));
        fn = occupied.select(
            -(kB_kcal * temperature) * (bs / safe_rc).log(),
            Eigen::ArrayXd::Zero(n_bins));

        // Shift so minimum = 0
        fn -= fn.minCoeff();

        // Check convergence
        {
            double maxdiff = (fn - fo).abs().maxCoeff();
            fo = fn;
            if (maxdiff < tolerance) break;
        }
    }

    // Build output
    std::vector<WHAMBin> result(static_cast<std::size_t>(n_bins));
    for (int b = 0; b < n_bins; ++b) {
        result[static_cast<std::size_t>(b)].coord_center = cmin + (b + 0.5) * bin_w;
        result[static_cast<std::size_t>(b)].count        = raw_count[static_cast<std::size_t>(b)];
        result[static_cast<std::size_t>(b)].free_energy  = f_new[static_cast<std::size_t>(b)];
    }
    return result;
}

// ─── thermodynamic integration ───────────────────────────────────────────────
// ΔG = ∫₀¹ ⟨∂V/∂λ⟩_λ dλ   (trapezoidal rule)

double StatMechEngine::thermodynamic_integration(std::span<const TIPoint> points) {
    if (points.size() < 2)
        throw std::invalid_argument("TI requires at least 2 points");

    double integral = 0.0;
    for (std::size_t i = 1; i < points.size(); ++i) {
        double dl = points[i].lambda - points[i-1].lambda;
        integral += 0.5 * dl * (points[i].dV_dlambda + points[i-1].dV_dlambda);
    }
    return integral;
}

// ─── BoltzmannLUT ────────────────────────────────────────────────────────────

BoltzmannLUT::BoltzmannLUT(double beta, double e_min, double e_max, int n_bins)
    : beta_(beta)
    , e_min_(e_min)
    , n_bins_(n_bins)
    , table_(static_cast<std::size_t>(n_bins))
{
    double range = e_max - e_min;
    if (range <= 0.0) range = 1.0;
    inv_bin_width_ = n_bins / range;

    // Eigen-vectorised LUT initialisation
    {
        Eigen::ArrayXd idx = Eigen::ArrayXd::LinSpaced(n_bins, 0.5, n_bins - 0.5);
        Eigen::ArrayXd E = e_min + idx * (range / n_bins);
        Eigen::Map<Eigen::ArrayXd>(table_.data(), n_bins) = (-beta * E).exp();
    }
}

double BoltzmannLUT::operator()(double energy) const noexcept {
    int idx = static_cast<int>((energy - e_min_) * inv_bin_width_);
    if (idx < 0) idx = 0;
    if (idx >= n_bins_) idx = n_bins_ - 1;
    return table_[static_cast<std::size_t>(idx)];
}

// ─── ensemble merging (parallel grid-decomposed docking) ─────────────────────

void StatMechEngine::merge(const StatMechEngine& other) {
    if (std::fabs(other.T_ - T_) > 1e-6)
        throw std::invalid_argument("Cannot merge engines at different temperatures");
    ensemble_.insert(ensemble_.end(),
                     other.ensemble_.begin(), other.ensemble_.end());
}

void StatMechEngine::merge_samples(std::span<const double> energies,
                                    std::span<const double> multiplicities) {
    if (energies.size() != multiplicities.size())
        throw std::invalid_argument("energies and multiplicities must have same size");
    for (size_t i = 0; i < energies.size(); ++i)
        ensemble_.push_back({energies[i], multiplicities[i]});
}

std::vector<double> StatMechEngine::serialize_energies() const {
    std::vector<double> out(ensemble_.size());
    for (size_t i = 0; i < ensemble_.size(); ++i)
        out[i] = ensemble_[i].energy;
    return out;
}

std::vector<double> StatMechEngine::serialize_multiplicities() const {
    std::vector<double> out(ensemble_.size());
    for (size_t i = 0; i < ensemble_.size(); ++i)
        out[i] = ensemble_[i].count;
    return out;
}

}  // namespace statmech
