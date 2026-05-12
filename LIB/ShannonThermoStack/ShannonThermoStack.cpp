// ShannonThermoStack.cpp — multi-path implementation
//
// Hardware dispatch priority (runtime):
//   1. CUDA GPU   (FLEXAIDS_USE_CUDA)
//   2. Metal GPU  (FLEXAIDS_HAS_METAL_SHANNON, Apple Silicon)
//   3. AVX-512    (__AVX512F__)  — 8 doubles/cycle histogram binning
//   4. OpenMP     (_OPENMP)
//   5. Scalar     (always available)
//
// Eigen is used for vectorised log() / probability array ops on all CPU paths.
#include "ShannonThermoStack.h"

#ifdef FLEXAIDS_HAS_METAL_SHANNON
#  include "ShannonMetalBridge.h"
#endif

#ifdef FLEXAIDS_USE_CUDA
#  include "shannon_cuda.cuh"
static ShannonCudaCtx s_cuda_ctx;
static bool           s_cuda_ready = false;
#endif

#ifdef __AVX512F__
#  include <immintrin.h>
#endif

#ifdef _OPENMP
#  include <omp.h>
#endif

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace shannon_thermo {

// ─── ShannonEnergyMatrix ─────────────────────────────────────────────────────

ShannonEnergyMatrix& ShannonEnergyMatrix::instance() {
    static ShannonEnergyMatrix inst;
    return inst;
}

void ShannonEnergyMatrix::initialise() {
    // Thread-safe one-shot initialisation. Multiple threads may call this
    // concurrently before the matrix is built; std::call_once ensures the
    // body executes exactly once, while readers see initialised_ = true
    // only after the matrix has been fully populated.
    std::call_once(init_once_, [this]() {
        // No mutex needed: call_once guarantees single-thread execution,
        // and lookup() does not acquire mtx_ so the guard provides no
        // additional protection against concurrent readers.
        matrix_.resize(SHANNON_BINS * SHANNON_BINS);

        std::mt19937 rng(42);
        std::normal_distribution<double> perturb(0.0, 0.05);
        const double base = 1.0 / SHANNON_BINS;
        std::vector<double> p_i(SHANNON_BINS), p_j(SHANNON_BINS);
        for (int k = 0; k < SHANNON_BINS; ++k) {
            p_i[k] = std::max(1e-9, base + perturb(rng));
            p_j[k] = std::max(1e-9, base + perturb(rng));
        }
        double si = 0, sj = 0;
        for (int k = 0; k < SHANNON_BINS; ++k) { si += p_i[k]; sj += p_j[k]; }
        for (int k = 0; k < SHANNON_BINS; ++k) { p_i[k] /= si; p_j[k] /= sj; }

        const double kT = kB_kcal * TEMPERATURE_K;
        // Fill the entropy matrix using natural log (nats)
        for (int i = 0; i < SHANNON_BINS; ++i)
            for (int j = 0; j < SHANNON_BINS; ++j)
                matrix_[i * SHANNON_BINS + j] = -kT * p_i[i] * std::log(p_j[j]);
        initialised_ = true;
    });
}

bool ShannonEnergyMatrix::initialise_from_file(const std::string& path) {
    FILE* fp = fopen(path.c_str(), "rb");
    if (!fp) return false;

    // Read and verify SHNN magic header
    char magic[4];
    if (fread(magic, 1, 4, fp) != 4 ||
        magic[0] != 'S' || magic[1] != 'H' ||
        magic[2] != 'N' || magic[3] != 'N') {
        fclose(fp);
        return false;
    }

    uint32_t version = 0, dim = 0;
    if (fread(&version, sizeof(uint32_t), 1, fp) != 1 ||
        fread(&dim, sizeof(uint32_t), 1, fp) != 1) {
        fclose(fp);
        return false;
    }

    if (static_cast<int>(dim) != SHANNON_BINS) {
        fclose(fp);
        return false;
    }

    std::vector<float> buf(SHANNON_BINS * SHANNON_BINS);
    size_t nread = fread(buf.data(), sizeof(float),
                         SHANNON_BINS * SHANNON_BINS, fp);
    fclose(fp);
    if (static_cast<int>(nread) != SHANNON_BINS * SHANNON_BINS)
        return false;

    // Lock against races with initialise() / initialise_from_data()
    std::lock_guard<std::mutex> lk(mtx_);
    matrix_.resize(SHANNON_BINS * SHANNON_BINS);
    for (int k = 0; k < SHANNON_BINS * SHANNON_BINS; ++k)
        matrix_[k] = static_cast<double>(buf[k]);

    initialised_ = true;
    // Mark the once_flag as fired so initialise() becomes a no-op.
    std::call_once(init_once_, [](){});
    return true;
}

void ShannonEnergyMatrix::initialise_from_data(const float* data, int count) {
    std::lock_guard<std::mutex> lk(mtx_);
    int expected = SHANNON_BINS * SHANNON_BINS;
    matrix_.resize(expected);
    int safe_count = std::min(count, expected);
    for (int k = 0; k < safe_count; ++k)
        matrix_[k] = static_cast<double>(data[k]);
    for (int k = safe_count; k < expected; ++k)
        matrix_[k] = 0.0;
    initialised_ = true;
    std::call_once(init_once_, [](){});
}


// ─── entropy from bin counts (Eigen-vectorised) ───────────────────────────────
static double entropy_from_counts(const int* counts, int num_bins, int total) {
    if (total == 0) return 0.0;

    Eigen::ArrayXd prob(num_bins);
    for (int b = 0; b < num_bins; ++b)
        prob(b) = static_cast<double>(counts[b]);
    prob /= static_cast<double>(total);
    // Mask zeros before log to avoid -inf; Eigen evaluates log vectorised
    Eigen::ArrayXd safe_p = (prob > 1e-15).select(prob, Eigen::ArrayXd::Constant(num_bins, 1.0));
    Eigen::ArrayXd lp     = (prob > 1e-15).select(safe_p.log(), Eigen::ArrayXd::Zero(num_bins));
    return -(prob * lp).sum();
}

// ─── AVX-512 private histogram ────────────────────────────────────────────────
#ifdef __AVX512F__
static void histogram_avx512(const double* values, int n,
                               double min_v, double inv_bw, int num_bins,
                               std::vector<int>& priv)
{
    int i = 0;
    __m512d vmin   = _mm512_set1_pd(min_v);
    __m512d vinvbw = _mm512_set1_pd(inv_bw);

    for (; i + 7 < n; i += 8) {
        __m512d ve   = _mm512_loadu_pd(values + i);
        __m512d vrel = _mm512_fmadd_pd(ve, vinvbw,
                           _mm512_mul_pd(_mm512_set1_pd(-min_v), vinvbw));
        // Convert 8 doubles → 8 int32 (truncate)
        __m256i v32 = _mm512_cvttpd_epi32(vrel);
        alignas(32) int tmp[8];
        _mm256_storeu_si256((__m256i*)tmp, v32);
        for (int k = 0; k < 8; ++k) {
            int b = std::min(std::max(tmp[k], 0), num_bins - 1);
            priv[b]++;
        }
    }
    for (; i < n; ++i) {
        int b = static_cast<int>((values[i] - min_v) * inv_bw);
        b = std::min(std::max(b, 0), num_bins - 1);
        priv[b]++;
    }
}
#endif // __AVX512F__

// ─── compute_shannon_entropy ─────────────────────────────────────────────────

double compute_shannon_entropy(const std::vector<double>& values, int num_bins) {
    if (values.empty()) return 0.0;
    if (num_bins <= 0)  num_bins = DEFAULT_HIST_BINS;

    auto [it_min, it_max] = std::minmax_element(values.begin(), values.end());
    double min_v = *it_min;
    double max_v = *it_max;
    if (max_v - min_v < 1e-12) return 0.0;
    double bin_width = (max_v - min_v) / num_bins + 1e-10;
    double inv_bw    = 1.0 / bin_width;
    int    n         = static_cast<int>(values.size());

    std::vector<int> bins(num_bins, 0);

// GPU dispatch only for large datasets (N > GPU_DISPATCH_THRESHOLD).
// For typical docking populations (100–10K), scalar/OpenMP is faster
// than GPU kernel launch + memory transfer overhead.
if (n > GPU_DISPATCH_THRESHOLD) {
// ── 1. CUDA ───────────────────────────────────────────────────────────────────
#ifdef FLEXAIDS_USE_CUDA
    {
        if (!s_cuda_ready) {
            shannon_cuda_init(s_cuda_ctx, 1 << 20, num_bins);
            s_cuda_ready = true;
        }
        if (n <= s_cuda_ctx.capacity) {
            shannon_cuda_histogram(s_cuda_ctx, values.data(), n,
                                   min_v, bin_width, bins.data());
            return entropy_from_counts(bins.data(), num_bins, n);
        }
    }
#endif

// ── 2. Metal ──────────────────────────────────────────────────────────────────
#ifdef FLEXAIDS_HAS_METAL_SHANNON
    return ShannonMetalBridge::compute_shannon_entropy_metal(values, num_bins);
#endif
} // GPU_DISPATCH_THRESHOLD

// ── 3. AVX-512 (+ optional OpenMP for multi-threaded private histograms) ──────
#ifdef __AVX512F__
    {
#  ifdef _OPENMP
        int n_threads = omp_get_max_threads();
        // Flat Eigen matrix for cache-friendly thread-local histograms
        Eigen::MatrixXi t_bins = Eigen::MatrixXi::Zero(n_threads, num_bins);
        #pragma omp parallel
        {
            int tid    = omp_get_thread_num();
            int chunk  = (n + n_threads - 1) / n_threads;
            int start  = tid * chunk;
            int end_i  = std::min(start + chunk, n);
            if (start < end_i) {
                // histogram_avx512 needs std::vector<int>& — use row view
                std::vector<int> priv(num_bins, 0);
                histogram_avx512(values.data() + start, end_i - start,
                                 min_v, inv_bw, num_bins, priv);
                for (int b = 0; b < num_bins; ++b)
                    t_bins(tid, b) = priv[b];
            }
        }
        Eigen::VectorXi col_sums = t_bins.colwise().sum();
        for (int b = 0; b < num_bins; ++b) bins[b] = col_sums(b);
#  else
        histogram_avx512(values.data(), n, min_v, inv_bw, num_bins, bins);
#  endif
    }

// ── 4. OpenMP scalar ──────────────────────────────────────────────────────────
#elif defined(_OPENMP)
    {
        int n_threads = omp_get_max_threads();
        Eigen::MatrixXi t_bins = Eigen::MatrixXi::Zero(n_threads, num_bins);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i) {
            int tid = omp_get_thread_num();
            int b   = static_cast<int>((values[i] - min_v) * inv_bw);
            t_bins(tid, std::min(std::max(b, 0), num_bins - 1))++;
        }
        Eigen::VectorXi col_sums = t_bins.colwise().sum();
        for (int b = 0; b < num_bins; ++b) bins[b] = col_sums(b);
    }

// ── 5. Scalar ─────────────────────────────────────────────────────────────────
#else
    for (int i = 0; i < n; ++i) {
        int b = static_cast<int>((values[i] - min_v) * inv_bw);
        bins[std::min(std::max(b, 0), num_bins - 1)]++;
    }
#endif

    return entropy_from_counts(bins.data(), num_bins, n);
}

double compute_shannon_entropy_discrete(const std::vector<int>& counts) {
    int total = std::accumulate(counts.begin(), counts.end(), 0);
    return entropy_from_counts(counts.data(), static_cast<int>(counts.size()), total);
}

// ─── compute_torsional_vibrational_entropy (Eigen-vectorised) ────────────────

double compute_torsional_vibrational_entropy(
    const std::vector<tencm::NormalMode>& modes,
    double temperature_K)
{
    if (modes.empty()) return 0.0;
    const double kT = kB_kcal * temperature_K;

    // Threshold floors. Prefilter on eigenvalue λ ≥ kEigenvalueFloor (1e-6) skips
    // rigid-body modes; equivalent frequency floor is ω ≥ sqrt(1e-6) = 1e-3.
    // We keep both layers of protection (prefilter + post-sqrt floor) as
    // defense in depth: zero/negative eigenvalues from external sources or
    // numerical noise must not produce ∞ in kT/ω.
    constexpr double kEigenvalueFloor = 1e-6;
    constexpr double kFrequencyFloor  = 1e-3;  // sqrt(kEigenvalueFloor)

    // Collect valid eigenvalues into Eigen array, then vectorise
    std::vector<double> ev_buf;
    ev_buf.reserve(modes.size());
    for (size_t m = 6; m < modes.size(); ++m)
        if (modes[m].eigenvalue > kEigenvalueFloor) ev_buf.push_back(modes[m].eigenvalue);
    if (ev_buf.empty()) return 0.0;

    Eigen::Map<Eigen::ArrayXd> evals(ev_buf.data(), (int)ev_buf.size());
    // eigenvalue λ = ω²; need frequency ω = sqrt(λ) for the HO entropy formula
    Eigen::ArrayXd freqs  = evals.sqrt();
    // Defense in depth: clamp frequencies to floor so kT/ω is always finite,
    // and mask out below-floor modes from the entropy sum.
    Eigen::ArrayXd valid_mask = (freqs > kFrequencyFloor).cast<double>();
    Eigen::ArrayXd safe_freqs = freqs.max(kFrequencyFloor);
    Eigen::ArrayXd ln_arg = kT / safe_freqs;  // element-wise: kT/ω, always finite
    // S_mode = kB*(1 + ln(kBT/ω)) for valid modes only
    return kB_kcal * (valid_mask * (1.0 + ln_arg.log())).sum();
}

// ─── run_shannon_thermo_stack ────────────────────────────────────────────────

FullThermoResult run_shannon_thermo_stack(
    const statmech::StatMechEngine& stat_engine,
    const tencm::TorsionalENM&      tencm_model,
    double                          base_deltaG,
    double                          temperature_K)
{
    // Conformational Shannon entropy of the Boltzmann distribution:
    //     H = -Σ_i w_i · ln(w_i)        (nats)
    //     S_conf = k_B · H              (kcal/mol·K)
    //
    // Previously this code histogrammed the distribution of −log(w_i) values
    // and called compute_shannon_entropy() on that — which gives the entropy
    // of the binned log-weight distribution, NOT the conformational entropy
    // of the underlying Boltzmann distribution. The two differ by an
    // arbitrary binning factor and have no thermodynamic interpretation.
    auto weights = stat_engine.boltzmann_weights();
    double S_conf_nats = 0.0;
    for (double w : weights)
        if (w > 0.0) S_conf_nats -= w * std::log(w);

    double S_vib        = tencm_model.is_built()
                          ? compute_torsional_vibrational_entropy(tencm_model.modes(), temperature_K)
                          : 0.0;
    // Shannon H is in nats (natural log). Convert to physical units:
    //     S_conf [kcal/(mol·K)] = k_B · H [nats]
    double S_conf_phys  = S_conf_nats * kB_kcal;

    // Additive decomposition: S_total = S_conf + S_vib
    // Valid for independent conformational and vibrational DOFs
    // (standard assumption in rigid-body docking + normal-mode analysis).
    double total_S      = S_conf_phys + S_vib;
    double S_contrib    = -temperature_K * total_S;
    double final_dG     = base_deltaG + S_contrib;

    const char* hw =
#if defined(FLEXAIDS_USE_CUDA)
        "CUDA";
#elif defined(FLEXAIDS_HAS_METAL_SHANNON)
        "Metal";
#elif defined(__AVX512F__)
        "AVX-512";
#elif defined(_OPENMP)
        "OpenMP";
#else
        "scalar";
#endif

    std::string report =
        std::string("ShannonThermoStack[") + hw +
        "+Eigen"
        "]: S_conf=" + std::to_string(S_conf_nats) +
        " nats, S_vib=" + std::to_string(S_vib) +
        " kcal/mol/K, ΔG=" + std::to_string(final_dG) + " kcal/mol";

    return { final_dG, S_conf_nats, S_vib, S_contrib, report };
}

// ─── detect_entropy_plateau ──────────────────────────────────────────────────

bool detect_entropy_plateau(const std::vector<double>& history,
                            int window, double rel_threshold) {
    if (window <= 0 || static_cast<int>(history.size()) < window)
        return false;

    size_t start = history.size() - window;
    double H_ref = history[start];

    for (size_t k = start + 1; k < history.size(); ++k) {
        double rel_change = (H_ref > 1e-12)
            ? std::abs(history[k] - H_ref) / H_ref
            : std::abs(history[k] - H_ref);
        if (rel_change > rel_threshold)
            return false;
    }
    return true;
}

// ─── EntropyEventDetector ───────────────────────────────────────────────────

EntropyEventDetector::EntropyEventDetector(
    int window_size, double collapse_threshold,
    double expansion_threshold, int oscillation_window)
    : window_size_(window_size > 0 ? window_size : 8)
    , collapse_threshold_(collapse_threshold)
    , expansion_threshold_(expansion_threshold > 0 ? expansion_threshold : -collapse_threshold)
    , oscillation_window_(oscillation_window > 0 ? oscillation_window : 5)
    , event_history_(oscillation_window_, EntropyEventType::None) {}

void EntropyEventDetector::reset() {
    history_.clear();
    event_history_.assign(oscillation_window_, EntropyEventType::None);
}

bool EntropyEventDetector::detect_oscillation() const {
    int alternations = 0;
    for (size_t i = 1; i < event_history_.size(); ++i) {
        auto prev = event_history_[i - 1];
        auto curr = event_history_[i];
        if ((prev == EntropyEventType::Collapse && curr == EntropyEventType::Expansion) ||
            (prev == EntropyEventType::Expansion && curr == EntropyEventType::Collapse)) {
            ++alternations;
        }
    }
    return alternations >= 2;
}

EntropyEventResult EntropyEventDetector::push(double entropy) {
    history_.push_back(entropy);

    int count = static_cast<int>(history_.size());
    if (count > window_size_) {
        history_.erase(history_.begin(), history_.begin() + (count - window_size_));
        count = window_size_;
    }

    EntropyEventResult result;
    result.entropy = entropy;

    if (count < 2) {
        event_history_.push_back(EntropyEventType::None);
        if (static_cast<int>(event_history_.size()) > oscillation_window_)
            event_history_.erase(event_history_.begin());
        return result;
    }

    double mean = 0.0;
    for (double v : history_) mean += v;
    mean /= static_cast<double>(count);

    double var = 0.0;
    for (double v : history_) { double d = v - mean; var += d * d; }
    var /= static_cast<double>(count);
    double std_ = std::sqrt(std::max(0.0, var));

    double delta = entropy - mean;
    double z = (std_ > 1e-12) ? delta / std_ : 0.0;
    result.delta = delta;
    result.z_score = z;

    bool window_ready = (count >= window_size_);
    EntropyEventType event = EntropyEventType::None;

    if (window_ready && delta < collapse_threshold_) {
        event = EntropyEventType::Collapse;
    } else if (window_ready && delta > expansion_threshold_) {
        event = EntropyEventType::Expansion;
    }

    event_history_.push_back(event);
    if (static_cast<int>(event_history_.size()) > oscillation_window_)
        event_history_.erase(event_history_.begin());

    if (window_ready && event != EntropyEventType::None && detect_oscillation()) {
        event = EntropyEventType::Oscillation;
    }

    result.event = event;
    return result;
}

} // namespace shannon_thermo
