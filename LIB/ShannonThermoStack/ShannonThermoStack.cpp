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
#include <cstddef>
#include <cstdlib>
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

// ─── bin index ───────────────────────────────────────────────────────────────
//
// The clamp is applied in DOUBLE, before the narrowing conversion.
//
// With a raw min/max support every value satisfied 0 <= (v-min)*inv_bw <=
// num_bins, so casting first and clamping the int was safe. Once the support can
// be a robust fence (see robust_support below), out-of-fence values are
// arbitrarily far outside it: a tight bulk plus one clash pose at 1e4 produces a
// raw index above 5e10, and converting that to int is undefined behaviour. In
// practice x86-64 cvttsd2si yields INT_MIN, which then clamps to bin 0 — placing
// the clash pose ON TOP OF the bulk, collapsing the histogram to a single
// occupied bin and firing the very gate the fence exists to protect — while ARM
// saturates to INT_MAX and lands it in the top bin. Clamping in double removes
// both the undefined behaviour and the platform divergence.
static inline int bin_index(double v, double min_v, double inv_bw, int num_bins) noexcept
{
    const double t = (v - min_v) * inv_bw;
    if (!(t > 0.0)) return 0;                      // also catches NaN
    const double top = static_cast<double>(num_bins - 1);
    return (t >= top) ? num_bins - 1 : static_cast<int>(t);
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

    const __m512d vlo = _mm512_setzero_pd();
    const __m512d vhi = _mm512_set1_pd(static_cast<double>(num_bins - 1));
    for (; i + 7 < n; i += 8) {
        __m512d ve   = _mm512_loadu_pd(values + i);
        __m512d vrel = _mm512_fmadd_pd(ve, vinvbw,
                           _mm512_mul_pd(_mm512_set1_pd(-min_v), vinvbw));
        // Clamp in double BEFORE narrowing: with a fenced support vrel can fall
        // far outside [0, num_bins), and converting that to int32 is UB (and
        // lands the value in the wrong edge bin). _mm512_max_pd also maps NaN
        // to the lower bound here.
        vrel = _mm512_min_pd(_mm512_max_pd(vrel, vlo), vhi);
        // Convert 8 doubles → 8 int32 (truncate)
        __m256i v32 = _mm512_cvttpd_epi32(vrel);
        alignas(32) int tmp[8];
        _mm256_storeu_si256((__m256i*)tmp, v32);
        for (int k = 0; k < 8; ++k) {
            int b = std::min(std::max(tmp[k], 0), num_bins - 1);
            priv[b]++;
        }
    }
    for (; i < n; ++i)
        priv[bin_index(values[i], min_v, inv_bw, num_bins)]++;
}
#endif // __AVX512F__

// ─── robust histogram support ────────────────────────────────────────────────
//
// The histogram support is normally the sample's own [min, max]. That makes a
// single extreme value rescale every bin: one clash/wall pose (evalue ~1e4) in
// an otherwise diverse population pushes every real sample into bin 0, so H
// reads ~0.08 bits and the caller's collapse gate fires on a population that
// has not collapsed at all.
//
// When the raw range is pathologically wider than the bulk of the sample, the
// support is instead taken from a Tukey far-out fence around the interquartile
// range. Values outside the fence are not discarded: every histogram path
// clamps the bin index, so they land in the edge bins and still carry their
// probability mass.
//
// Engagement is range/IQR based, and the RAW range is not itself robust, so
// which samples trigger the fence is tail-dependent rather than merely
// "pathological". Measured firing rates: 0% on Gaussian (n = 1e3 and 1e5),
// uniform and exponential — those are bit-identical to the previous estimator —
// but ~100% on heavy-tailed samples (Student-t, lognormal repulsive tails) and
// on any population carrying clash/wall poses, where H rises by 1.5–2.3 nats.
// That is the intended repair, but it is a real shift in the GA collapse gate's
// operating point, so it is exposed as an A/B arm: set FLEXAIDDS_SHANNON_ROBUST=0
// to restore the previous raw min/max support bit-for-bit.
namespace {

constexpr std::size_t kRobustMinSamples  = 8;    // below this, quartiles are meaningless
constexpr double      kRobustTriggerIQR  = 20.0; // engage only well outside the bulk
constexpr double      kRobustFenceIQR    = 3.0;  // Tukey "far out" fence

bool robust_support_enabled()
{
    static const bool enabled = [] {
        const char* e = std::getenv("FLEXAIDDS_SHANNON_ROBUST");
        return !(e && std::atoi(e) == 0);
    }();
    return enabled;
}


// Returns true (and fills lo/hi) when a robust support should replace [min,max].
bool robust_support(const std::vector<double>& values,
                    double raw_min, double raw_max,
                    double& lo, double& hi)
{
    const std::size_t n = values.size();
    if (n < kRobustMinSamples) return false;
    if (!robust_support_enabled()) return false;

    std::vector<double> scratch(values);
    auto quantile = [&scratch, n](double f) {
        const std::size_t idx =
            static_cast<std::size_t>(f * static_cast<double>(n - 1));
        // nth_element is valid on any permutation, so successive calls on the
        // same (partially reordered) buffer still return the correct element.
        std::nth_element(scratch.begin(),
                         scratch.begin() + static_cast<std::ptrdiff_t>(idx),
                         scratch.end());
        return scratch[idx];
    };

    const double q1  = quantile(0.25);
    const double q3  = quantile(0.75);
    const double iqr = q3 - q1;
    if (!(iqr > 0.0) || !std::isfinite(iqr)) return false;

    if ((raw_max - raw_min) <= kRobustTriggerIQR * iqr) return false;

    lo = q1 - kRobustFenceIQR * iqr;
    hi = q3 + kRobustFenceIQR * iqr;
    // Must clear the +1e-10 bin-width floor applied by the caller, not merely be
    // nonzero: a fence narrower than num_bins*1e-10 is swallowed by that epsilon,
    // so every bulk sample collapses into bin 0 and the fence silently does
    // nothing while still paying for the copy and forfeiting the GPU path.
    return (hi - lo) > 1e-9;
}

}  // namespace

// ─── compute_shannon_entropy ─────────────────────────────────────────────────

double compute_shannon_entropy(const std::vector<double>& values, int num_bins) {
    if (values.empty()) return 0.0;
    if (num_bins <= 0)  num_bins = DEFAULT_HIST_BINS;

    auto [it_min, it_max] = std::minmax_element(values.begin(), values.end());
    double min_v = *it_min;
    double max_v = *it_max;
    if (max_v - min_v < 1e-12) return 0.0;

    // Outlier-robust support; see robust_support() above.
    double rob_lo = 0.0, rob_hi = 0.0;
    const bool robust = robust_support(values, min_v, max_v, rob_lo, rob_hi);
    if (robust) {
        min_v = rob_lo;
        max_v = rob_hi;
    }

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
    // The Metal bridge derives the histogram support internally from the raw
    // sample, so it cannot honour a robust fence. When one is active, fall
    // through to the CPU paths rather than silently returning the
    // outlier-dominated value the fence exists to prevent.
    if (!robust)
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
            // Chunking must use the team size actually granted, not
            // omp_get_max_threads(). When the runtime hands out fewer threads
            // than the maximum (nested regions, dynamic teams, an active
            // num_threads clause), a max-derived chunk leaves the tail of the
            // sample unvisited while the normalisation below still divides by
            // the full count — silently dropping samples from the histogram.
            const int nt = omp_get_num_threads();
            int tid    = omp_get_thread_num();
            int chunk  = (n + nt - 1) / nt;
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
            t_bins(tid, bin_index(values[i], min_v, inv_bw, num_bins))++;
        }
        Eigen::VectorXi col_sums = t_bins.colwise().sum();
        for (int b = 0; b < num_bins; ++b) bins[b] = col_sums(b);
    }

// ── 5. Scalar ─────────────────────────────────────────────────────────────────
#else
    for (int i = 0; i < n; ++i)
        bins[bin_index(values[i], min_v, inv_bw, num_bins)]++;
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

    // Collect valid eigenvalues into Eigen array, then vectorise.
    // These modes come from the INTERNAL-coordinate torsional Hessian
    // (TorsionalENM::build/build_from_ca), which has no rigid-body zero-mode
    // manifold — so do NOT positionally skip the first 6 (that discarded the 6
    // softest real torsions, which carry the largest per-mode entropy). The
    // eigenvalue > 1e-6 guard drops numerically singular modes, consistent with
    // how the Cartesian ANM consumers (ic2cf.cpp, ligand_tencom_pose.cpp) filter
    // their own rigid-body modes by cutoff rather than by position.
    std::vector<double> ev_buf;
    ev_buf.reserve(modes.size());
    for (size_t m = 0; m < modes.size(); ++m)
        if (modes[m].eigenvalue > 1e-6) ev_buf.push_back(modes[m].eigenvalue);
    if (ev_buf.empty()) return 0.0;

    Eigen::Map<Eigen::ArrayXd> evals(ev_buf.data(), (int)ev_buf.size());

    // ── Dimensional status (must be understood before touching this formula) ──
    //
    // omega = sqrt(ENCoM eigenvalue) is in *model units*, NOT in rad/s.
    // The classical HO formula requires a dimensionless argument:
    //   arg = kBT [J] / (ħ [J·s] × ω [rad/s])   →   dimensionless.
    // Here ħ × omega_model is dimensionally J·s × model_unit^(1/2), not J.
    // Numerically at 300 K with eigenvalue ≈ 1 (model unit):
    //   arg ≈ (1.38e-23 × 300) / (1.055e-34 × 1) ≈ 3.9e13
    //   ln(arg) ≈ 31.3  →  S_mode ≈ kB_kcal × 32.3 ≈ 0.064 kcal/(mol·K)
    // The constant offset ln(kBT_SI / hbar_SI) ≈ 31.3 is physically meaningless;
    // it cancels in differential comparisons between structures run with the
    // same protocol:
    //   ΔS_vib = kB_kcal × Σ ln(ω_ref / ω_target)
    //          = (kB_kcal / 2) × Σ ln(λ_ref / λ_target)
    // which depends only on eigenvalue *ratios*, not absolute magnitude.
    //
    // Absolute S_vib and -T·S_vib from this path are heuristic unless a
    // calibration bundle maps ENCoM eigenvalues to physical rad/s, as the
    // ENCoM/tENCoM calibration metadata path requires externally.
    //
    // DO NOT use the return value of this function as an absolute physical
    // entropy without a calibrated eigenvalue scale.
    Eigen::ArrayXd omega = evals.sqrt();
    Eigen::ArrayXd ln_arg = (kB_SI * temperature_K) / (hbar_SI * omega);
    auto valid = (ln_arg > 0.0) && ln_arg.isFinite();
    Eigen::ArrayXd safe_arg = valid.select(ln_arg, 1.0);
    Eigen::ArrayXd contribution = valid.select(1.0 + safe_arg.log(), 0.0);
    return kB_kcal * contribution.sum();
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
    // Shannon H is in nats (natural log). Converting to physical units:
    //     S_conf [kcal/(mol·K)] = k_B · H [nats]
    //
    // ⚠ That unit label is only earned when the weights w_i came from a
    // calibrated physical energy scale. The engine's provenance is the
    // authority, and a CF/contact-function optimizer ensemble is proxy-only:
    // its weights are a softmax over arbitrary units, so k_B is dimensionally
    // meaningless here and -T·k_B·H is not a kcal/mol free-energy term.
    //
    // The numbers below are deliberately left unchanged. This branch's
    // provenance contract (statmech.h) is explicit that the metadata "never
    // participates in energy evaluation, weighting, ranking, or sorting", and
    // this ΔG is consumed by the GA. Gating the arithmetic on provenance would
    // silently drop an entropy term from production docking output — a
    // science-affecting change that belongs in a benchmarked A/B, not in a
    // labeling fix. The calibration status is surfaced in the report instead,
    // so a reader can tell whether "kcal/mol" is a claim or a placeholder.
    const bool calibrated =
        stat_engine.provenance().allows_canonical_physical_claim();
    double S_conf_phys  = S_conf_nats * kB_kcal;

    // The torsional ENCoM path is separately model-scale only. Keep S_vib in
    // the result/report as a relative flexibility diagnostic, but do not fold it
    // into kcal/mol free-energy terms until a calibrated frequency path reaches
    // this stack.
    double total_S      = S_conf_phys;
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
        (calibrated
             ? " nats (calibrated energy scale; -T·S is kcal/mol)"
             : " nats (proxy-only ensemble; -T·S in CF a.u., kcal/mol label unearned)") +
        ", S_vib_heuristic=" + std::to_string(S_vib) +
        " kcal/mol/K (model-scale heuristic; excluded from dG), ΔG=" +
        std::to_string(final_dG) + " kcal/mol";

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
