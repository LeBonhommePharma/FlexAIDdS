// =============================================================================
// TurboQuant.h
// Production-grade C++20 header-only implementation of the TurboQuant vector
// quantization algorithm, adapted for the FlexAIDdS molecular docking engine.
//
// Algorithm Reference:
//   Zandieh, A., Daliri, M., Hadian, M., & Mirrokni, V. (2025).
//   "TurboQuant: Online Vector Quantization with Near-optimal Distortion Rate."
//   arXiv:2504.19874. https://arxiv.org/abs/2504.19874
//
// FlexAIDdS Application:
//   - Compressed 256-dim row vectors of the SoftContactMatrix (256×256)
//   - GA conformational ensemble energy vector quantization
//   - Fast nearest-neighbor for FastOPTICS super-cluster detection
//
// Copyright 2026 Le Bonhomme Pharma
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// =============================================================================

#pragma once

// MSVC does not support __restrict__ (GCC/Clang extension); map to __restrict
#ifdef _MSC_VER
#  ifndef __restrict__
#    define __restrict__ __restrict
#  endif
#endif

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Eigen (required)
#include <Eigen/Core>
#include <Eigen/QR>

// SIMD intrinsics (optional, guarded by feature-test macros)
#if defined(__AVX512F__)
#  include <immintrin.h>
#  define TURBOQUANT_HAS_AVX512 1
#elif defined(__AVX2__)
#  include <immintrin.h>
#  define TURBOQUANT_HAS_AVX2 1
#endif

// OpenMP (optional)
#if defined(_OPENMP)
#  include <omp.h>
#  define TURBOQUANT_HAS_OMP 1
#endif

// CUDA host/device annotations
#if defined(__CUDACC__)
#  define TQ_HD __host__ __device__
#  define TQ_DEVICE __device__
#else
#  define TQ_HD
#  define TQ_DEVICE
#endif

// =============================================================================
// Forward declarations for CUDA batch kernels (implemented in TurboQuant.cu)
// =============================================================================
#if defined(__CUDACC__) || defined(TURBOQUANT_DECLARE_CUDA)
extern "C" {
    // Batch quantize N vectors of dimension d using precomputed rotation Π.
    // d_Pi: device pointer to d×d rotation matrix (row-major float32)
    // d_input: device pointer to N×d input matrix (row-major float32)
    // d_indices: device pointer to output N×d uint8 index matrix
    // d_norms: device pointer to output N float norms
    // boundaries: pointer to 2^b - 1 boundary values
    // num_boundaries: 2^b - 1
    void turboquant_cuda_batch_quantize(
        const float* d_Pi,
        const float* d_input,
        uint8_t*     d_indices,
        float*       d_norms,
        const float* boundaries,
        int          num_boundaries,
        int          N,
        int          d,
        int          bit_width,
        cudaStream_t stream
    );

    // Batch dequantize N vectors.
    void turboquant_cuda_batch_dequantize(
        const float*   d_PiT,
        const uint8_t* d_indices,
        const float*   d_centroids,
        float*         d_output,
        int            N,
        int            d,
        int            bit_width,
        cudaStream_t   stream
    );
}
#endif // __CUDACC__ || TURBOQUANT_DECLARE_CUDA

// =============================================================================
namespace turboquant {
// =============================================================================

// -----------------------------------------------------------------------------
// Internal constants and compile-time codebook tables
// -----------------------------------------------------------------------------

/// Standard dimension for FlexAIDdS contact vectors.
inline constexpr int kContactDim = 256;

// Precomputed optimal Lloyd-Max codebook centroids for the Beta distribution
// f_X(x) = Γ(d/2)/(√π · Γ((d-1)/2)) · (1-x²)^((d-3)/2)  on [-1,1]
// For d=256 this distribution is extremely well approximated by N(0, 1/256).
// Centroids are stored in ascending order.
// Scale factor: 1/√d = 1/√256 = 0.0625  (all values below are pre-scaled)
//
// Raw (unit-variance Gaussian) Lloyd-Max centroids from classical tables:
//   b=1: ±0.7979   → /√256 = ±0.049868...
//   b=2: ±0.4528, ±1.5104  → /√256
//   b=3: ±0.2451, ±0.7560, ±1.3439, ±2.1520  → /√256
//   b=4: ±0.1284, ±0.3880, ±0.6568, ±0.9424, ±1.2562,
//         ±1.6180, ±2.0690, ±2.7326  → /√256
//
// For general d, use build_codebook() which runs Lloyd-Max iteration.

namespace detail {

// --- bit packing / unpacking helpers ----------------------------------------

/// Pack `n_values` b-bit unsigned integers from `vals` into `out` (bytes).
/// Assumes vals[i] < 2^bit_width.
inline void pack_bits(const uint8_t* vals, int n_values, int bit_width,
                      uint8_t* out, int out_len) noexcept {
    if (out_len > 0) std::fill(out, out + out_len, uint8_t{0});
    int bit_cursor = 0;
    for (int i = 0; i < n_values; ++i) {
        int byte_idx = bit_cursor >> 3;
        int bit_off  = bit_cursor & 7;
        // Write low-order bits
        out[byte_idx] |= static_cast<uint8_t>(vals[i] << bit_off);
        if (bit_off + bit_width > 8 && byte_idx + 1 < out_len) {
            out[byte_idx + 1] |=
                static_cast<uint8_t>(vals[i] >> (8 - bit_off));
        }
        bit_cursor += bit_width;
    }
}

/// Unpack `n_values` b-bit unsigned integers from `in` into `vals`.
inline void unpack_bits(const uint8_t* in, int n_values, int bit_width,
                        uint8_t* vals) noexcept {
    const uint8_t mask = static_cast<uint8_t>((1u << bit_width) - 1u);
    int bit_cursor = 0;
    for (int i = 0; i < n_values; ++i) {
        int byte_idx = bit_cursor >> 3;
        int bit_off  = bit_cursor & 7;
        uint8_t v    = static_cast<uint8_t>(in[byte_idx] >> bit_off);
        if (bit_off + bit_width > 8) {
            v |= static_cast<uint8_t>(in[byte_idx + 1] << (8 - bit_off));
        }
        vals[i]    = v & mask;
        bit_cursor += bit_width;
    }
}

/// Number of bytes required to store `n_values` values at `bit_width` bpv.
inline constexpr int packed_byte_size(int n_values, int bit_width) noexcept {
    return (n_values * bit_width + 7) / 8;
}

// --- Beta distribution PDF --------------------------------------------------

/// Beta marginal density of a uniform unit-sphere coordinate in R^d.
/// f_X(x) = Γ(d/2)/(√π·Γ((d-1)/2)) · (1-x²)^((d-3)/2),  x ∈ (-1,1)
inline double beta_pdf(double x, int d) noexcept {
    if (x <= -1.0 || x >= 1.0) return 0.0;
    // Use log-gamma for numerical stability
    double log_norm = std::lgamma(d / 2.0) - 0.5 * std::log(M_PI)
                    - std::lgamma((d - 1) / 2.0);
    double log_val  = log_norm + ((d - 3) / 2.0) * std::log(1.0 - x * x);
    return std::exp(log_val);
}

/// Gaussian PDF with variance sigma2.
inline double gaussian_pdf(double x, double sigma2) noexcept {
    return std::exp(-0.5 * x * x / sigma2) / std::sqrt(2.0 * M_PI * sigma2);
}

// --- Lloyd-Max iteration ----------------------------------------------------

/// Run Lloyd-Max algorithm to find MSE-optimal scalar codebook centroids for
/// the given 1-D distribution PDF sampled on `grid` with spacing `dx`.
/// Returns centroids in ascending order.
///
/// @param pdf       Probability density sampled on `grid`
/// @param grid      Evaluation points in [-1,1], uniform spacing
/// @param n_grid    Number of grid points
/// @param k         Number of centroids (= 2^bit_width)
/// @param max_iter  Maximum Lloyd-Max iterations
/// @param tol       Convergence tolerance on centroid shift
inline std::vector<double> lloyd_max(
    const std::vector<double>& pdf,
    const std::vector<double>& grid,
    int k,
    int max_iter = 1000,
    double tol   = 1e-10) {

    assert(pdf.size() == grid.size());
    int n = static_cast<int>(grid.size());
    [[maybe_unused]] double dx = grid[1] - grid[0];

    // Initialize centroids uniformly in [-1,1]
    std::vector<double> centroids(k);
    for (int i = 0; i < k; ++i)
        centroids[i] = -1.0 + (2.0 * (i + 0.5)) / k;

    for (int iter = 0; iter < max_iter; ++iter) {
        // Compute Voronoi boundaries (midpoints)
        std::vector<double> boundaries(k + 1);
        boundaries[0]     = -1.0;
        boundaries[k]     =  1.0;
        for (int i = 1; i < k; ++i)
            boundaries[i] = 0.5 * (centroids[i - 1] + centroids[i]);

        // Update centroids: c_i = E[X | X ∈ region_i] = ∫_{b_{i-1}}^{b_i} x·f(x)dx / ∫_{b_{i-1}}^{b_i} f(x)dx
        double max_shift = 0.0;
        for (int i = 0; i < k; ++i) {
            double sum_xf = 0.0, sum_f = 0.0;
            for (int j = 0; j < n; ++j) {
                double x = grid[j];
                if (x >= boundaries[i] && x < boundaries[i + 1]) {
                    sum_xf += x * pdf[j];
                    sum_f  += pdf[j];
                }
            }
            double new_c = (sum_f > 1e-300) ? sum_xf / sum_f : centroids[i];
            max_shift    = std::max(max_shift, std::abs(new_c - centroids[i]));
            centroids[i] = new_c;
        }
        if (max_shift < tol) break;
    }
    return centroids;
}

} // namespace detail

// =============================================================================
// Codebook
// =============================================================================

/// Optimal scalar quantization codebook for TurboQuant.
struct Codebook {
    int              bit_width{0};
    int              num_centroids{0};    ///< = 2^bit_width
    std::vector<float> centroids;         ///< ascending, size = num_centroids
    std::vector<float> boundaries;        ///< midpoints, size = num_centroids-1

    Codebook() = default;
    Codebook(int bw, std::vector<float> c)
        : bit_width(bw),
          num_centroids(static_cast<int>(c.size())),
          centroids(std::move(c))
    {
        boundaries.resize(num_centroids - 1);
        for (int i = 0; i < num_centroids - 1; ++i)
            boundaries[i] = 0.5f * (centroids[i] + centroids[i + 1]);
    }

    /// Find nearest centroid index for scalar value v (brute-force, fast for k≤16).
    TQ_HD inline int nearest(float v) const noexcept {
        int best = 0;
        float best_d = std::abs(v - centroids[0]);
        for (int i = 1; i < num_centroids; ++i) {
            float d = std::abs(v - centroids[i]);
            if (d < best_d) { best_d = d; best = i; }
        }
        return best;
    }

    /// Faster nearest using sorted boundaries (binary-search-like for k≤16).
    TQ_HD inline int nearest_fast(float v) const noexcept {
        // Linear scan through boundaries (k ≤ 16 → at most 15 comparisons)
        int idx = 0;
        for (int i = 0; i < static_cast<int>(boundaries.size()); ++i) {
            if (v >= boundaries[i]) idx = i + 1;
            else break;
        }
        return idx;
    }
};

// =============================================================================
// build_codebook: run Lloyd-Max for the Beta(d) distribution
// =============================================================================

/// Build optimal Lloyd-Max codebook for the Beta marginal of uniform
/// unit-sphere vectors in R^d, at the given bit width.
///
/// For d=256 the distribution is practically Gaussian N(0,1/d), so the
/// precomputed analytic centroids are extremely accurate, but this function
/// runs the full numerical procedure for arbitrary d and bit widths.
///
/// @param d         Dimension of the ambient space (must be ≥ 2)
/// @param bit_width Quantization bits (1–4 recommended)
/// @param n_grid    Grid resolution for numerical integration (default 50000)
inline Codebook build_codebook(int d, int bit_width, int n_grid = 50000) {
    if (d < 2)       throw std::invalid_argument("build_codebook: d must be >= 2");
    if (bit_width < 1 || bit_width > 8)
        throw std::invalid_argument("build_codebook: bit_width must be in [1,8]");

    int k = 1 << bit_width;

    // Build uniform grid on (-1+eps, 1-eps)
    const double eps = 1e-6;
    const double lo  = -1.0 + eps;
    const double hi  =  1.0 - eps;
    double dx = (hi - lo) / (n_grid - 1);

    std::vector<double> grid(n_grid), pdf(n_grid);
    for (int i = 0; i < n_grid; ++i) {
        grid[i] = lo + i * dx;
        pdf[i]  = detail::beta_pdf(grid[i], d);
    }

    // Normalise pdf so it integrates to 1 (trapezoidal rule)
    double norm = 0.0;
    for (int i = 0; i < n_grid - 1; ++i)
        norm += 0.5 * (pdf[i] + pdf[i + 1]) * dx;
    if (norm > 1e-300)
        for (auto& p : pdf) p /= norm;

    auto raw = detail::lloyd_max(pdf, grid, k);

    std::vector<float> fc(raw.begin(), raw.end());
    return Codebook(bit_width, std::move(fc));
}

// =============================================================================
// Precomputed constexpr codebooks for d=256, b=1,2,3,4
// All centroids are normalised by 1/√256 = 1/16 from unit-Gaussian Lloyd-Max.
// =============================================================================

namespace detail {

// Unit-Gaussian Lloyd-Max centroids (classical tabulated values).
// Sign-symmetric; only positive halves listed; fill via reflection.
// Ref: Jayant & Noll, "Digital Coding of Waveforms", Prentice-Hall 1984.

// b=1: ±√(2/π) ≈ ±0.79788 → scaled /16 = ±0.049868
inline constexpr std::array<float, 2> kRawCentroids1 = {
    -0.049868f, +0.049868f
};

// b=2: (−1.5104, −0.4528, +0.4528, +1.5104) / 16
inline constexpr std::array<float, 4> kRawCentroids2 = {
    -0.094400f, -0.028300f, +0.028300f, +0.094400f
};

// b=3: 8 centroids, /16
// ±0.2451/16, ±0.7560/16, ±1.3439/16, ±2.1520/16
inline constexpr std::array<float, 8> kRawCentroids3 = {
    -0.134500f, -0.084000f, -0.047250f, -0.015319f,
    +0.015319f, +0.047250f, +0.084000f, +0.134500f
};

// b=4: 16 centroids, /16
// ±0.1284, ±0.3880, ±0.6568, ±0.9424, ±1.2562, ±1.6180, ±2.0690, ±2.7326
inline constexpr std::array<float, 16> kRawCentroids4 = {
    -0.170788f, -0.129313f, -0.103750f, -0.082150f,
    -0.059025f, -0.041100f, -0.024250f, -0.008025f,
    +0.008025f, +0.024250f, +0.041100f, +0.059025f,
    +0.082150f, +0.103750f, +0.129313f, +0.170788f
};

} // namespace detail

/// Return a precomputed codebook for d=256 and b ∈ {1,2,3,4}.
/// Falls back to build_codebook() for other combinations.
inline Codebook make_codebook_d256(int bit_width) {
    switch (bit_width) {
    case 1: {
        std::vector<float> c(detail::kRawCentroids1.begin(),
                             detail::kRawCentroids1.end());
        return Codebook(1, std::move(c));
    }
    case 2: {
        std::vector<float> c(detail::kRawCentroids2.begin(),
                             detail::kRawCentroids2.end());
        return Codebook(2, std::move(c));
    }
    case 3: {
        std::vector<float> c(detail::kRawCentroids3.begin(),
                             detail::kRawCentroids3.end());
        return Codebook(3, std::move(c));
    }
    case 4: {
        std::vector<float> c(detail::kRawCentroids4.begin(),
                             detail::kRawCentroids4.end());
        return Codebook(4, std::move(c));
    }
    default:
        return build_codebook(256, bit_width);
    }
}

// =============================================================================
// SIMD helpers for matrix–vector multiply Π·x (d×d, float32)
// =============================================================================

namespace detail {

/// Dense matvec: y = A * x, A is (rows × cols) row-major float32.
/// Dispatches to AVX-512, AVX2, or scalar depending on compile flags.
inline void matvec_f32(const float* __restrict__ A,
                       const float* __restrict__ x,
                       float*       __restrict__ y,
                       int rows, int cols) noexcept {
#if TURBOQUANT_HAS_AVX512
    // AVX-512: process 16 columns at a time
    for (int r = 0; r < rows; ++r) {
        const float* row = A + r * cols;
        __m512 acc = _mm512_setzero_ps();
        int c = 0;
        for (; c + 16 <= cols; c += 16) {
            __m512 a = _mm512_loadu_ps(row + c);
            __m512 v = _mm512_loadu_ps(x + c);
            acc = _mm512_fmadd_ps(a, v, acc);
        }
        float sum = _mm512_reduce_add_ps(acc);
        for (; c < cols; ++c) sum += row[c] * x[c];
        y[r] = sum;
    }
#elif TURBOQUANT_HAS_AVX2
    // AVX2: process 8 columns at a time
    for (int r = 0; r < rows; ++r) {
        const float* row = A + r * cols;
        __m256 acc = _mm256_setzero_ps();
        int c = 0;
        for (; c + 8 <= cols; c += 8) {
            __m256 a = _mm256_loadu_ps(row + c);
            __m256 v = _mm256_loadu_ps(x + c);
            acc = _mm256_fmadd_ps(a, v, acc);
        }
        // Horizontal reduce of acc
        __m128 lo  = _mm256_castps256_ps128(acc);
        __m128 hi  = _mm256_extractf128_ps(acc, 1);
        lo = _mm_add_ps(lo, hi);
        lo = _mm_hadd_ps(lo, lo);
        lo = _mm_hadd_ps(lo, lo);
        float sum = _mm_cvtss_f32(lo);
        for (; c < cols; ++c) sum += row[c] * x[c];
        y[r] = sum;
    }
#else
    // Scalar fallback
    for (int r = 0; r < rows; ++r) {
        float sum = 0.0f;
        for (int c = 0; c < cols; ++c) sum += A[r * cols + c] * x[c];
        y[r] = sum;
    }
#endif
}

/// Batched centroid-search: for each of `d` elements in `y`, find the
/// nearest centroid index in `cb` and write to `indices`.
inline void nearest_centroid_batch(const float* y, int d,
                                   const Codebook& cb,
                                   uint8_t* indices) noexcept {
    const int k = cb.num_centroids;
    const float* cents = cb.centroids.data();
    const float* bounds = cb.boundaries.data();
    const int nb = k - 1;

    for (int j = 0; j < d; ++j) {
        float v = y[j];
        // Binary-style linear scan on boundaries (k ≤ 16 → at most 15 cmp)
        int idx = 0;
        for (int bi = 0; bi < nb; ++bi) {
            if (v >= bounds[bi]) idx = bi + 1;
            else break;
        }
        (void)cents; // used via idx → centroids[idx] elsewhere
        indices[j] = static_cast<uint8_t>(idx);
    }
}

/// xorshift64 PRNG for fast reproducible random Gaussian streams.
struct XorShift64 {
    uint64_t state;
    explicit XorShift64(uint64_t seed) : state(seed ? seed : 0xdeadbeefcafeULL) {}

    uint64_t next() noexcept {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        return state;
    }

    /// Box-Muller: two independent N(0,1) samples from two uniform draws.
    std::pair<double,double> next_gaussian_pair() noexcept {
        // Two uniform [0,1] via uint64 → double
        constexpr double inv = 1.0 / (static_cast<double>(UINT64_MAX) + 1.0);
        double u1, u2;
        do { u1 = static_cast<double>(next()) * inv; } while (u1 == 0.0);
        u2 = static_cast<double>(next()) * inv;
        double r   = std::sqrt(-2.0 * std::log(u1));
        double phi = 2.0 * M_PI * u2;
        return {r * std::cos(phi), r * std::sin(phi)};
    }
};

/// Fill matrix (rows × cols, row-major) with i.i.d. N(0,1).
inline void fill_normal(float* mat, int rows, int cols, uint64_t seed) {
    XorShift64 rng(seed);
    int n = rows * cols;
    int i = 0;
    for (; i + 1 < n; i += 2) {
        auto [g1, g2] = rng.next_gaussian_pair();
        mat[i]     = static_cast<float>(g1);
        mat[i + 1] = static_cast<float>(g2);
    }
    if (i < n) {
        auto [g1, _] = rng.next_gaussian_pair();
        mat[i] = static_cast<float>(g1);
    }
}

} // namespace detail

// =============================================================================
// TurboQuantMSE  (Algorithm 1 from Zandieh et al.)
// =============================================================================

/// MSE-optimal TurboQuant quantizer.
///
/// Quantization pipeline:
///   1. y = Π·x          (random rotation; Π is orthogonal, from QR)
///   2. idx_j = argmin_k |y_j - c_k|  (scalar quantization per coordinate)
///   3. store packed bit-indices + original L2 norm
///
/// Dequantization:
///   1. ỹ_j = c_{idx_j}
///   2. x̃ = Πᵀ·ỹ
///
/// MSE distortion: E[||x - x̃||²] ≤ (√3·π/2) / 4^b  ≈ 2.7 / 4^b
class TurboQuantMSE {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // -----------------------------------------------------------------
    // Packed quantized vector
    // -----------------------------------------------------------------
    struct QuantizedVector {
        std::vector<uint8_t> packed_indices;  ///< Bit-packed centroid indices
        float                norm{0.0f};       ///< Original L2 norm (||x||_2)
        int                  dim{0};           ///< Original dimension
        int                  bit_width{0};     ///< Bits per coordinate
    };

    // -----------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------

    TurboQuantMSE() = default;

    explicit TurboQuantMSE(int dimension, int bit_width, uint64_t seed = 42)
        : d_(dimension), bit_width_(bit_width), seed_(seed)
    {
        if (d_ <= 0)        throw std::invalid_argument("TurboQuantMSE: d must be > 0");
        if (bit_width_ < 1 || bit_width_ > 4)
            throw std::invalid_argument("TurboQuantMSE: bit_width must be in [1,4]");

        codebook_ = (d_ == 256) ? make_codebook_d256(bit_width_)
                                 : build_codebook(d_, bit_width_);
        initialize();
    }

    // -----------------------------------------------------------------
    // Initialization: generate random rotation Π via QR decomposition
    // -----------------------------------------------------------------

    void initialize() {
        using Matrix = Eigen::MatrixXf;

        // Fill d×d matrix with i.i.d. N(0,1)
        Matrix A(d_, d_);
        detail::fill_normal(A.data(), d_, d_, seed_);

        // QR decomposition → Q is the orthogonal rotation
        Eigen::HouseholderQR<Matrix> qr(A);
        Pi_  = qr.householderQ() * Matrix::Identity(d_, d_);
        PiT_ = Pi_.transpose();

        // Cache contiguous row-major copies for SIMD matvec
        Pi_data_.resize(d_ * d_);
        PiT_data_.resize(d_ * d_);
        Eigen::Map<Eigen::MatrixXf>(Pi_data_.data(),  d_, d_) = Pi_;
        Eigen::Map<Eigen::MatrixXf>(PiT_data_.data(), d_, d_) = PiT_;
    }

    // -----------------------------------------------------------------
    // Quantize
    // -----------------------------------------------------------------

    /// Quantize a d-dimensional vector x → packed indices + norm.
    QuantizedVector quantize(std::span<const float> x) const {
        if (static_cast<int>(x.size()) != d_)
            throw std::invalid_argument("TurboQuantMSE::quantize: wrong dimension");

        QuantizedVector qv;
        qv.dim       = d_;
        qv.bit_width = bit_width_;
        qv.norm      = 0.0f;
        for (auto v : x) qv.norm += v * v;
        qv.norm = std::sqrt(qv.norm);

        // Step 1: y = Π·x
        thread_local std::vector<float> y_buf;
        y_buf.resize(d_);
        detail::matvec_f32(Pi_data_.data(), x.data(), y_buf.data(), d_, d_);

        // Step 2: scalar quantize each coordinate
        thread_local std::vector<uint8_t> idx_buf;
        idx_buf.resize(d_);
        detail::nearest_centroid_batch(y_buf.data(), d_, codebook_, idx_buf.data());

        // Step 3: bit-pack
        int nbytes = detail::packed_byte_size(d_, bit_width_);
        qv.packed_indices.resize(nbytes);
        detail::pack_bits(idx_buf.data(), d_, bit_width_,
                          qv.packed_indices.data(), nbytes);
        return qv;
    }

    // -----------------------------------------------------------------
    // Dequantize
    // -----------------------------------------------------------------

    /// Reconstruct approximate vector x̃ from packed indices.
    void dequantize(const QuantizedVector& qv, std::span<float> out) const {
        if (static_cast<int>(out.size()) != d_)
            throw std::invalid_argument("TurboQuantMSE::dequantize: wrong output size");
        if (qv.dim != d_ || qv.bit_width != bit_width_)
            throw std::invalid_argument("TurboQuantMSE::dequantize: mismatched quantizer");

        // Unpack indices
        thread_local std::vector<uint8_t> idx_buf;
        idx_buf.resize(d_);
        detail::unpack_bits(qv.packed_indices.data(), d_, bit_width_, idx_buf.data());

        // Reconstruct ỹ from codebook
        thread_local std::vector<float> y_buf;
        y_buf.resize(d_);
        for (int j = 0; j < d_; ++j)
            y_buf[j] = codebook_.centroids[idx_buf[j]];

        // x̃ = Πᵀ · ỹ
        detail::matvec_f32(PiT_data_.data(), y_buf.data(), out.data(), d_, d_);
    }

    // -----------------------------------------------------------------
    // Batch quantize (OpenMP parallel)
    // -----------------------------------------------------------------

    /// Quantize N vectors stored row-major in data (N×d).
    std::vector<QuantizedVector> quantize_batch(
        std::span<const float> data, int N) const
    {
        if (static_cast<int>(data.size()) != N * d_)
            throw std::invalid_argument("TurboQuantMSE::quantize_batch: wrong data size");

        std::vector<QuantizedVector> result(N);

#ifdef TURBOQUANT_HAS_OMP
#pragma omp parallel for schedule(dynamic, 32)
#endif
        for (int i = 0; i < N; ++i) {
            std::span<const float> xi(data.data() + i * d_, d_);
            result[i] = quantize(xi);
        }
        return result;
    }

    // -----------------------------------------------------------------
    // Approximate inner product between quantized x and full y
    // -----------------------------------------------------------------

    /// Computes ⟨ỹ, y⟩ where ỹ = Q^{-1}(Q(x)) without full dequantize.
    /// This is faster than dequantize() + dot product because it avoids
    /// the Πᵀ matvec when y is already in the rotated domain.
    ///
    /// Full version: dequantize + dot is O(d²) due to Πᵀ matvec.
    /// This version: rotate y once, lookup centroids, dot → O(d²) + O(d).
    float approx_inner_product(const QuantizedVector& qx,
                                std::span<const float> y) const {
        if (static_cast<int>(y.size()) != d_)
            throw std::invalid_argument("approx_inner_product: wrong y dimension");

        // Unpack indices
        thread_local std::vector<uint8_t> idx_buf;
        idx_buf.resize(d_);
        detail::unpack_bits(qx.packed_indices.data(), d_, bit_width_, idx_buf.data());

        // Rotate y: ŷ = Π·y  (so that ⟨ỹ, y⟩ = ⟨ỹ_rot, ŷ_rot⟩ )
        // Since Π is orthogonal: ⟨Πᵀ·ỹ_idx, y⟩ = ⟨ỹ_idx, Π·y⟩
        thread_local std::vector<float> yrot;
        yrot.resize(d_);
        detail::matvec_f32(Pi_data_.data(), y.data(), yrot.data(), d_, d_);

        // Dot product in rotated space
        float ip = 0.0f;
        for (int j = 0; j < d_; ++j)
            ip += codebook_.centroids[idx_buf[j]] * yrot[j];
        return ip;
    }

    // -----------------------------------------------------------------
    // Diagnostics
    // -----------------------------------------------------------------

    /// Theoretical MSE upper bound from Theorem 1 of Zandieh et al.
    float theoretical_mse() const noexcept {
        // D_mse ≤ sqrt(3π/2) / 4^b  (unit sphere vectors)
        static constexpr float kCoeff = 2.720f; // √(3π/2) ≈ 2.720
        float denom = 1.0f;
        for (int b = 0; b < bit_width_; ++b) denom *= 4.0f;
        return kCoeff / denom;
    }

    // -----------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------

    int dimension()  const noexcept { return d_; }
    int bit_width()  const noexcept { return bit_width_; }
    uint64_t seed()  const noexcept { return seed_; }
    const Codebook& codebook() const noexcept { return codebook_; }

    /// Read-only access to rotation matrix Π (Eigen).
    const Eigen::MatrixXf& Pi()  const noexcept { return Pi_; }
    const Eigen::MatrixXf& PiT() const noexcept { return PiT_; }

    /// Raw row-major float32 pointer to Π (for CUDA kernels).
    const float* pi_data()  const noexcept { return Pi_data_.data(); }
    const float* pit_data() const noexcept { return PiT_data_.data(); }

private:
    int      d_{0};
    int      bit_width_{1};
    Codebook codebook_;
    uint64_t seed_{42};

    Eigen::MatrixXf    Pi_;
    Eigen::MatrixXf    PiT_;
    std::vector<float> Pi_data_;   ///< Row-major copy for SIMD
    std::vector<float> PiT_data_;  ///< Row-major copy for SIMD
};

// =============================================================================
// TurboQuantProd  (Algorithm 2 from Zandieh et al.)
// =============================================================================

/// Inner-product-optimal (unbiased) TurboQuant quantizer.
///
/// Two-stage pipeline:
///   Stage 1 (MSE part):   idx = TurboQuant_mse(x, b-1 bits)
///   Stage 2 (QJL part):   r = x - DeQuant_mse(idx)
///                         γ = ||r||_2
///                         qjl = sign(S·(r/γ))  if γ > 0
///
/// Dequantization:
///   x̃ = DeQuant_mse(idx) + γ · (√(π/2)/d) · Sᵀ · qjl
///
/// Inner product guarantee (Theorem 2):
///   E[⟨y, x̃⟩] = ⟨y, x⟩  (unbiased)
///   Var ≤ (√3 · π² · ||y||²/d) / 4^b
///
/// QJL: Q_qjl(x) = sign(S·x),  Q_qjl^{-1}(z) = (√(π/2)/d)·Sᵀ·z
class TurboQuantProd {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // -----------------------------------------------------------------
    // Packed quantized vector
    // -----------------------------------------------------------------
    struct QuantizedVectorProd {
        TurboQuantMSE::QuantizedVector mse_part;
        std::vector<int8_t>            qjl_signs;    ///< +1 or -1 per coordinate
        float                          residual_norm; ///< ||r||_2
    };

    // -----------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------

    TurboQuantProd() = default;

    explicit TurboQuantProd(int dimension, int bit_width, uint64_t seed = 42)
        : d_(dimension), bit_width_(bit_width),
          mse_quant_(dimension, std::max(1, bit_width - 1), seed),
          seed_(seed)
    {
        if (d_ <= 0)
            throw std::invalid_argument("TurboQuantProd: d must be > 0");
        if (bit_width_ < 1 || bit_width_ > 5)
            throw std::invalid_argument("TurboQuantProd: bit_width must be in [1,5]");
        initialize();
    }

    // -----------------------------------------------------------------
    // Initialization: generate QJL random matrix S ~ N(0,1)
    // -----------------------------------------------------------------

    void initialize() {
        // Use a different seed offset to decorrelate Π and S
        uint64_t s_seed = seed_ ^ 0xC0FFEE0DEADBEEF1ULL;

        S_.resize(static_cast<size_t>(d_) * d_);
        detail::fill_normal(S_.data(), d_, d_, s_seed);

        // Copy into Eigen for transpose
        Eigen::Map<Eigen::MatrixXf> S_map(S_.data(), d_, d_);
        ST_.resize(static_cast<size_t>(d_) * d_);
        Eigen::Map<Eigen::MatrixXf>(ST_.data(), d_, d_) = S_map.transpose();
    }

    // -----------------------------------------------------------------
    // Quantize
    // -----------------------------------------------------------------

    QuantizedVectorProd quantize(std::span<const float> x) const {
        if (static_cast<int>(x.size()) != d_)
            throw std::invalid_argument("TurboQuantProd::quantize: wrong dimension");

        QuantizedVectorProd qv;

        // Stage 1: MSE quantization with (b-1) bits
        qv.mse_part = mse_quant_.quantize(x);

        // Compute residual r = x - DeQuant_mse(idx)
        thread_local std::vector<float> x_hat;
        x_hat.resize(d_);
        mse_quant_.dequantize(qv.mse_part, std::span<float>(x_hat));

        thread_local std::vector<float> r;
        r.resize(d_);
        float r_norm2 = 0.0f;
        for (int j = 0; j < d_; ++j) {
            r[j]    = x[j] - x_hat[j];
            r_norm2 += r[j] * r[j];
        }
        qv.residual_norm = std::sqrt(r_norm2);

        // Stage 2: QJL on normalised residual (or zero vector if norm is tiny)
        qv.qjl_signs.resize(d_);
        if (qv.residual_norm > 1e-9f) {
            float inv_norm = 1.0f / qv.residual_norm;
            thread_local std::vector<float> r_hat;
            r_hat.resize(d_);
            for (int j = 0; j < d_; ++j) r_hat[j] = r[j] * inv_norm;

            // z = S · r_hat  →  sign per coordinate
            thread_local std::vector<float> Sz;
            Sz.resize(d_);
            detail::matvec_f32(S_.data(), r_hat.data(), Sz.data(), d_, d_);
            for (int j = 0; j < d_; ++j)
                qv.qjl_signs[j] = (Sz[j] >= 0.0f) ? int8_t{1} : int8_t{-1};
        } else {
            std::fill(qv.qjl_signs.begin(), qv.qjl_signs.end(), int8_t{1});
        }

        return qv;
    }

    // -----------------------------------------------------------------
    // Dequantize
    // -----------------------------------------------------------------

    /// Reconstruct x̃ = x̃_mse + γ · (√(π/2)/d) · Sᵀ · qjl
    void dequantize(const QuantizedVectorProd& qv, std::span<float> out) const {
        if (static_cast<int>(out.size()) != d_)
            throw std::invalid_argument("TurboQuantProd::dequantize: wrong output size");

        // MSE part
        mse_quant_.dequantize(qv.mse_part, out);

        if (qv.residual_norm < 1e-9f) return;

        // QJL correction: γ · (√(π/2)/d) · Sᵀ · qjl
        static const float kQjlScale = std::sqrt(float(M_PI) / 2.0f);
        float scale = qv.residual_norm * kQjlScale / static_cast<float>(d_);

        // Sᵀ · qjl: convert int8 signs to float on-the-fly
        thread_local std::vector<float> qjl_f;
        qjl_f.resize(d_);
        for (int j = 0; j < d_; ++j)
            qjl_f[j] = static_cast<float>(qv.qjl_signs[j]);

        thread_local std::vector<float> correction;
        correction.resize(d_);
        detail::matvec_f32(ST_.data(), qjl_f.data(), correction.data(), d_, d_);

        for (int j = 0; j < d_; ++j)
            out[j] += scale * correction[j];
    }

    // -----------------------------------------------------------------
    // Unbiased inner product estimator
    // -----------------------------------------------------------------

    /// Compute ⟨y, x̃⟩ where x̃ is the TurboQuantProd approximation of x.
    /// By Theorem 2: E[result] = ⟨y, x⟩ (unbiased).
    float inner_product(const QuantizedVectorProd& qx,
                         std::span<const float> y) const {
        if (static_cast<int>(y.size()) != d_)
            throw std::invalid_argument("TurboQuantProd::inner_product: wrong dimension");

        // MSE part: ⟨y, x̃_mse⟩
        float ip_mse = mse_quant_.approx_inner_product(qx.mse_part, y);

        if (qx.residual_norm < 1e-9f) return ip_mse;

        // QJL correction: γ · (√(π/2)/d) · y^T Sᵀ qjl  =  γ · (√(π/2)/d) · (Sy)^T qjl
        static const float kQjlScale = std::sqrt(float(M_PI) / 2.0f);
        float scale = qx.residual_norm * kQjlScale / static_cast<float>(d_);

        // Sy = S · y
        thread_local std::vector<float> Sy;
        Sy.resize(d_);
        detail::matvec_f32(S_.data(), y.data(), Sy.data(), d_, d_);

        // dot(Sy, qjl)
        float dot_Sy_qjl = 0.0f;
        for (int j = 0; j < d_; ++j)
            dot_Sy_qjl += Sy[j] * static_cast<float>(qx.qjl_signs[j]);

        return ip_mse + scale * dot_Sy_qjl;
    }

    // -----------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------

    int dimension() const noexcept { return d_; }
    int bit_width() const noexcept { return bit_width_; }

    const TurboQuantMSE& mse_quantizer() const noexcept { return mse_quant_; }

    /// Theoretical distortion bound (Theorem 2): sqrt(3)·π²·||y||²/(d·4^b)
    /// Returned as coefficient; caller multiplies by ||y||²/d.
    float theoretical_distortion_coeff() const noexcept {
        static const float kCoeff = float(std::sqrt(3.0) * M_PI * M_PI); // ≈17.3
        float denom = 1.0f;
        for (int b = 0; b < bit_width_; ++b) denom *= 4.0f;
        return kCoeff / denom;
    }

private:
    int           d_{0};
    int           bit_width_{1};
    TurboQuantMSE mse_quant_;
    uint64_t      seed_{42};

    std::vector<float> S_;    ///< QJL matrix d×d, row-major, N(0,1)
    std::vector<float> ST_;   ///< Transpose of S, row-major
};

// =============================================================================
// QuantizedContactMatrix
// FlexAIDdS-specific: compresses the 256×256 SoftContactMatrix
// =============================================================================

/// Quantized representation of the SoftContactMatrix.
///
/// The 256×256 soft_contact_matrix represents pairwise interaction potentials
/// between the 256 atom types. Each row i is a 256-dim interaction profile
/// vector. TurboQuantMSE compresses each row to b bits/coordinate.
///
/// Compression factor (b=2): 256×256×4 bytes = 256 KB → 256×64 B = 16 KB
/// → 16× compression with theoretically bounded MSE distortion.
///
/// Usage in gaboom.cpp:
///   QuantizedContactMatrix qcm(/*bit_width=*/2);
///   qcm.build(soft_contact_matrix.data());
///   float score = qcm.approximate_score(atom_type_i, atom_type_j);
class QuantizedContactMatrix {
public:
    static constexpr int kNumAtomTypes = 256;
    static constexpr int kDim          = 256;

    explicit QuantizedContactMatrix(int bit_width = 2, uint64_t seed = 42)
        : quantizer_(kDim, bit_width, seed), bit_width_(bit_width)
    {}

    // -----------------------------------------------------------------
    // Build from raw 256×256 float matrix
    // -----------------------------------------------------------------

    /// Build compressed representation from flat row-major float array.
    /// matrix_data must point to at least 256*256 floats.
    void build(const float* matrix_data) {
        std::span<const float> mat(matrix_data, kNumAtomTypes * kDim);
        auto batch = quantizer_.quantize_batch(mat, kNumAtomTypes);
        for (int i = 0; i < kNumAtomTypes; ++i)
            rows_[i] = std::move(batch[i]);
    }

    // -----------------------------------------------------------------
    // Approximate scoring
    // -----------------------------------------------------------------

    /// Approximate the (type_i, type_j) interaction score.
    ///
    /// In the original matrix: score = row_i[type_j]
    /// Here: we dequantize row_i and return element type_j.
    /// For the full inner-product approximation (Boltzmann weighting),
    /// use batch_score with one-hot y vectors.
    float approximate_score(int type_i, int type_j) const {
        if (type_i < 0 || type_i >= kNumAtomTypes ||
            type_j < 0 || type_j >= kNumAtomTypes)
            throw std::out_of_range("QuantizedContactMatrix::approximate_score: index out of range");

        // Dequantize row type_i
        thread_local std::vector<float> buf(kDim);
        quantizer_.dequantize(rows_[type_i], std::span<float>(buf));
        return buf[type_j];
    }

    /// Batch score: compute approximate scores for a list of (type_i, type_j) pairs.
    void batch_score(std::span<const std::pair<int,int>> pairs,
                     std::span<float> scores_out) const
    {
        if (pairs.size() != scores_out.size())
            throw std::invalid_argument("batch_score: size mismatch");

        int n = static_cast<int>(pairs.size());

#ifdef TURBOQUANT_HAS_OMP
#pragma omp parallel for schedule(dynamic, 64)
#endif
        for (int k = 0; k < n; ++k) {
            scores_out[k] = approximate_score(pairs[k].first, pairs[k].second);
        }
    }

    /// Fast approximate inner product ⟨row_i, y⟩ without full dequantize.
    float row_inner_product(int type_i, std::span<const float> y) const {
        if (type_i < 0 || type_i >= kNumAtomTypes)
            throw std::out_of_range("row_inner_product: type_i out of range");
        return quantizer_.approx_inner_product(rows_[type_i], y);
    }

    // -----------------------------------------------------------------
    // Memory & compression stats
    // -----------------------------------------------------------------

    size_t memory_bytes() const noexcept {
        size_t total = 0;
        for (const auto& r : rows_)
            total += r.packed_indices.size() + sizeof(r.norm)
                   + sizeof(r.dim) + sizeof(r.bit_width);
        // Rotation matrix: d×d float32
        total += static_cast<size_t>(kDim * kDim) * sizeof(float) * 2; // Π + Πᵀ
        return total;
    }

    float compression_ratio() const noexcept {
        float original = static_cast<float>(kNumAtomTypes * kDim) * sizeof(float);
        float compressed = 0.0f;
        for (const auto& r : rows_)
            compressed += static_cast<float>(r.packed_indices.size());
        return (compressed > 0.0f) ? original / compressed : 0.0f;
    }

    int bit_width() const noexcept { return bit_width_; }
    const TurboQuantMSE& quantizer() const noexcept { return quantizer_; }

private:
    TurboQuantMSE                              quantizer_;
    std::array<TurboQuantMSE::QuantizedVector, kNumAtomTypes> rows_;
    int                                        bit_width_;
};

// =============================================================================
// QuantizedEnsemble
// FlexAIDdS-specific: GA population energy vector quantization for
// StatMechEngine Boltzmann-weighted thermodynamics
// =============================================================================

/// Quantized representation of a GA conformational ensemble.
///
/// Each chromosome/conformer in the GA population is associated with a
/// d-dimensional energy descriptor (e.g., component-decomposed CF values).
/// TurboQuantProd compresses these vectors while preserving inner product
/// structure, which is essential for Boltzmann weight computations:
///   w_i = exp(-β·E_i) / Z,   Z = Σ_i exp(-β·E_i)
/// and partition-function contractions over the ensemble.
///
/// Usage in StatMechEngine:
///   QuantizedEnsemble qens(energy_dim, /*bit_width=*/3);
///   for (auto& chr : population) qens.add_state(chr.energy_descriptor);
///   float bw = qens.approximate_boltzmann_weight(i, beta_vector);
class QuantizedEnsemble {
public:
    explicit QuantizedEnsemble(int energy_dim, int bit_width = 3)
        : quantizer_(energy_dim, bit_width), energy_dim_(energy_dim)
    {
        if (energy_dim_ <= 0)
            throw std::invalid_argument("QuantizedEnsemble: energy_dim must be > 0");
    }

    // -----------------------------------------------------------------
    // Add states
    // -----------------------------------------------------------------

    void add_state(std::span<const float> energy_descriptor) {
        if (static_cast<int>(energy_descriptor.size()) != energy_dim_)
            throw std::invalid_argument("QuantizedEnsemble::add_state: dimension mismatch");
        states_.push_back(quantizer_.quantize(energy_descriptor));
    }

    void reserve(int n) { states_.reserve(n); }
    void clear()        { states_.clear(); }

    // -----------------------------------------------------------------
    // Boltzmann-weighted inner product
    // -----------------------------------------------------------------

    /// Compute approximate ⟨beta_E, x̃_i⟩ for state i.
    ///
    /// In StatMechEngine: beta_E is β times the energy gradient/profile vector.
    /// The Boltzmann weight for state i is proportional to exp(-⟨beta_E, x_i⟩).
    /// TurboQuantProd provides an unbiased estimate of this inner product.
    float approximate_boltzmann_weight(int state_i,
                                        std::span<const float> beta_E) const
    {
        if (state_i < 0 || state_i >= static_cast<int>(states_.size()))
            throw std::out_of_range("QuantizedEnsemble::approximate_boltzmann_weight: index out of range");
        if (static_cast<int>(beta_E.size()) != energy_dim_)
            throw std::invalid_argument("QuantizedEnsemble::approximate_boltzmann_weight: dimension mismatch");

        return quantizer_.inner_product(states_[state_i], beta_E);
    }

    /// Compute all Boltzmann weights exp(-ip_i) and return partition function.
    /// weights_out must have size() == states_count().
    float compute_partition_function(std::span<const float> beta_E,
                                     std::span<float> weights_out) const
    {
        if (static_cast<int>(weights_out.size()) != static_cast<int>(states_.size()))
            throw std::invalid_argument("compute_partition_function: size mismatch");

        float log_Z = 0.0f;
        (void)log_Z;
        int n = static_cast<int>(states_.size());

        // Compute unnormalised log-weights
        thread_local std::vector<float> log_w;
        log_w.resize(n);

#ifdef TURBOQUANT_HAS_OMP
#pragma omp parallel for schedule(dynamic, 64)
#endif
        for (int i = 0; i < n; ++i)
            log_w[i] = -approximate_boltzmann_weight(i, beta_E);

        // Numerically stable softmax-style
        float log_max = *std::max_element(log_w.begin(), log_w.end());
        float sum = 0.0f;
        for (int i = 0; i < n; ++i) {
            weights_out[i] = std::exp(log_w[i] - log_max);
            sum += weights_out[i];
        }
        // Normalise
        if (sum > 0.0f)
            for (auto& w : weights_out) w /= sum;

        return std::log(sum) + log_max;  // log(Z)
    }

    // -----------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------

    size_t size()       const noexcept { return states_.size(); }
    int    energy_dim() const noexcept { return energy_dim_; }
    const TurboQuantProd& quantizer() const noexcept { return quantizer_; }

    const TurboQuantProd::QuantizedVectorProd& state(int i) const {
        return states_.at(i);
    }

    size_t memory_bytes() const noexcept {
        size_t total = 0;
        for (const auto& s : states_) {
            total += s.mse_part.packed_indices.size()
                   + s.qjl_signs.size() * sizeof(int8_t)
                   + sizeof(s.residual_norm)
                   + sizeof(s.mse_part.norm)
                   + sizeof(s.mse_part.dim)
                   + sizeof(s.mse_part.bit_width);
        }
        // Rotation matrices (Π, Πᵀ, S, Sᵀ) for both sub-quantizers
        size_t mat = static_cast<size_t>(energy_dim_) * energy_dim_ * sizeof(float);
        total += 4 * mat;
        return total;
    }

private:
    TurboQuantProd                                          quantizer_;
    std::vector<TurboQuantProd::QuantizedVectorProd>        states_;
    int                                                     energy_dim_;
};

// =============================================================================
// NearestNeighborIndex
// FlexAIDdS-specific: compressed NN index for FastOPTICS super-cluster detection
// =============================================================================

/// Approximate nearest-neighbor index backed by TurboQuantMSE.
///
/// FastOPTICS in gaboom.cpp computes pairwise distances over the GA population
/// conformations. This index stores quantized conformer descriptors and
/// approximates squared L2 distances via:
///   ||x - y||² ≈ ||x̃ - ỹ||²  (preserved up to 2.7/4^b per unit sphere)
///
/// For unnormalised vectors: scale bound by ||x||² + ||y||².
class NearestNeighborIndex {
public:
    explicit NearestNeighborIndex(int dim, int bit_width = 2, uint64_t seed = 42)
        : quantizer_(dim, bit_width, seed), dim_(dim)
    {}

    void add(std::span<const float> v) {
        vectors_.push_back(quantizer_.quantize(v));
        // Dequantize to get x̃ and store its norm for consistent distance formula
        thread_local std::vector<float> xhat;
        xhat.resize(dim_);
        quantizer_.dequantize(vectors_.back(), std::span<float>(xhat));
        float norm2 = 0.0f;
        for (auto x : xhat) norm2 += x * x;
        norms_.push_back(norm2);  // store ||x̃||²
    }

    void reserve(int n) { vectors_.reserve(n); norms_.reserve(n); }
    void clear()        { vectors_.clear(); norms_.clear(); }
    size_t size() const { return vectors_.size(); }

    /// Approximate squared L2 distance between stored vector i and query q.
    /// Uses the identity ||x̃_i - q||^2 ≈ ||x_i - q||^2.
    /// Formula: ||x̃_i||^2 - 2⟨x̃_i, q⟩ + ||q||^2  (all in original space)
    float approx_sq_distance(int i, std::span<const float> q) const {
        // approx_inner_product returns ⟨x̃_i, Π·q⟩ = ⟨Πᵀx̃_rot_i, q⟩ = ⟨x̃_i, q⟩
        float qi_dot = quantizer_.approx_inner_product(vectors_[i], q);
        float norm_q2 = 0.0f;
        for (auto v : q) norm_q2 += v * v;
        // norms_[i] already stores ||x̃_i||^2
        return norms_[i] - 2.0f * qi_dot + norm_q2;
    }

    /// Find k approximate nearest neighbors of query q.
    /// Returns (index, approx_sq_dist) pairs, sorted by distance.
    std::vector<std::pair<int,float>> knn(std::span<const float> q, int k) const {
        int n = static_cast<int>(vectors_.size());
        std::vector<std::pair<float,int>> dists(n);

#ifdef TURBOQUANT_HAS_OMP
#pragma omp parallel for schedule(dynamic, 64)
#endif
        for (int i = 0; i < n; ++i)
            dists[i] = {approx_sq_distance(i, q), i};

        int top_k = std::min(k, n);
        std::partial_sort(dists.begin(), dists.begin() + top_k, dists.end());

        std::vector<std::pair<int,float>> result(top_k);
        for (int i = 0; i < top_k; ++i)
            result[i] = {dists[i].second, dists[i].first};
        return result;
    }

    const TurboQuantMSE& quantizer() const noexcept { return quantizer_; }

private:
    TurboQuantMSE                              quantizer_;
    std::vector<TurboQuantMSE::QuantizedVector> vectors_;
    std::vector<float>                          norms_;
    int                                         dim_;
};

// =============================================================================
// Utility: distortion diagnostics
// =============================================================================

namespace util {

/// Measure empirical MSE of TurboQuantMSE on a random sample.
/// Returns empirical MSE (mean over n_samples of ||x - x̃||²).
inline float empirical_mse(const TurboQuantMSE& q, int n_samples = 1024,
                            uint64_t seed = 0xFEEDFACE) {
    int d = q.dimension();
    detail::XorShift64 rng(seed);

    std::vector<float> x(d);
    std::vector<float> xhat(d);
    double total_mse = 0.0;

    for (int s = 0; s < n_samples; ++s) {
        // Sample unit-sphere vector
        double norm = 0.0;
        for (int j = 0; j < d; j += 2) {
            auto [g1, g2] = rng.next_gaussian_pair();
            x[j]     = static_cast<float>(g1);
            if (j + 1 < d) x[j + 1] = static_cast<float>(g2);
            norm += g1 * g1;
            if (j + 1 < d) norm += g2 * g2;
        }
        float inv_norm = 1.0f / static_cast<float>(std::sqrt(norm));
        for (auto& v : x) v *= inv_norm;

        // Quantize → dequantize
        auto qv = q.quantize(std::span<const float>(x));
        q.dequantize(qv, std::span<float>(xhat));

        // MSE
        double mse = 0.0;
        for (int j = 0; j < d; ++j) {
            double diff = x[j] - xhat[j];
            mse += diff * diff;
        }
        total_mse += mse;
    }

    return static_cast<float>(total_mse / n_samples);
}

/// Measure empirical inner-product bias of TurboQuantProd.
/// Returns: (mean_error, std_error)
inline std::pair<float,float> empirical_ip_bias(const TurboQuantProd& q,
                                                  int n_samples = 512,
                                                  uint64_t seed = 0xDEADBEEF) {
    int d = q.dimension();
    detail::XorShift64 rng(seed);

    std::vector<float> x(d), y(d);
    double sum_err = 0.0, sum_err2 = 0.0;

    for (int s = 0; s < n_samples; ++s) {
        // Random unit sphere x
        double nx = 0.0;
        for (int j = 0; j < d; j += 2) {
            auto [g1, g2] = rng.next_gaussian_pair();
            x[j] = static_cast<float>(g1);
            if (j+1 < d) x[j+1] = static_cast<float>(g2);
            nx += g1*g1; if (j+1 < d) nx += g2*g2;
        }
        float ix = 1.0f / static_cast<float>(std::sqrt(nx));
        for (auto& v : x) v *= ix;

        // Random y (not normalised for generality)
        for (int j = 0; j < d; j += 2) {
            auto [g1, g2] = rng.next_gaussian_pair();
            y[j] = static_cast<float>(g1);
            if (j+1 < d) y[j+1] = static_cast<float>(g2);
        }

        float true_ip = 0.0f;
        for (int j = 0; j < d; ++j) true_ip += x[j] * y[j];

        auto qv   = q.quantize(std::span<const float>(x));
        float est  = q.inner_product(qv, std::span<const float>(y));

        double err = est - true_ip;
        sum_err  += err;
        sum_err2 += err * err;
    }

    float mean_err = static_cast<float>(sum_err  / n_samples);
    float var      = static_cast<float>(sum_err2 / n_samples - (sum_err/n_samples)*(sum_err/n_samples));
    return {mean_err, std::sqrt(std::max(0.0f, var))};
}

/// Return a human-readable summary string for a Codebook.
inline std::string codebook_summary(const Codebook& cb) {
    std::string s = "Codebook(b=" + std::to_string(cb.bit_width)
                  + ", k=" + std::to_string(cb.num_centroids) + "): [";
    for (int i = 0; i < cb.num_centroids; ++i) {
        if (i > 0) s += ", ";
        // Format to 5 decimal places
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%.5f", cb.centroids[i]);
        s += buf;
    }
    s += "]";
    return s;
}

} // namespace util

// =============================================================================
} // namespace turboquant
// =============================================================================

// Restore any platform-specific state
#undef TQ_HD
#undef TQ_DEVICE
