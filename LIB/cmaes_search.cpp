// cmaes_search.cpp — Clean-room CMA-ES (Hansen 2006 / Tutorial CMA-ES style)
// Apache-2.0. No GPL code. No modifications to gaboom.cpp / ic2cf.cpp.
//
// Algorithm: rank-μ update with cumulative step-size adaptation (CSA).
//   dim <= 40 → full covariance (eigen decomposition each generation as needed)
//   dim >  40 → diagonal CMA (O(n) covariance storage/update)
//
// Engine seam (extern via gaboom.h / flexaid.h):
//   set_gene_lim, set_bins, eval_chromosome, get_cf_evalue, get_apparent_cf_evalue

#include "cmaes_search.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <utility>

// gaboom.h already included via cmaes_search.h (set_gene_lim, set_bins,
// eval_chromosome, gene/chromosome types). flexaid.h via gaboom.h
// (get_cf_evalue, get_apparent_cf_evalue, FA_Global, etc.).

namespace {

// Boltzmann constant kcal/mol/K (matches kB_kcal convention in engine)
constexpr double kB_kcal = 0.001987204258;
constexpr double kT_300  = kB_kcal * 300.0;  // ≈ 0.596

// ── Small dense linear algebra (no Eigen dependency for adapter isolation) ───

using Vec = std::vector<double>;
using Mat = std::vector<double>;  // row-major n×n

inline double& mat_at(Mat& M, int n, int i, int j) {
    return M[static_cast<std::size_t>(i * n + j)];
}
inline double mat_at(const Mat& M, int n, int i, int j) {
    return M[static_cast<std::size_t>(i * n + j)];
}

// Symmetric Jacobi eigen-decomposition: A (destroyed) → V D V^T.
// Eigenvalues written to d (length n). Sufficient for CMA (n ≤ 40 full path).
void jacobi_eigen_sym(Mat& A, int n, Vec& d, Mat& V, int max_sweeps = 64) {
    V.assign(static_cast<std::size_t>(n * n), 0.0);
    d.assign(static_cast<std::size_t>(n), 0.0);
    for (int i = 0; i < n; ++i) {
        mat_at(V, n, i, i) = 1.0;
        d[static_cast<std::size_t>(i)] = mat_at(A, n, i, i);
    }
    if (n == 1) return;

    Vec b = d;
    Vec z(static_cast<std::size_t>(n), 0.0);

    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
        double off = 0.0;
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                off += std::fabs(mat_at(A, n, i, j));
        if (off < 1e-15 * static_cast<double>(n * n)) break;

        const double thresh = (sweep < 3) ? (0.2 * off / (n * n)) : 0.0;
        for (int p = 0; p < n; ++p) {
            for (int q = p + 1; q < n; ++q) {
                double& apq = mat_at(A, n, p, q);
                const double g = 100.0 * std::fabs(apq);
                if (sweep > 3 &&
                    (std::fabs(d[static_cast<std::size_t>(p)]) + g) ==
                        std::fabs(d[static_cast<std::size_t>(p)]) &&
                    (std::fabs(d[static_cast<std::size_t>(q)]) + g) ==
                        std::fabs(d[static_cast<std::size_t>(q)])) {
                    apq = 0.0;
                    continue;
                }
                if (std::fabs(apq) <= thresh) continue;

                double h = d[static_cast<std::size_t>(q)] - d[static_cast<std::size_t>(p)];
                double t;
                if ((std::fabs(h) + g) == std::fabs(h)) {
                    t = apq / h;
                } else {
                    const double theta = 0.5 * h / apq;
                    t = 1.0 / (std::fabs(theta) + std::sqrt(1.0 + theta * theta));
                    if (theta < 0.0) t = -t;
                }
                const double c = 1.0 / std::sqrt(1.0 + t * t);
                const double s = t * c;
                const double tau = s / (1.0 + c);
                h = t * apq;
                z[static_cast<std::size_t>(p)] -= h;
                z[static_cast<std::size_t>(q)] += h;
                d[static_cast<std::size_t>(p)] -= h;
                d[static_cast<std::size_t>(q)] += h;
                apq = 0.0;

                for (int j = 0; j < p; ++j) {
                    double& ajp = mat_at(A, n, j, p);
                    double& ajq = mat_at(A, n, j, q);
                    const double x = ajp, y = ajq;
                    ajp = x - s * (y + x * tau);
                    ajq = y + s * (x - y * tau);
                }
                for (int j = p + 1; j < q; ++j) {
                    double& apj = mat_at(A, n, p, j);
                    double& ajq = mat_at(A, n, j, q);
                    const double x = apj, y = ajq;
                    apj = x - s * (y + x * tau);
                    ajq = y + s * (x - y * tau);
                }
                for (int j = q + 1; j < n; ++j) {
                    double& apj = mat_at(A, n, p, j);
                    double& aqj = mat_at(A, n, q, j);
                    const double x = apj, y = aqj;
                    apj = x - s * (y + x * tau);
                    aqj = y + s * (x - y * tau);
                }
                for (int j = 0; j < n; ++j) {
                    double& vjp = mat_at(V, n, j, p);
                    double& vjq = mat_at(V, n, j, q);
                    const double x = vjp, y = vjq;
                    vjp = x - s * (y + x * tau);
                    vjq = y + s * (x - y * tau);
                }
            }
        }
        for (int i = 0; i < n; ++i) {
            b[static_cast<std::size_t>(i)] += z[static_cast<std::size_t>(i)];
            d[static_cast<std::size_t>(i)] = b[static_cast<std::size_t>(i)];
            z[static_cast<std::size_t>(i)] = 0.0;
        }
    }
    for (int i = 0; i < n; ++i)
        if (d[static_cast<std::size_t>(i)] < 0.0) d[static_cast<std::size_t>(i)] = 0.0;
}

// ── CMA-ES core state ────────────────────────────────────────────────────────

struct CmaState {
    int    n        = 0;
    int    lambda   = 0;
    int    mu       = 0;
    bool   diagonal = false;
    double sigma    = 0.3;
    double chiN     = 0.0;

    double c_sigma = 0.0, d_sigma = 0.0, c_c = 0.0, c_1 = 0.0, c_mu = 0.0;
    double mu_eff  = 0.0;

    Vec mean;     // n
    Vec pc, ps;   // evolution paths
    Vec weights;  // μ
    Vec diag_C;   // n
    Mat C;        // n×n when full
    Mat B;        // eigenvectors
    Vec D;        // sqrt(eigenvalues)
    int  count_eval = 0;
    int  count_eval_since_eigen = 0;
    int  generation = 0;
};

void cma_set_weights(CmaState& st, int mu) {
    st.mu = mu;
    st.weights.resize(static_cast<std::size_t>(mu));
    double wsum = 0.0;
    for (int i = 0; i < mu; ++i) {
        st.weights[static_cast<std::size_t>(i)] =
            std::log(static_cast<double>(mu) + 0.5) - std::log(static_cast<double>(i) + 1.0);
        wsum += st.weights[static_cast<std::size_t>(i)];
    }
    for (double& w : st.weights) w /= wsum;
    double w2 = 0.0;
    for (double w : st.weights) w2 += w * w;
    st.mu_eff = 1.0 / w2;
}

void cma_init(CmaState& st, int n, int lambda, double sigma0) {
    st.n = n;
    st.lambda = std::max(4, lambda);
    st.diagonal = (n > 40);
    st.sigma = (sigma0 > 0.0) ? sigma0 : 0.3;
    st.mean.assign(static_cast<std::size_t>(n), 0.0);
    st.pc.assign(static_cast<std::size_t>(n), 0.0);
    st.ps.assign(static_cast<std::size_t>(n), 0.0);
    st.diag_C.assign(static_cast<std::size_t>(n), 1.0);
    st.count_eval = 0;
    st.count_eval_since_eigen = 0;
    st.generation = 0;

    cma_set_weights(st, st.lambda / 2);
    if (st.mu < 1) cma_set_weights(st, 1);

    const double nn = static_cast<double>(n);
    st.c_sigma = (st.mu_eff + 2.0) / (nn + st.mu_eff + 5.0);
    st.d_sigma = 1.0 + 2.0 * std::max(0.0, std::sqrt((st.mu_eff - 1.0) / (nn + 1.0)) - 1.0) +
                 st.c_sigma;
    st.c_c = (4.0 + st.mu_eff / nn) / (nn + 4.0 + 2.0 * st.mu_eff / nn);
    st.c_1 = 2.0 / ((nn + 1.3) * (nn + 1.3) + st.mu_eff);
    st.c_mu = std::min(1.0 - st.c_1,
                       2.0 * (st.mu_eff - 2.0 + 1.0 / st.mu_eff) /
                           ((nn + 2.0) * (nn + 2.0) + st.mu_eff));

    st.chiN = std::sqrt(nn) * (1.0 - 1.0 / (4.0 * nn) + 1.0 / (21.0 * nn * nn));

    if (!st.diagonal) {
        st.C.assign(static_cast<std::size_t>(n * n), 0.0);
        for (int i = 0; i < n; ++i) mat_at(st.C, n, i, i) = 1.0;
        st.B = st.C;
        st.D.assign(static_cast<std::size_t>(n), 1.0);
    }
}

void cma_update_eigensystem(CmaState& st) {
    if (st.diagonal) return;
    Mat A = st.C;
    for (int i = 0; i < st.n; ++i)
        for (int j = i + 1; j < st.n; ++j) {
            const double v = 0.5 * (mat_at(A, st.n, i, j) + mat_at(A, st.n, j, i));
            mat_at(A, st.n, i, j) = v;
            mat_at(A, st.n, j, i) = v;
        }
    Vec evals;
    Mat V;
    jacobi_eigen_sym(A, st.n, evals, V);
    st.B = std::move(V);
    st.D.resize(static_cast<std::size_t>(st.n));
    for (int i = 0; i < st.n; ++i) {
        const double e = std::max(evals[static_cast<std::size_t>(i)], 1e-30);
        st.D[static_cast<std::size_t>(i)] = std::sqrt(e);
        st.diag_C[static_cast<std::size_t>(i)] = e;
    }
    st.count_eval_since_eigen = 0;
}

// Sample one offspring: x = m + σ * y, with y ~ N(0, C).
void cma_sample(const CmaState& st, std::mt19937& rng, Vec& x, Vec& y) {
    std::normal_distribution<double> N01(0.0, 1.0);
    const int n = st.n;
    x.resize(static_cast<std::size_t>(n));
    y.resize(static_cast<std::size_t>(n));

    if (st.diagonal) {
        for (int i = 0; i < n; ++i) {
            const double z = N01(rng);
            y[static_cast<std::size_t>(i)] =
                std::sqrt(std::max(st.diag_C[static_cast<std::size_t>(i)], 1e-30)) * z;
            x[static_cast<std::size_t>(i)] =
                st.mean[static_cast<std::size_t>(i)] + st.sigma * y[static_cast<std::size_t>(i)];
        }
        return;
    }

    Vec z(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) z[static_cast<std::size_t>(i)] = N01(rng);
    for (int i = 0; i < n; ++i) {
        double yi = 0.0;
        for (int j = 0; j < n; ++j)
            yi += mat_at(st.B, n, i, j) * st.D[static_cast<std::size_t>(j)] *
                  z[static_cast<std::size_t>(j)];
        y[static_cast<std::size_t>(i)] = yi;
        x[static_cast<std::size_t>(i)] = st.mean[static_cast<std::size_t>(i)] + st.sigma * yi;
    }
}

void cma_tell(CmaState& st, const std::vector<Vec>& Y_sorted) {
    const int n = st.n;
    const int mu = st.mu;

    Vec y_w(static_cast<std::size_t>(n), 0.0);
    for (int k = 0; k < mu; ++k) {
        const double w = st.weights[static_cast<std::size_t>(k)];
        for (int i = 0; i < n; ++i)
            y_w[static_cast<std::size_t>(i)] +=
                w * Y_sorted[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)];
    }

    for (int i = 0; i < n; ++i)
        st.mean[static_cast<std::size_t>(i)] += st.sigma * y_w[static_cast<std::size_t>(i)];

    // C^{-1/2} y_w for ps update
    Vec invsqrt_y(static_cast<std::size_t>(n), 0.0);
    if (st.diagonal) {
        for (int i = 0; i < n; ++i)
            invsqrt_y[static_cast<std::size_t>(i)] =
                y_w[static_cast<std::size_t>(i)] /
                std::sqrt(std::max(st.diag_C[static_cast<std::size_t>(i)], 1e-30));
    } else {
        Vec tmp(static_cast<std::size_t>(n), 0.0);
        for (int j = 0; j < n; ++j) {
            double s = 0.0;
            for (int i = 0; i < n; ++i)
                s += mat_at(st.B, n, i, j) * y_w[static_cast<std::size_t>(i)];
            tmp[static_cast<std::size_t>(j)] =
                s / std::max(st.D[static_cast<std::size_t>(j)], 1e-30);
        }
        for (int i = 0; i < n; ++i) {
            double s = 0.0;
            for (int j = 0; j < n; ++j)
                s += mat_at(st.B, n, i, j) * tmp[static_cast<std::size_t>(j)];
            invsqrt_y[static_cast<std::size_t>(i)] = s;
        }
    }

    const double cs = st.c_sigma;
    const double fac_ps = std::sqrt(cs * (2.0 - cs) * st.mu_eff);
    for (int i = 0; i < n; ++i)
        st.ps[static_cast<std::size_t>(i)] =
            (1.0 - cs) * st.ps[static_cast<std::size_t>(i)] +
            fac_ps * invsqrt_y[static_cast<std::size_t>(i)];

    double ps_norm2 = 0.0;
    for (double v : st.ps) ps_norm2 += v * v;
    const double ps_norm = std::sqrt(ps_norm2);

    const double gen = static_cast<double>(st.generation + 1);
    const double ps_expect =
        std::sqrt(1.0 - std::pow(1.0 - cs, 2.0 * gen)) * st.chiN;
    const double hsig_thresh = 1.4 + 2.0 / (static_cast<double>(n) + 1.0);
    const int hsig =
        (ps_norm / std::max(ps_expect, 1e-30) < hsig_thresh) ? 1 : 0;

    const double cc = st.c_c;
    const double fac_pc = static_cast<double>(hsig) * std::sqrt(cc * (2.0 - cc) * st.mu_eff);
    for (int i = 0; i < n; ++i)
        st.pc[static_cast<std::size_t>(i)] =
            (1.0 - cc) * st.pc[static_cast<std::size_t>(i)] +
            fac_pc * y_w[static_cast<std::size_t>(i)];

    if (st.diagonal) {
        for (int i = 0; i < n; ++i) {
            const double pci = st.pc[static_cast<std::size_t>(i)];
            double rank_mu = 0.0;
            for (int k = 0; k < mu; ++k) {
                const double yi =
                    Y_sorted[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)];
                rank_mu += st.weights[static_cast<std::size_t>(k)] * yi * yi;
            }
            st.diag_C[static_cast<std::size_t>(i)] =
                (1.0 - st.c_1 - st.c_mu) * st.diag_C[static_cast<std::size_t>(i)] +
                st.c_1 * (pci * pci +
                          (1 - hsig) * cc * (2.0 - cc) *
                              st.diag_C[static_cast<std::size_t>(i)]) +
                st.c_mu * rank_mu;
            if (st.diag_C[static_cast<std::size_t>(i)] < 1e-30)
                st.diag_C[static_cast<std::size_t>(i)] = 1e-30;
        }
    } else {
        Mat Cnew = st.C;
        const double decay = 1.0 - st.c_1 - st.c_mu;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
                double v = decay * mat_at(st.C, n, i, j);
                v += st.c_1 *
                     (st.pc[static_cast<std::size_t>(i)] * st.pc[static_cast<std::size_t>(j)] +
                      (1 - hsig) * cc * (2.0 - cc) * mat_at(st.C, n, i, j));
                double rank_mu = 0.0;
                for (int k = 0; k < mu; ++k) {
                    rank_mu +=
                        st.weights[static_cast<std::size_t>(k)] *
                        Y_sorted[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)] *
                        Y_sorted[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)];
                }
                v += st.c_mu * rank_mu;
                mat_at(Cnew, n, i, j) = v;
                mat_at(Cnew, n, j, i) = v;
            }
        }
        st.C = std::move(Cnew);
        for (int i = 0; i < n; ++i)
            st.diag_C[static_cast<std::size_t>(i)] = mat_at(st.C, n, i, i);
    }

    st.sigma *= std::exp((cs / st.d_sigma) * (ps_norm / st.chiN - 1.0));
    if (st.sigma < 1e-20) st.sigma = 1e-20;
    if (st.sigma > 1e6) st.sigma = 1e6;

    st.count_eval += st.lambda;
    st.count_eval_since_eigen += st.lambda;
    ++st.generation;

    const int period = std::max(
        1, static_cast<int>(static_cast<double>(st.lambda) /
                            (st.c_1 * static_cast<double>(n) * 10.0 + 1.0)));
    if (!st.diagonal && st.count_eval_since_eigen >= period) cma_update_eigensystem(st);
}

// ── Helpers ──────────────────────────────────────────────────────────────────

void clamp_vec(Vec& x, const Vec& lo, const Vec& hi) {
    for (std::size_t i = 0; i < x.size(); ++i) {
        if (x[i] < lo[i]) x[i] = lo[i];
        if (x[i] > hi[i]) x[i] = hi[i];
    }
}

double shannon_from_probs(const Vec& p) {
    double H = 0.0;
    for (double v : p) {
        if (v > 0.0) H -= v * std::log(v);
    }
    return H;
}

// Shannon of normalized rank-μ selection weights over the λ sample.
double H_search_from_fitness(const Vec& f, int mu) {
    const int lambda = static_cast<int>(f.size());
    if (lambda <= 0) return 0.0;
    std::vector<int> order(static_cast<std::size_t>(lambda));
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return f[static_cast<std::size_t>(a)] < f[static_cast<std::size_t>(b)];
    });

    Vec w(static_cast<std::size_t>(lambda), 0.0);
    double wsum = 0.0;
    const int m = std::min(mu, lambda);
    for (int r = 0; r < m; ++r) {
        const double wi =
            std::log(static_cast<double>(m) + 0.5) - std::log(static_cast<double>(r) + 1.0);
        w[static_cast<std::size_t>(order[static_cast<std::size_t>(r)])] = wi;
        wsum += wi;
    }
    if (wsum <= 0.0) return 0.0;
    for (double& v : w) v /= wsum;
    return shannon_from_probs(w);
}

// Boltzmann Shannon over sample energies (CF). T in energy units (kcal/mol).
double H_energy_boltzmann(const Vec& f, double T) {
    if (f.empty() || T <= 0.0) return 0.0;
    const double fmin = *std::min_element(f.begin(), f.end());
    Vec p(f.size());
    double Z = 0.0;
    for (std::size_t i = 0; i < f.size(); ++i) {
        p[i] = std::exp(-(f[i] - fmin) / T);
        Z += p[i];
    }
    if (Z <= 0.0) return 0.0;
    for (double& v : p) v /= Z;
    return shannon_from_probs(p);
}

double temperature_energy(FA_Global* FA) {
    // T (kcal/mol) for Boltzmann / F-proxy: kB * T_K.
    // FA->temperature is Kelvin (uint). Fallback 300 K → kT ≈ 0.596.
    if (FA != nullptr && FA->temperature > 0)
        return kB_kcal * static_cast<double>(FA->temperature);
    return kT_300;
}

void archive_insert(CmaesResult& res, const Vec& genes, double cf, double app_cf, int cap) {
    if (cap <= 0) cap = 1;
    for (const auto& g : res.archive_genes) {
        if (g.size() != genes.size()) continue;
        double d2 = 0.0;
        for (std::size_t i = 0; i < g.size(); ++i) {
            const double d = g[i] - genes[i];
            d2 += d * d;
        }
        if (d2 < 1e-12) return;
    }
    res.archive_genes.push_back(genes);
    res.archive_cfs.push_back(cf);
    res.archive_app_cfs.push_back(app_cf);

    std::vector<std::size_t> idx(res.archive_cfs.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
              [&](std::size_t a, std::size_t b) { return res.archive_cfs[a] < res.archive_cfs[b]; });

    std::vector<Vec> ng;
    Vec nc, na;
    const std::size_t keep = std::min(static_cast<std::size_t>(cap), idx.size());
    ng.reserve(keep);
    nc.reserve(keep);
    na.reserve(keep);
    for (std::size_t k = 0; k < keep; ++k) {
        ng.push_back(std::move(res.archive_genes[idx[k]]));
        nc.push_back(res.archive_cfs[idx[k]]);
        na.push_back(res.archive_app_cfs[idx[k]]);
    }
    res.archive_genes = std::move(ng);
    res.archive_cfs = std::move(nc);
    res.archive_app_cfs = std::move(na);
}

// Generic CMA loop over an objective f(x) → double, with box bounds.
int cma_optimize(int n, const Vec& lo, const Vec& hi, const CmaesConfig& config,
                 const std::function<double(const double*)>& objective, CmaesResult* result,
                 std::vector<EntropyTraceSample>* optional_trace, double T_energy,
                 const std::function<double()>& last_app_cf = {}) {
    if (!result || n <= 0) return -1;

    try {
        result->best_cf = std::numeric_limits<double>::infinity();
        result->best_app_cf = std::numeric_limits<double>::infinity();
        result->best_genes.assign(static_cast<std::size_t>(n), 0.0);
        result->n_evals = 0;
        result->n_gens = 0;
        result->status = 0;
        result->archive_genes.clear();
        result->archive_cfs.clear();
        result->archive_app_cfs.clear();
        if (optional_trace) optional_trace->clear();

        const int lambda = std::max(4, config.population);
        CmaState st;
        cma_init(st, n, lambda, config.sigma0);

        for (int i = 0; i < n; ++i) {
            st.mean[static_cast<std::size_t>(i)] =
                0.5 * (lo[static_cast<std::size_t>(i)] + hi[static_cast<std::size_t>(i)]);
        }
        if (config.sigma0 <= 0.0) {
            double span = 0.0;
            for (int i = 0; i < n; ++i)
                span += (hi[static_cast<std::size_t>(i)] - lo[static_cast<std::size_t>(i)]);
            span /= static_cast<double>(n);
            st.sigma = 0.3 * std::max(span, 1e-6);
        }

        std::mt19937 rng(config.seed);
        if (!st.diagonal) cma_update_eigensystem(st);

        const std::int64_t max_evals = std::max<std::int64_t>(lambda, config.max_evals);
        const int mu_nominal = st.mu;

        while (result->n_evals < max_evals) {
            const int remaining = static_cast<int>(
                std::min<std::int64_t>(lambda, max_evals - result->n_evals));
            if (remaining <= 0) break;

            std::vector<Vec> X(static_cast<std::size_t>(remaining));
            std::vector<Vec> Y(static_cast<std::size_t>(remaining));
            Vec fit(static_cast<std::size_t>(remaining));
            Vec app(static_cast<std::size_t>(remaining), 0.0);

            for (int k = 0; k < remaining; ++k) {
                cma_sample(st, rng, X[static_cast<std::size_t>(k)],
                           Y[static_cast<std::size_t>(k)]);
                clamp_vec(X[static_cast<std::size_t>(k)], lo, hi);
                for (int i = 0; i < n; ++i) {
                    Y[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)] =
                        (X[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)] -
                         st.mean[static_cast<std::size_t>(i)]) /
                        st.sigma;
                }
                fit[static_cast<std::size_t>(k)] =
                    objective(X[static_cast<std::size_t>(k)].data());
                if (last_app_cf) app[static_cast<std::size_t>(k)] = last_app_cf();
                ++result->n_evals;

                if (fit[static_cast<std::size_t>(k)] < result->best_cf) {
                    result->best_cf = fit[static_cast<std::size_t>(k)];
                    result->best_app_cf = app[static_cast<std::size_t>(k)];
                    result->best_genes = X[static_cast<std::size_t>(k)];
                }
                archive_insert(*result, X[static_cast<std::size_t>(k)],
                               fit[static_cast<std::size_t>(k)],
                               app[static_cast<std::size_t>(k)], config.archive_size);
            }

            std::vector<int> order(static_cast<std::size_t>(remaining));
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(), [&](int a, int b) {
                return fit[static_cast<std::size_t>(a)] < fit[static_cast<std::size_t>(b)];
            });

            const int mu_use = std::min(mu_nominal, remaining);
            cma_set_weights(st, mu_use);
            // Refresh strategy constants that depend on mu_eff (partial last gen)
            {
                const double nn = static_cast<double>(n);
                st.c_sigma = (st.mu_eff + 2.0) / (nn + st.mu_eff + 5.0);
                st.d_sigma =
                    1.0 +
                    2.0 * std::max(0.0, std::sqrt((st.mu_eff - 1.0) / (nn + 1.0)) - 1.0) +
                    st.c_sigma;
                st.c_c = (4.0 + st.mu_eff / nn) / (nn + 4.0 + 2.0 * st.mu_eff / nn);
                st.c_1 = 2.0 / ((nn + 1.3) * (nn + 1.3) + st.mu_eff);
                st.c_mu = std::min(1.0 - st.c_1,
                                   2.0 * (st.mu_eff - 2.0 + 1.0 / st.mu_eff) /
                                       ((nn + 2.0) * (nn + 2.0) + st.mu_eff));
            }

            std::vector<Vec> Y_sorted(static_cast<std::size_t>(mu_use));
            for (int r = 0; r < mu_use; ++r)
                Y_sorted[static_cast<std::size_t>(r)] =
                    Y[static_cast<std::size_t>(order[static_cast<std::size_t>(r)])];

            cma_tell(st, Y_sorted);
            result->n_gens = st.generation;

            if (config.enable_entropy_trace && optional_trace) {
                EntropyTraceSample s;
                s.gen = st.generation;
                s.H_search = H_search_from_fitness(fit, mu_use);
                s.H_energy = H_energy_boltzmann(fit, T_energy);
                s.best_cf = result->best_cf;
                s.F = result->best_cf - T_energy * s.H_energy;
                s.n_evals = result->n_evals;
                optional_trace->push_back(s);
            }
        }

        if (!config.write_trace.empty() && optional_trace && !optional_trace->empty()) {
            cmaes_write_trace_csv(config.write_trace, *optional_trace);
        }

        result->status = 0;
        return 0;
    } catch (const std::exception& ex) {
        std::fprintf(stderr, "[cmaes] exception: %s\n", ex.what());
        if (result) result->status = -2;
        return -2;
    } catch (...) {
        std::fprintf(stderr, "[cmaes] unknown exception\n");
        if (result) result->status = -2;
        return -2;
    }
}

}  // namespace

// ── Mock objective ───────────────────────────────────────────────────────────

double cmaes_mock_objective(const double* x, int n) {
    if (!x || n <= 0) return 0.0;
    // Smooth multi-dim well: weighted quadratic + mild nearest-neighbor coupling.
    // Global min at x=0 with f=0.
    double f = 0.0;
    for (int i = 0; i < n; ++i) {
        const double xi = x[i];
        f += static_cast<double>(i + 1) * xi * xi;
        if (i + 1 < n) f += 0.15 * xi * x[i + 1];
    }
    return f;
}

int cmaes_run_mock(int dim, const CmaesConfig& config, CmaesResult* result,
                   std::vector<EntropyTraceSample>* optional_trace) {
    if (dim <= 0 || !result) {
        if (result) result->status = -1;
        return -1;
    }
    Vec lo(static_cast<std::size_t>(dim), -5.0);
    Vec hi(static_cast<std::size_t>(dim), 5.0);
    const double T = kT_300;
    return cma_optimize(
        dim, lo, hi, config, [dim](const double* x) { return cmaes_mock_objective(x, dim); },
        result, optional_trace, T);
}

// ── Live dock path ───────────────────────────────────────────────────────────

int cmaes_run_dock(FA_Global* FA, GB_Global* GB, VC_Global* VC, genlim* gene_lim, atom* atoms,
                   resid* residue, gridpoint* cleftgrid, CmaesTargetFn target,
                   const CmaesConfig& config, CmaesResult* result,
                   std::vector<EntropyTraceSample>* optional_trace) {
    if (!FA || !GB || !gene_lim || !result || !target || GB->num_genes <= 0) {
        if (result) result->status = -1;
        return -1;
    }

    try {
        set_gene_lim(FA, GB, gene_lim);
        set_bins(gene_lim, GB->num_genes);

        const int n = GB->num_genes;
        Vec lo(static_cast<std::size_t>(n));
        Vec hi(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) {
            lo[static_cast<std::size_t>(i)] = gene_lim[i].min;
            hi[static_cast<std::size_t>(i)] = gene_lim[i].max;
            if (hi[static_cast<std::size_t>(i)] < lo[static_cast<std::size_t>(i)])
                std::swap(hi[static_cast<std::size_t>(i)], lo[static_cast<std::size_t>(i)]);
        }

        std::vector<gene> john(static_cast<std::size_t>(n));
        double last_app = 0.0;
        gene* john_ptr = john.data();

        auto objective = [&](const double* x) -> double {
            for (int i = 0; i < n; ++i) {
                john_ptr[i].to_ic = x[i];
                john_ptr[i].to_int32 = 0;
            }
            cfstr cf = eval_chromosome(FA, GB, VC, gene_lim, atoms, residue, cleftgrid, john_ptr,
                                       target);
            last_app = get_apparent_cf_evalue(&cf);
            return get_cf_evalue(&cf, FA);
        };

        const double T = temperature_energy(FA);
        return cma_optimize(n, lo, hi, config, objective, result, optional_trace, T,
                            [&]() { return last_app; });
    } catch (const std::exception& ex) {
        std::fprintf(stderr, "[cmaes_run_dock] exception: %s\n", ex.what());
        if (result) result->status = -2;
        return -2;
    } catch (...) {
        std::fprintf(stderr, "[cmaes_run_dock] unknown exception\n");
        if (result) result->status = -2;
        return -2;
    }
}

// ── Snapshot helper ──────────────────────────────────────────────────────────

int cmaes_fill_chromosomes(const CmaesResult& result, int num_genes, chromosome* chrom_out,
                           int max_chrom, gene* gene_storage) {
    if (!chrom_out || !gene_storage || num_genes <= 0 || max_chrom <= 0) return 0;

    int k = 0;

    auto fill_one = [&](const std::vector<double>& genes, double cf, double app) {
        if (k >= max_chrom) return;
        if (static_cast<int>(genes.size()) != num_genes) return;
        gene* slot =
            gene_storage + static_cast<std::size_t>(k) * static_cast<std::size_t>(num_genes);
        for (int g = 0; g < num_genes; ++g) {
            slot[g].to_ic = genes[static_cast<std::size_t>(g)];
            slot[g].to_int32 = 0;
        }
        chrom_out[k].genes = slot;
        chrom_out[k].evalue = cf;
        chrom_out[k].app_evalue = app;
        chrom_out[k].fitnes = 0.0;
        chrom_out[k].boltzmann_weight = 0.0;
        chrom_out[k].free_energy = 0.0;
        chrom_out[k].status = 'n';
        chrom_out[k].cf = cfstr{};
        for (int r = 0; r < MAX_RING_FLEX; ++r) {
            chrom_out[k].ring_phases[r] = 0.0f;
            chrom_out[k].ring_six[r] = 0;
            chrom_out[k].ring_five[r] = 0;
        }
        ++k;
    };

    const int n_arch = static_cast<int>(result.archive_genes.size());
    if (n_arch > 0) {
        for (int i = 0; i < n_arch && k < max_chrom; ++i) {
            fill_one(result.archive_genes[static_cast<std::size_t>(i)],
                     result.archive_cfs[static_cast<std::size_t>(i)],
                     result.archive_app_cfs[static_cast<std::size_t>(i)]);
        }
    } else if (!result.best_genes.empty()) {
        fill_one(result.best_genes, result.best_cf, result.best_app_cf);
    }
    return k;
}

// ── Trace CSV ────────────────────────────────────────────────────────────────

void cmaes_write_trace_csv(const std::string& path,
                           const std::vector<EntropyTraceSample>& samples) {
    if (path.empty()) return;
    std::ofstream out(path);
    if (!out) {
        std::fprintf(stderr, "[cmaes] failed to open trace CSV: %s\n", path.c_str());
        return;
    }
    out << "gen,H_search,H_energy,F,best_cf,n_evals\n";
    out.setf(std::ios::fmtflags(0), std::ios::floatfield);
    out.precision(10);
    for (const auto& s : samples) {
        out << s.gen << ',' << s.H_search << ',' << s.H_energy << ',' << s.F << ',' << s.best_cf
            << ',' << s.n_evals << '\n';
    }
}
