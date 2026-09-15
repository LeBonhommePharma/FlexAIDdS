// tencm_cpu_fallback.cpp — CPU fallback for TENCoM contact discovery & Hessian
//
// Uses Eigen for matrix operations and OpenMP for parallelism.
// Produces results equivalent to tencm_cuda.cu and tencm_metal.mm.
//
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include "tencm_cpu_fallback.h"
#include <Eigen/Dense>
#include <cmath>
#include <cstring>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace tencm { namespace cpu_fallback {

// ─── contact discovery ──────────────────────────────────────────────────────

int build_contacts_cpu(const float* ca_xyz, int N,
                       float cutoff, float k0,
                       int*   contacts_ij,
                       float* contacts_k,
                       float* contacts_r0)
{
    if (N < 3) return 0;

    const float cutoff2 = cutoff * cutoff;

    // Collect contacts in thread-local buffers, then merge
#ifdef _OPENMP
    const int n_threads = omp_get_max_threads();
#else
    const int n_threads = 1;
#endif

    struct Contact { int i, j; float k, r0; };
    std::vector<std::vector<Contact>> thread_contacts(n_threads);

#ifdef _OPENMP
    // Same unchecked-index defect as tencm.cpp: thread_contacts is sized from
    // omp_get_max_threads() outside the region. Clamp the team.
    #pragma omp parallel num_threads(n_threads)
#endif
    {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        auto& local = thread_contacts[tid];

#ifdef _OPENMP
        #pragma omp for schedule(dynamic, 16)
#endif
        for (int i = 0; i < N - 1; ++i) {
            for (int j = i + 2; j < N; ++j) {
                const float dx = ca_xyz[i*3+0] - ca_xyz[j*3+0];
                const float dy = ca_xyz[i*3+1] - ca_xyz[j*3+1];
                const float dz = ca_xyz[i*3+2] - ca_xyz[j*3+2];
                const float r2 = dx*dx + dy*dy + dz*dz;

                if (r2 <= cutoff2) {
                    const float r0_val = std::sqrt(r2);
                    const float ratio  = cutoff / r0_val;
                    const float r3     = ratio * ratio * ratio;
                    const float k_val  = k0 * (r3 * r3);  // k0 * (cutoff/r0)^6

                    local.push_back({i, j, k_val, r0_val});
                }
            }
        }
    }

    // Merge thread-local results
    int total = 0;
    for (const auto& v : thread_contacts) {
        for (const auto& c : v) {
            contacts_ij[total*2+0] = c.i;
            contacts_ij[total*2+1] = c.j;
            contacts_k[total]      = c.k;
            contacts_r0[total]     = c.r0;
            ++total;
        }
    }
    return total;
}

// ─── Hessian assembly ───────────────────────────────────────────────────────

// Torsional Jacobian: J_k(atom_i) = axis_k × (r_i - pivot_k)  if i > k, else 0
static inline void jac(const float* ca_xyz,
                        const float* bond_axis, const float* bond_pivot,
                        int bond_k, int atom_i,
                        float* out)
{
    if (atom_i <= bond_k) {
        out[0] = out[1] = out[2] = 0.0f;
        return;
    }
    const float dx = ca_xyz[atom_i*3+0] - bond_pivot[bond_k*3+0];
    const float dy = ca_xyz[atom_i*3+1] - bond_pivot[bond_k*3+1];
    const float dz = ca_xyz[atom_i*3+2] - bond_pivot[bond_k*3+2];

    const float ax = bond_axis[bond_k*3+0];
    const float ay = bond_axis[bond_k*3+1];
    const float az = bond_axis[bond_k*3+2];

    // cross product: axis × d
    out[0] = ay*dz - az*dy;
    out[1] = az*dx - ax*dz;
    out[2] = ax*dy - ay*dx;
}

void assemble_hessian_cpu(const float* ca_xyz, int N,
                          const int*   contacts_ij,
                          const float* contacts_k,
                          int M, int C,
                          double* H_out)
{
    if (C == 0 || M == 0) return;

    // Build bond axes and pivots
    std::vector<float> bond_axis(M * 3);
    std::vector<float> bond_pivot(M * 3);
    for (int k = 0; k < M; ++k) {
        float ax = ca_xyz[(k+1)*3+0] - ca_xyz[k*3+0];
        float ay = ca_xyz[(k+1)*3+1] - ca_xyz[k*3+1];
        float az = ca_xyz[(k+1)*3+2] - ca_xyz[k*3+2];
        float inv = 1.0f / std::sqrt(ax*ax + ay*ay + az*az + 1e-20f);
        bond_axis[k*3+0] = ax*inv;
        bond_axis[k*3+1] = ay*inv;
        bond_axis[k*3+2] = az*inv;
        bond_pivot[k*3+0] = 0.5f*(ca_xyz[k*3+0] + ca_xyz[(k+1)*3+0]);
        bond_pivot[k*3+1] = 0.5f*(ca_xyz[k*3+1] + ca_xyz[(k+1)*3+1]);
        bond_pivot[k*3+2] = 0.5f*(ca_xyz[k*3+2] + ca_xyz[(k+1)*3+2]);
    }

    // Zero the output Hessian
    std::memset(H_out, 0, static_cast<size_t>(M) * M * sizeof(double));

    // Use Eigen Map for accumulation
    Eigen::Map<Eigen::MatrixXd> H(H_out, M, M);

    // Assemble: one contact at a time, all (k,l) pairs
    // With OpenMP, each thread accumulates into a private Hessian, then reduce
#ifdef _OPENMP
    const int n_threads = omp_get_max_threads();
    std::vector<Eigen::MatrixXd> thread_H(n_threads, Eigen::MatrixXd::Zero(M, M));

    // Same unchecked-index defect as tencm.cpp: thread_H is sized from
    // omp_get_max_threads() outside the region. Clamp the team.
    #pragma omp parallel num_threads(n_threads)
    {
        const int tid = omp_get_thread_num();
        auto& local_H = thread_H[tid];

        #pragma omp for schedule(dynamic, 4)
        for (int ci = 0; ci < C; ++ci) {
            const int ci_i = contacts_ij[ci*2+0];
            const int ci_j = contacts_ij[ci*2+1];
            const float kij = contacts_k[ci];

            for (int k = 0; k < M; ++k) {
                float jki[3], jkj[3];
                jac(ca_xyz, bond_axis.data(), bond_pivot.data(), k, ci_i, jki);
                jac(ca_xyz, bond_axis.data(), bond_pivot.data(), k, ci_j, jkj);
                float djk[3] = { jki[0]-jkj[0], jki[1]-jkj[1], jki[2]-jkj[2] };

                for (int l = k; l < M; ++l) {
                    float jli[3], jlj[3];
                    jac(ca_xyz, bond_axis.data(), bond_pivot.data(), l, ci_i, jli);
                    jac(ca_xyz, bond_axis.data(), bond_pivot.data(), l, ci_j, jlj);
                    float djl[3] = { jli[0]-jlj[0], jli[1]-jlj[1], jli[2]-jlj[2] };

                    double contrib = static_cast<double>(kij) *
                        static_cast<double>(djk[0]*djl[0] + djk[1]*djl[1] + djk[2]*djl[2]);

                    local_H(k, l) += contrib;
                    if (l != k) local_H(l, k) += contrib;
                }
            }
        }
    }

    // Reduce thread-local Hessians
    for (int t = 0; t < n_threads; ++t) {
        H += thread_H[t];
    }
#else
    // Single-threaded path
    for (int ci = 0; ci < C; ++ci) {
        const int ci_i = contacts_ij[ci*2+0];
        const int ci_j = contacts_ij[ci*2+1];
        const float kij = contacts_k[ci];

        for (int k = 0; k < M; ++k) {
            float jki[3], jkj[3];
            jac(ca_xyz, bond_axis.data(), bond_pivot.data(), k, ci_i, jki);
            jac(ca_xyz, bond_axis.data(), bond_pivot.data(), k, ci_j, jkj);
            float djk[3] = { jki[0]-jkj[0], jki[1]-jkj[1], jki[2]-jkj[2] };

            for (int l = k; l < M; ++l) {
                float jli[3], jlj[3];
                jac(ca_xyz, bond_axis.data(), bond_pivot.data(), l, ci_i, jli);
                jac(ca_xyz, bond_axis.data(), bond_pivot.data(), l, ci_j, jlj);
                float djl[3] = { jli[0]-jlj[0], jli[1]-jlj[1], jli[2]-jlj[2] };

                double contrib = static_cast<double>(kij) *
                    static_cast<double>(djk[0]*djl[0] + djk[1]*djl[1] + djk[2]*djl[2]);

                H(k, l) += contrib;
                if (l != k) H(l, k) += contrib;
            }
        }
    }
#endif
}

}}  // namespace tencm::cpu_fallback
