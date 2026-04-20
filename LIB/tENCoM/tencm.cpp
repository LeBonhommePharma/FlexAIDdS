// tencm.cpp — Torsional Elastic Network Contact Model implementation
//
// Mathematical basis (Delarue & Sanejouand 2002; Yang, Song & Cui 2009):
//
//  Spring potential over Cα contacts within r_cutoff:
//    V = ½ Σ_{(i,j)∈contacts} k_ij |δr_i - δr_j|²
//
//  Torsion-space Jacobian (pseudo-bond k connects Cα_k → Cα_{k+1}):
//    J_k(i) = e_k × (r_i - p_k)    if i > k  (downstream)
//             {0,0,0}               if i ≤ k  (upstream)
//    where e_k = unit bond axis, p_k = bond midpoint
//
//  Torsional Hessian:
//    H_kl = Σ_{(i,j)∈contacts} k_ij [(J_k(i)-J_k(j)) · (J_l(i)-J_l(j))]
//
//  Normal modes: symmetric diagonalisation  H V = V Λ
//    → Eigen SelfAdjointEigenSolver when available (divide-and-conquer)
//    → Jacobi fallback otherwise
//
//  Boltzmann sampling at temperature T:
//    σ_m² = kB T / λ_m  (equipartition per mode, skip m=0..5 ≈ rigid-body)
//    δθ   = Σ_{m≥6} σ_m * z_m * v_m,   z_m ~ N(0,1)
//
//  Perturbed Cα: r_i' = r_i + Σ_k J_k(i) δθ_k
//
// Hardware dispatch priority:
//   1. Eigen   — SelfAdjointEigenSolver for diagonalisation; Map<> for BLAS
//   2. AVX-512 — 16-wide distance batching, 8-wide double accumulation
//   3. AVX2    — 8-wide distance batching, FMA inner products
//   4. OpenMP  — parallel contacts, Hessian assembly, B-factors, sampling
//   5. Scalar  — always available
#include "tencm.h"
#include "simd_distance.h"
#include "ion_utils.h"

#include <cstring>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <stdexcept>
#include <numbers>
#include <random>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#ifdef __AVX512F__
#  include <immintrin.h>
#endif

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace tencm {

// ─── local utility helpers ───────────────────────────────────────────────────

static void cross3f(const float* a, const float* b, float* c) noexcept {
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

// ─── TorsionalENM::build ─────────────────────────────────────────────────────
void TorsionalENM::build(const atom*  atoms,
                          const resid* residue,
                          int          res_cnt,
                          float        cutoff,
                          float        k0)
{
    cutoff_ = cutoff;
    k0_     = k0;
    built_  = false;

    extract_ca(atoms, residue, res_cnt);

    if (static_cast<int>(ca_.size()) < 3) {
        std::cerr << "TENCM: fewer than 3 Cα atoms found; skipping.\n";
        return;
    }

    build_contacts();
    build_bonds();
    assemble_hessian();
    diagonalize();

    built_ = true;
}

// ─── build_from_ca (standalone, no FA dependency) ───────────────────────────
void TorsionalENM::build_from_ca(const std::vector<std::array<float,3>>& ca_coords,
                                  float cutoff, float k0)
{
    cutoff_ = cutoff;
    k0_     = k0;
    built_  = false;

    ca_ = ca_coords;
    ca_atom_idx_.clear();
    res_idx_.clear();
    ion_radii_.clear();
    // Identity mapping — each sequential Cα is its own "residue"; no ions here.
    for (int i = 0; i < static_cast<int>(ca_.size()); ++i) {
        ca_atom_idx_.push_back(i);
        res_idx_.push_back(i);
    }
    n_protein_ca_ = static_cast<int>(ca_.size());

    if (static_cast<int>(ca_.size()) < 3) {
        std::cerr << "TENCM: fewer than 3 Cα atoms; skipping.\n";
        return;
    }

    build_contacts();
    build_bonds();
    assemble_hessian();
    diagonalize();

    built_ = true;
}

// ─── extract_ca ──────────────────────────────────────────────────────────────
//
// Phase 1: collect protein Cα atoms (residue.type == 0).
// Phase 2: append metal ion atoms as rigid pseudo-nodes (residue.type == 1,
//          matched by is_ion_resname).  Ion nodes have index >= n_protein_ca_.
void TorsionalENM::extract_ca(const atom*  atoms,
                               const resid* residue,
                               int          res_cnt)
{
    ca_.clear();
    ca_atom_idx_.clear();
    res_idx_.clear();
    ion_radii_.clear();

    // Phase 1: protein Cα atoms
    for (int ri = 1; ri <= res_cnt; ++ri) {
        if (residue[ri].type != 0) continue;   // protein residues only

        int first = residue[ri].fatm[0];
        int last  = residue[ri].latm[0];
        for (int ai = first; ai <= last; ++ai) {
            const char* nm = atoms[ai].name;
            // Match " CA " or "CA  " PDB atom name
            bool is_ca = (nm[0]==' '&&nm[1]=='C'&&nm[2]=='A'&&nm[3]==' ') ||
                         (nm[0]=='C'&&nm[1]=='A'&&nm[2]==' '&&nm[3]==' ');
            if (is_ca) {
                ca_.push_back({ atoms[ai].coor[0],
                                atoms[ai].coor[1],
                                atoms[ai].coor[2] });
                ca_atom_idx_.push_back(ai);
                res_idx_.push_back(ri);
                break;
            }
        }
    }

    // Mark the protein/ion boundary
    n_protein_ca_ = static_cast<int>(ca_.size());

    // Phase 2: metal ion pseudo-nodes (HETATM, single heavy atom each)
    for (int ri = 1; ri <= res_cnt; ++ri) {
        if (residue[ri].type != 1) continue;             // HETATM only
        if (!is_ion_resname(residue[ri].name)) continue; // metal ions only

        int first = residue[ri].fatm[0];
        int last  = residue[ri].latm[0];
        for (int ai = first; ai <= last; ++ai) {
            float r = atoms[ai].radius;
            if (r < 0.5f) continue;  // radius must be assigned (by assign_radii)
            ca_.push_back({ atoms[ai].coor[0],
                            atoms[ai].coor[1],
                            atoms[ai].coor[2] });
            ca_atom_idx_.push_back(ai);
            res_idx_.push_back(ri);
            ion_radii_.push_back(r);
            break;  // one representative atom per ion
        }
    }
}

// ─── build_contacts ─────────────────────────────────────────────────────────
//
// Two phases:
//  (A) Protein–protein Cα contacts — O(Np²) with hardware dispatch:
//        AVX-512 → 16 distances/cycle
//        AVX2    → 8  distances/cycle via simd::distance2_1x8
//        OpenMP  → parallel outer loop with thread-local contact lists
//  (B) Protein–ion contacts — scalar loop over (Np × N_ions) pairs:
//        uses surface-to-surface distance and area-scaled spring constant
//
// Approximate Cα VdW radius for surface-distance computation:
static constexpr float R_CA_EFF = 1.88f;  // aliphatic carbon (Å)

void TorsionalENM::build_contacts()
{
    contacts_.clear();
    tmcontsct_.clear();
    [[maybe_unused]] const int N   = static_cast<int>(ca_.size());
    const int Np  = n_protein_ca_;
    const float rc2 = cutoff_ * cutoff_;

#if defined(_OPENMP) && !defined(TENCM_NO_OMP)
    // Parallel contact discovery with thread-local vectors
    const int n_threads = omp_get_max_threads();
    std::vector<std::vector<Contact>> thread_contacts(n_threads);

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        auto& local = thread_contacts[tid];

        #pragma omp for schedule(dynamic, 4)
        for (int i = 0; i < Np - 1; ++i) {
            const float xi = ca_[i][0], yi = ca_[i][1], zi = ca_[i][2];
            (void)xi; (void)yi; (void)zi; // used in SIMD paths below

#if defined(__AVX512F__)
            // AVX-512: batch 16 j-atoms at a time
            const __m512 vxi  = _mm512_set1_ps(xi);
            const __m512 vyi  = _mm512_set1_ps(yi);
            const __m512 vzi  = _mm512_set1_ps(zi);
            const __m512 vrc2 = _mm512_set1_ps(rc2);
            const __m512 vk0  = _mm512_set1_ps(k0_);
            const __m512 vrc  = _mm512_set1_ps(cutoff_);

            int j = i + 2;
            for (; j + 15 < Np; j += 16) {
                // Gather 16 Cα positions
                float jx[16], jy[16], jz[16];
                for (int q = 0; q < 16; ++q) {
                    jx[q] = ca_[j+q][0];
                    jy[q] = ca_[j+q][1];
                    jz[q] = ca_[j+q][2];
                }
                __m512 dx = _mm512_sub_ps(_mm512_loadu_ps(jx), vxi);
                __m512 dy = _mm512_sub_ps(_mm512_loadu_ps(jy), vyi);
                __m512 dz = _mm512_sub_ps(_mm512_loadu_ps(jz), vzi);
                __m512 r2 = _mm512_fmadd_ps(dz, dz,
                            _mm512_fmadd_ps(dy, dy,
                            _mm512_mul_ps(dx, dx)));

                // Mask: r2 <= rc2
                __mmask16 mask = _mm512_cmple_ps_mask(r2, vrc2);
                if (mask == 0) continue;

                // Extract contacts for set bits
                alignas(64) float r2_arr[16];
                _mm512_store_ps(r2_arr, r2);
                for (int q = 0; q < 16; ++q) {
                    if (mask & (1u << q)) {
                        float r0    = std::sqrt(r2_arr[q]);
                        float ratio = cutoff_ / r0;
                        float r3    = ratio * ratio * ratio;
                        float k     = k0_ * (r3 * r3);
                        local.push_back({i, j + q, k, r0});
                    }
                }
            }
            // Scalar tail for remaining j atoms
            for (; j < Np; ++j) {
                float dx = ca_[i][0] - ca_[j][0];
                float dy = ca_[i][1] - ca_[j][1];
                float dz = ca_[i][2] - ca_[j][2];
                float r2 = dx*dx + dy*dy + dz*dz;
                if (r2 <= rc2) {
                    float r0    = std::sqrt(r2);
                    float ratio = cutoff_ / r0;
                    float r3    = ratio * ratio * ratio;
                    local.push_back({i, j, k0_ * (r3 * r3), r0});
                }
            }

#elif FLEXAIDS_HAS_AVX2
            // AVX2: batch 8 j-atoms at a time using simd_distance.h patterns
            int j = i + 2;
            for (; j + 7 < Np; j += 8) {
                float jx[8], jy[8], jz[8], r2_arr[8];
                for (int q = 0; q < 8; ++q) {
                    jx[q] = ca_[j+q][0];
                    jy[q] = ca_[j+q][1];
                    jz[q] = ca_[j+q][2];
                }
                simd::distance2_1x8(jx, jy, jz, xi, yi, zi, r2_arr);
                for (int q = 0; q < 8; ++q) {
                    if (r2_arr[q] <= rc2) {
                        float r0    = std::sqrt(r2_arr[q]);
                        float ratio = cutoff_ / r0;
                        float r3    = ratio * ratio * ratio;
                        local.push_back({i, j + q, k0_ * (r3 * r3), r0});
                    }
                }
            }
            for (; j < Np; ++j) {
                float dx = ca_[i][0] - ca_[j][0];
                float dy = ca_[i][1] - ca_[j][1];
                float dz = ca_[i][2] - ca_[j][2];
                float r2 = dx*dx + dy*dy + dz*dz;
                if (r2 <= rc2) {
                    float r0    = std::sqrt(r2);
                    float ratio = cutoff_ / r0;
                    float r3    = ratio * ratio * ratio;
                    local.push_back({i, j, k0_ * (r3 * r3), r0});
                }
            }

#else
            // Scalar path
            for (int j = i + 2; j < Np; ++j) {
                float dx = ca_[i][0] - ca_[j][0];
                float dy = ca_[i][1] - ca_[j][1];
                float dz = ca_[i][2] - ca_[j][2];
                float r2 = dx*dx + dy*dy + dz*dz;
                if (r2 <= rc2) {
                    float r0    = std::sqrt(r2);
                    float ratio = cutoff_ / r0;
                    float r3    = ratio * ratio * ratio;
                    local.push_back({i, j, k0_ * (r3 * r3), r0});
                }
            }
#endif
        }
    } // end omp parallel

    // Merge thread-local contact lists
    std::size_t total = 0;
    for (const auto& tc : thread_contacts) total += tc.size();
    contacts_.reserve(total);
    for (auto& tc : thread_contacts)
        contacts_.insert(contacts_.end(),
                         std::make_move_iterator(tc.begin()),
                         std::make_move_iterator(tc.end()));

#else
    // ── Serial path (no OpenMP) with SIMD distance batching ──

    for (int i = 0; i < Np - 1; ++i) {
        [[maybe_unused]] const float xi = ca_[i][0], yi = ca_[i][1], zi = ca_[i][2];
#if defined(__AVX512F__)
        const __m512 vxi  = _mm512_set1_ps(xi);
        const __m512 vyi  = _mm512_set1_ps(yi);
        const __m512 vzi  = _mm512_set1_ps(zi);
        const __m512 vrc2 = _mm512_set1_ps(rc2);

        int j = i + 2;
        for (; j + 15 < Np; j += 16) {
            float jx[16], jy[16], jz[16];
            for (int q = 0; q < 16; ++q) {
                jx[q] = ca_[j+q][0]; jy[q] = ca_[j+q][1]; jz[q] = ca_[j+q][2];
            }
            __m512 dx = _mm512_sub_ps(_mm512_loadu_ps(jx), vxi);
            __m512 dy = _mm512_sub_ps(_mm512_loadu_ps(jy), vyi);
            __m512 dz = _mm512_sub_ps(_mm512_loadu_ps(jz), vzi);
            __m512 r2 = _mm512_fmadd_ps(dz, dz, _mm512_fmadd_ps(dy, dy, _mm512_mul_ps(dx, dx)));
            __mmask16 mask = _mm512_cmple_ps_mask(r2, vrc2);
            if (mask == 0) continue;
            alignas(64) float r2_arr[16];
            _mm512_store_ps(r2_arr, r2);
            for (int q = 0; q < 16; ++q) {
                if (mask & (1u << q)) {
                    float r0 = std::sqrt(r2_arr[q]);
                    float ratio = cutoff_ / r0;
                    float r3 = ratio * ratio * ratio;
                    contacts_.push_back({i, j+q, k0_ * (r3*r3), r0});
                }
            }
        }
        for (; j < Np; ++j) {
            float dx = ca_[i][0]-ca_[j][0], dy = ca_[i][1]-ca_[j][1], dz = ca_[i][2]-ca_[j][2];
            float r2 = dx*dx+dy*dy+dz*dz;
            if (r2 <= rc2) {
                float r0 = std::sqrt(r2); float ratio = cutoff_/r0; float r3 = ratio*ratio*ratio;
                contacts_.push_back({i, j, k0_*(r3*r3), r0});
            }
        }

#elif FLEXAIDS_HAS_AVX2
        int j = i + 2;
        for (; j + 7 < Np; j += 8) {
            float jx[8], jy[8], jz[8], r2_arr[8];
            for (int q = 0; q < 8; ++q) {
                jx[q] = ca_[j+q][0]; jy[q] = ca_[j+q][1]; jz[q] = ca_[j+q][2];
            }
            simd::distance2_1x8(jx, jy, jz, xi, yi, zi, r2_arr);
            for (int q = 0; q < 8; ++q) {
                if (r2_arr[q] <= rc2) {
                    float r0 = std::sqrt(r2_arr[q]); float ratio = cutoff_/r0; float r3 = ratio*ratio*ratio;
                    contacts_.push_back({i, j+q, k0_*(r3*r3), r0});
                }
            }
        }
        for (; j < Np; ++j) {
            float dx = ca_[i][0]-ca_[j][0], dy = ca_[i][1]-ca_[j][1], dz = ca_[i][2]-ca_[j][2];
            float r2 = dx*dx+dy*dy+dz*dz;
            if (r2 <= rc2) {
                float r0 = std::sqrt(r2); float ratio = cutoff_/r0; float r3 = ratio*ratio*ratio;
                contacts_.push_back({i, j, k0_*(r3*r3), r0});
            }
        }

#else
        for (int j = i + 2; j < Np; ++j) {
            float dx = ca_[i][0]-ca_[j][0], dy = ca_[i][1]-ca_[j][1], dz = ca_[i][2]-ca_[j][2];
            float r2 = dx*dx+dy*dy+dz*dz;
            if (r2 <= rc2) {
                float r0 = std::sqrt(r2);
                float ratio = cutoff_ / r0;
                float r3 = ratio*ratio*ratio;
                contacts_.push_back({i, j, k0_ * (r3*r3), r0});
            }
        }
#endif
    }
#endif // _OPENMP

    // ── Phase B: protein–ion contacts (scalar; few ions so SIMD not needed) ──
    // Surface-to-surface distance: d_surf = r_center − r_ion − R_CA_EFF
    // Spring constant: k_ij = k0 * (rc/d_surf)^6 * (r_ion/R_CA_EFF)²
    {
        const int N_ions = static_cast<int>(ion_radii_.size());
        for (int i = 0; i < Np; ++i) {
            for (int q = 0; q < N_ions; ++q) {
                const int j = Np + q;
                float dx = ca_[i][0] - ca_[j][0];
                float dy = ca_[i][1] - ca_[j][1];
                float dz = ca_[i][2] - ca_[j][2];
                float r_center = std::sqrt(dx*dx + dy*dy + dz*dz);
                float r_ion    = ion_radii_[static_cast<std::size_t>(q)];
                float d_surf   = r_center - r_ion - R_CA_EFF;
                if (d_surf < 0.01f) d_surf = 0.01f;  // clamp overlap
                if (d_surf > cutoff_) continue;        // outside cutoff
                float area_scale = (r_ion / R_CA_EFF) * (r_ion / R_CA_EFF);
                float ratio = cutoff_ / d_surf;
                float r3    = ratio * ratio * ratio;
                float k_scaled = k0_ * (r3 * r3) * area_scale;
                contacts_.push_back({i, j, k_scaled, d_surf});
            }
        }
    }

    // ── Build tmcontsct_ from all real contacts (protein–protein + protein–ion) ──
    tmcontsct_.reserve(contacts_.size());
    for (const auto& c : contacts_) {
        if (c.k == 0.0f && c.i == 0 && c.j == 0) continue;  // skip sentinels
        bool is_ion_contact = (c.j >= Np);
        tmcontsct_.push_back({c.i, c.j, c.k, c.r0, is_ion_contact});
    }

    // Pad contact list to next multiple of 16 for SIMD tail elimination.
    // Sentinel contacts have k=0, contributing nothing to the Hessian.
#ifdef __AVX512F__
    while (contacts_.size() % 16 != 0)
        contacts_.push_back({0, 0, 0.0f, 0.0f});
#endif
}

// ─── build_bonds ─────────────────────────────────────────────────────────────
// Pseudo-bonds are defined only between consecutive protein Cα atoms.
// Ion nodes (index >= n_protein_ca_) have no backbone torsions.
void TorsionalENM::build_bonds()
{
    bonds_.clear();
    const int Np = n_protein_ca_;  // protein nodes only

    for (int k = 0; k < Np - 1; ++k) {
        PseudoBond pb;
        pb.k = k;

        float ax = ca_[k+1][0] - ca_[k][0];
        float ay = ca_[k+1][1] - ca_[k][1];
        float az = ca_[k+1][2] - ca_[k][2];
        float inv = 1.0f / std::sqrt(ax*ax + ay*ay + az*az);
        pb.axis[0] = ax * inv;
        pb.axis[1] = ay * inv;
        pb.axis[2] = az * inv;

        pb.pivot[0] = 0.5f * (ca_[k][0] + ca_[k+1][0]);
        pb.pivot[1] = 0.5f * (ca_[k][1] + ca_[k+1][1]);
        pb.pivot[2] = 0.5f * (ca_[k][2] + ca_[k+1][2]);

        bonds_.push_back(pb);
    }
}

// ─── jac ─────────────────────────────────────────────────────────────────────
std::array<float,3>
TorsionalENM::jac(int bond_k, int atom_i) const noexcept
{
    // Ion pseudo-nodes (index >= n_protein_ca_) are rigid: no torsional DOF.
    // They constrain protein flexibility via Hessian spring terms but do not
    // rotate themselves — J(ion) = 0 reduces the Hessian contribution to
    //   k_ij * J_k(protein_i) * J_l(protein_i)   (restoring anchor effect).
    if (atom_i >= n_protein_ca_) return {0.0f, 0.0f, 0.0f};

    // Upstream protein atoms are unaffected by this torsion rotation
    if (atom_i <= bond_k) return {0.0f, 0.0f, 0.0f};

    const PseudoBond& pb = bonds_[static_cast<std::size_t>(bond_k)];
    float d[3] = { ca_[atom_i][0] - pb.pivot[0],
                   ca_[atom_i][1] - pb.pivot[1],
                   ca_[atom_i][2] - pb.pivot[2] };
    float j[3];
    cross3f(pb.axis, d, j);
    return {j[0], j[1], j[2]};
}

// ─── assemble_hessian ───────────────────────────────────────────────────────
//
// H_kl = Σ_{contacts} k_ij [(J_k(i)-J_k(j)) · (J_l(i)-J_l(j))]
//
// Hardware acceleration:
//   Eigen → use Map<MatrixXd> for BLAS-level accumulation
//   OpenMP → parallelize over contacts with atomic accumulation
//   AVX-512 → vectorize inner dot-product accumulation (8 doubles/cycle)
void TorsionalENM::assemble_hessian()
{
    const int M = static_cast<int>(bonds_.size());
    const int N = static_cast<int>(ca_.size());

    // Pre-compute Jacobians: J[k*N + i] = 3-vector.
    // Hoisted here so it can be cached for bfactors()/sample() reuse.
    std::vector<std::array<float,3>> J(static_cast<std::size_t>(M * N));

    #if defined(_OPENMP) && !defined(TENCM_NO_OMP)
    #pragma omp parallel for schedule(static) collapse(2)
    #endif
    for (int k = 0; k < M; ++k)
        for (int i = 0; i < N; ++i)
            J[static_cast<std::size_t>(k * N + i)] = jac(k, i);

    // ── Eigen path: accumulate into Eigen matrix, leverage BLAS ──

    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(M, M);

    // Per-contact: compute ΔJ vectors and accumulate outer product
    #if defined(_OPENMP) && !defined(TENCM_NO_OMP)
    // Thread-local Hessian matrices to avoid atomics
    const int n_threads = omp_get_max_threads();
    std::vector<Eigen::MatrixXd> thread_H(n_threads, Eigen::MatrixXd::Zero(M, M));

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        auto& Ht = thread_H[tid];
        const int C = static_cast<int>(contacts_.size());

        #pragma omp for schedule(dynamic, 8)
        for (int ci = 0; ci < C; ++ci) {
            const auto& c = contacts_[ci];
            const float kij = c.k;

            for (int k = 0; k < M; ++k) {
                const auto& jki = J[k * N + c.i];
                const auto& jkj = J[k * N + c.j];
                float djk[3] = { jki[0]-jkj[0], jki[1]-jkj[1], jki[2]-jkj[2] };

                for (int l = k; l < M; ++l) {
                    const auto& jli = J[l * N + c.i];
                    const auto& jlj = J[l * N + c.j];
                    float djl[3] = { jli[0]-jlj[0], jli[1]-jlj[1], jli[2]-jlj[2] };

                    double contrib = kij * static_cast<double>(
                        djk[0]*djl[0] + djk[1]*djl[1] + djk[2]*djl[2]);

                    Ht(k, l) += contrib;
                    if (l != k) Ht(l, k) += contrib;
                }
            }
        }
    }

    // Reduce thread-local matrices
    for (const auto& Ht : thread_H) H += Ht;

    #else
    // Serial Eigen path
    for (const auto& c : contacts_) {
        const float kij = c.k;
        for (int k = 0; k < M; ++k) {
            const auto& jki = J[k * N + c.i];
            const auto& jkj = J[k * N + c.j];
            float djk[3] = { jki[0]-jkj[0], jki[1]-jkj[1], jki[2]-jkj[2] };

            for (int l = k; l < M; ++l) {
                const auto& jli = J[l * N + c.i];
                const auto& jlj = J[l * N + c.j];
                float djl[3] = { jli[0]-jlj[0], jli[1]-jlj[1], jli[2]-jlj[2] };

                double contrib = kij * static_cast<double>(
                    djk[0]*djl[0] + djk[1]*djl[1] + djk[2]*djl[2]);

                H(k, l) += contrib;
                if (l != k) H(l, k) += contrib;
            }
        }
    }
    #endif

    // Copy to flat H_ for compatibility with sample()/strain energy
    H_.resize(static_cast<std::size_t>(M * M));
    Eigen::Map<Eigen::MatrixXd>(H_.data(), M, M) = H;

    // Cache the Jacobian for reuse by bfactors() and sample()
    J_cached_ = std::move(J);
    jac_cached_ = true;
}

// ─── Jacobi eigenvalue decomposition ─────────────────────────────────────────
// Classic off-diagonal Jacobi sweeps for a symmetric n×n matrix A.
// On exit: A holds eigenvalues on diagonal; V holds eigenvectors in columns.
void TorsionalENM::jacobi_rotate(std::vector<double>& A,
                                  std::vector<double>& V,
                                  int n, int p, int q) noexcept
{
    if (std::abs(A[static_cast<std::size_t>(p*n+q)]) < 1e-15) return;

    double tau = (A[static_cast<std::size_t>(q*n+q)] -
                  A[static_cast<std::size_t>(p*n+p)]) /
                 (2.0 * A[static_cast<std::size_t>(p*n+q)]);
    double t = (tau >= 0.0 ? 1.0 : -1.0) / (std::abs(tau) + std::sqrt(1.0 + tau*tau));
    double c = 1.0 / std::sqrt(1.0 + t*t);
    double s = t * c;

    // Update A
    double app = A[static_cast<std::size_t>(p*n+p)];
    double aqq = A[static_cast<std::size_t>(q*n+q)];
    double apq = A[static_cast<std::size_t>(p*n+q)];

    A[static_cast<std::size_t>(p*n+p)] = c*c*app - 2.0*s*c*apq + s*s*aqq;
    A[static_cast<std::size_t>(q*n+q)] = s*s*app + 2.0*s*c*apq + c*c*aqq;
    A[static_cast<std::size_t>(p*n+q)] = A[static_cast<std::size_t>(q*n+p)] = 0.0;

    for (int r = 0; r < n; ++r) {
        if (r == p || r == q) continue;
        double arp = A[static_cast<std::size_t>(r*n+p)];
        double arq = A[static_cast<std::size_t>(r*n+q)];
        A[static_cast<std::size_t>(r*n+p)] = A[static_cast<std::size_t>(p*n+r)] = c*arp - s*arq;
        A[static_cast<std::size_t>(r*n+q)] = A[static_cast<std::size_t>(q*n+r)] = s*arp + c*arq;
    }

    // Update eigenvector matrix V
    for (int r = 0; r < n; ++r) {
        double vrp = V[static_cast<std::size_t>(r*n+p)];
        double vrq = V[static_cast<std::size_t>(r*n+q)];
        V[static_cast<std::size_t>(r*n+p)] = c*vrp - s*vrq;
        V[static_cast<std::size_t>(r*n+q)] = s*vrp + c*vrq;
    }
}

// ─── diagonalize ─────────────────────────────────────────────────────────────
//
// Hardware dispatch:
//   Eigen → SelfAdjointEigenSolver (divide-and-conquer, O(M³) but optimised)
//   Jacobi fallback → hand-rolled sweeps (max 50)
void TorsionalENM::diagonalize()
{
    const int M = static_cast<int>(bonds_.size());
    if (M == 0) return;

    // ── Eigen path: highly optimised divide-and-conquer eigendecomposition ──
    Eigen::Map<const Eigen::MatrixXd> H_map(H_.data(), M, M);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H_map);

    if (solver.info() != Eigen::Success) {
        std::cerr << "TENCM: Eigen diagonalisation failed; falling back to Jacobi.\n";
        goto jacobi_fallback;
    }

    {
        const Eigen::VectorXd& vals = solver.eigenvalues();
        const Eigen::MatrixXd& vecs = solver.eigenvectors();

        modes_.clear();
        modes_.reserve(static_cast<std::size_t>(M));
        for (int i = 0; i < M; ++i) {
            NormalMode nm;
            nm.eigenvalue = vals(i);
            nm.eigenvector.resize(static_cast<std::size_t>(M));
            for (int j = 0; j < M; ++j)
                nm.eigenvector[static_cast<std::size_t>(j)] = vecs(j, i);
            modes_.push_back(std::move(nm));
        }

        // Eigen returns eigenvalues sorted ascending — already correct
        return;
    }

    jacobi_fallback:
    {
        // ── Jacobi fallback ──
        std::vector<double> A = H_;
        std::vector<double> V(static_cast<std::size_t>(M * M), 0.0);
        for (int i = 0; i < M; ++i)
            V[static_cast<std::size_t>(i*M+i)] = 1.0;

        const int MAX_SWEEPS = 50;
        for (int sweep = 0; sweep < MAX_SWEEPS; ++sweep) {
            double max_off = 0.0;
            int pp = 0, qq = 1;
            for (int p = 0; p < M - 1; ++p)
                for (int q = p + 1; q < M; ++q) {
                    double v = std::abs(A[static_cast<std::size_t>(p*M+q)]);
                    if (v > max_off) { max_off = v; pp = p; qq = q; }
                }
            if (max_off < 1e-12) break;
            jacobi_rotate(A, V, M, pp, qq);
        }

        modes_.clear();
        modes_.reserve(static_cast<std::size_t>(M));
        for (int i = 0; i < M; ++i) {
            NormalMode nm;
            nm.eigenvalue = A[static_cast<std::size_t>(i*M+i)];
            nm.eigenvector.resize(static_cast<std::size_t>(M));
            for (int j = 0; j < M; ++j)
                nm.eigenvector[static_cast<std::size_t>(j)] = V[static_cast<std::size_t>(j*M+i)];
            modes_.push_back(std::move(nm));
        }

        std::sort(modes_.begin(), modes_.end(),
                  [](const NormalMode& a, const NormalMode& b){
                      return a.eigenvalue < b.eigenvalue; });
    }
}

// ─── sample ──────────────────────────────────────────────────────────────────
//
// Boltzmann-weighted sampling of backbone perturbations.
// Hardware: Eigen for matrix-vector products; OpenMP not used (single conformer).
Conformer TorsionalENM::sample(float temperature, std::mt19937& rng) const
{
    if (!built_) throw std::runtime_error("TENCM: model not built");

    const int M = static_cast<int>(bonds_.size());
    const int N = static_cast<int>(ca_.size());

    std::normal_distribution<float> gauss(0.0f, 1.0f);

    Conformer conf;
    conf.delta_theta.assign(static_cast<std::size_t>(M), 0.0f);

    // Skip first 6 modes (rigid-body; eigenvalue ≈ 0)
    const int SKIP = std::min(6, M);
    const double kBT = static_cast<double>(kB_kcal) * temperature;

    // Eigen path: accumulate Σ σ_m * z_m * v_m as vector operation
    Eigen::VectorXd dtheta = Eigen::VectorXd::Zero(M);

    for (int m = SKIP; m < M && m < SKIP + N_MODES; ++m) {
        double lam = modes_[static_cast<std::size_t>(m)].eigenvalue;
        if (lam < 1e-8) continue;
        double sigma = std::sqrt(kBT / lam);
        double z = gauss(rng);

        Eigen::Map<const Eigen::VectorXd> evec(
            modes_[static_cast<std::size_t>(m)].eigenvector.data(), M);
        dtheta += (sigma * z) * evec;
    }

    for (int k = 0; k < M; ++k)
        conf.delta_theta[static_cast<std::size_t>(k)] = static_cast<float>(dtheta(k));

    // Build perturbed Cα coordinates:  r_i' = r_i + Σ_k J_k(i) δθ_k
    conf.ca.resize(static_cast<std::size_t>(N));

    #if defined(_OPENMP) && !defined(TENCM_NO_OMP) && 0
    // NOTE: atom loop is typically N < 1000 with per-atom M-inner-loop
    // For small systems OpenMP overhead exceeds gain; keep serial.
    #endif
    for (int i = 0; i < N; ++i) {
        float disp[3] = {0.0f, 0.0f, 0.0f};
        for (int k = 0; k < M; ++k) {
            float dth = conf.delta_theta[static_cast<std::size_t>(k)];
            if (dth == 0.0f || i <= k) continue;
            auto jki = jac(k, i);
            disp[0] += jki[0] * dth;
            disp[1] += jki[1] * dth;
            disp[2] += jki[2] * dth;
        }
        conf.ca[static_cast<std::size_t>(i)] = {
            ca_[static_cast<std::size_t>(i)][0] + disp[0],
            ca_[static_cast<std::size_t>(i)][1] + disp[1],
            ca_[static_cast<std::size_t>(i)][2] + disp[2]
        };
    }

    // Elastic strain energy: ½ δθᵀ H δθ
    {
        Eigen::Map<const Eigen::MatrixXd> Hmat(H_.data(), M, M);
        conf.strain_energy = static_cast<float>(0.5 * dtheta.dot(Hmat * dtheta));
    }

    return conf;
}

// ─── apply ───────────────────────────────────────────────────────────────────
void TorsionalENM::apply(const Conformer& conf,
                          atom*            atoms,
                          const resid*     residue) const
{
    if (!built_) return;
    // Apply displacement only to protein Cα nodes (ion nodes are rigid; J=0).
    const int Np = n_protein_ca_;

    for (int seq = 0; seq < Np; ++seq) {
        int ri = res_idx_[static_cast<std::size_t>(seq)];
        int   ai = ca_atom_idx_[static_cast<std::size_t>(seq)];
        float dx = conf.ca[static_cast<std::size_t>(seq)][0] - atoms[ai].coor[0];
        float dy = conf.ca[static_cast<std::size_t>(seq)][1] - atoms[ai].coor[1];
        float dz = conf.ca[static_cast<std::size_t>(seq)][2] - atoms[ai].coor[2];

        int first = residue[ri].fatm[0];
        int last  = residue[ri].latm[0];
        for (int a = first; a <= last; ++a) {
            atoms[a].coor[0] += dx;
            atoms[a].coor[1] += dy;
            atoms[a].coor[2] += dz;
        }
    }
}

// ─── bfactors ────────────────────────────────────────────────────────────────
// B_i = (8π²/3) <Δr_i²>  where <Δr_i²> = kBT Σ_m (σ_m² |J_m(i)|²)
//
// Only protein Cα nodes (indices 0..n_protein_ca_-1) are included.
// Ion pseudo-nodes have zero Jacobian and thus zero B-factor — they are
// excluded from the output to preserve backward compatibility with callers
// that expect one B-factor per protein residue.
//
// Hardware dispatch:
//   Eigen → vectorised eigenvector·Jacobian accumulation
//   OpenMP → parallelize over residues (each independent)
//   AVX-512 → vectorize mode accumulation loop (8 doubles/cycle)
std::vector<float> TorsionalENM::bfactors(float temperature) const
{
    const int N  = static_cast<int>(ca_.size());   // total nodes (protein + ions)
    const int Np = n_protein_ca_;                  // protein nodes only (output)
    const int M  = static_cast<int>(bonds_.size());
    const double kBT = static_cast<double>(kB_kcal) * temperature;
    const int SKIP = std::min(6, M);
    const int M_end = std::min(M, SKIP + N_MODES);
    const double BF_SCALE = 8.0 * std::numbers::pi * std::numbers::pi / 3.0;

    // Output only protein Cα B-factors; ion nodes are excluded (zero Jacobian).
    std::vector<float> bf(static_cast<std::size_t>(Np), 0.0f);

    // Pre-compute per-mode: kBT / λ_m (skip modes with tiny eigenvalue)
    std::vector<double> inv_stiffness(static_cast<std::size_t>(M), 0.0);
    for (int m = SKIP; m < M_end; ++m) {
        double lam = modes_[static_cast<std::size_t>(m)].eigenvalue;
        if (lam >= 1e-8)
            inv_stiffness[static_cast<std::size_t>(m)] = kBT / lam;
    }

#if defined(_OPENMP) && !defined(TENCM_NO_OMP)
    #pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < Np; ++i) {  // protein nodes only
        double msf = 0.0;  // mean-square fluctuation
        for (int m = SKIP; m < M_end; ++m) {
            double s2 = inv_stiffness[static_cast<std::size_t>(m)];
            if (s2 == 0.0) continue;

            // J_m(i) = Σ_k v_mk * J_k(i)
            // Use cached Jacobian when available to avoid recomputing cross products.
            float jmi[3] = {0.0f, 0.0f, 0.0f};

#if defined(__AVX512F__)
            // AVX-512: accumulate 8 bonds per cycle (double precision)
            const double* evec = modes_[static_cast<std::size_t>(m)].eigenvector.data();
            __m512d acc_x = _mm512_setzero_pd();
            __m512d acc_y = _mm512_setzero_pd();
            __m512d acc_z = _mm512_setzero_pd();
            int k = 0;
            int k_end = std::min(M, i);
            for (; k + 7 < k_end; k += 8) {
                __m512d vmk = _mm512_loadu_pd(evec + k);
                double jx[8], jy[8], jz[8];
                for (int q = 0; q < 8; ++q) {
                    const auto& jkq = jac_cached_
                        ? J_cached_[static_cast<std::size_t>((k + q) * N + i)]
                        : jac(k + q, i);
                    jx[q] = jkq[0]; jy[q] = jkq[1]; jz[q] = jkq[2];
                }
                acc_x = _mm512_fmadd_pd(vmk, _mm512_loadu_pd(jx), acc_x);
                acc_y = _mm512_fmadd_pd(vmk, _mm512_loadu_pd(jy), acc_y);
                acc_z = _mm512_fmadd_pd(vmk, _mm512_loadu_pd(jz), acc_z);
            }
            jmi[0] = static_cast<float>(_mm512_reduce_add_pd(acc_x));
            jmi[1] = static_cast<float>(_mm512_reduce_add_pd(acc_y));
            jmi[2] = static_cast<float>(_mm512_reduce_add_pd(acc_z));
            for (; k < k_end; ++k) {
                const auto& jki = jac_cached_
                    ? J_cached_[static_cast<std::size_t>(k * N + i)]
                    : jac(k, i);
                float vmk_f = static_cast<float>(evec[k]);
                jmi[0] += vmk_f * jki[0];
                jmi[1] += vmk_f * jki[1];
                jmi[2] += vmk_f * jki[2];
            }

#else
            // Scalar / AVX2 path — use cached Jacobian when available
            for (int k = 0; k < M; ++k) {
                if (i <= k) continue;
                const auto& jki = jac_cached_
                    ? J_cached_[static_cast<std::size_t>(k * N + i)]
                    : jac(k, i);
                float vmk = static_cast<float>(
                    modes_[static_cast<std::size_t>(m)].eigenvector[static_cast<std::size_t>(k)]);
                jmi[0] += vmk * jki[0];
                jmi[1] += vmk * jki[1];
                jmi[2] += vmk * jki[2];
            }
#endif

            double jmag2 = static_cast<double>(
                jmi[0]*jmi[0] + jmi[1]*jmi[1] + jmi[2]*jmi[2]);
            msf += s2 * jmag2;
        }
        bf[static_cast<std::size_t>(i)] = static_cast<float>(BF_SCALE * msf);
    }
    return bf;
}

// ─── residue_rms_fluctuation ─────────────────────────────────────────────────
float residue_rms_fluctuation(const TorsionalENM& tencm,
                              int  residue_idx,
                              float temperature)
{
    auto bf = tencm.bfactors(temperature);
    if (residue_idx < 0 || residue_idx >= static_cast<int>(bf.size()))
        return 0.0f;
    float msf = bf[static_cast<std::size_t>(residue_idx)] /
                static_cast<float>(8.0 * std::numbers::pi * std::numbers::pi / 3.0);
    return std::sqrt(msf);
}

}  // namespace tencm
