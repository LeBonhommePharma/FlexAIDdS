// CavityDetect.cpp — Task 1.3: FULL AVX-512 + OpenMP SURFNET port
// In-memory, zero-deps, Apache-2.0 © 2026 Le Bonhomme Pharma
// LP's Get_Cleft probe-sphere logic — native, parallel, zero GPL contact

#include "CavityDetect.h"
#include "CavityDetectMetalBridge.h"
#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <numeric>
#include <functional>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#endif

#if defined(FLEXAIDS_USE_AVX512) && defined(__AVX512F__)
#  include <immintrin.h>
#  define USE_AVX512 1
#elif defined(FLEXAIDS_USE_AVX2) && defined(__AVX2__)
#  include <immintrin.h>
#  define USE_AVX2 1
#endif

namespace cavity_detect {

// ─── helpers ─────────────────────────────────────────────────────────────────

float CavityDetector::distance(const float* a, const float* b) const {
    const float dx = a[0] - b[0];
    const float dy = a[1] - b[1];
    const float dz = a[2] - b[2];
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

// Compute the maximum clash-free radius for a probe centred at `mid`.
// Returns the minimum remaining gap between `mid` and any atom surface:
//   r_free = min_over_k ( dist(mid, atoms[k].coor) - atoms[k].radius )
// This is the inner SURFNET reduction: scalar fallback for correctness.
[[maybe_unused]]
static float max_free_radius_scalar(
        const float* mid,
        const std::vector<atom>& atoms,
        const std::size_t* neighbor_indices, std::size_t neighbor_count,
        std::size_t skip_i, std::size_t skip_j,
        float r_current)
{
    float r = r_current;
    for (std::size_t n = 0; n < neighbor_count; ++n) {
        const std::size_t k = neighbor_indices[n];
        if (k == skip_i || k == skip_j) continue;
        const float* ck = atoms[k].coor;
        const float dx = mid[0] - ck[0];
        const float dy = mid[1] - ck[1];
        const float dz = mid[2] - ck[2];
        const float gap = std::sqrt(dx*dx + dy*dy + dz*dz) - atoms[k].radius;
        if (gap < r) r = gap;
    }
    return r;
}

#ifdef USE_AVX512
// AVX-512 inner loop: process 16 neighbor atoms at once.
static float max_free_radius_avx512(
        const float* mid,
        const std::vector<atom>& atoms,
        const std::size_t* neighbor_indices, std::size_t neighbor_count,
        std::size_t skip_i, std::size_t skip_j,
        float r_init)
{
    __m512 vr = _mm512_set1_ps(r_init);
    const __m512 vmx = _mm512_set1_ps(mid[0]);
    const __m512 vmy = _mm512_set1_ps(mid[1]);
    const __m512 vmz = _mm512_set1_ps(mid[2]);

    const __mmask16 full = 0xFFFF;
    alignas(64) float cx[16], cy[16], cz[16], cr[16];

    std::size_t n = 0;
    for (; n + 16 <= neighbor_count; n += 16) {
        __mmask16 mask = full;
        for (int lane = 0; lane < 16; ++lane) {
            std::size_t idx = neighbor_indices[n + static_cast<std::size_t>(lane)];
            cx[lane] = atoms[idx].coor[0];
            cy[lane] = atoms[idx].coor[1];
            cz[lane] = atoms[idx].coor[2];
            cr[lane] = atoms[idx].radius;
            if (idx == skip_i || idx == skip_j)
                mask &= ~(static_cast<__mmask16>(1) << lane);
        }
        __m512 vx = _mm512_load_ps(cx);
        __m512 vy = _mm512_load_ps(cy);
        __m512 vz = _mm512_load_ps(cz);
        __m512 vrad = _mm512_load_ps(cr);

        __m512 ddx = _mm512_sub_ps(vmx, vx);
        __m512 ddy = _mm512_sub_ps(vmy, vy);
        __m512 ddz = _mm512_sub_ps(vmz, vz);
        __m512 dist2 = _mm512_fmadd_ps(ddx, ddx,
                        _mm512_fmadd_ps(ddy, ddy,
                        _mm512_mul_ps(ddz, ddz)));
        __m512 dist  = _mm512_sqrt_ps(dist2);
        __m512 gap   = _mm512_sub_ps(dist, vrad);

        vr = _mm512_mask_min_ps(vr, mask, vr, gap);
    }

    // Reduce 16-lane vector to scalar
    float r = _mm512_reduce_min_ps(vr);

    // Scalar tail for remaining neighbors
    for (; n < neighbor_count; ++n) {
        const std::size_t k = neighbor_indices[n];
        if (k == skip_i || k == skip_j) continue;
        const float* ck = atoms[k].coor;
        const float dx = mid[0] - ck[0];
        const float dy = mid[1] - ck[1];
        const float dz = mid[2] - ck[2];
        const float gap = std::sqrt(dx*dx + dy*dy + dz*dz) - atoms[k].radius;
        if (gap < r) r = gap;
    }
    return r;
}
#endif // USE_AVX512

#ifdef USE_AVX2
// AVX2 inner loop: 8 neighbor atoms at once.
static float max_free_radius_avx2(
        const float* mid,
        const std::vector<atom>& atoms,
        const std::size_t* neighbor_indices, std::size_t neighbor_count,
        std::size_t skip_i, std::size_t skip_j,
        float r_init)
{
    __m256 vr = _mm256_set1_ps(r_init);
    const __m256 vmx = _mm256_set1_ps(mid[0]);
    const __m256 vmy = _mm256_set1_ps(mid[1]);
    const __m256 vmz = _mm256_set1_ps(mid[2]);

    alignas(32) float cx[8], cy[8], cz[8], cr[8];

    std::size_t n = 0;
    for (; n + 8 <= neighbor_count; n += 8) {
        for (int lane = 0; lane < 8; ++lane) {
            std::size_t idx = neighbor_indices[n + static_cast<std::size_t>(lane)];
            cx[lane] = atoms[idx].coor[0];
            cy[lane] = atoms[idx].coor[1];
            cz[lane] = atoms[idx].coor[2];
            cr[lane] = atoms[idx].radius;
        }
        __m256 vx = _mm256_load_ps(cx);
        __m256 vy = _mm256_load_ps(cy);
        __m256 vz = _mm256_load_ps(cz);
        __m256 vrad = _mm256_load_ps(cr);

        __m256 ddx = _mm256_sub_ps(vmx, vx);
        __m256 ddy = _mm256_sub_ps(vmy, vy);
        __m256 ddz = _mm256_sub_ps(vmz, vz);
        __m256 dist2 = _mm256_fmadd_ps(ddx, ddx,
                        _mm256_fmadd_ps(ddy, ddy,
                        _mm256_mul_ps(ddz, ddz)));
        __m256 dist  = _mm256_sqrt_ps(dist2);
        __m256 gap   = _mm256_sub_ps(dist, vrad);

        // Mask skip lanes with large value so they don't affect min
        for (int lane = 0; lane < 8; ++lane) {
            std::size_t idx = neighbor_indices[n + static_cast<std::size_t>(lane)];
            if (idx == skip_i || idx == skip_j) {
                alignas(32) float tmp[8];
                _mm256_store_ps(tmp, gap);
                tmp[lane] = 1e38f;
                gap = _mm256_load_ps(tmp);
            }
        }
        vr = _mm256_min_ps(vr, gap);
    }

    // Reduce 8-lane to scalar
    alignas(32) float lane_vals[8];
    _mm256_store_ps(lane_vals, vr);
    float r = r_init;
    for (int i = 0; i < 8; ++i) r = std::min(r, lane_vals[i]);

    // Scalar tail for remaining neighbors
    for (; n < neighbor_count; ++n) {
        const std::size_t k = neighbor_indices[n];
        if (k == skip_i || k == skip_j) continue;
        const float* ck = atoms[k].coor;
        const float dx = mid[0] - ck[0];
        const float dy = mid[1] - ck[1];
        const float dz = mid[2] - ck[2];
        const float gap = std::sqrt(dx*dx + dy*dy + dz*dz) - atoms[k].radius;
        if (gap < r) r = gap;
    }
    return r;
}
#endif // USE_AVX2

// Dispatch to best available SIMD path.
static inline float max_free_radius(
        const float* mid,
        const std::vector<atom>& atoms,
        const std::size_t* neighbor_indices, std::size_t neighbor_count,
        std::size_t skip_i, std::size_t skip_j,
        float r_init)
{
#ifdef USE_AVX512
    return max_free_radius_avx512(mid, atoms, neighbor_indices, neighbor_count, skip_i, skip_j, r_init);
#elif defined(USE_AVX2)
    return max_free_radius_avx2(mid, atoms, neighbor_indices, neighbor_count, skip_i, skip_j, r_init);
#else
    return max_free_radius_scalar(mid, atoms, neighbor_indices, neighbor_count, skip_i, skip_j, r_init);
#endif
}


// ─── data loading ─────────────────────────────────────────────────────────────

void CavityDetector::load_from_fa(const atom* atoms, const resid* residues, int res_cnt) {
    m_atoms.clear();
    m_clefts.clear();

    // Walk each residue's atom range; residues[r].fatm[0]..latm[0] give the
    // 0-based atom indices into the `atoms` array (FlexAID convention).
    for (int r = 0; r < res_cnt; ++r) {
        if (!residues[r].fatm || !residues[r].latm) continue;
        const int first = residues[r].fatm[0];
        const int last  = residues[r].latm[0];
        for (int a = first; a <= last; ++a) {
            // Shallow copy — only coor and radius are used by the detector.
            // Pointer fields are explicitly nulled to prevent use-after-free.
            atom copy{};
            copy.coor[0]  = atoms[a].coor[0];
            copy.coor[1]  = atoms[a].coor[1];
            copy.coor[2]  = atoms[a].coor[2];
            copy.radius   = atoms[a].radius;
            copy.number   = atoms[a].number;
            copy.ofres    = atoms[a].ofres;
            copy.coor_ref = nullptr;
            copy.par      = nullptr;
            copy.cons     = nullptr;
            copy.optres   = nullptr;
            copy.eigen    = nullptr;
            std::memcpy(copy.name, atoms[a].name, sizeof(copy.name));
            m_atoms.push_back(copy);
        }
    }
}

void CavityDetector::load_from_pdb(const std::string& pdb_file) {
    m_atoms.clear();
    m_clefts.clear();

    std::ifstream in(pdb_file);
    if (!in) return;

    int atom_idx = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (line.size() < 54) continue;
        const std::string rec = line.substr(0, 6);
        if (rec != "ATOM  " && rec != "HETATM") continue;

        atom a{};
        a.number = std::stoi(line.substr(6, 5));
        std::string name_str = line.substr(12, 4);
        // Left-trim name
        std::size_t ns = name_str.find_first_not_of(' ');
        if (ns != std::string::npos) name_str = name_str.substr(ns);
        std::strncpy(a.name, name_str.c_str(), 4);
        a.name[4] = '\0';

        a.coor[0] = std::stof(line.substr(30, 8));
        a.coor[1] = std::stof(line.substr(38, 8));
        a.coor[2] = std::stof(line.substr(46, 8));

        // Van der Waals radius from B-factor column if present; default 1.7 Å
        if (line.size() >= 66) {
            try { a.radius = std::stof(line.substr(60, 6)); }
            catch (...) { a.radius = 1.7f; }
        } else {
            a.radius = 1.7f;
        }
        if (a.radius <= 0.0f) a.radius = 1.7f;

        a.ofres  = atom_idx;
        a.coor_ref = nullptr; a.par = nullptr;
        a.cons = nullptr; a.optres = nullptr; a.eigen = nullptr;
        m_atoms.push_back(a);
        ++atom_idx;
    }
}


// ─── probe-sphere placement (SURFNET core) ───────────────────────────────────

void CavityDetector::detect(float min_radius, float max_radius) {
    m_sphere_lwb = min_radius;
    m_sphere_upb = max_radius;
    m_clefts.clear();

    const std::size_t N = m_atoms.size();
    if (N < 2) return;

    // Build spatial grid for O(K) neighbor lookups instead of O(N)
    m_grid.build(m_atoms);

    // ── Metal primary path (Apple Silicon, FLEXAIDS_USE_METAL) ────────────
#ifdef FLEXAIDS_USE_METAL
    {
        std::vector<MetalAtom> gpu_atoms(N);
        for (std::size_t i = 0; i < N; ++i) {
            gpu_atoms[i].pos[0] = m_atoms[i].coor[0];
            gpu_atoms[i].pos[1] = m_atoms[i].coor[1];
            gpu_atoms[i].pos[2] = m_atoms[i].coor[2];
            gpu_atoms[i].radius = m_atoms[i].radius;
        }
        std::vector<MetalSphereResult> metal_spheres;
        if (cavity_detect_metal_dispatch(
                gpu_atoms.data(), static_cast<int>(N),
                min_radius, max_radius, metal_spheres))
        {
            // Convert Metal results to DetectedSphere and cluster on CPU
            std::vector<DetectedSphere> all_spheres;
            all_spheres.reserve(metal_spheres.size());
            for (const auto& ms : metal_spheres) {
                DetectedSphere s;
                s.center[0] = ms.center[0];
                s.center[1] = ms.center[1];
                s.center[2] = ms.center[2];
                s.radius    = ms.radius;
                s.cleft_id  = 0;
                all_spheres.push_back(s);
            }
            _cluster_and_finalize(std::move(all_spheres));
            return;
        }
        // Metal unavailable or failed — fall through to CPU path
    }
#endif

    // ── CPU path: AVX-512 / AVX2 / scalar + OpenMP ───────────────────────
    // Collect accepted probe spheres from all threads.
    // Each thread writes to its own local vector; we merge afterwards to avoid
    // false sharing and to keep the critical section out of the hot loop.
    const int n_threads =
#ifdef _OPENMP
        omp_get_max_threads();
#else
        1;
#endif
    std::vector<std::vector<DetectedSphere>> thread_spheres(
            static_cast<std::size_t>(n_threads));

    // Per-thread scratch buffers for grid neighbor queries (zero-allocation hot path)
    const std::size_t max_neighbors = 512;
    std::vector<std::vector<std::size_t>> thread_scratch(
            static_cast<std::size_t>(n_threads),
            std::vector<std::size_t>(max_neighbors));

    // NOTE: MSVC's OpenMP 2.0 requires a *signed* integer loop variable.
    // std::size_t (unsigned) is rejected; use int with a cast inside.
    const int N_signed = static_cast<int>(N);
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
#endif
    for (int _i = 0; _i < N_signed; ++_i) {
        const std::size_t i = static_cast<std::size_t>(_i);
        const int tid =
#ifdef _OPENMP
            omp_get_thread_num();
#else
            0;
#endif
        auto& scratch = thread_scratch[static_cast<std::size_t>(tid)];

        for (std::size_t j = i + 1; j < N; ++j) {
            // Midpoint between atoms i and j
            const float mid[3] = {
                (m_atoms[i].coor[0] + m_atoms[j].coor[0]) * 0.5f,
                (m_atoms[i].coor[1] + m_atoms[j].coor[1]) * 0.5f,
                (m_atoms[i].coor[2] + m_atoms[j].coor[2]) * 0.5f,
            };

            // Initial probe radius: half the gap between the two atom surfaces
            const float d_ij   = distance(m_atoms[i].coor, m_atoms[j].coor);
            float r = d_ij * 0.5f
                    - 0.5f * (m_atoms[i].radius + m_atoms[j].radius);

            if (r < m_sphere_lwb) continue;
            if (r > m_sphere_upb) r = m_sphere_upb;

            // Query spatial grid for nearby atoms, then shrink r (SURFNET clash test)
            std::size_t n_near = m_grid.query_neighbors(mid, scratch.data(), max_neighbors);
            r = max_free_radius(mid, m_atoms, scratch.data(), n_near, i, j, r);

            if (r < m_sphere_lwb) continue;
            if (r > m_sphere_upb) r = m_sphere_upb;

            DetectedSphere s;
            s.center[0] = mid[0];
            s.center[1] = mid[1];
            s.center[2] = mid[2];
            s.radius    = r;
            s.cleft_id  = 0;  // assigned by merge_clefts()
            thread_spheres[static_cast<std::size_t>(tid)].push_back(s);
        }
    }

    // Merge thread-local vectors
    std::vector<DetectedSphere> all_spheres;
    for (auto& v : thread_spheres)
        all_spheres.insert(all_spheres.end(), v.begin(), v.end());

    if (all_spheres.empty()) return;
    _cluster_and_finalize(std::move(all_spheres));
}

void CavityDetector::_cluster_and_finalize(std::vector<DetectedSphere> all_spheres) {
    // Cluster accepted spheres into clefts by overlap graph (union-find).
    // Two spheres belong to the same cleft when their surfaces overlap:
    //   dist(c_i, c_j) < r_i + r_j
    const std::size_t S = all_spheres.size();
    std::vector<int> parent(S);
    std::iota(parent.begin(), parent.end(), 0);

    std::function<int(int)> find = [&](int x) -> int {
        while (parent[static_cast<std::size_t>(x)] != x) {
            parent[static_cast<std::size_t>(x)] =
                parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(x)])];
            x = parent[static_cast<std::size_t>(x)];
        }
        return x;
    };
    auto unite = [&](int a, int b) {
        a = find(a); b = find(b);
        if (a != b) parent[static_cast<std::size_t>(a)] = b;
    };

    for (std::size_t si = 0; si < S; ++si) {
        for (std::size_t sj = si + 1; sj < S; ++sj) {
            const float dx = all_spheres[si].center[0] - all_spheres[sj].center[0];
            const float dy = all_spheres[si].center[1] - all_spheres[sj].center[1];
            const float dz = all_spheres[si].center[2] - all_spheres[sj].center[2];
            const float d  = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (d < all_spheres[si].radius + all_spheres[sj].radius)
                unite(static_cast<int>(si), static_cast<int>(sj));
        }
    }

    // Build cleft list from connected components
    std::unordered_map<int, int> root_to_cleft;
    int next_id = 1;
    for (std::size_t si = 0; si < S; ++si) {
        const int root = find(static_cast<int>(si));
        auto it = root_to_cleft.find(root);
        int cid;
        if (it == root_to_cleft.end()) {
            cid = next_id++;
            root_to_cleft[root] = cid;
            DetectedCleft cleft{};
            cleft.id    = cid;
            cleft.label = cid;
            m_clefts.push_back(std::move(cleft));
        } else {
            cid = it->second;
        }
        all_spheres[si].cleft_id = cid;
        m_clefts[static_cast<std::size_t>(cid - 1)].spheres.push_back(
            all_spheres[si]);
    }

    // Compute per-cleft center of mass and approximate volume
    for (auto& cleft : m_clefts) {
        float cx = 0, cy = 0, cz = 0;
        for (const auto& s : cleft.spheres) {
            cx += s.center[0]; cy += s.center[1]; cz += s.center[2];
        }
        const float inv = 1.0f / static_cast<float>(cleft.spheres.size());
        cleft.center[0] = cx * inv;
        cleft.center[1] = cy * inv;
        cleft.center[2] = cz * inv;
        float vol = 0.0f;
        for (const auto& s : cleft.spheres)
            vol += (4.0f / 3.0f) * 3.14159265f * s.radius * s.radius * s.radius;
        cleft.volume = vol;
        cleft.effrad = std::cbrt(vol * (3.0f / (4.0f * 3.14159265f)));
    }

    sort_clefts();
    assign_atoms_to_clefts();
}


// ─── post-processing ──────────────────────────────────────────────────────────

void CavityDetector::merge_clefts() {
    // Iteratively merge pairs of clefts whose centers are within
    // (effrad_i + effrad_j) of each other (i.e. overlapping effective volumes).
    bool merged = true;
    while (merged) {
        merged = false;
        for (std::size_t i = 0; i < m_clefts.size() && !merged; ++i) {
            for (std::size_t j = i + 1; j < m_clefts.size() && !merged; ++j) {
                const float dx = m_clefts[i].center[0] - m_clefts[j].center[0];
                const float dy = m_clefts[i].center[1] - m_clefts[j].center[1];
                const float dz = m_clefts[i].center[2] - m_clefts[j].center[2];
                const float d  = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (d < m_clefts[i].effrad + m_clefts[j].effrad) {
                    // Absorb j into i
                    for (auto& s : m_clefts[j].spheres) {
                        s.cleft_id = m_clefts[i].id;
                        m_clefts[i].spheres.push_back(s);
                    }
                    m_clefts[i].volume += m_clefts[j].volume;
                    m_clefts[i].effrad =
                        std::cbrt(m_clefts[i].volume * (3.0f / (4.0f * 3.14159265f)));
                    // Recompute center
                    float cx = 0, cy = 0, cz = 0;
                    for (const auto& s : m_clefts[i].spheres) {
                        cx += s.center[0]; cy += s.center[1]; cz += s.center[2];
                    }
                    const float inv = 1.0f / static_cast<float>(m_clefts[i].spheres.size());
                    m_clefts[i].center[0] = cx * inv;
                    m_clefts[i].center[1] = cy * inv;
                    m_clefts[i].center[2] = cz * inv;
                    m_clefts.erase(m_clefts.begin() + static_cast<std::ptrdiff_t>(j));
                    merged = true;
                }
            }
        }
    }
    sort_clefts();
}

void CavityDetector::sort_clefts() {
    std::sort(m_clefts.begin(), m_clefts.end(),
        [](const DetectedCleft& a, const DetectedCleft& b) {
            return a.spheres.size() > b.spheres.size();
        });
    // Re-assign sequential IDs after sort
    for (int i = 0; i < static_cast<int>(m_clefts.size()); ++i)
        m_clefts[static_cast<std::size_t>(i)].id = i + 1;
}

void CavityDetector::assign_atoms_to_clefts(float contact_threshold) {
    // Assign each protein atom to the cleft (if any) whose nearest sphere
    // surface is within contact_threshold Å.
    for (auto& a : m_atoms) {
        float best_dist = 1e38f;  // Use large value instead of infinity for -ffast-math safety
        int   best_cid  = 0;
        for (const auto& cleft : m_clefts) {
            for (const auto& s : cleft.spheres) {
                const float dx = a.coor[0] - s.center[0];
                const float dy = a.coor[1] - s.center[1];
                const float dz = a.coor[2] - s.center[2];
                const float d  = std::sqrt(dx*dx + dy*dy + dz*dz) - s.radius;
                if (d < best_dist) { best_dist = d; best_cid = cleft.id; }
            }
        }
        // Store cleft assignment in ofres field (repurposed as cavity tag)
        if (best_dist <= contact_threshold) a.ofres = best_cid;
    }
}

void CavityDetector::filter_anchor_residues(const std::string& anchor_residues) {
    // Parse comma-separated residue numbers (e.g. "45,72,103").
    // Keep only clefts that contain at least one sphere within 5 Å of any
    // anchor residue's Cα atom.
    if (anchor_residues.empty()) return;

    std::vector<int> anchor_nums;
    std::istringstream ss(anchor_residues);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        tok.erase(std::remove_if(tok.begin(), tok.end(), ::isspace), tok.end());
        if (!tok.empty()) {
            try { anchor_nums.push_back(std::stoi(tok)); }
            catch (...) {}
        }
    }
    if (anchor_nums.empty()) return;

    // Collect Cα coordinates for anchor residues
    std::vector<std::array<float, 3>> anchor_pos;
    for (const auto& a : m_atoms) {
        for (int rn : anchor_nums) {
            if (a.ofres == rn && std::strncmp(a.name, "CA", 2) == 0) {
                // MSVC fix: explicit std::array<float,3> — bare braced-init-list
                // is ambiguous to MSVC when deducing the push_back overload.
                anchor_pos.push_back(std::array<float, 3>{a.coor[0], a.coor[1], a.coor[2]});
            }
        }
    }
    if (anchor_pos.empty()) return;

    m_clefts.erase(std::remove_if(m_clefts.begin(), m_clefts.end(),
        [&](const DetectedCleft& cleft) {
            for (const auto& s : cleft.spheres) {
                for (const auto& ap : anchor_pos) {
                    const float dx = s.center[0] - ap[0];
                    const float dy = s.center[1] - ap[1];
                    const float dz = s.center[2] - ap[2];
                    if (std::sqrt(dx*dx + dy*dy + dz*dz) < 5.0f)
                        return false;  // keep
                }
            }
            return true;  // no anchor nearby → remove
        }),
        m_clefts.end());

    sort_clefts();
}


// ─── FlexAID compatibility ────────────────────────────────────────────────────

sphere* CavityDetector::to_flexaid_spheres(int cleft_id) const {
    // Find the requested cleft
    const DetectedCleft* target = nullptr;
    for (const auto& c : m_clefts) {
        if (c.id == cleft_id) { target = &c; break; }
    }
    if (!target || target->spheres.empty()) return nullptr;

    // Build a singly-linked list (prev-chained, as FlexAID expects)
    // Ownership is transferred to the caller (free with delete / FlexAID cleanup).
    sphere* head = nullptr;
    for (const auto& s : target->spheres) {
        sphere* node  = new sphere{};
        node->center[0] = s.center[0];
        node->center[1] = s.center[1];
        node->center[2] = s.center[2];
        node->radius    = s.radius;
        node->prev      = head;
        head            = node;
    }
    return head;  // caller passes to generate_grid()
}


// ─── output ───────────────────────────────────────────────────────────────────

void CavityDetector::write_sphere_pdb(const std::string& filename, int cleft_id) const {
    std::ofstream out(filename);
    if (!out) return;

    out << "REMARK  Cleft spheres — FlexAIDΔS CavityDetector (Task 1.3)\n";
    int atom_id = 1;
    for (const auto& cleft : m_clefts) {
        if (cleft_id > 0 && cleft.id != cleft_id) continue;
        for (const auto& s : cleft.spheres) {
            out << "HETATM"
                << std::setw(5) << atom_id++
                << "  SPH SURF "
                << std::setw(4) << cleft.id
                << "    "
                << std::fixed << std::setprecision(3)
                << std::setw(8) << s.center[0]
                << std::setw(8) << s.center[1]
                << std::setw(8) << s.center[2]
                << "  1.00"
                << std::setw(6) << std::setprecision(2) << s.radius
                << "           S\n";
        }
    }
    out << "END\n";
}

} // namespace cavity_detect
