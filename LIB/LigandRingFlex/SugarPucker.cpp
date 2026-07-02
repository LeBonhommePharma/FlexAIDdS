// SugarPucker.cpp — Cremer-Pople pseudorotation implementation
// Eigen is used for 3D vector arithmetic in apply_sugar_puckers().
// AVX-512 is used for batch phase-to-torsion computation in the population path.
#define _USE_MATH_DEFINES
#include <cmath>

#include "SugarPucker.h"
#include "../RngSeed.h"

#include <Eigen/Dense>
#ifdef __AVX512F__
#  include <immintrin.h>
#endif
#include <cstring>
#include <random>
#include <algorithm>
#include <numeric>

// Bring in FA atom struct for atom name access
#include "../flexaid.h"

namespace sugar_pucker {

static constexpr float PI_F = 3.14159265f;
static constexpr float DEG2RAD = PI_F / 180.0f;
static constexpr float RAD2DEG = 180.0f / PI_F;

// ─── detect_sugar_type ───────────────────────────────────────────────────────
SugarType detect_sugar_type(const atom* atoms,
                             const int*  ring_atom_indices,
                             int         ring_size)
{
    if (!atoms || !ring_atom_indices) return SugarType::Unknown;

    for (int k = 0; k < ring_size; ++k) {
        int idx = ring_atom_indices[k];
        // Check for O2' atom name (ribose indicator)
        const char* name = atoms[idx].name;
        if (name && (strstr(name, "O2'") || strstr(name, "O2*") ||
                     strstr(name, "O2P") || strcmp(name, "O2'") == 0))
            return SugarType::Ribose;
    }
    return SugarType::Deoxyribose;
}

// ─── compute_ring_torsions ───────────────────────────────────────────────────
// Altona-Sundaralingam (1972) equations for furanose pseudorotation:
//   νk = nu_max * cos(P + (4π/5)*(k-2))  for k = 0..4
void compute_ring_torsions(const PuckerParams& p, float torsions[5]) noexcept {
    const float P_rad      = p.P * DEG2RAD;
    const float delta_rad  = 4.0f * PI_F / 5.0f; // 144°

    for (int k = 0; k < 5; ++k) {
        torsions[k] = p.nu_max *
                      std::cos(P_rad + delta_rad * static_cast<float>(k - 2));
    }
}

// ─── compute_pucker_energy ───────────────────────────────────────────────────
// Cosine-based energy landscape:
//   Ribose:      E(P) = A * [1 - cos(P - 18°)] * [1 - cos(P - 162°)]  kcal/mol
//   Deoxyribose: E(P) = A * [1 - 0.8*cos(P - 162°)] with secondary min at 18°
double compute_pucker_energy(float phase_deg, SugarType stype) noexcept {
    const double p = static_cast<double>(phase_deg) * M_PI / 180.0;

    if (stype == SugarType::Ribose) {
        // Two minima: C3'-endo (18°) and C2'-endo (162°)
        double A = 1.5; // kcal/mol barrier
        return A * (1.0 - std::cos(p - 18.0 * M_PI / 180.0))
                 * (1.0 - std::cos(p - 162.0 * M_PI / 180.0)) * 0.25;
    } else {
        // Single dominant minimum at C2'-endo (162°), small shoulder at 18°
        double A = 1.2;
        double B = 0.4;
        return A * (1.0 - std::cos(p - 162.0 * M_PI / 180.0))
             - B * std::cos(p - 18.0 * M_PI / 180.0);
    }
}

// ─── apply_sugar_puckers ─────────────────────────────────────────────────────
void apply_sugar_puckers(
    atom*                                atoms,
    const std::vector<std::vector<int>>& ring_indices,
    const std::vector<float>&            phases_deg,
    const std::vector<SugarType>&        sugar_types)
{
    if (!atoms) return;
    size_t n = std::min({ring_indices.size(), phases_deg.size(), sugar_types.size()});

#ifdef __AVX512F__
    // AVX-512 batch processing disabled due to SVML dependency and type issues.
    // Falls through to Eigen-optimized or scalar path below.
#endif

    // Use Eigen for batched cosine evaluation across all rings
    for (size_t i = 0; i < n; ++i) {
        const auto& ring = ring_indices[i];
        if (ring.size() != 5) continue;

        float P_rad = phases_deg[i] * DEG2RAD;
        const float nu_max = 38.0f;
        const float delta  = 4.0f * PI_F / 5.0f;

        // Build Eigen array of offsets k-2 for k=0..4
        Eigen::Array<float,5,1> k_offsets;
        k_offsets << -2.0f, -1.0f, 0.0f, 1.0f, 2.0f;
        Eigen::Array<float,5,1> angles = P_rad + delta * k_offsets;
        Eigen::Array<float,5,1> torsions = nu_max * angles.cos() * RAD2DEG;

        for (int k = 0; k < 5; ++k)
            atoms[ring[k]].dih = torsions(k);
    }
}

// ─── mutate_phase ────────────────────────────────────────────────────────────
float mutate_phase(float current_phase_deg, float sigma_deg) {
    auto& rng = flexaids_rng::lazy_thread_rng(0x5A6A9ULL);
    std::normal_distribution<float> dist(0.0f, sigma_deg);
    float new_phase = current_phase_deg + dist(rng);
    // Wrap to [0, 360)
    while (new_phase < 0.0f)   new_phase += 360.0f;
    while (new_phase >= 360.0f) new_phase -= 360.0f;
    return new_phase;
}

// ─── compute_pucker_entropy ──────────────────────────────────────────────────
double compute_pucker_entropy(const std::vector<float>& phase_ensemble_deg) {
    if (phase_ensemble_deg.empty()) return 0.0;

    constexpr int BINS = 36; // 10° resolution
    std::vector<int> counts(BINS, 0);
    for (float p : phase_ensemble_deg) {
        float pp = p;
        while (pp < 0.0f)    pp += 360.0f;
        while (pp >= 360.0f) pp -= 360.0f;
        int b = std::min(static_cast<int>(pp / 10.0f), BINS - 1);
        counts[b]++;
    }

    int total = static_cast<int>(phase_ensemble_deg.size());
    double H = 0.0;
    const double log2_inv = 1.0 / std::log(2.0);
    for (int c : counts) {
        if (c > 0) {
            double prob = static_cast<double>(c) / total;
            H -= prob * std::log(prob) * log2_inv;
        }
    }
    return H;
}

} // namespace sugar_pucker
