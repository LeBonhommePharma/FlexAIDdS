// new_search_arch.h — residual Phase-4 search architecture (S4 options A+B)
//
// Default: all OFF (no behavior change). Product gates via env only.
//
// A) FLEXAIDDS_PHENOTYPE_UNIQUE=1
//    Classic bit-flip mutation often only flips dead low bits (phenotype
//    to_ic unchanged). When enabled, after a classic mutate pass, if the
//    phenotype bin signature is unchanged, apply a ±1-bin phenotype-live
//    step so uniqueness / search actually moves IC space.
//
// B) FLEXAIDDS_BASIN_REINJECT=1
//    On diversity collapse, prefer reinjecting worst individuals that land
//    outside a Cartesian ligand RMSD basin around the current best, instead
//    of pure random gene reinit keyed only on allele entropy.
//
// Convenience: FLEXAIDDS_NEW_SEARCH_ARCH=1 enables both A and B.
//             FLEXAIDDS_NEW_SEARCH_ARCH=phenotype_unique,basin_reinject
//             (comma/space/semicolon-separated tokens).
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>

namespace flexaids {
namespace new_search {

inline bool env_truthy(const char* e) {
    if (!e || e[0] == '\0') return false;
    if (e[0] == '0' || e[0] == 'n' || e[0] == 'N' || e[0] == 'f' || e[0] == 'F')
        return false;
    return true;
}

/// Case-insensitive token search in comma/space/semicolon-separated list.
inline bool arch_token_present(const char* list, const char* token) {
    if (!list || !token || !*token) return false;
    // Exact "1" enables all tokens.
    if (list[0] == '1' && list[1] == '\0') return true;
    const size_t tlen = std::strlen(token);
    const char* p = list;
    while (*p) {
        while (*p == ' ' || *p == ',' || *p == ';' || *p == '|') ++p;
        if (!*p) break;
        const char* start = p;
        while (*p && *p != ' ' && *p != ',' && *p != ';' && *p != '|') ++p;
        const size_t n = static_cast<size_t>(p - start);
        if (n == tlen) {
            bool match = true;
            for (size_t i = 0; i < n; ++i) {
                char a = start[i], b = token[i];
                if (a >= 'A' && a <= 'Z') a = static_cast<char>(a - 'A' + 'a');
                if (b >= 'A' && b <= 'Z') b = static_cast<char>(b - 'A' + 'a');
                if (a != b) {
                    match = false;
                    break;
                }
            }
            if (match) return true;
        }
    }
    return false;
}

inline bool phenotype_unique_enabled() {
    if (env_truthy(std::getenv("FLEXAIDDS_PHENOTYPE_UNIQUE"))) return true;
    return arch_token_present(std::getenv("FLEXAIDDS_NEW_SEARCH_ARCH"),
                              "phenotype_unique");
}

inline bool basin_reinject_enabled() {
    if (env_truthy(std::getenv("FLEXAIDDS_BASIN_REINJECT"))) return true;
    return arch_token_present(std::getenv("FLEXAIDDS_NEW_SEARCH_ARCH"),
                              "basin_reinject");
}

/// Bin index for one gene's IC phenotype (stable for uniqueness).
inline int phenotype_bin_index(double to_ic, double gmin, double gmax, double nbin) {
    double range = gmax - gmin;
    if (range <= 0.0) range = 1.0;
    int nb = static_cast<int>(std::lround(nbin > 1.0 ? nbin : 2.0));
    if (nb < 2) nb = 2;
    double u = (to_ic - gmin) / range;
    if (u < 0.0) u = 0.0;
    if (u > 1.0) u = 1.0;
    int b = static_cast<int>(std::floor(u * static_cast<double>(nb)));
    if (b >= nb) b = nb - 1;
    if (b < 0) b = 0;
    return b;
}

/// Hash phenotype bins only (no dead gene-bit dependence).
template <typename GeneT, typename GenlimT>
inline std::size_t hash_phenotype_bins(const GeneT* genes, int n_genes,
                                       const GenlimT* lim) {
    std::size_t h = 0;
    for (int i = 0; i < n_genes; ++i) {
        const int b = phenotype_bin_index(genes[i].to_ic, lim[i].min, lim[i].max,
                                          lim[i].nbin);
        h ^= std::hash<int>{}(b) + 0x9e3779b9 + (h << 6) + (h >> 2);
    }
    return h;
}

/// Integer gene step for ±k phenotype bins (matches G4.3 granular step size).
inline int32_t phenotype_bin_step_int(double nbin, int max_random_value = 0x7fffffff) {
    const double nb = nbin > 1.0 ? nbin : 2.0;
    const double step = std::floor((static_cast<double>(max_random_value) + 1.0) / nb);
    const int32_t s = static_cast<int32_t>(std::max(1.0, step));
    return s;
}

/// Apply ±1 (or ±2) bin step to gene j; clamps to [0, MAX]. Returns true if
/// to_int32 changed. Does not recompute to_ic (caller should genetoic).
template <typename GeneT>
inline bool apply_phenotype_bin_step(GeneT* gene, double nbin, int sign, int k_bins,
                                     int32_t max_random_value = 0x7fffffff) {
    if (!gene || sign == 0) return false;
    if (k_bins < 1) k_bins = 1;
    const int32_t step = phenotype_bin_step_int(nbin, max_random_value);
    const int64_t delta =
        static_cast<int64_t>(sign) * static_cast<int64_t>(step) * static_cast<int64_t>(k_bins);
    int64_t ng = static_cast<int64_t>(gene->to_int32) + delta;
    if (ng < 0) ng = 0;
    if (ng > static_cast<int64_t>(max_random_value))
        ng = static_cast<int64_t>(max_random_value);
    const int32_t before = gene->to_int32;
    gene->to_int32 = static_cast<int32_t>(ng);
    return gene->to_int32 != before;
}

/// True if Cartesian RMSD is outside the basin (strictly greater than sigma).
inline bool outside_basin(double rmsd_ang, double sigma_ang) {
    return rmsd_ang > sigma_ang;
}

inline double basin_sigma_ang(double default_ang = 2.0) {
    if (const char* s = std::getenv("FLEXAIDDS_BASIN_SIGMA_ANG")) {
        const double v = std::atof(s);
        if (v > 0.0) return v;
    }
    return default_ang;
}

}  // namespace new_search
}  // namespace flexaids
