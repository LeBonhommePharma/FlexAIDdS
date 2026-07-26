// G4.2 niche pair distance — pure, header-only, unit-testable.
//
// Default (FLEXAIDDS_NICHE_CARTESIAN unset/0): gene-space RMSP over to_ic values
// (historical FlexAID calc_rmsp formula). Gene 0 is a cleft-grid ORDINAL mixed
// with angular genes — structural defect documented in PHASE4_GATES_ACTUALIZED.
//
// Cartesian (FLEXAIDDS_NICHE_CARTESIAN=1): RMSD over ligand heavy-atom XYZ in
// Angstroms. Does NOT use gene vectors; no ordinal gene0 contribution.
//
// Product path: gaboom.cpp calculate_fitness PSHARE/SMFREE call these helpers.
#pragma once

#include <cmath>
#include <cstddef>
#include <cstdlib>

namespace flexaids {

/// Gene-space RMSP: sqrt(mean((ic_a[i]-ic_b[i])^2)) — same formula as calc_rmsp
/// over to_ic. Pure; no FA/grid dependency. Used when Cartesian niche is OFF.
inline double niche_gene_rmsp(const double* ic_a, const double* ic_b, int npar) {
    if (npar <= 0 || ic_a == nullptr || ic_b == nullptr) {
        return 0.0;
    }
    double sum_sq = 0.0;
    for (int i = 0; i < npar; ++i) {
        const double d = ic_a[i] - ic_b[i];
        sum_sq += d * d;
    }
    return std::sqrt(sum_sq / static_cast<double>(npar));
}

/// Cartesian ligand heavy-atom RMSD (Å). xyz_a/xyz_b are interleaved x,y,z
/// (length 3*n_atoms). No gene / ordinal input.
inline double niche_cartesian_rmsd(const float* xyz_a, const float* xyz_b,
                                   int n_atoms) {
    if (n_atoms <= 0 || xyz_a == nullptr || xyz_b == nullptr) {
        return 0.0;
    }
    double sum_sq = 0.0;
    for (int t = 0; t < n_atoms; ++t) {
        const double dx = static_cast<double>(xyz_a[t * 3 + 0] - xyz_b[t * 3 + 0]);
        const double dy = static_cast<double>(xyz_a[t * 3 + 1] - xyz_b[t * 3 + 1]);
        const double dz = static_cast<double>(xyz_a[t * 3 + 2] - xyz_b[t * 3 + 2]);
        sum_sq += dx * dx + dy * dy + dz * dz;
    }
    return std::sqrt(sum_sq / static_cast<double>(n_atoms));
}

/// True iff FLEXAIDDS_NICHE_CARTESIAN is a non-zero integer (product env gate).
inline bool niche_cartesian_env_enabled() {
    const char* e = std::getenv("FLEXAIDDS_NICHE_CARTESIAN");
    return e != nullptr && e[0] != '\0' && std::atoi(e) != 0;
}

/// Default sigma (Å) when Cartesian niche is on; override via FLEXAIDDS_NICHE_SIGMA_ANG.
inline double niche_cartesian_sigma_ang(double default_ang = 2.0) {
    if (const char* s = std::getenv("FLEXAIDDS_NICHE_SIGMA_ANG")) {
        const double v = std::atof(s);
        if (v > 0.0) {
            return v;
        }
    }
    return default_ang;
}

/// Pair distance for niching: Cartesian RMSD if use_cartesian, else gene RMSP.
/// When use_cartesian, ic_* are ignored; when not, xyz_* are ignored.
inline double niche_pair_distance(bool use_cartesian,
                                  const double* ic_a, const double* ic_b, int npar,
                                  const float* xyz_a, const float* xyz_b, int n_atoms) {
    if (use_cartesian) {
        return niche_cartesian_rmsd(xyz_a, xyz_b, n_atoms);
    }
    return niche_gene_rmsp(ic_a, ic_b, npar);
}

}  // namespace flexaids
