// =============================================================================
// local_orient.h — offline LOCAL ORIENTATION refinement of an existing pose
//
// WHY THIS EXISTS (WO-LOCALORIENT-1, see LOCAL_ORIENT_PREREGISTRATION.md)
//   Fifteen benchmark arms tested GLOBAL search knobs (selection, weighting,
//   temperature, diversity, matrix cells, crossover, restarts, budget, fitness
//   model). Every matched-and-powered contrast was null, and the two surviving
//   leads were refuted by paired tests with matched controls.
//
//   Measured on 1800 emitted poses: the median pose sits 6.20 A off the native
//   ligand centroid while the best sits at 1.42 A — budget is spent OUTSIDE the
//   pocket. But conditioning on the 29 poses already within 1.5 A of the native
//   centroid inverts the error attribution:
//
//       squared-error share   translation 0.143   ROTATION 0.495   torsion 0.349
//       sub-2 A in place                                          3 / 29  (10.3%)
//       ATTAINABLE by a centroid-fixed rotation, sqrt(off^2+conf^2) 8 / 29  (27.6%)
//
//   That 3 -> 8 gap is the headroom this stage exists to convert. It is a
//   REFINEMENT hypothesis: convert a pose the search ALREADY found. It is not a
//   sampling claim and cannot become one — nothing here searches.
//
//   CORRECTION, and it bounds what this code can deliver. An earlier version of
//   this comment claimed 20/29 "if orientation alone were corrected". That was
//   WRONG: the 20 came from an OPTIMAL RIGID SUPERPOSITION, which fits
//   TRANSLATION as well as rotation. This stage holds the centroid FIXED, so it
//   cannot use the translational part of that fit. The decomposition is exact
//   and orthogonal —
//       total^2 = offset^2 + centred^2,   centred^2 = rotation^2 + conformer^2
//   — so driving the rotation residual to zero at a fixed centroid leaves offset
//   and conformer, and the attainable in-place RMSD is sqrt(offset^2+conf^2).
//   By that measure 8 of 29, not 20. The lever is real and 2.4x smaller than
//   first claimed. It also SATURATES at 8 for every placement cutoff >= 1.5 A:
//   no pose further than 1.5 A from the native centroid is reachable by rotation
//   alone, which is exactly why the eventual arm selects targets on placement.
//
//   NOTE the unconditional picture is the opposite (translation share 0.631 vs
//   rotation 0.289). An earlier claim that "rotation separates a hit from a
//   miss" was WITHDRAWN as circular. The conditional claim above is the one this
//   code tests, and the conditioning is load-bearing.
//
// WHAT IT DOES
//   For one pose, optimise the 3 rigid-body ORIENTATION degrees of freedom with
//   the ligand CENTROID HELD FIXED and TORSIONS FROZEN, against the production
//   CF, under a bounded evaluation budget.
//
//   Centroid-fixed and torsion-frozen are structural, not policed: a rigid
//   rotation about the ligand centroid cannot move the centroid and cannot
//   change any internal coordinate. There is no code path here that could.
//
// WHY NOT AN EXISTING KNOB — verified in source, and one name is a trap
//   FLEXAIDDS_MEDOID_REFINE       DEPRECATED (ClusterRepMode.h:53-58); selects
//                                 cluster REPRESENTATIVES, refines no geometry.
//   FLEXAIDDS_BASIN_SIGMA_ANG     NOT AN ANGLE. gaboom.cpp:1036-1047 — a
//                                 catastrophic re-inject that re-randomises the
//                                 WORST fraction until Cartesian ligand RMSD vs
//                                 best exceeds sigma (default 2.0 ANGSTROM).
//                                 A GLOBAL DIVERSITY knob, the opposite of local
//                                 refinement. "_ANG" means Angstrom, not angle.
//   FLEXAIDDS_NICHE_SIGMA_ANG     niche_distance.h:58 — niching, same units.
//   FLEXAIDDS_COARSE_ORIENTATIONS config_parser.cpp:389 — orientation COUNT at
//                                 search start. Global, not local.
//   No --refine / --minimize / --polish / --local CLI flag exists.
//
// ENV SURFACE (all default OFF/inert; unset reproduces prior behaviour exactly)
//   FLEXAIDDS_LOCAL_ORIENT=1        enable the stage (default OFF)
//   FLEXAIDDS_LOCAL_ORIENT_MODE     "orient" (treatment, default) | "jitter"
//                                   (MATCHED CONTROL: identical search, identical
//                                   budget, perturbing the 3 TRANSLATIONAL DoF
//                                   instead). The jitter arm exists because a
//                                   gate-OFF control confounds refinement with
//                                   "extra CF evaluations help"; gate-OFF is a
//                                   third reference arm, not the primary control.
//   FLEXAIDDS_LOCAL_ORIENT_STEPS    CF-evaluation budget per pose (default 200)
//   FLEXAIDDS_LOCAL_ORIENT_STEP0    initial step: degrees (orient) or Angstrom
//                                   (jitter). Default 15.0 / 0.5.
//   FLEXAIDDS_LOCAL_ORIENT_MINSTEP  convergence floor, same units (default 0.5 /
//                                   0.02)
//
// DETERMINISM
//   Pattern search, no RNG: three axes x two signs per sweep, accept the best
//   improving move, halve the step when a sweep finds none. Every candidate is
//   built from the ORIGINAL coordinates rather than composed incrementally, so
//   there is no orthonormality drift and no dependence on evaluation order.
//   Same pose in, same pose out, on any machine.
//
// Apache-2.0 (c) 2026 Le Bonhomme Pharma
// =============================================================================

#pragma once

#include "flexaid.h"
#include "Vcontacts.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

// Defined in vcfunction.cpp / ic2cf.cpp; declared in flexaid.h.
double vcfunction(FA_Global*, VC_Global*, atom*, resid*,
                  std::vector<std::pair<int,int>>&, bool*);
double get_cf_evalue(cfstr*, FA_Global*);

namespace flexaids {
namespace local_orient {

enum class Mode { Orient, Jitter };

struct Config {
    bool   enabled  = false;
    Mode   mode     = Mode::Orient;
    int    budget   = 200;
    double step0    = 0.0;   // resolved per mode below
    double min_step = 0.0;
};

struct Result {
    bool   ran            = false;
    double cf_before      = 0.0;
    double cf_after       = 0.0;
    int    evals_used     = 0;
    double applied[3]     = {0.0, 0.0, 0.0};  // deg (orient) or A (jitter)
    double magnitude      = 0.0;              // deg of net rotation, or |dt| in A
    double centroid_shift = 0.0;              // MUST be ~0 in orient mode
    double max_bond_drift = 0.0;              // MUST be ~0 in both modes
};

/// Read the env surface. Absent or "0" leaves `enabled` false and the caller
/// must then not alter a single coordinate.
inline Config config_from_env()
{
    Config c;
    const char* g = std::getenv("FLEXAIDDS_LOCAL_ORIENT");
    if (!g || g[0] == '\0' || std::strcmp(g, "0") == 0) return c;
    c.enabled = true;

    if (const char* m = std::getenv("FLEXAIDDS_LOCAL_ORIENT_MODE")) {
        if (std::strcmp(m, "jitter") == 0)      c.mode = Mode::Jitter;
        else if (std::strcmp(m, "orient") == 0) c.mode = Mode::Orient;
        else {
            // Fail loud: a typo'd mode must not silently run the treatment and
            // be reported as the control.
            std::fprintf(stderr,
                "[LOCAL-ORIENT] FATAL: unknown MODE '%s' (expected "
                "'orient' or 'jitter')\n", m);
            std::exit(4);
        }
    }
    c.step0    = (c.mode == Mode::Orient) ? 15.0 : 0.5;
    c.min_step = (c.mode == Mode::Orient) ? 0.5  : 0.02;

    if (const char* s = std::getenv("FLEXAIDDS_LOCAL_ORIENT_STEPS")) {
        const int v = std::atoi(s);
        if (v <= 0) {
            std::fprintf(stderr, "[LOCAL-ORIENT] FATAL: STEPS must be > 0\n");
            std::exit(4);
        }
        c.budget = v;
    }
    if (const char* s = std::getenv("FLEXAIDDS_LOCAL_ORIENT_STEP0")) {
        const double v = std::atof(s);
        if (v > 0.0) c.step0 = v;
    }
    if (const char* s = std::getenv("FLEXAIDDS_LOCAL_ORIENT_MINSTEP")) {
        const double v = std::atof(s);
        if (v > 0.0) c.min_step = v;
    }
    if (c.min_step >= c.step0) {
        std::fprintf(stderr,
            "[LOCAL-ORIENT] FATAL: MINSTEP (%.4f) >= STEP0 (%.4f); the search "
            "would terminate before evaluating anything\n", c.min_step, c.step0);
        std::exit(4);
    }
    return c;
}

/// Rotation matrix for intrinsic X->Y->Z angles (radians). Built fresh from the
/// accumulated angles on every candidate, never composed, so R stays exactly
/// orthonormal to floating-point precision.
inline void rot_matrix(const double rx, const double ry, const double rz,
                       double R[3][3])
{
    const double cx = std::cos(rx), sx = std::sin(rx);
    const double cy = std::cos(ry), sy = std::sin(ry);
    const double cz = std::cos(rz), sz = std::sin(rz);
    R[0][0] =  cy * cz;
    R[0][1] = -cy * sz;
    R[0][2] =  sy;
    R[1][0] =  sx * sy * cz + cx * sz;
    R[1][1] = -sx * sy * sz + cx * cz;
    R[1][2] = -sx * cy;
    R[2][0] = -cx * sy * cz + sx * sz;
    R[2][1] =  cx * sy * sz + sx * cz;
    R[2][2] =  cx * cy;
}

/// Net rotation angle of R, in degrees: theta = acos((trace(R) - 1) / 2).
inline double net_angle_deg(const double R[3][3])
{
    double t = (R[0][0] + R[1][1] + R[2][2] - 1.0) * 0.5;
    if (t >  1.0) t =  1.0;
    if (t < -1.0) t = -1.0;
    return std::acos(t) * 180.0 / M_PI;
}

// -----------------------------------------------------------------------------
// The stage. `fa`..`la` are the inclusive ligand slot bounds in atoms[].
//
// Contract: on return, atoms[] holds the BEST pose found (which is the input
// pose if nothing improved). Receptor slots are never touched. If cfg.enabled
// is false this function writes nothing and returns {ran = false}.
// -----------------------------------------------------------------------------
inline Result refine(const Config& cfg,
                     FA_Global* FA, VC_Global* VC, atom* atoms, resid* residue,
                     const int fa, const int la)
{
    Result out;
    if (!cfg.enabled) return out;          // gate OFF: not one coordinate moves
    const int n = la - fa + 1;
    if (n < 2) return out;                 // a single atom has no orientation

    // Snapshot the input ligand coordinates. Every candidate transform is
    // applied to THIS, never to the working coordinates.
    std::vector<double> orig(static_cast<size_t>(n) * 3);
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k)
            orig[static_cast<size_t>(i) * 3 + k] = atoms[fa + i].coor[k];

    double c0[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k) c0[k] += orig[static_cast<size_t>(i) * 3 + k];
    for (int k = 0; k < 3; ++k) c0[k] /= static_cast<double>(n);

    // Reference intra-ligand distances, for the torsion-frozen assertion below.
    // Consecutive pairs suffice: a rigid transform preserves ALL distances, so
    // any internal change shows up here, and this stays O(n).
    std::vector<double> ref_d(static_cast<size_t>(n > 1 ? n - 1 : 0));
    for (int i = 0; i + 1 < n; ++i) {
        double s = 0.0;
        for (int k = 0; k < 3; ++k) {
            const double d = orig[static_cast<size_t>(i) * 3 + k]
                           - orig[static_cast<size_t>(i + 1) * 3 + k];
            s += d * d;
        }
        ref_d[static_cast<size_t>(i)] = std::sqrt(s);
    }

    // Write candidate (p[0],p[1],p[2]) into atoms[] and score it.
    // Orient: p = XYZ rotation angles in RADIANS about the fixed centroid.
    // Jitter: p = translation in ANGSTROM (the matched control).
    int evals = 0;
    auto score = [&](const double p[3]) -> double {
        if (cfg.mode == Mode::Orient) {
            double R[3][3];
            rot_matrix(p[0], p[1], p[2], R);
            for (int i = 0; i < n; ++i) {
                double v[3];
                for (int k = 0; k < 3; ++k)
                    v[k] = orig[static_cast<size_t>(i) * 3 + k] - c0[k];
                for (int k = 0; k < 3; ++k)
                    atoms[fa + i].coor[k] = static_cast<float>(
                        c0[k] + R[k][0] * v[0] + R[k][1] * v[1] + R[k][2] * v[2]);
            }
        } else {
            for (int i = 0; i < n; ++i)
                for (int k = 0; k < 3; ++k)
                    atoms[fa + i].coor[k] = static_cast<float>(
                        orig[static_cast<size_t>(i) * 3 + k] + p[k]);
        }
        std::vector<std::pair<int,int>> intra;
        bool err = false;
        const double penalty = vcfunction(FA, VC, atoms, residue, intra, &err);
        ++evals;
        cfstr cf{};
        if (err) { cf.wal = penalty; cf.rclash = 1; }
        else {
            for (int j = 0; j < FA->num_optres; ++j) {
                cf.com         += FA->optres[j].cf.com;
                cf.wal         += FA->optres[j].cf.wal;
                cf.sas         += FA->optres[j].cf.sas;
                cf.con         += FA->optres[j].cf.con;
                cf.elec        += FA->optres[j].cf.elec;
                cf.gist_desolv += FA->optres[j].cf.gist_desolv;
                cf.metal_coord += FA->optres[j].cf.metal_coord;
                cf.hbond       += FA->optres[j].cf.hbond;
                cf.entropy     += FA->optres[j].cf.entropy;
                cf.pb_clash    += FA->optres[j].cf.pb_clash;
            }
        }
        return get_cf_evalue(&cf, FA);
    };

    double best_p[3] = {0.0, 0.0, 0.0};
    const double cf0 = score(best_p);      // identity: reproduces the input pose
    double best_cf = cf0;

    const double unit = (cfg.mode == Mode::Orient) ? (M_PI / 180.0) : 1.0;
    double step = cfg.step0;

    while (step >= cfg.min_step && evals < cfg.budget) {
        bool improved = false;
        double sweep_p[3] = {best_p[0], best_p[1], best_p[2]};
        double sweep_cf = best_cf;
        for (int axis = 0; axis < 3 && evals < cfg.budget; ++axis) {
            for (const int sign : {+1, -1}) {
                if (evals >= cfg.budget) break;
                double cand[3] = {best_p[0], best_p[1], best_p[2]};
                cand[axis] += sign * step * unit;
                const double cf = score(cand);
                if (cf < sweep_cf) {
                    sweep_cf = cf;
                    for (int k = 0; k < 3; ++k) sweep_p[k] = cand[k];
                    improved = true;
                }
            }
        }
        if (improved) {
            best_cf = sweep_cf;
            for (int k = 0; k < 3; ++k) best_p[k] = sweep_p[k];
        } else {
            step *= 0.5;                    // no improving move at this scale
        }
    }

    // Leave atoms[] holding the BEST pose, not the last candidate probed.
    score(best_p);
    --evals;                                // the restore is bookkeeping, not search

    out.ran        = true;
    out.cf_before  = cf0;
    out.cf_after   = best_cf;
    out.evals_used = evals;
    for (int k = 0; k < 3; ++k)
        out.applied[k] = (cfg.mode == Mode::Orient)
                       ? best_p[k] * 180.0 / M_PI : best_p[k];

    if (cfg.mode == Mode::Orient) {
        double R[3][3];
        rot_matrix(best_p[0], best_p[1], best_p[2], R);
        out.magnitude = net_angle_deg(R);
    } else {
        out.magnitude = std::sqrt(best_p[0]*best_p[0] + best_p[1]*best_p[1]
                                + best_p[2]*best_p[2]);
    }

    // Post-conditions, measured rather than asserted in comments. In orient mode
    // centroid_shift must be ~0; in BOTH modes max_bond_drift must be ~0, since
    // both transforms are rigid. A nonzero value here means the stage is not
    // doing what it claims and the receipt will show it.
    double c1[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k) c1[k] += atoms[fa + i].coor[k];
    for (int k = 0; k < 3; ++k) c1[k] /= static_cast<double>(n);
    double cs = 0.0;
    for (int k = 0; k < 3; ++k) { const double d = c1[k] - c0[k]; cs += d * d; }
    out.centroid_shift = std::sqrt(cs);

    double drift = 0.0;
    for (int i = 0; i + 1 < n; ++i) {
        double s = 0.0;
        for (int k = 0; k < 3; ++k) {
            const double d = atoms[fa + i].coor[k] - atoms[fa + i + 1].coor[k];
            s += d * d;
        }
        drift = std::max(drift, std::fabs(std::sqrt(s) - ref_d[static_cast<size_t>(i)]));
    }
    out.max_bond_drift = drift;
    return out;
}

/// Write `src` with the ligand coordinate columns replaced by atoms[fa..la],
/// to `dst`. Never overwrites `src`. Receptor records are copied byte-for-byte.
inline bool write_refined_pdb(const char* src, const char* dst,
                              const atom* atoms, const int fa, const int la)
{
    FILE* in = std::fopen(src, "r");
    if (!in) return false;
    FILE* o = std::fopen(dst, "w");
    if (!o) { std::fclose(in); return false; }

    char buf[512];
    while (std::fgets(buf, sizeof(buf), in)) {
        const bool is_atom = std::strncmp(buf, "ATOM  ", 6) == 0 ||
                             std::strncmp(buf, "HETATM", 6) == 0;
        if (!is_atom || std::strlen(buf) < 54) { std::fputs(buf, o); continue; }
        char sb[8];
        std::memcpy(sb, buf + 6, 5);
        sb[5] = '\0';
        char* endp = nullptr;
        const long serial = std::strtol(sb, &endp, 10);
        int slot = -1;
        if (endp != sb)
            for (int i = fa; i <= la; ++i)
                if (atoms[i].number == static_cast<int>(serial)) { slot = i; break; }
        if (slot < 0) { std::fputs(buf, o); continue; }
        char line[512];
        std::snprintf(line, sizeof(line), "%.30s%8.3f%8.3f%8.3f%s",
                      buf, atoms[slot].coor[0], atoms[slot].coor[1],
                      atoms[slot].coor[2],
                      (std::strlen(buf) > 54) ? buf + 54 : "\n");
        std::fputs(line, o);
    }
    std::fclose(in);
    std::fclose(o);
    return true;
}

} // namespace local_orient
} // namespace flexaids
