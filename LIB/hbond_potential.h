// hbond_potential.h — Angular-dependent hydrogen bond potential
//
// Implements a Gaussian bell potential for H-bond scoring that accounts for
// both donor-acceptor distance and D-H...A angle. Also differentiates
// standard H-bonds from salt bridges based on atom charges.
//
// Integration: called from vcfunction.cpp during the pairwise contact loop.
// assign_virtual_h_geometry() called from top.cpp after type256 population.
//
// Virtual-H architecture: for heavy-atom-only PDB inputs, H positions are
// computed from standard bond geometry (DonorGeom recipe stored per atom in
// atom_struct::vH_kind/vH_nbr). The reconstruction reads live coords so the
// vH tracks GA moves of flexible donors without desync.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cmath>
#include <cstdint>
#include "atom_typing_256.h"

// Forward-declare atom_struct to avoid circular include with flexaid.h.
// Member access in inline functions below works because every TU that
// includes this header also includes flexaid.h (via gaboom.h) first.
struct atom_struct;

namespace hbond {

// ── angle_deg: angle at vertex `a` between vectors a→b and a→c ──────────────
inline double angle_deg(const float* a, const float* b, const float* c) {
    double ab[3] = { b[0]-a[0], b[1]-a[1], b[2]-a[2] };
    double ac[3] = { c[0]-a[0], c[1]-a[1], c[2]-a[2] };
    double dot    = ab[0]*ac[0] + ab[1]*ac[1] + ab[2]*ac[2];
    double mag_ab = std::sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
    double mag_ac = std::sqrt(ac[0]*ac[0] + ac[1]*ac[1] + ac[2]*ac[2]);
    if (mag_ab < 1e-8 || mag_ac < 1e-8) return 0.0;
    double cos_theta = dot / (mag_ab * mag_ac);
    if (cos_theta >  1.0) cos_theta =  1.0;
    if (cos_theta < -1.0) cos_theta = -1.0;
    return std::acos(cos_theta) * 180.0 / 3.14159265358979323846;
}

// ── DonorGeom: virtual-H placement recipe ─────────────────────────────────
// Encodes which standard-geometry rule reconstructs the virtual H from live
// heavy-neighbor coordinates. Stored as atom_struct::vH_kind (uint8_t).
// Set once per complex by assign_virtual_h_geometry() in top.cpp.
// Reconstructed by build_virtual_H() at every scoring call.
enum DonorGeom : uint8_t {
    VHG_NONE     = 0, // no directional donor (or explicit H handled separately)
    VHG_AMIDE    = 1, // N.am: external bisector of 2 heavy nbrs (in-plane)
    VHG_SP2_1NBR = 2, // N.2/N.pl3 with 1 heavy nbr: anti to single bond
    VHG_SP2_2NBR = 3, // N.2/N.pl3 with 2 heavy nbrs: external bisector
    VHG_SP3_2NBR = 4, // N.3 secondary: tetrahedral out-of-plane (1 H)
    VHG_SP3_1NBR = 5, // N.3 primary with 1 heavy nbr: 2 H on tetrahedral cone
    VHG_HYDROXYL = 6, // O.3/S.3 with 1 heavy nbr: canonical 104.5° bend
};

// ── Inline 3D vector math (private implementation detail) ─────────────────
namespace detail {
    inline void sub3(double out[3], const float* a, const float* b) {
        out[0]=a[0]-b[0]; out[1]=a[1]-b[1]; out[2]=a[2]-b[2];
    }
    inline double dot3(const double a[3], const double b[3]) {
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }
    inline void cross3(double out[3], const double a[3], const double b[3]) {
        out[0]=a[1]*b[2]-a[2]*b[1];
        out[1]=a[2]*b[0]-a[0]*b[2];
        out[2]=a[0]*b[1]-a[1]*b[0];
    }
    // Returns false if vector is degenerate (near-zero magnitude).
    inline bool normalize3(double v[3]) {
        double m = std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        if (m < 1e-8) return false;
        v[0]/=m; v[1]/=m; v[2]/=m; return true;
    }
    // Arbitrary unit vector perpendicular to unit vector v.
    inline void perp3(double out[3], const double v[3]) {
        double ref[3] = {0.0,0.0,0.0};
        if (std::fabs(v[0]) < 0.9) ref[0]=1.0; else ref[2]=1.0;
        double d = dot3(ref,v);
        out[0]=ref[0]-d*v[0]; out[1]=ref[1]-d*v[1]; out[2]=ref[2]-d*v[2];
        normalize3(out);
    }
} // namespace detail

// ── build_virtual_H: reconstruct up to 2 vH positions from live coords ───────
// Reads atoms[vH_nbr[k]].coor live — correct even for GA-moved flexible donors.
// Returns H count (0, 1, or 2). out[k][3] receives the k-th vH position.
// Bond lengths: N-H=1.01Å, O-H=0.97Å (type 14), S-H=1.34Å (type 18).
inline int build_virtual_H(const atom_struct* atoms,
                           const atom_struct& donor,
                           float out[2][3])
{
    if (donor.vH_kind == VHG_NONE) return 0;

    const int k0 = donor.vH_nbr[0];
    const int k1 = donor.vH_nbr[1];

    double d0[3]={0,0,0}, d1[3]={0,0,0};
    bool has0=false, has1=false;
    if (k0>=0) { detail::sub3(d0, atoms[k0].coor, donor.coor); has0=detail::normalize3(d0); }
    if (k1>=0) { detail::sub3(d1, atoms[k1].coor, donor.coor); has1=detail::normalize3(d1); }

    float blen = 1.01f;
    if (donor.type==14) blen=0.97f;
    if (donor.type==18) blen=1.34f;

    auto place = [&](const double dir[3], int slot) {
        out[slot][0]=donor.coor[0]+blen*(float)dir[0];
        out[slot][1]=donor.coor[1]+blen*(float)dir[1];
        out[slot][2]=donor.coor[2]+blen*(float)dir[2];
    };

    switch (donor.vH_kind) {
    case VHG_AMIDE:    // N.am planar: external bisector of C-N-Cα
    case VHG_SP2_2NBR: // same formula — external bisector of 2 heavy bonds
    {
        if (!has0||!has1) return 0;
        double dir[3]={-(d0[0]+d1[0]),-(d0[1]+d1[1]),-(d0[2]+d1[2])};
        if (!detail::normalize3(dir)) return 0;
        place(dir,0); return 1;
    }
    case VHG_SP2_1NBR: // anti to single bond (imine, terminal amide-like)
    {
        if (!has0) return 0;
        double dir[3]={-d0[0],-d0[1],-d0[2]};
        place(dir,0); return 1;
    }
    case VHG_SP3_2NBR: // sp3 secondary amine, 2 heavy nbrs → 1 H out-of-plane
    {
        // Derivation: sp3 N with bonds to d0,d1; H satisfies H·d0=H·d1=-1/3
        // → H = normalize(-0.5*(d0+d1) + sqrt(2/3)*normalize(d0×d1))
        if (!has0||!has1) return 0;
        double nrm[3]; detail::cross3(nrm,d0,d1);
        if (!detail::normalize3(nrm)) return 0;
        double dir[3]={
            -0.5*(d0[0]+d1[0])+0.8165*nrm[0],
            -0.5*(d0[1]+d1[1])+0.8165*nrm[1],
            -0.5*(d0[2]+d1[2])+0.8165*nrm[2]
        };
        if (!detail::normalize3(dir)) return 0;
        place(dir,0); return 1;
    }
    case VHG_SP3_1NBR: // primary amine, 1 heavy nbr → 2 H on tetrahedral cone
    {
        if (!has0) return 0;
        double a[3]={-d0[0],-d0[1],-d0[2]}; // unit axis away from heavy nbr
        double p[3],q[3];
        detail::perp3(p,a);
        detail::cross3(q,a,p); detail::normalize3(q);
        // cos(109.47°)≈-0.3333, sin(109.47°)≈0.9428; H0 and H1 at 120° apart
        const double base=-0.3333, rad=0.9428;
        double dir0[3]={base*a[0]+rad*p[0],base*a[1]+rad*p[1],base*a[2]+rad*p[2]};
        detail::normalize3(dir0); place(dir0,0);
        double dir1[3]={
            base*a[0]+rad*(-0.5*p[0]+0.8660*q[0]),
            base*a[1]+rad*(-0.5*p[1]+0.8660*q[1]),
            base*a[2]+rad*(-0.5*p[2]+0.8660*q[2])
        };
        detail::normalize3(dir1); place(dir1,1); return 2;
    }
    case VHG_HYDROXYL: // O.3/S.3: H at 104.5° from C-O axis, arbitrary torsion
    {
        // d0=O→C; H direction: cos(104.5°)*d0 + sin(104.5°)*perp
        // cos(104.5°)≈-0.2504, sin(104.5°)≈0.9682
        if (!has0) return 0;
        double p[3]; detail::perp3(p,d0);
        double dir[3]={
            -0.2504*d0[0]+0.9682*p[0],
            -0.2504*d0[1]+0.9682*p[1],
            -0.2504*d0[2]+0.9682*p[2]
        };
        detail::normalize3(dir); place(dir,0); return 1;
    }
    default: return 0;
    }
}

inline bool is_hydrogen_atom(const atom_struct& atom) {
    return atom.element[0] == 'H' ||
           (atom.element[0] == ' ' && atom.element[1] == 'H') ||
           atom.name[0] == 'H';
}

// Find the index of a bonded hydrogen atom for a given donor atom.
// Returns internal atom index of the H, or -1 if none found.
inline int find_bonded_hydrogen(const atom_struct* atoms, const atom_struct& donor) {
    int nbonds = donor.bond[0];
    for (int b = 1; b <= nbonds && b <= 6; ++b) {
        int idx = donor.bond[b];
        if (idx < 0) continue;
        if (is_hydrogen_atom(atoms[idx])) return idx;
    }
    return -1;
}

// ── assign_virtual_h_geometry: populate vH recipe per donor atom ─────────────
// Called once per atom in the type256 loop (top.cpp) after encode_from_sybyl
// and bond topology are final. explicit_h and heavy_bonds are pre-computed at
// the call site — pass them directly to avoid redundant bond iteration.
inline void assign_virtual_h_geometry(atom_struct* atoms, int i,
                                      int explicit_h, int heavy_bonds,
                                      bool is_pro = false)
{
    atom_struct& a = atoms[i];
    a.vH_kind   = VHG_NONE;
    a.vH_n      = 0;
    a.vH_nbr[0] = -1;
    a.vH_nbr[1] = -1;

    // Explicit H present: find_bonded_hydrogen handles the angle term.
    if (explicit_h > 0) return;

    // Collect up to 2 heavy-neighbor internal indices (bond[] order).
    int nheavy = 0, hidx[2] = {-1,-1};
    for (int b=1; b<=a.bond[0] && b<=6 && nheavy<2; ++b) {
        int nb = a.bond[b];
        if (nb < 0) continue;
        if (is_hydrogen_atom(atoms[nb])) continue;
        hidx[nheavy++] = nb;
    }

    // N.3 is aliased onto N.am (row 11) for scoring because matrix row 8 is
    // all-zero, so a.type never holds 8 at this point. Geometry must still
    // follow the real sp3 chemistry: dispatch on the recorded original row so
    // an aliphatic amine gets pyramidal SP3 geometry (and, for a primary
    // amine, two virtual H) rather than the planar amide bisector.
    const int geom_type = (a.sybyl_orig != 0) ? a.sybyl_orig : a.type;

    switch (geom_type) {
    case 7: // N.2 — sp2 imine/enamine; donor when 1 heavy bond
        if (heavy_bonds<=1 && nheavy>=1) {
            a.vH_kind=VHG_SP2_1NBR; a.vH_n=1; a.vH_nbr[0]=hidx[0];
        }
        break;
    case 8: // N.3 — sp3 amine
    {
        const int valence = (a.charge>=0.3f) ? 4 : 3;
        const int n_h = valence - heavy_bonds;
        if (n_h<=0) break;
        if (heavy_bonds>=2 && nheavy>=2) {
            a.vH_kind=VHG_SP3_2NBR; a.vH_n=1;
            a.vH_nbr[0]=hidx[0]; a.vH_nbr[1]=hidx[1];
        } else if (heavy_bonds==1 && nheavy>=1) {
            a.vH_kind=VHG_SP3_1NBR; a.vH_n=2; a.vH_nbr[0]=hidx[0];
        }
        break;
    }
    case 10: // N.ar: VHG only when type256 confirms N-H (pyrrole/indole/benzimidazole-NH).
        // Pyridine-like N.ar (acceptor-only, n_hydrogens=0) must stay VHG_NONE —
        // type256 donor bit is the single source of truth from conservative_implicit_h_count.
        // Without this guard, pyridine N.ar would get ghost VHG_SP2_2NBR, which the
        // D/A gate in compute_hbond_energy blocks, but it wastes vH reconstruction work
        // and obscures type256 intent.
        if (!atom256::get_hbond_donor(a.type256)) break;
        if (heavy_bonds == 2 && nheavy >= 2) {
            a.vH_kind=VHG_SP2_2NBR; a.vH_n=1;
            a.vH_nbr[0]=hidx[0]; a.vH_nbr[1]=hidx[1];
        } else if (heavy_bonds == 1 && nheavy >= 1) {
            a.vH_kind=VHG_SP2_1NBR; a.vH_n=1; a.vH_nbr[0]=hidx[0];
        }
        break;
    case 11: // N.am — amide N; donor with planar VHG_AMIDE geometry
        // PRO backbone N is tertiary (no labile H): skip all donor assignment.
        if (is_pro) break;
        // External bisector of C-N-Cα gives correct in-plane H direction.
        // Angular discrimination (not blanket suppression) prevents false minima.
        // Restrict VHG_AMIDE to primary/secondary amide (heavy_bonds==2);
        // tertiary N.am (ring junction, N-methyl amide) has no labile H.
        if (heavy_bonds==2 && nheavy>=2) {
            a.vH_kind=VHG_AMIDE; a.vH_n=1;
            a.vH_nbr[0]=hidx[0]; a.vH_nbr[1]=hidx[1];
        } else if (heavy_bonds==1 && nheavy>=1) {
            a.vH_kind=VHG_SP2_1NBR; a.vH_n=1; a.vH_nbr[0]=hidx[0];
        }
        break;
    case 12: // N.pl3 — planar N (guanidinium, urea, etc.)
    {
        const int n_h = 3 - heavy_bonds;
        if (n_h<=0) break;
        if (heavy_bonds>=2 && nheavy>=2) {
            a.vH_kind=VHG_SP2_2NBR; a.vH_n=1;
            a.vH_nbr[0]=hidx[0]; a.vH_nbr[1]=hidx[1];
        } else if (heavy_bonds==1 && nheavy>=1) {
            a.vH_kind=VHG_SP2_1NBR; a.vH_n=1; a.vH_nbr[0]=hidx[0];
        }
        break;
    }
    case 14: // O.3 — hydroxyl (blen=0.97Å applied in build_virtual_H)
    case 18: // S.3 — thiol   (blen=1.34Å applied in build_virtual_H)
        if (heavy_bonds<=1 && nheavy>=1) {
            a.vH_kind=VHG_HYDROXYL; a.vH_n=1; a.vH_nbr[0]=hidx[0];
        }
        break;
    default:
        break;
    }
}

// donor_angle_term: D-H...A angle Gaussian term.
// Angle convention: vertex = H position, so D-H...A = 180° is ideal.
// Priority: (1) explicit bonded H from find_bonded_hydrogen,
//           (2) virtual H from build_virtual_H (reads live coords),
//           (3) returns 0.0 (caller applies 0.3 fallback).
inline double donor_angle_term(const atom_struct* atoms,
                               const atom_struct& donor,
                               const atom_struct& acceptor,
                               double optimal_angle,
                               double sigma_angle) {
    // 1) Explicit bonded H
    int h_idx = find_bonded_hydrogen(atoms, donor);
    if (h_idx >= 0) {
        double theta = angle_deg(atoms[h_idx].coor, donor.coor, acceptor.coor);
        double da = (theta - optimal_angle) / sigma_angle;
        return std::exp(-0.5 * da * da);
    }
    // 2) Virtual H reconstructed from live heavy-neighbor coords
    float vH[2][3];
    int n = build_virtual_H(atoms, donor, vH);
    if (n <= 0) return 0.0;
    double best = 0.0;
    for (int k = 0; k < n; ++k) {
        double theta = angle_deg(vH[k], donor.coor, acceptor.coor);
        double da = (theta - optimal_angle) / sigma_angle;
        double t = std::exp(-0.5 * da * da);
        if (t > best) best = t;
    }
    return best; // multi-H donors (NH2): best-satisfied H wins
}

inline bool is_charge_supported_salt_bridge(const atom_struct& a,
                                            const atom_struct& b) {
    constexpr float Q_SALT = 0.30f;
    return (a.charge <= -Q_SALT && b.charge >= Q_SALT) ||
           (b.charge <= -Q_SALT && a.charge >= Q_SALT);
}

// Compute the angular-dependent H-bond energy between two contacting atoms.
//
// E_hb = weight * exp(-0.5*((d-d0)/σd)²) * angle_term
//
// Directionality: donor_angle_term() uses explicit H first, then virtual-H
// reconstruction from live heavy-neighbor coords (set by
// assign_virtual_h_geometry in top.cpp). N.am backbone donors use VHG_AMIDE
// (planar external-bisector geometry) — correctly rejecting acceptors outside
// the amide plane without the v57d false-minima that arose from implicit-H
// without angular discrimination.
//
// Salt bridge: ionic pairs bypass the angular gate and use salt_bridge_weight.
// Per-pair cap: -2.0 kcal/mol units.
inline double compute_hbond_energy(
    const atom_struct* atoms,
    int idx_a, int idx_b,
    double dist,
    double optimal_dist,
    double optimal_angle,
    double sigma_dist,
    double sigma_angle,
    double weight,
    double salt_bridge_weight)
{
    const atom_struct& a = atoms[idx_a];
    const atom_struct& b = atoms[idx_b];

    const bool donor_a    = atom256::get_hbond_donor(a.type256);
    const bool acceptor_a = atom256::get_hbond_acceptor(a.type256);
    const bool donor_b    = atom256::get_hbond_donor(b.type256);
    const bool acceptor_b = atom256::get_hbond_acceptor(b.type256);

    const bool a_to_b    = donor_a && acceptor_b;
    const bool b_to_a    = donor_b && acceptor_a;
    const bool salt_bridge = is_charge_supported_salt_bridge(a, b);
    if (!a_to_b && !b_to_a && !salt_bridge) return 0.0;

    double dd = (dist - optimal_dist) / sigma_dist;
    double E_dist = std::exp(-0.5 * dd * dd);

    // Angular term: explicit H → virtual H → 0.0 (handled per direction)
    double best_angle_term = 0.0;
    if (a_to_b) {
        double term = donor_angle_term(atoms, a, b, optimal_angle, sigma_angle);
        if (term > best_angle_term) best_angle_term = term;
    }
    if (b_to_a) {
        double term = donor_angle_term(atoms, b, a, optimal_angle, sigma_angle);
        if (term > best_angle_term) best_angle_term = term;
    }

    // Reduced fallback only when both explicit and virtual H are unavailable
    // (rare: donor with underdetermined geometry — no recorded heavy neighbors).
    if (best_angle_term == 0.0) best_angle_term = 0.3;

    double w = salt_bridge ? salt_bridge_weight : weight;
    double energy = w * E_dist * best_angle_term;
    constexpr double HBOND_PAIR_MIN = -2.0;
    if (energy < HBOND_PAIR_MIN) energy = HBOND_PAIR_MIN;

    // Per-pair H-bond debug: activated by FLEXAIDS_VH_DEBUG=1.
    // Prints every non-zero H-bond pair with angle term, distance term, and
    // whether the angle came from explicit-H, virtual-H, or 0.3 fallback.
    // Used to diagnose cf_native collapse (e.g. 1JD0: expected ~-23, got -1.23).
    static const bool s_vh_debug = (std::getenv("FLEXAIDS_VH_DEBUG") != nullptr);
    if (s_vh_debug && energy != 0.0) {
        const char* angle_src = "fallback(0.3)";
        if (a_to_b || b_to_a) {
            // Re-determine source for logging (lightweight, debug-only path)
            const atom_struct& donor_atom = a_to_b ? a : b;
            int h_idx = find_bonded_hydrogen(atoms, donor_atom);
            float vH[2][3]; int nv = build_virtual_H(atoms, donor_atom, vH);
            if      (h_idx >= 0) angle_src = "explicit-H";
            else if (nv > 0)     angle_src = "virtual-H";
        }
        printf("[hbdbg] %5s[%d]->%5s[%d] dist=%.3f E_dist=%.4f "
               "angle_term=%.4f src=%-14s E=%.4f%s\n",
               a.name, idx_a, b.name, idx_b,
               dist, E_dist, best_angle_term, angle_src, energy,
               salt_bridge ? " [salt]" : "");
    }

    return energy;
}

} // namespace hbond
