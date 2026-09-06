// LIB/receptor_strain.h — E_strain(r), the price of evicting a receptor side chain
//
// WHY THIS EXISTS
// ---------------
// FlexAID uses implicit solvation: the absence of an atom means the space is
// water, and the contact function carries a ligand-solvent term eps(i,w)*S(i,w).
// That is coherent for a RIGID receptor, where empty space really is bulk
// solvent. It breaks for a FLEXIBLE receptor: when the GA rotates a side chain
// away, the vacated space reads as "water", water is cheap to displace, and the
// flexible arm is rewarded for manufacturing a pocket it invented. Nothing in
// the CF charges for the eviction.
//
// `Pose::receptor_strain` (BindingMode.h:51) was declared for exactly this term
//     E_strain(r) = E_conformer(r) - E_conformer(r_ref)
// with `E_total = CF + receptor_strain` (BindingMode.h:54), and has been
// identically 0.0 since it was introduced: initialised at BindingMode.cpp:1154
// and statmech.h:164, only ever copied (BindingMode.cpp:582) or averaged
// (statmech.cpp:507/920/933), never computed. This header computes it.
//
// THE ESTIMATOR (and what it is allowed to claim)
// -----------------------------------------------
// Two measured probability sources already present in this repository are used;
// no constant is invented here.
//
//   (1) EVICTION / CONFORMATIONAL-CHANGE PRIOR — set_intprob.cpp.
//       p_HAP2(restype) is the measured probability that a side chain of that
//       residue type changes conformation between apo and holo, from the
//       Holo-Apo Protein Pairs analysis (HAP2db). Cost of moving AT ALL:
//           dE_switch = -kT ln p_HAP2
//       This is the term the implicit-solvent model is missing. A VAL
//       (p=0.070) is expensive to evict; a LYS (p=0.532) is cheap. It is zero
//       for a residue that stays on its crystal rotamer.
//
//   (2) ROTAMER-LIBRARY TERM — read_rotlib.cpp:89, `rot.pro = obs/total`,
//       the Lovell/Dunbrack-style population of each library rotamer.
//           dE_lib = +kT ln( p_ref / p_rot )
//
// LIMITATION, STATED UP FRONT — THE REFERENCE PROBABILITY IS A PROXY.
//   r_ref is the CRYSTAL side-chain conformation (rotamer gene index 0). The
//   crystal conformer is NOT a library entry, so it HAS no `pro`. p_ref is
//   therefore taken as the MAXIMUM `pro` over the rotamers accepted for that
//   residue, i.e. "the crystal side chain is assumed to sit in the best-
//   populated well available to it". Consequences, in the open:
//     * dE_lib >= 0 always, and dE_lib = 0 for the most-populated rotamer.
//     * A residue whose crystal conformer is genuinely a rare rotamer is
//       OVER-charged for leaving it by up to kT ln(p_max/p_crystal).
//   This is why term (1) exists and is not optional: without it, switching to
//   the modal rotamer would be free, which is precisely the eviction loophole.
//   The honest alternative — recomputing E_conformer from a force field — is
//   not reachable from any code currently in this tree.
//
// SECOND LIMITATION — ROTAMER-OBSERVATION MODE.
//   read_rotobs.cpp never assigns `rot.pro` (it sets only obs/tot/nid), so with
//   FA->rotobs the library term is UNAVAILABLE and only term (1) contributes.
//   Every such residue is counted in StrainResult::n_unresolved and reported in
//   the REMARK. A zero must never masquerade as a computed value.
//
// GATE — FLEXAIDDS_RECEPTOR_STRAIN, DEFAULT OFF.
//   Unset => record_*() are no-ops, evaluate_*() return computed=false and
//   0.0, and every caller's behaviour is bit-identical to HEAD.
//   FLEXAIDDS_RECEPTOR_STRAIN_T sets the temperature in kelvin (default 298.15,
//   accepted range 1..1000 K); it only scales kT.
//
// UNITS. kB is in kcal/mol/K and the probabilities are dimensionless, so the
// result is in kcal/mol. CF is NOT calibrated to kcal/mol (see the proxy-only
// provenance notes in BindingMode.cpp), so any sum of the two is a proxy and
// must be reported as one.
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <map>
#include <vector>

#include "flexaid.h"    // FA_Global, atom, resid, optmap
#include "EnvFlags.h"   // flexaids::env_bool

namespace flexaids {
namespace receptor_strain {

/// Boltzmann constant in kcal/(mol*K) — same value as statmech's kB_kcal.
inline constexpr double kB_kcal_mol_K = 0.0019872041;

/// FLEXAIDDS_RECEPTOR_STRAIN — DEFAULT OFF. Snapshotted once (magic static):
/// the table writers are called from build_rotamers() and the readers from the
/// pose-emission path, and both must agree for the whole process.
inline bool enabled() noexcept
{
    static const bool on = flexaids::env_bool("FLEXAIDDS_RECEPTOR_STRAIN", false);
    return on;
}

/// FLEXAIDDS_RECEPTOR_STRAIN_T — temperature in kelvin for kT. Default 298.15.
/// Out-of-range / unparseable values fall back to the default rather than
/// silently scaling the whole term by garbage.
inline double temperature_K() noexcept
{
    static const double T = []() noexcept -> double {
        const char* s = std::getenv("FLEXAIDDS_RECEPTOR_STRAIN_T");
        if (s == nullptr || s[0] == '\0') return 298.15;
        const double v = std::atof(s);
        if (!(v > 1.0) || !(v < 1000.0)) return 298.15;
        return v;
    }();
    return T;
}

/// Probability side-tables, keyed by INTERNAL residue index (`resid` array
/// index, i.e. FA->flex_res[].inum / atoms[].ofres).
///
/// Kept out of `resid` on purpose: `resid` is realloc'd and memset in several
/// places (read_lig, build_rotamers, read_input) and adding an owning member
/// would need a matching allocation at every one of them.
struct RotamerTable {
    /// ires -> per-rotamer library probability, indexed by the SAME rotamer
    /// index the GA gene resolves to. Slot 0 is the crystal conformer and is
    /// always -1.0 (no library probability exists for it). A slot left at
    /// -1.0 means "no usable probability", never "probability zero".
    std::map<int, std::vector<double> > rot_pro;
    /// ires -> p_HAP2 for that residue type (set_intprob.cpp), 0.0 if unknown.
    std::map<int, double> hap2_prob;
};

inline RotamerTable& table()
{
    static RotamerTable t;
    return t;
}

/// Drop everything. Called at the head of build_rotamers() so a second call
/// (read_input.cpp has two call sites, top.cpp a third) cannot double-append.
inline void reset()
{
    table().rot_pro.clear();
    table().hap2_prob.clear();
}

/// Record the library probability backing ACCEPTED rotamer `rot_index` of
/// internal residue `ires`. MUST be called only after the rigid-clash test
/// accepted the rotamer, because build_rotamers() rolls `trot` back on
/// rejection and the surviving indices are the ones the GA gene addresses.
inline void record_rotamer(int ires, int rot_index, double pro)
{
    if (!enabled()) return;
    if (ires < 0 || rot_index < 1) return;
    std::vector<double>& v = table().rot_pro[ires];
    if (static_cast<std::size_t>(rot_index) >= v.size())
        v.resize(static_cast<std::size_t>(rot_index) + 1u, -1.0);
    v[static_cast<std::size_t>(rot_index)] = pro;
}

/// Record the HAP2 conformational-change prior for internal residue `ires`
/// (FA->flex_res[i].prob, filled by set_intprob()).
inline void record_residue_prior(int ires, double hap2_prob)
{
    if (!enabled()) return;
    if (ires < 0) return;
    table().hap2_prob[ires] = hap2_prob;
}

struct StrainResult {
    double total_kcal_mol = 0.0;  ///< sum of per-residue E_strain
    int    n_flexible     = 0;    ///< side-chain genes seen
    int    n_moved        = 0;    ///< residues off their crystal rotamer
    int    n_unresolved   = 0;    ///< moved residues with NO usable probability
    bool   computed       = false;///< false => gate off / no data; NOT "0.0 strain"
};

/// E_strain for one residue at rotamer index `rot_idx`.
/// rot_idx == 0 is the crystal reference by construction => exactly 0.0.
inline double residue_term(int ires, int rot_idx, double kT, int* n_unresolved)
{
    if (rot_idx <= 0) return 0.0;

    const RotamerTable& t = table();
    double e = 0.0;
    bool resolved = false;

    // (1) eviction prior — HAP2 apo/holo conformational-change probability.
    const std::map<int, double>::const_iterator ph = t.hap2_prob.find(ires);
    if (ph != t.hap2_prob.end() && ph->second > 0.0 && ph->second < 1.0) {
        e += -kT * std::log(ph->second);
        resolved = true;
    }

    // (2) rotamer-library term, relative to the best-populated accepted rotamer.
    const std::map<int, std::vector<double> >::const_iterator pr = t.rot_pro.find(ires);
    if (pr != t.rot_pro.end()) {
        const std::vector<double>& v = pr->second;
        const std::size_t k = static_cast<std::size_t>(rot_idx);
        if (k < v.size() && v[k] > 0.0) {
            double p_ref = 0.0;
            for (std::size_t m = 1; m < v.size(); ++m)
                if (v[m] > p_ref) p_ref = v[m];
            if (p_ref > 0.0) {
                e += kT * std::log(p_ref / v[k]);
                resolved = true;
            }
        }
    }

    if (!resolved && n_unresolved != nullptr) ++(*n_unresolved);
    return e;
}

/// Evaluate against the LIVE receptor state, i.e. residue[].rot as already set
/// by ic2cf() for the pose whose coordinates are about to be written. This is
/// the same field write_pdb.cpp:40 uses to pick the side-chain atoms it emits,
/// so it is by construction the conformer that lands in the PDB.
inline StrainResult evaluate_live(const FA_Global* FA, const atom* atoms, const resid* residue)
{
    StrainResult r;
    if (!enabled()) return r;
    if (FA == nullptr || atoms == nullptr || residue == nullptr) return r;
    if (FA->map_par == nullptr) return r;

    const double kT = kB_kcal_mol_K * temperature_K();
    for (int i = 0; i < FA->npar; ++i) {
        if (FA->map_par[i].typ != 4) continue;          // 4 == side-chain rotamer gene
        const int ires = atoms[FA->map_par[i].atm].ofres;
        if (ires < 0) continue;
        ++r.n_flexible;
        const int rot_idx = residue[ires].rot;
        if (rot_idx > 0) ++r.n_moved;
        r.total_kcal_mol += residue_term(ires, rot_idx, kT, &r.n_unresolved);
    }
    r.computed = true;
    return r;
}

/// Evaluate straight from a chromosome's gene vector, without disturbing any
/// Local, self-contained copy of the rot_gene_index() clamp (rot_gene_index.cpp:19-29).
/// Deliberately NOT a call to that function: twelve GoogleTest targets compile
/// BindingMode.cpp without linking LIB/rot_gene_index.cpp and without flexaid_core,
/// and a header-only term must not add a link dependency to them. Semantics are
/// identical; only the ROT-GUARD stderr defect report is omitted, which is correct
/// here because this evaluator runs once per pose and the guard is already reported
/// from the five real derivation sites.
inline int rot_index_clamped(double gene, const resid* res) noexcept
{
    const int trot = (res != nullptr) ? res->trot : 0;
    const double r = gene + 0.5;
    int idx = static_cast<int>(r);   // truncation toward zero; r may be negative
    if (r < 0.0) idx = 0;            // any negative draw floors at the rigid rotamer
    if (idx < 0)    idx = 0;
    if (idx > trot) idx = trot;
    return idx;
}

/// global state. `genes[i]` must be the `to_ic` value for optimisation
/// parameter i (the same array ic2cf() consumes). Uses the shared
/// rot_index_clamped() clamp so an out-of-range gene cannot index fatm/latm.
inline StrainResult evaluate_genes(const FA_Global* FA, const atom* atoms, const resid* residue,
                                   const double* genes, int n_genes)
{
    StrainResult r;
    if (!enabled()) return r;
    if (FA == nullptr || atoms == nullptr || residue == nullptr || genes == nullptr) return r;
    if (FA->map_par == nullptr) return r;

    const double kT = kB_kcal_mol_K * temperature_K();
    const int n = (n_genes < FA->npar) ? n_genes : FA->npar;
    for (int i = 0; i < n; ++i) {
        if (FA->map_par[i].typ != 4) continue;
        const int ires = atoms[FA->map_par[i].atm].ofres;
        if (ires < 0) continue;
        ++r.n_flexible;
        const int rot_idx = rot_index_clamped(genes[i], &residue[ires]);
        if (rot_idx > 0) ++r.n_moved;
        r.total_kcal_mol += residue_term(ires, rot_idx, kT, &r.n_unresolved);
    }
    r.computed = true;
    return r;
}

}  // namespace receptor_strain
}  // namespace flexaids
