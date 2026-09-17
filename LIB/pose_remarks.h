#ifndef FLEXAIDDS_POSE_REMARKS_H
#define FLEXAIDDS_POSE_REMARKS_H

// ─── SHARED POSE REMARK WRITER ──────────────────────────────────────────────
//
// WHY THIS FILE EXISTS
// --------------------
// The pose REMARK block was composed by FOUR hand-maintained code paths that
// drifted apart. Measured over LIB/ on 2026-09-16:
//
//     76 distinct fields, 197 REMARK format literals
//     14 fields emitted by all four writers   (the common core)
//     56 fields emitted by some but not all
//     34 fields duplicated at >1 site inside a single file
//        (BindingMode.cpp alone carries two full copies of the CF block:
//         CF.app at :856 and :1058, the RMSD pair at :984/:990 and :1154/:1160)
//
// Three concrete failures came out of that drift, none of them cosmetic:
//
//   1. DSVIB.* -- the dS_vib receipt, the quantity this project exists to
//      measure -- is emitted by cluster.cpp ONLY. DensityPeak_Cluster.cpp,
//      FOPTICS.cpp, BindingMode.cpp and FastOPTICS_cluster.cpp contain zero
//      DSVIB tokens. Selecting a non-CF backend therefore drops the entropy
//      provenance from every pose: exit code 0, correct file count, nothing
//      missing to notice, and unrecoverable once the pool is archived because
//      the entropy is a function of the pose AND the calculation at emit time.
//
//   2. The STICKY HUNGARIAN FLAG. `bool Hungarian` is declared once OUTSIDE
//      the rank loop and set true INSIDE it, never reset:
//          legacy    cluster.c:6    -> :249
//          FlexAIDdS cluster.cpp:206 -> :897      (still live at time of writing)
//      So only rank 0's "(no symmetry correction)" line is genuinely serial.
//      For ranks >= 1 BOTH lines carry the same Hungarian value -- confirmed in
//      campaign data: nosym != sym in 560/560 rank-0 heads, nosym == sym in
//      5040/5040 higher ranks. The serial value then exists NOWHERE in the
//      outputs, and the "(no symmetry correction)" label is false on ~90% of
//      poses. FOPTICS.cpp:440/:444 and BindingMode.cpp:982/:988,:1152/:1158
//      happen to declare adjacent to use and are correct; cluster.cpp and
//      DensityPeak_Cluster.cpp:39/:685 are not. Same field, different
//      correctness, decided by which backend ran.
//
//   3. Two DSVIB emitters with different mechanisms -- cluster.cpp:800-822 uses
//      "%.8f" format literals, DatasetRunner.cpp:9214-9238 uses string
//      concatenation -- so the same number can reach disk at different
//      precision depending on the path that wrote it.
//
// WHAT THIS FIXES BY CONSTRUCTION, NOT BY DISCIPLINE
// --------------------------------------------------
// Defects 1 and 2 stop being possible rather than stopping by convention:
//
//   * emit_rmsd_pair() declares its own `Hungarian` local at each call. There
//     is no enclosing scope for the flag to leak from, so a rank-N pose cannot
//     inherit rank-(N-1)'s flag state. The sticky bug has nowhere to live.
//   * append_dsvib() is part of the standard composition order. A backend that
//     forgets it does not silently omit it -- it fails
//     tests/test_pose_remark_contract.py, which pins the field union.
//
// WHAT THIS DELIBERATELY DOES NOT TOUCH
// -------------------------------------
// Clustering SEMANTICS. The pose-limit cap means a structurally different thing
// in each backend BY DESIGN: CF (cluster.cpp:399) bounds cluster CONSTRUCTION,
// DensityPeak (DensityPeak_Cluster.cpp:434-437) bounds EMISSION only with its
// creation-side cap intentionally disabled. That is LP's design decision, it
// traces to original FlexAID, and it is guarded by
// tests/test_dp_partition_invariant.py. Unify the WRITER. Never the semantics.
//
// GATING AND REPRODUCIBILITY
// --------------------------
// Behind FLEXAIDDS_UNIFIED_REMARKS, DEFAULT OFF. This is not timidity: turning
// it on CHANGES BYTES for DP and FOPTICS, because those backends gain the DSVIB
// block and the rmsd_raw/rmsd_sym debug lines they currently lack. Every stored
// arm must stay reproducible through the gate, so the default reproduces each
// backend's existing output exactly and the new path is opt-in.
//
// FLIP CONDITION, STATED AT LANDING TIME so this does not become another gate
// that is present but never reached: the default flips to ON once (a) a CF-arm
// A/B shows byte-identical poses with the gate on and off, and (b) a DP arm and
// an FO arm each show the DSVIB block present. Not "once it shows an
// improvement" -- there is no improvement to show. This is a correctness fix,
// and a null is a pass.
//
// THE NON-PURITY CONTRACT -- READ BEFORE REFACTORING
// --------------------------------------------------
// calc_rmsd() IS NOT A PURE FUNCTION. gaboom.cpp:1452 documents that it writes
// state. Every existing RMSD-REMARK site calls it TWICE IN IMMEDIATE SUCCESSION
// with the flag flipped between calls, and the raw-then-sym ORDER is
// load-bearing. Do not memoize it, do not hoist it out of a loop, do not
// compute sym before raw, and do not "optimize away" the second call because
// the arguments look identical -- they are not, the flag differs, and the
// function's internal state carries between them.
// tests/test_pose_remark_contract.py asserts this ordering.
//
// A NOTE ON WHAT THE "(symmetry corrected)" LINE ACTUALLY IS
// ----------------------------------------------------------
// It is calc_rmsd(..., Hungarian=true) -> calc_Hungarian_RMSD: a per-FlexAID-
// atom-type assignment. It is NOT a graph-automorphism symmetry correction and
// it is NOT serial. In the legacy f766a14 engine its cost matrix was also
// malformed (atoms[k+l] walked past the ligand; unfilled cells stayed 0.0 and
// the assignment could select them, dropping atoms from the sum while the
// divisor stayed num_het_atm), which made the emitted value read low -- below
// the exact value on 54/56 audited targets, and not even a bound, since a
// correct Hungarian sat below it on 3. FlexAIDdS fixed that walk at
// calc_rmsd.cpp:179-183. The label here describes the mechanism honestly so no
// downstream consumer mistakes it for a claim-grade metric; a claim-grade RMSD
// is computed post-hoc from coordinates, never read from this REMARK.

#include <cstdio>
#include <cstring>
#include <string>

#include "flexaid.h"

namespace flexaidds {
namespace remarks {

// Which pose-writing path is composing this block. Used ONLY to select the
// one discriminator line each backend legitimately owns -- never to decide
// which scientific fields are present.
enum class Emitter {
    CF,            // cluster.cpp
    DensityPeak,   // DensityPeak_Cluster.cpp
    FastOPTICS,    // FOPTICS.cpp
    BindingMode,   // BindingMode.cpp
};

inline const char* emitter_line(Emitter e)
{
    switch (e) {
        case Emitter::CF:          return NULL;  // CF writes its cluster/rank line separately
        case Emitter::DensityPeak: return "REMARK Density Peak clustering algorithm used to cluster Poses\n";
        case Emitter::FastOPTICS:  return "REMARK Fast OPTICS clustering algorithm used to order Poses in OPTICS\n";
        case Emitter::BindingMode: return NULL;
    }
    return NULL;
}

// Accumulates the block. Deliberately a plain std::string rather than the
// historical fixed char[MAX_REMARK] + strcat pattern: the old pattern is how a
// long `inputs:` path silently truncated a block mid-field.
class PoseRemarkBuilder {
public:
    explicit PoseRemarkBuilder(Emitter e) : emitter_(e) { buf_.reserve(4096); }

    void line(const char* fmt, ...)
#if defined(__GNUC__) || defined(__clang__)
        __attribute__((format(printf, 2, 3)))
#endif
        ;

    void raw(const char* s) { if (s) buf_ += s; }

    // ── the fixed composition order ──────────────────────────────────────
    // Order is part of the contract: parsers in scripts/ key on first-match,
    // so moving a field can change which value a downstream table reports.
    void append_header();                       // generated-by / optimized structure
    // CF= / CF.app come from the chromosome or the ic2cf() return, not from
    // FA_Global. The provided FA-only append_cf_terms() cannot see those
    // totals; this method is the stated amendment so backends can pass the
    // values they already computed without inventing a FA_Global field.
    void append_cf_totals(double cf, double cf_app);
    void append_cf_terms(const FA_Global* FA);  // CF.com, CF.sas, CF.wal, CF.con, CF.gist, CF.hbond
    // Residue names live on resid[], not FA_Global. Overload so the
    // "optimizable residue" line can be composed without a dangling FA lookup.
    void append_cf_terms(const FA_Global* FA, const resid* residue);
    void append_residue_sas(const FA_Global* FA);
    void append_torsions(const FA_Global* FA);  // the REMARK [%8.3f] block

    // Emits BOTH RMSD lines, in the mandatory raw-then-sym order, each with its
    // OWN local flag. This signature is the fix for the sticky-flag defect:
    // callers cannot pass a flag in, so they cannot leak one between poses.
    void emit_rmsd_pair(FA_Global* FA, atom* atoms, resid* residue,
                        gridpoint* cleftgrid);

    // The dS_vib receipt. Composed for EVERY emitter, not just CF. If the
    // entropy is unavailable for this pose, emits DSVIB.status with the reason
    // rather than emitting nothing -- "not computed" and "not recorded" must
    // never be indistinguishable downstream.
    void append_dsvib(const FA_Global* FA);
    // ic2cf() writes dS_vib onto the RETURNED cfstr, not onto FA->optres[].cf
    // (cluster.cpp:776-780). FA-only append_dsvib() therefore cannot recover
    // the numbers; this overload is the stated amendment.
    void append_dsvib(const cfstr* scored);

    void append_emitter_line() { raw(emitter_line(emitter_)); }
    void append_inputs(const char* config_path, const char* ga_path);
    void append_interaction_contributions(const FA_Global* FA, atom* atoms);

    const std::string& str() const { return buf_; }
    const char* c_str() const { return buf_.c_str(); }

private:
    Emitter emitter_;
    std::string buf_;
};

// Compose the standard block in the standard order. Every pose-writing backend
// calls THIS, then hands the result to write_pdb(..., char remark[]) at
// LIB/write_pdb.cpp:119 -- the common entry point all four already use, so this
// needs no new plumbing.
std::string compose_pose_remarks(Emitter e,
                                 FA_Global* FA,
                                 atom* atoms,
                                 resid* residue,
                                 gridpoint* cleftgrid,
                                 const char* config_path,
                                 const char* ga_path);

// Runtime gate. Default FALSE -> each backend keeps its existing hand-rolled
// block and every stored arm stays byte-reproducible.
// Uses flexaids::env_bool rather than a raw getenv first-char test: the house
// idiom `if (v && v[0] != '0')` treats an EMPTY STRING as ENABLED, which is the
// opposite of what anyone typing FLEXAIDDS_UNIFIED_REMARKS= intends. That
// defect has already cost this project one voided verification run.
bool unified_remarks_enabled();

}  // namespace remarks
}  // namespace flexaidds

#endif  // FLEXAIDDS_POSE_REMARKS_H
