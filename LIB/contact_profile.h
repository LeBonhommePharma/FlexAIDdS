// =============================================================================
// contact_profile.h — atom-type-pair contact-surface vector
//                     (gate FLEXAIDDS_CONTACT_PROFILE, DEFAULT OFF)
//
// Apache-2.0 (c) 2026 Le Bonhomme Pharma
//
// WHAT THIS IS
// ------------
// The FlexAID contact function contracts a rich per-type-pair profile into a
// single scalar:
//
//     CF.com = sum over contacts  eps(type_i, type_j) * S(i, j)
//
// The contraction is lossy and it is lossy in a way that matters for the
// flexible-receptor experiment: once CF is a number you can no longer ask
// WHICH surface produced it, so you cannot separate "the ligand made good
// contacts" from "the search evicted a side chain and the ligand then filled
// the hole". This header keeps S(i, j) BEFORE the contraction: the total
// Voronoi contact surface accumulated per UNORDERED atom-type pair.
//
// The resulting vector is a pose descriptor that is RECEPTOR-FRAME AGNOSTIC.
// It does not reference receptor coordinates, receptor atom indices, receptor
// residue identity or a rotamer index — only atom TYPES and areas. That is
// precisely the property a flexible-receptor comparison needs, because the
// receptor frame is exactly what differs between the rigid and the flexible
// arm and is therefore exactly what makes a coordinate-space overlap metric
// unattributable.
//
// IT ADDS NOTHING TO CF. This module only reads areas that vcfunction has
// already computed and stores them. It has no term, no weight, no constant,
// and no effect on ranking. With the gate unset it does not even allocate.
//
// THE ZERO-EPSILON BLIND SPOT (the reason this instrument exists)
// --------------------------------------------------------------
// A large fraction of the MC_st0r5.2_6 interaction-matrix entries are exactly
// 0.0, and several whole atom-type rows are entirely zero — an atom on such a
// row is INVISIBLE to CF no matter how much surface it buries. That is a known
// defect of the scoring matrix, and it is undetectable from CF alone because
// the missing contribution is, by construction, zero.
//
// This vector sees it. Accumulation is deliberately NOT filtered on
// energy_matrix->weight or on any matrix value: a pair contributes to
// area_total the moment vcfunction visits its contact. The sidecar also
// carries, per pair, the CF contribution the engine actually credited
// (FA->contributions, which vcfunction already accumulates per type pair), so
// a row with
//
//     area_total > 0  AND  cf_pair == 0
//
// is a directly observable instance of the blind spot: real buried surface
// that the contact function priced at nothing.
//
// WHAT IS COUNTED
// ---------------
// Exactly the contact set that vcfunction's per-contact loop reaches, i.e.
// every Voronoi contact that survives:
//   * the intra-residue BONDED exclusion (bond/angle partners are not
//     contacts and never were), and
//   * the per-evaluation already-visited dedup (FA->contacts[]), which is what
//     makes each unordered pair counted exactly ONCE.
// Note that vcfunction's outer loop `continue`s on any atom whose
// atoms[].optres is NULL, so the counted contacts are only those seen from an
// OPTIMISABLE atom: the ligand, plus any flexed receptor side chain. Bulk
// receptor/receptor contacts are never enumerated, so this vector is bounded by
// the interface, not by the size of the protein.
//
// Contacts are split into two disjoint channels using the same
// `intramolecular` flag vcfunction itself computes at that point, so the split
// is by construction the same one the CF uses:
//   * area_inter — the two atoms are in different molecules. For a docking run
//     this is the ligand/receptor INTERFACE, and it INCLUDES ligand contacts
//     against a flexed side chain (that side chain's residue is type 0 while
//     the ligand residue is not, so the pair is intermolecular — which is
//     correct, it is part of the interface). This is the column to use as the
//     pose descriptor.
//   * area_intra — same molecule. Be precise about what lands here: it is
//     ligand-internal non-bonded contacts AND flexed-side-chain-versus-receptor
//     contacts, because both endpoints of the latter are residue type 0. It is
//     therefore a mixed channel describing conformers (ligand and receptor),
//     not the interface. Do not read it as "ligand strain".
//   * area_total = area_inter + area_intra, exactly.
//
// ORACLE STATUS
// -------------
// Comparing the profiles of TWO DOCKED POSES is oracle-free and is legitimate
// inside a production run. Comparing a docked pose against the profile of the
// NATIVE complex is an ORACLE metric: it consumes the answer. It is valid for
// benchmark analysis and for diagnosis, and it must NEVER be turned into a
// production scoring or selection term. See scripts/contact_profile_tanimoto.py.
//
// DEFAULT-OFF CONTRACT
// --------------------
// enabled() is false unless FLEXAIDDS_CONTACT_PROFILE is set to a value other
// than "" or "0". Callers hoist that bool into a file-scope static so the hot
// loop tests a plain bool. With the gate off:
//   * begin() is never called, so live().ntypes stays 0 and all four vectors
//     stay default-constructed and empty — no allocation, ever;
//   * add() is never called;
//   * snapshot()/write_csv() are never called, and write_csv() would refuse
//     anyway (ntypes <= 0);
//   * no CF channel, no REMARK and no existing artifact is touched.
// The engine is therefore bit-identical with the gate unset.
// =============================================================================

#pragma once

#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace flexaids {
namespace contact_profile {

/// Gate: FLEXAIDDS_CONTACT_PROFILE. Unset, empty or exactly "0" -> disabled.
/// Read once (thread-safe magic static). Hot callers must hoist the result
/// into a file-scope static so the per-contact loop never touches the guard.
inline bool enabled()
{
    static const bool v = [](){
        const char* s = std::getenv("FLEXAIDDS_CONTACT_PROFILE");
        if (s == nullptr || s[0] == '\0') return false;
        if (s[0] == '0' && s[1] == '\0')  return false;
        return true;
    }();
    return v;
}

/// Number of unordered pairs over `ntypes` 1-based atom types.
/// 40 types -> 820, matching the interaction-matrix dimensionality.
inline std::size_t packed_size(int ntypes)
{
    if (ntypes <= 0) return 0;
    const long long n = static_cast<long long>(ntypes);
    return static_cast<std::size_t>(n * (n + 1) / 2);
}

/// Packed upper-triangular index of the unordered pair {type_i, type_j}.
/// Types are 1-based ([1, ntypes]) exactly as stored in atom::type. Callers
/// MUST have range-checked both types first; add() does.
inline std::size_t pair_index(int ntypes, int type_i, int type_j)
{
    const long long lo = (type_i < type_j ? type_i : type_j) - 1;  // 0-based row
    const long long hi = (type_i < type_j ? type_j : type_i) - 1;  // 0-based col
    const long long n  = static_cast<long long>(ntypes);
    // rows 0..lo-1 hold (n) + (n-1) + ... = lo*n - lo*(lo-1)/2 entries.
    return static_cast<std::size_t>(lo * n - (lo * (lo - 1)) / 2 + (hi - lo));
}

/// One evaluation's profile. All four vectors are packed upper-triangular and
/// have length packed_size(ntypes); ntypes == 0 means "never populated".
struct Accumulator {
    int       ntypes    = 0;
    long long ncontacts = 0;           ///< number of contacts folded in
    std::vector<double> area_total;    ///< A^2, all counted contacts
    std::vector<double> area_inter;    ///< A^2, different-molecule contacts
    std::vector<double> area_intra;    ///< A^2, same-molecule contacts
    std::vector<double> cf_pair;       ///< CF contribution, filled by snapshot()

    void reset(int n)
    {
        const std::size_t sz = packed_size(n);
        ntypes    = n;
        ncontacts = 0;
        area_total.assign(sz, 0.0);
        area_inter.assign(sz, 0.0);
        area_intra.assign(sz, 0.0);
        cf_pair.assign(sz, 0.0);
    }

    void clear_all()
    {
        ntypes    = 0;
        ncontacts = 0;
        area_total.clear();
        area_inter.clear();
        area_intra.clear();
        cf_pair.clear();
    }
};

/// Accumulator for the evaluation currently in progress on THIS thread.
/// vcfunction runs under OpenMP; thread_local keeps concurrent GA fitness
/// evaluations from interleaving into one another's profile.
inline Accumulator& live()
{
    static thread_local Accumulator a;
    return a;
}

/// Frozen copy of the last evaluation this thread explicitly snapshotted.
/// Emitters read from here, never from live(), so an intervening evaluation
/// cannot silently replace the profile that is about to be written.
inline Accumulator& last()
{
    static thread_local Accumulator a;
    return a;
}

/// Start a new evaluation. Called once at the top of vcfunction, gate-ON only.
inline void begin(int ntypes)
{
    if (ntypes <= 0) { live().clear_all(); return; }
    live().reset(ntypes);
}

/// Fold one Voronoi contact into the live profile. Hot path: one bounds check,
/// one index computation, two adds. Types are 1-based; out-of-range types are
/// dropped rather than indexed (a corrupt type must never write out of bounds).
inline void add(int type_i, int type_j, double area, bool intramolecular)
{
    Accumulator& a = live();
    if (a.ntypes <= 0) return;
    if (type_i < 1 || type_j < 1 || type_i > a.ntypes || type_j > a.ntypes) return;
    const std::size_t k = pair_index(a.ntypes, type_i, type_j);
    if (k >= a.area_total.size()) return;
    a.area_total[k] += area;
    if (intramolecular) a.area_intra[k] += area;
    else                a.area_inter[k] += area;
    ++a.ncontacts;
}

/// Freeze the live profile and fold in the per-type-pair CF contributions the
/// engine credited for the same evaluation.
///
/// @param contributions FA->contributions — a dense ntypes*ntypes float matrix
///        that vcfunction memsets at entry and accumulates symmetrically, so
///        cell (a,b) already holds the full unordered-pair sum. May be nullptr,
///        in which case cf_pair stays all-zero.
/// @param ntypes must equal the ntypes passed to begin(); mismatch leaves
///        cf_pair all-zero rather than mis-indexing.
inline void snapshot(const float* contributions, int ntypes)
{
    last() = live();                      // deep copy of all four vectors
    Accumulator& s = last();
    if (s.ntypes <= 0) return;
    if (contributions == nullptr || ntypes != s.ntypes) return;
    const std::size_t n = static_cast<std::size_t>(ntypes);
    for (int a = 0; a < s.ntypes; ++a) {
        for (int b = a; b < s.ntypes; ++b) {
            const std::size_t k = pair_index(s.ntypes, a + 1, b + 1);
            if (k >= s.cf_pair.size()) continue;
            s.cf_pair[k] = static_cast<double>(
                contributions[static_cast<std::size_t>(a) * n +
                              static_cast<std::size_t>(b)]);
        }
    }
}

/// Write the frozen profile as a commented CSV sidecar.
///
/// Layout — `#`-prefixed metadata header, then a CSV header line, then one row
/// per pair that has ANY signal (area_total != 0 or cf_pair != 0). Pairs with
/// no signal are omitted; `ntypes` in the header is sufficient to rebuild the
/// dense packed_size(ntypes) vector with zeros, which is what the Tanimoto
/// script does. Nothing is filtered on the interaction matrix.
///
/// @param path       destination file
/// @param source_tag which emitter wrote this ("cluster_emitted_pose", ...)
/// @param pose_file  the artifact this profile belongs to (join key)
/// @param pose_index rank / binding-mode index, or -1 when not applicable
/// @param cf_kind    label for the cf value ("apparent", "total", ...)
/// @param cf_value   that value, recorded for the join, never used here
/// @return true on success; false if there is nothing to write or fopen failed
inline bool write_csv(const std::string& path,
                      const char*        source_tag,
                      const std::string& pose_file,
                      int                pose_index,
                      const char*        cf_kind,
                      double             cf_value)
{
    const Accumulator& s = last();
    if (s.ntypes <= 0) return false;
    if (s.area_total.size() != packed_size(s.ntypes)) return false;

    // Staleness check for the reader. FA->contributions is accumulated with the
    // very same `contribution` that vcfunction adds to cfs->com, so the sum of
    // cf_pair over all unordered pairs reproduces the pose's total CF.com to
    // float32 precision. If a caller ever snapshots a profile that does NOT
    // belong to the pose it is written beside — e.g. because the scoring call
    // short-circuited before reaching vcfunction — this sum will disagree with
    // the pose's REMARK CF.com while everything else still looks plausible.
    // Compare the two before trusting a profile.
    double cf_pair_sum = 0.0;
    for (std::size_t k = 0; k < s.cf_pair.size(); ++k) cf_pair_sum += s.cf_pair[k];
    double area_inter_sum = 0.0, area_intra_sum = 0.0;
    for (std::size_t k = 0; k < s.area_inter.size(); ++k) area_inter_sum += s.area_inter[k];
    for (std::size_t k = 0; k < s.area_intra.size(); ++k) area_intra_sum += s.area_intra[k];

    std::FILE* f = std::fopen(path.c_str(), "w");
    if (f == nullptr) return false;

    std::fprintf(f, "# flexaidds_contact_profile v1\n");
    std::fprintf(f, "# gate = FLEXAIDDS_CONTACT_PROFILE\n");
    std::fprintf(f, "# source = %s\n", source_tag ? source_tag : "unknown");
    std::fprintf(f, "# pose_file = %s\n", pose_file.c_str());
    std::fprintf(f, "# pose_index = %d\n", pose_index);
    std::fprintf(f, "# cf_kind = %s\n", cf_kind ? cf_kind : "unknown");
    std::fprintf(f, "# cf_value = %.6f\n", cf_value);
    std::fprintf(f, "# ntypes = %d\n", s.ntypes);
    std::fprintf(f, "# npairs = %llu\n",
                 static_cast<unsigned long long>(packed_size(s.ntypes)));
    std::fprintf(f, "# ncontacts = %lld\n", s.ncontacts);
    std::fprintf(f, "# area_inter_sum = %.6f\n", area_inter_sum);
    std::fprintf(f, "# area_intra_sum = %.6f\n", area_intra_sum);
    // Must match the pose's REMARK CF.com to float32 precision; see above.
    std::fprintf(f, "# cf_pair_sum = %.6f\n", cf_pair_sum);
    std::fprintf(f, "# units = angstrom^2 (areas); cf_pair in CF units\n");
    std::fprintf(f, "# types are 1-based; rows with no signal are omitted\n");
    std::fprintf(f, "# area_total = area_inter + area_intra\n");
    std::fprintf(f, "# NOT filtered on the interaction matrix: a row with\n");
    std::fprintf(f, "#   area_total > 0 and cf_pair == 0 is surface the CF priced at nothing\n");
    std::fprintf(f, "type_i,type_j,area_total,area_inter,area_intra,cf_pair\n");

    for (int a = 1; a <= s.ntypes; ++a) {
        for (int b = a; b <= s.ntypes; ++b) {
            const std::size_t k = pair_index(s.ntypes, a, b);
            if (k >= s.area_total.size()) continue;
            if (s.area_total[k] == 0.0 && s.cf_pair[k] == 0.0) continue;
            std::fprintf(f, "%d,%d,%.6f,%.6f,%.6f,%.6f\n",
                         a, b,
                         s.area_total[k], s.area_inter[k],
                         s.area_intra[k], s.cf_pair[k]);
        }
    }

    std::fclose(f);
    return true;
}

} // namespace contact_profile
} // namespace flexaids
