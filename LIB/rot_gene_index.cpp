#include "flexaid.h"
#include <cstdio>
#include <cstdlib>

// Clamped rotamer-gene -> array index. See the contract in flexaid.h.
//
// Returns an index guaranteed to satisfy 0 <= idx <= res->trot, which is the
// exact range for which build_rotamers() wrote residue[].fatm[]/latm[].
//
// A clamp firing is a DEFECT REPORT, not a normal event: the gene left the
// bounds add2_optimiz_vec.cpp:44 declared for it. We print the first few
// occurrences per site to stderr and keep a global tally so a run that relied
// on clamping is never mistaken for a clean one. Set FLEXAIDDS_ROT_GUARD_FATAL=1
// to abort instead, for bisecting the upstream cause.
static long rot_guard_hits = 0;

int rot_gene_index(double gene, const resid* res, const char* site)
{
    const int trot = (res != NULL) ? res->trot : 0;

    // Round half-up in DOUBLE, never through an unsigned cast: (uint) of a
    // negative float is undefined behaviour and was the original fault.
    double r = gene + 0.5;
    int idx = (int)r;          // truncation toward zero; r may be negative
    if (r < 0.0) idx = 0;      // any negative draw floors at the rigid rotamer

    int clamped = idx;
    if (clamped < 0)     clamped = 0;
    if (clamped > trot)  clamped = trot;

    if (clamped != idx) {
        long n = ++rot_guard_hits;
        if (n <= 8) {
            fprintf(stderr,
                "[ROT-GUARD] %s: rotamer gene %.6f -> index %d out of range "
                "[0,%d]; clamped to %d (occurrence %ld)\n",
                (site != NULL ? site : "?"), gene, idx, trot, clamped, n);
        } else if (n == 9) {
            fprintf(stderr, "[ROT-GUARD] further occurrences suppressed; "
                            "total reported at exit\n");
        }
        const char* fatal = getenv("FLEXAIDDS_ROT_GUARD_FATAL");
        if (fatal && fatal[0] == '1') {
            fprintf(stderr, "[ROT-GUARD] FLEXAIDDS_ROT_GUARD_FATAL=1 -> abort\n");
            abort();
        }
    }
    return clamped;
}

long rot_gene_guard_hits(void) { return rot_guard_hits; }
