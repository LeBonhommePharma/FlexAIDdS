// tests/cf_aggregator_stubs.cpp — slim stubs for test_cf_aggregator.
//
// This test target links the REAL LIB/ic2cf.cpp so that get_cf_evalue() is
// the production aggregator (not the reduced reimplementation in
// tests/stubs.cpp). ic2cf.cpp's ic2cf()/compute_ligand_h_rep() pull in a
// handful of heavy symbols that get_cf_evalue() itself never touches; we stub
// only those so the aggregator can be exercised in isolation. Deliberately
// does NOT define get_cf_evalue / get_apparent_cf_evalue — those come from the
// real ic2cf.cpp under test.

#include "gaboom.h"
#include <vector>
#include <utility>
#include "tENCoM/tencm.h"

namespace vibentropy {
double compute_vib_entropy_collapse(const std::vector<std::vector<double>>&) { return 0.0; }
}
namespace sugar_pucker {
enum class SugarType;
void apply_sugar_puckers(atom*, const std::vector<std::vector<int>>&,
                         const std::vector<float>&, const std::vector<SugarType>&) {}
}
namespace tencm {
void TorsionalENM::build_from_ligand(const atom*, int, int, float, float) {}
}

// Controllable fail points for serial contamination tests (ic2cf restore).
namespace ic2cf_test_hooks {
int g_buildcc_fail_next = 0;
int g_vcfunction_error_next = 0;
int g_buildcc_calls = 0;
int g_vcfunction_calls = 0;
// When buildcc runs and succeeds, optionally scribble moved atom coords to
// prove restore rolls them back on a subsequent vcfunction error.
bool g_buildcc_scribble_coords = false;
float g_scribble_value = 9999.0f;
}

void alter_mode(atom*, resid*, float*, int, int) {}
double vcfunction(FA_Global* FA, VC_Global*, atom* atoms, resid*,
                  std::vector<std::pair<int,int>>&, bool* error) {
    ++ic2cf_test_hooks::g_vcfunction_calls;
    if (error) *error = false;
    if (ic2cf_test_hooks::g_vcfunction_error_next > 0) {
        --ic2cf_test_hooks::g_vcfunction_error_next;
        if (error) *error = true;
        return 1.0e6;
    }
    (void)FA;
    (void)atoms;
    return 0.0;
}
bool buildcc(FA_Global* FA, atom* atoms, int nmov, int mov[]) {
    ++ic2cf_test_hooks::g_buildcc_calls;
    if (ic2cf_test_hooks::g_buildcc_fail_next > 0) {
        --ic2cf_test_hooks::g_buildcc_fail_next;
        return false;
    }
    if (ic2cf_test_hooks::g_buildcc_scribble_coords && atoms && mov && nmov > 0) {
        for (int m = 0; m < nmov; ++m) {
            const int ai = mov[m];
            atoms[ai].coor[0] = ic2cf_test_hooks::g_scribble_value;
            atoms[ai].coor[1] = ic2cf_test_hooks::g_scribble_value;
            atoms[ai].coor[2] = ic2cf_test_hooks::g_scribble_value;
        }
    }
    (void)FA;
    return true;
}
void dee_first(psFlexDEE_Node, psFlexDEE_Node) {}
void dee_last(psFlexDEE_Node, psFlexDEE_Node) {}
int  dee_pivot(psFlexDEE_Node, psFlexDEE_Node*, int, int, int, int, int) { return 0; }
