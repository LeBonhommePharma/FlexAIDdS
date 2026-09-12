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
#include <stdexcept>
#include "tENCoM/tencm.h"
#include "encom.h"        // canonical declaration of ENCoMEngine/VibrationalEntropy
#include "VibEntropy.h"   // canonical declaration; see the ODR note below

namespace vibentropy {
// SIGNATURE MUST MATCH VibEntropy.h EXACTLY. This stub returned `double` while
// the canonical declaration at VibEntropy.h:83 returns VibEntropyResult, and the
// production caller (ic2cf.cpp:50) reads `.H_pop` off it. Because the two live in
// separate translation units the compiler could not see the conflict: it is a
// silent ODR violation that links, and would return garbage the moment this path
// became reachable in the fixture. The canonical header is now included so any
// future drift is a COMPILE error instead.
//
// Throws rather than returning a zeroed result: this fixture deliberately does
// not exercise the vibrational-entropy path (its tENCoM model builder is a
// no-op), so reaching this is a test-construction error that should be loud, not
// a silently plausible 0.0. Found by external audit.
VibEntropyResult compute_vib_entropy_collapse(const std::vector<std::vector<double>>&) {
    throw std::logic_error(
        "cf_aggregator fixture: compute_vib_entropy_collapse is not exercised by "
        "these tests; if you reached it, wire a real model or extend the stub");
}
}
namespace sugar_pucker {
enum class SugarType;
void apply_sugar_puckers(atom*, const std::vector<std::vector<int>>&,
                         const std::vector<float>&, const std::vector<SugarType>&) {}
}
namespace tencm {
void TorsionalENM::build_from_ligand(const atom*, int, int, float, float) {}
std::vector<double> TorsionalENM::vibrational_eigenvalues(int*, int*) const {
    throw std::logic_error("CF aggregator fixture must not execute tENCoM");
}
// The thermodynamic dS_vib cycle added to ic2cf.cpp in 4e1043b6 made these two
// torsional builders new collaborators of this target. Stubbed here rather than
// linked, per this target's hermetic contract (CMakeLists.txt: "If ic2cf.cpp
// acquires a NEW dependency, stub it") -- linking LIB/tENCoM/tencm.cpp is the
// mistake the CMake comment records twice, because tencm.cpp needs Eigen that
// this target deliberately has no include path for.
//
// No-ops, matching build_from_ligand above: builders only populate internal
// state. The EVALUATORS stay loud -- vibrational_eigenvalues() above and
// ENCoMEngine::compute_vibrational_entropy() below both throw, so a fixture that
// actually reaches the vibrational path fails audibly instead of scoring a
// silently plausible 0.0.
void TorsionalENM::build_from_ligand_torsional(const std::array<float,3>*, int,
                                               const LigandTopology&,
                                               float, float) {}
void TorsionalENM::build_from_ligand_torsional_in_field(const std::array<float,3>*, int,
                                                        const LigandTopology&,
                                                        const std::array<float,3>*, int,
                                                        float, float) {}
}

namespace encom {
// Throws for the same reason vibentropy::compute_vib_entropy_collapse does: this
// fixture does not exercise the vibrational-entropy path (its tENCoM builders are
// no-ops, so any mode set reaching here is empty and meaningless). A zeroed
// VibrationalEntropy would be indistinguishable from a real result and would let
// get_cf_evalue() be pinned against a term that silently contributed nothing.
VibrationalEntropy ENCoMEngine::compute_vibrational_entropy(
        const std::vector<NormalMode>&, double, double) {
    throw std::logic_error(
        "cf_aggregator fixture: ENCoMEngine::compute_vibrational_entropy is not "
        "exercised by these tests; if you reached it, wire a real model or extend "
        "the stub");
}
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

// ic2cf.cpp calls rot_gene_index on the flexible side-chain path
// (94d2d0ba). This target links the real ic2cf.cpp, not rot_gene_index.cpp;
// get_cf_evalue() never reaches the clamp, so a no-op stub is enough to
// link. Identity on in-range values; clamp to [0, res->trot] otherwise.
int rot_gene_index(double gene, const resid* res, const char* site) {
    (void)site;
    const int trot = (res != NULL) ? res->trot : 0;
    int idx = (int)(gene + 0.5);
    if (idx < 0) idx = 0;
    if (idx > trot) idx = trot;
    return idx;
}
long rot_gene_guard_hits(void) { return 0; }
