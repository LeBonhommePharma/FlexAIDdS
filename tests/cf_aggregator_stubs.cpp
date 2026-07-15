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

void alter_mode(atom*, resid*, float*, int, int) {}
double vcfunction(FA_Global*, VC_Global*, atom*, resid*,
                  std::vector<std::pair<int,int>>&, bool*) { return 0.0; }
void buildcc(FA_Global*, atom*, int, int[]) {}
void dee_first(psFlexDEE_Node, psFlexDEE_Node) {}
void dee_last(psFlexDEE_Node, psFlexDEE_Node) {}
int  dee_pivot(psFlexDEE_Node, psFlexDEE_Node*, int, int, int, int, int) { return 0; }
