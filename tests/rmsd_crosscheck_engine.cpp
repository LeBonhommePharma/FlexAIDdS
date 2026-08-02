// Engine-side half of the RMSD cross-check. Isolated TU: this is the only file
// that includes flexaid.h (which #defines the macro E), keeping it out of the
// gtest translation unit that also uses std headers named E-adjacent.
#include "rmsd_crosscheck_engine.h"

#include "flexaid.h"

#include <vector>

// Defined in LIB/calc_rmsd.cpp (linked via flexaid_core).
float calc_Hungarian_RMSD(FA_Global* FA, atom* atoms, resid* residue,
                          gridpoint* cleftgrid, int npar, const double* icv);

namespace crosscheck {

float engine_hungarian_rmsd(const std::vector<AtomSpec>& specs) {
    const int n = static_cast<int>(specs.size());
    if (n <= 0) return 0.0f;

    // Value-init: POD members zeroed, and the two std::vector members
    // (model_coords / model_strain) properly constructed empty. Do NOT memset
    // over this — that overwrites the vectors' internal pointers and is UB
    // (Honey, #371 review).
    FA_Global FA{};
    FA.num_het = 1;
    FA.num_het_atm = n;
    FA.res_cnt = 1;
    FA.het_res[1] = 1;  // ligand is residue index 1

    // 1-based atom indexing (FlexAID convention): atoms[1..n] are the ligand.
    std::vector<atom> atoms(static_cast<size_t>(n) + 1);
    std::vector<std::array<float, 3>> refstore(specs.size());
    int fatm_arr[1] = {1};
    int latm_arr[1] = {n};

    std::vector<resid> residue(2);  // 0 unused, 1 ligand
    residue[1].fatm = fatm_arr;
    residue[1].latm = latm_arr;
    residue[1].type = 1;

    for (int i = 1; i <= n; ++i) {
        const AtomSpec& s = specs[static_cast<size_t>(i - 1)];
        atoms[i].type = s.type;
        atoms[i].coor[0] = s.coor[0];
        atoms[i].coor[1] = s.coor[1];
        atoms[i].coor[2] = s.coor[2];
        refstore[i - 1] = s.coor_ref;
        atoms[i].coor_ref = refstore[i - 1].data();
    }

    return calc_Hungarian_RMSD(&FA, atoms.data(), residue.data(), nullptr, 0, nullptr);
}

}  // namespace crosscheck
