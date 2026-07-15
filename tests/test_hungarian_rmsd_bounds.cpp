// Regression: calc_Hungarian_RMSD must not walk past ligand atom bounds.
//
// Historical bug (cad_only_no_pose_pdb): after CF clustering wrote *.cad, pose
// emission called calc_rmsd(..., Hungarian=true). The old cost-matrix fill used
// atoms[k + l] for l in [0, num_het_atm), reading protein/unmapped memory and
// often NULL coor_ref → SIGSEGV before any ranked *_N.pdb was written.
//
// This unit test builds a minimal FA/atom/residue layout that would have
// crashed the old walker (ligand followed by protein atoms without coor_ref)
// and asserts the fixed Hungarian path returns a finite RMSD.

#include <gtest/gtest.h>

// flexaid.h defines E (Euler) as a macro — must come after gtest headers.
#include "flexaid.h"

#include <cmath>
#include <cstring>
#include <vector>

// Defined in LIB/calc_rmsd.cpp
float calc_Hungarian_RMSD(FA_Global* FA, atom* atoms, resid* residue,
                          gridpoint* cleftgrid, int npar, const double* icv);

namespace {

struct MiniSystem {
    FA_Global FA{};
    std::vector<atom> atoms;
    std::vector<resid> residue;
    // fatm/latm are pointers in resid — own storage here.
    int prot_fatm[1]{1};
    int prot_latm[1]{4};
    int lig_fatm[1]{5};
    int lig_latm[1]{8};
    float ref_store[4][3]{};

    MiniSystem() {
        // indices: 0 unused (1-based FlexAID convention), 1..4 protein, 5..8 ligand
        atoms.assign(9, atom{});
        residue.assign(3, resid{});  // 0 unused, 1 protein, 2 ligand

        std::memset(static_cast<void*>(&FA), 0, sizeof(FA));
        FA.num_het = 1;
        FA.num_het_atm = 4;
        FA.res_cnt = 2;
        FA.het_res[1] = 2;

        // Protein residue (no coor_ref) — old k+l walker would enter here.
        residue[1].fatm = prot_fatm;
        residue[1].latm = prot_latm;
        residue[1].type = 0;
        for (int i = 1; i <= 4; ++i) {
            atoms[i].type = 4;  // same type as ligand carbons → would match
            atoms[i].coor[0] = 100.0f + static_cast<float>(i);
            atoms[i].coor[1] = 0.0f;
            atoms[i].coor[2] = 0.0f;
            atoms[i].coor_ref = nullptr;
        }

        // Ligand residue with reference coords
        residue[2].fatm = lig_fatm;
        residue[2].latm = lig_latm;
        residue[2].type = 1;
        for (int i = 0; i < 4; ++i) {
            const int ai = 5 + i;
            atoms[ai].type = 4;
            atoms[ai].coor[0] = static_cast<float>(i);
            atoms[ai].coor[1] = 0.0f;
            atoms[ai].coor[2] = 0.0f;
            ref_store[i][0] = static_cast<float>(i) + 0.1f;
            ref_store[i][1] = 0.0f;
            ref_store[i][2] = 0.0f;
            atoms[ai].coor_ref = ref_store[i];
        }
    }
};

}  // namespace

TEST(HungarianRmsdBounds, DoesNotSegfaultWhenProteinFollowsLigand) {
    MiniSystem sys;
    const float rmsd = calc_Hungarian_RMSD(&sys.FA, sys.atoms.data(),
                                           sys.residue.data(), nullptr, 0,
                                           nullptr);
    EXPECT_TRUE(std::isfinite(rmsd));
    EXPECT_GE(rmsd, 0.0f);
}

TEST(HungarianRmsdBounds, NullRefAtomsSkipped) {
    MiniSystem sys;
    // Drop coor_ref on one ligand atom — must not crash.
    sys.atoms[6].coor_ref = nullptr;
    const float rmsd = calc_Hungarian_RMSD(&sys.FA, sys.atoms.data(),
                                           sys.residue.data(), nullptr, 0,
                                           nullptr);
    EXPECT_TRUE(std::isfinite(rmsd));
    EXPECT_GE(rmsd, 0.0f);
}

TEST(HungarianRmsdBounds, EmptySystemSafe) {
    FA_Global FA{};
    FA.num_het_atm = 0;
    const float rmsd = calc_Hungarian_RMSD(&FA, nullptr, nullptr, nullptr, 0,
                                           nullptr);
    EXPECT_FLOAT_EQ(rmsd, 0.0f);
}
