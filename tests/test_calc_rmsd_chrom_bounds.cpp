// test_calc_rmsd_chrom_bounds.cpp — OOB gene indices must not crash coord cache
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>
#include <cstring>
#include <vector>

#include "gaboom.h"

namespace {

struct LigandFixtureStorage {
    int fatm = 0;
    int latm = 0;
};

void init_minimal_ligand_fixture(FA_Global& fa, GB_Global& gb, optmap& map0,
                                 std::vector<atom>& atoms, std::vector<resid>& residue,
                                 std::vector<gridpoint>& grid, genlim& lim,
                                 chromosome& chrom, LigandFixtureStorage& store) {
    std::memset(&fa, 0, sizeof(fa));
    std::memset(&gb, 0, sizeof(gb));

    fa.npar = 1;
    fa.num_grd = 2;
    fa.nors = 0;
    fa.atm_cnt = 1;
    fa.res_cnt = 0;
    fa.map_par = &map0;
    map0.typ = -1;
    map0.atm = 0;

    atoms.resize(2);
    std::memset(atoms.data(), 0, atoms.size() * sizeof(atom));
    atoms[0].coor[0] = 1.0f;
    atoms[0].coor[1] = 2.0f;
    atoms[0].coor[2] = 3.0f;
    atoms[0].ofres = 0;

    store.fatm = 0;
    store.latm = 0;
    residue.resize(2);
    std::memset(residue.data(), 0, residue.size() * sizeof(resid));
    residue[0].trot = 1;
    residue[0].fatm = &store.fatm;
    residue[0].latm = &store.latm;

    grid.resize(2);
    grid[0].dis = 0.0f;
    grid[0].ang = 0.0f;
    grid[0].dih = 0.0f;
    grid[1].dis = 1.0f;
    grid[1].ang = 0.0f;
    grid[1].dih = 0.0f;

    gb.num_genes = 1;
    lim.min = 0.0;
    lim.max = 1.0;
    lim.del = 1.0;
    lim.map = -1;

    chrom.genes = new gene[1];
    chrom.genes[0].to_ic = 1.0e9;  // far outside gene_lim → must clamp
    chrom.cf = {};
}

}  // namespace

TEST(CalcRmsdChromBounds, OutOfBoundsGridGeneDoesNotCrash) {
    FA_Global fa{};
    GB_Global gb{};
    optmap map0{};
    std::vector<atom> atoms;
    std::vector<resid> residue;
    std::vector<gridpoint> grid;
    genlim lim{};
    chromosome chrom{};
    LigandFixtureStorage store{};
    init_minimal_ligand_fixture(fa, gb, map0, atoms, residue, grid, lim, chrom, store);

    const int stride = 3;
    std::vector<float> coor(static_cast<size_t>(stride), 0.0f);
    int n_atoms = 0;

    EXPECT_NO_THROW({
        calc_rmsd_chrom(&fa, &gb, &chrom, &lim, atoms.data(), residue.data(),
                        grid.data(), gb.num_genes, 0, 0, coor.data(), nullptr,
                        false, &n_atoms);
    });
    EXPECT_GE(n_atoms, 1);
    EXPECT_LE(n_atoms, 1);

    delete[] chrom.genes;
}

TEST(CalcRmsdChromBounds, FuzzRandomGeneOverflow) {
    FA_Global fa{};
    GB_Global gb{};
    optmap map0{};
    std::vector<atom> atoms;
    std::vector<resid> residue;
    std::vector<gridpoint> grid;
    genlim lim{};
    chromosome chrom{};
    LigandFixtureStorage store{};
    init_minimal_ligand_fixture(fa, gb, map0, atoms, residue, grid, lim, chrom, store);

    const int stride = 3;
    std::vector<float> coor(static_cast<size_t>(stride), 0.0f);

    for (int trial = 0; trial < 64; ++trial) {
        chrom.genes[0].to_ic = (trial % 2 == 0) ? -1.0e6 : 1.0e6;
        int n_atoms = 0;
        EXPECT_NO_THROW({
            calc_rmsd_chrom(&fa, &gb, &chrom, &lim, atoms.data(), residue.data(),
                            grid.data(), gb.num_genes, 0, 0, coor.data(), nullptr,
                            false, &n_atoms);
        });
        EXPECT_GE(n_atoms, 0);
        EXPECT_LE(n_atoms, 1);
    }

    delete[] chrom.genes;
}