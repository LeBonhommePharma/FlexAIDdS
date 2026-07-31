// tests/test_read_lig_latm.cpp
// Regression: LIG.inp HETTYP load must set inclusive latm so the last ligand
// atom is in fatm..latm (write_pdb / buildlist / emission).
//
// Bug class (pilot8): latm = atm_cnt - 1 after atm_cnt++ dropped last atom
// (1P62 90017, 1T40 90027) from INI and ranked poses.
//
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>

#include "../LIB/flexaid.h"
#include "../LIB/fileio.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>

namespace fs = std::filesystem;

static void init_fa(FA_Global* FA, atom** atoms, resid** residue) {
    std::memset(FA, 0, sizeof(FA_Global));
    FA->MIN_NUM_ATOM = 100;
    FA->MIN_NUM_RESIDUE = 10;
    FA->MIN_FLEX_BONDS = 8;
    FA->MIN_OPTRES = 1;
    FA->MIN_ROTAMER = 2;
    FA->atm_cnt = 0;
    FA->atm_cnt_real = 0;
    FA->res_cnt = 0;
    FA->num_het = 0;
    FA->num_het_atm = 0;
    FA->ntypes = 40;
    FA->num_atm = (int*)calloc(100000, sizeof(int));
    *atoms = (atom*)calloc(FA->MIN_NUM_ATOM, sizeof(atom));
    *residue = (resid*)calloc(FA->MIN_NUM_RESIDUE, sizeof(resid));
}

static void cleanup_fa(FA_Global* FA, atom* atoms, resid* residue) {
    for (int r = 0; r <= FA->res_cnt; ++r) {
        // read_lig allocates bonded / shortpath / shortflex per residue via
        // update_bonded / shortest_path / assign_shortflex. Nothing in LIB
        // released them, so the process exited dirty and LeakSanitizer
        // reported 193896 B in 2360 allocations on linux-gcc-asan and
        // linux-clang-asan. Freed before fatm/latm because the dimension the
        // allocators used is read back out of them.
        if (residue[r].fatm != nullptr && residue[r].latm != nullptr) {
            const int natm = residue[r].latm[0] - residue[r].fatm[0] + 1;
            free_bonded(&residue[r], natm);
            free_shortpath(&residue[r], natm);
            free_shortflex(&residue[r], natm);
        }
        free(residue[r].fatm);
        free(residue[r].latm);
        free(residue[r].bond);
        free(residue[r].gpa);
    }
    free(FA->optres);
    free(FA->num_atm);
    free(atoms);
    free(residue);
}

// Match ProcessLigand fixed-width HETTYP layout used by live pilot LIG.inp:
// HETTYP90001 3 C    m 900029000390004    1
static std::string write_min_lig(const fs::path& dir, int n_atoms) {
    fs::create_directories(dir);
    const auto inp = dir / "LIG.inp";
    const auto ic = dir / "LIG.ic";
    {
        std::ofstream o(inp);
        o << "RESIDU LIG   9999\n";
        for (int i = 0; i < n_atoms; ++i) {
            const int ser = 90000 + i;
            const int r0 = (i > 0) ? (90000 + i - 1) : 0;
            const int r1 = (i + 1 < n_atoms) ? (90000 + i + 1) : 0;
            const int r2 = 0;
            char line[160];
            // type column: two chars at buffer[11..12] — " 3" for type 3
            std::snprintf(line, sizeof(line),
                          "HETTYP%05d 3 C    m %5d%5d%5d    1\n",
                          ser, r0, r1, r2);
            o << line;
        }
        const int g0 = 90000;
        const int g1 = 90000 + (n_atoms > 1 ? 1 : 0);
        const int g2 = 90000 + (n_atoms > 2 ? 2 : (n_atoms > 1 ? 1 : 0));
        o << "GPATOM " << g0 << " " << g1 << " " << g2 << "\n";
        for (int i = 0; i < n_atoms; ++i) {
            const int ser = 90000 + i;
            // CONECT: fixed 5-digit fields, no spaces (classic FlexAID)
            char line[80];
            if (i == 0 && n_atoms == 1) {
                std::snprintf(line, sizeof(line), "CONECT%05d\n", ser);
            } else if (i == 0) {
                std::snprintf(line, sizeof(line), "CONECT%05d%05d\n", ser, ser + 1);
            } else if (i + 1 == n_atoms) {
                std::snprintf(line, sizeof(line), "CONECT%05d%05d\n", ser, ser - 1);
            } else {
                std::snprintf(line, sizeof(line), "CONECT%05d%05d%05d\n",
                              ser, ser - 1, ser + 1);
            }
            o << line;
        }
        o << "ENDINP\n";
    }
    {
        std::ofstream o(ic);
        // read_lig IC: serial in cols 0-4, values from col 7 (col 5 often ':')
        o << "REFPCG  0.000  0.000  0.000\n";
        for (int i = 0; i < n_atoms; ++i) {
            const int ser = 90000 + i;
            char line[96];
            std::snprintf(line, sizeof(line),
                          "%05d:    1.500  109.500  180.000\n", ser);
            o << line;
        }
    }
    return inp.string();
}

TEST(ReadLigLatm, InclusiveLatmIncludesLastHettyp) {
    // N=18 matches 1P62 heavy count; last serial 90017 was dropped by the bug.
    constexpr int N = 18;
    const auto tmp = fs::temp_directory_path() / "flexaids_read_lig_latm";
    fs::remove_all(tmp);
    const std::string ligpath = write_min_lig(tmp, N);

    FA_Global FA{};
    atom* atoms = nullptr;
    resid* residue = nullptr;
    init_fa(&FA, &atoms, &residue);

    ASSERT_NO_THROW({
        read_lig(&FA, &atoms, &residue, const_cast<char*>(ligpath.c_str()));
    });

    ASSERT_GE(FA.res_cnt, 1);
    const int r = FA.res_cnt;
    ASSERT_NE(residue[r].fatm, nullptr);
    ASSERT_NE(residue[r].latm, nullptr);
    const int fatm = residue[r].fatm[0];
    const int latm = residue[r].latm[0];
    const int n_in_range = latm - fatm + 1;

    EXPECT_EQ(n_in_range, N) << "fatm=" << fatm << " latm=" << latm
                             << " (bug was latm=atm_cnt-1 → N-1 atoms)";
    EXPECT_EQ(FA.num_het_atm, N);

    const int last_serial = 90000 + N - 1;
    bool found_last = false;
    for (int i = fatm; i <= latm; ++i) {
        if (atoms[i].number == last_serial) {
            found_last = true;
            break;
        }
    }
    EXPECT_TRUE(found_last) << "last serial " << last_serial
                            << " missing from fatm..latm emission range";

    cleanup_fa(&FA, atoms, residue);
    fs::remove_all(tmp);
}

TEST(ReadLigLatm, OneT40Style28Atoms) {
    constexpr int N = 28;  // 1T40 heavy count; last serial 90027
    const auto tmp = fs::temp_directory_path() / "flexaids_read_lig_latm_28";
    fs::remove_all(tmp);
    const std::string ligpath = write_min_lig(tmp, N);

    FA_Global FA{};
    atom* atoms = nullptr;
    resid* residue = nullptr;
    init_fa(&FA, &atoms, &residue);

    ASSERT_NO_THROW({
        read_lig(&FA, &atoms, &residue, const_cast<char*>(ligpath.c_str()));
    });

    const int r = FA.res_cnt;
    const int n_in_range = residue[r].latm[0] - residue[r].fatm[0] + 1;
    EXPECT_EQ(n_in_range, N);
    EXPECT_EQ(FA.num_het_atm, N);

    cleanup_fa(&FA, atoms, residue);
    fs::remove_all(tmp);
}
