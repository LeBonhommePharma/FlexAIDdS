// tests/test_tencom_output.cpp
// tENCoM FlexPopulation mode↔structure pairing after sort_by_free_energy.
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>

#include "tENCoM/tencom_output.h"
#include "tENCoM/pdb_calpha.h"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace {

std::string write_synthetic_pdb(const std::string& prefix, int n_residues)
{
    std::string path = std::filesystem::temp_directory_path().string()
                       + "/" + prefix + ".pdb";
    std::ofstream ofs(path);
    const float turn = 100.0f * 3.14159265f / 180.0f;
    const float radius = 2.3f;
    const float rise = 1.5f;
    for (int r = 0; r < n_residues; ++r) {
        float x = radius * std::cos(r * turn);
        float y = radius * std::sin(r * turn);
        float z = r * rise;
        char line[82];
        std::snprintf(line, sizeof(line),
            "ATOM  %5d  CA  ALA A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C",
            r + 1, r + 1, x, y, z);
        ofs << line << "\n";
    }
    ofs << "END\n";
    return path;
}

}  // namespace

TEST(TencomOutput, SortPreservesStructurePairingWhenNgt2)
{
    using tencom_output::FlexMode;
    using tencom_output::FlexPopulation;

    FlexPopulation pop;
    std::vector<tencom_pdb::CalphaStructure> structures;
    const int nres[] = {8, 12, 16, 20};
    const double dF[] = {0.0, 3.0, 1.0, 2.0};
    std::vector<std::string> paths;

    for (int i = 0; i < 4; ++i) {
        auto path = write_synthetic_pdb("tencom_pair_" + std::to_string(i), nres[i]);
        paths.push_back(path);
        auto st = tencom_pdb::read_pdb_calpha(path);
        FlexMode m;
        m.mode_id = i;
        m.structure_index = static_cast<std::size_t>(i);
        m.pdb_path = path;
        m.n_residues = st.res_cnt;
        m.delta_F_vib = dF[i];
        pop.modes.push_back(std::move(m));
        structures.push_back(std::move(st));
    }

    for (std::size_t i = 0; i < pop.modes.size(); ++i) {
        EXPECT_EQ(structures[i].res_cnt, pop.modes[i].n_residues);
    }

    pop.sort_by_free_energy();

    ASSERT_EQ(pop.modes.size(), 4u);
    EXPECT_EQ(pop.modes[0].structure_index, 0u);
    EXPECT_EQ(pop.modes[1].structure_index, 2u);  // ΔF=1 was original index 2
    EXPECT_EQ(pop.modes[2].structure_index, 3u);  // ΔF=2 was original index 3
    EXPECT_EQ(pop.modes[3].structure_index, 1u);  // ΔF=3 was original index 1

    bool old_index_pairing_broken = false;
    for (std::size_t i = 0; i < pop.modes.size(); ++i) {
        const auto* st = FlexPopulation::paired_structure(pop.modes[i], structures);
        ASSERT_NE(st, nullptr);
        EXPECT_EQ(st->res_cnt, pop.modes[i].n_residues);
        EXPECT_EQ(st->filename, pop.modes[i].pdb_path);
        if (structures[i].res_cnt != pop.modes[i].n_residues) {
            old_index_pairing_broken = true;
        }
    }
    EXPECT_TRUE(old_index_pairing_broken)
        << "this test is vacuous unless sort permutes modes relative to structures";

    for (const auto& p : paths) {
        std::remove(p.c_str());
    }
}
