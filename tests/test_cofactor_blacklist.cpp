// =============================================================================
// test_cofactor_blacklist.cpp — GTest suite for DatasetRunner cofactor blacklist
//
// Regression guard for the 1HNN bug: the PDB co-crystallises the cofactor SAH
// (S-adenosyl-L-homocysteine) alongside the real PNMT inhibitor SKF.  The
// largest-HETATM heuristic used to extract SAH as the "ligand", producing a
// ~22.5 Å self-dock RMSD because the wrong molecule was docked.
//
// These tests assert that DatasetRunner::extract_ligand() skips blacklisted
// biochemical cofactors and falls through to the genuine docking target.
//
// This file is intentionally separate from test_dataset_runner.cpp so it owns
// its own executable and does not entangle with unrelated in-flight work there.
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
// =============================================================================

#include <gtest/gtest.h>
#include "DatasetRunner.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

namespace fs = std::filesystem;
using dataset::DatasetRunner;

// Path to the repository source tree, injected by CMake so the real 1HNN
// benchmark structure can be located regardless of the build directory.
#ifndef FLEXAIDS_SOURCE_DIR
#define FLEXAIDS_SOURCE_DIR "."
#endif

namespace {

// Read the first line (molecule title = residue name) of an SDF file.
std::string sdf_title(const std::string& sdf_path) {
    std::ifstream ifs(sdf_path);
    std::string line;
    std::getline(ifs, line);
    return line;
}

// Read the atom count from the SDF V2000 counts line (line 4, columns 0-2).
int sdf_atom_count(const std::string& sdf_path) {
    std::ifstream ifs(sdf_path);
    std::string line;
    for (int i = 0; i < 4 && std::getline(ifs, line); ++i) { /* advance */ }
    try { return std::stoi(line.substr(0, 3)); } catch (...) { return -1; }
}

}  // namespace

// =============================================================================
// Real-structure regression: 1HNN must extract SKF, never the SAH cofactor.
// =============================================================================

TEST(CofactorBlacklist, Hnn1ExtractsSkfNotSah) {
    const std::string cif =
        std::string(FLEXAIDS_SOURCE_DIR) +
        "/benchmarks/astex_diverse/astex_diverse/1HNN/1HNN.cif";

    if (!fs::exists(cif)) {
        GTEST_SKIP() << "1HNN.cif not present in source tree: " << cif;
    }

    const std::string out_sdf =
        (fs::temp_directory_path() / "flexaidds_1hnn_blacklist.sdf").string();
    fs::remove(out_sdf);

    DatasetRunner runner((fs::temp_directory_path() / "flexaidds_blacklist_cache").string());
    const bool ok = runner.extract_ligand(cif, out_sdf);

    ASSERT_TRUE(ok) << "extract_ligand failed on 1HNN.cif";
    ASSERT_TRUE(fs::exists(out_sdf));

    const std::string title = sdf_title(out_sdf);
    EXPECT_EQ(title, "SKF")
        << "1HNN should dock the PNMT inhibitor SKF, got '" << title << "'";
    EXPECT_NE(title, "SAH")
        << "SAH is a biochemical cofactor and must never be the docking target";

    // 1HNN is a homodimer: SKF appears in two chains at 14 heavy atoms each
    // (28 total), SAH at 26 atoms per copy. extract_ligand keeps one molecule
    // (largest connected component), so the SKF ligand is exactly 14 atoms.
    // Guard against any silent regression back to the SAH cofactor.
    const int natoms = sdf_atom_count(out_sdf);
    EXPECT_EQ(natoms, 14) << "expected one 14-atom SKF molecule";
    EXPECT_NE(natoms, 26) << "atom count matches SAH — cofactor leaked through";

    fs::remove(out_sdf);
}

// =============================================================================
// Synthetic: a blacklisted cofactor (SAH) larger than the true ligand (LIG)
// must be skipped in favour of LIG.
// =============================================================================

TEST(CofactorBlacklist, SkipsLargerCofactorForRealLigand) {
    const std::string dir = (fs::temp_directory_path() / "flexaidds_blacklist_syn").string();
    fs::create_directories(dir);
    const std::string pdb = dir + "/syn.pdb";
    const std::string sdf = dir + "/syn_ligand.sdf";
    fs::remove(sdf);

    {
        std::ofstream ofs(pdb);
        ofs << "HEADER    SYNTHETIC\n";
        // Protein (ignored by HETATM extraction)
        ofs << "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 10.00           N\n";
        // Blacklisted cofactor SAH — 7 atoms, the LARGEST HETATM group.
        ofs << "HETATM    2  N1  SAH B   1      10.000  10.000  10.000  1.00 20.00           N\n";
        ofs << "HETATM    3  C1  SAH B   1      11.000  10.000  10.000  1.00 20.00           C\n";
        ofs << "HETATM    4  C2  SAH B   1      12.000  10.000  10.000  1.00 20.00           C\n";
        ofs << "HETATM    5  O1  SAH B   1      13.000  10.000  10.000  1.00 20.00           O\n";
        ofs << "HETATM    6  C3  SAH B   1      14.000  10.000  10.000  1.00 20.00           C\n";
        ofs << "HETATM    7  S1  SAH B   1      15.000  10.000  10.000  1.00 20.00           S\n";
        ofs << "HETATM    8  N2  SAH B   1      16.000  10.000  10.000  1.00 20.00           N\n";
        // Real ligand LIG — 5 connected atoms (~1.5 Å spacing).
        ofs << "HETATM    9  C1  LIG C   1       0.000   0.000  20.000  1.00 20.00           C\n";
        ofs << "HETATM   10  C2  LIG C   1       1.500   0.000  20.000  1.00 20.00           C\n";
        ofs << "HETATM   11  N1  LIG C   1       3.000   0.000  20.000  1.00 20.00           N\n";
        ofs << "HETATM   12  O1  LIG C   1       4.500   0.000  20.000  1.00 20.00           O\n";
        ofs << "HETATM   13  C3  LIG C   1       6.000   0.000  20.000  1.00 20.00           C\n";
        // Water + ion (excluded)
        ofs << "HETATM   14  O   HOH D   1      30.000  30.000  30.000  1.00  5.00           O\n";
        ofs << "HETATM   15 ZN   ZN  E   1      35.000  35.000  35.000  1.00  5.00          ZN\n";
        ofs << "END\n";
    }

    DatasetRunner runner(dir + "/cache");
    const bool ok = runner.extract_ligand(pdb, sdf);

    ASSERT_TRUE(ok) << "extract_ligand should succeed with a non-cofactor ligand present";
    EXPECT_EQ(sdf_title(sdf), "LIG")
        << "the blacklisted SAH cofactor must be skipped in favour of LIG";

    fs::remove_all(dir);
}

// =============================================================================
// Synthetic: when EVERY non-water HETATM is a blacklisted cofactor, extraction
// must fail rather than return a cofactor.
// =============================================================================

TEST(CofactorBlacklist, FailsWhenOnlyCofactorsPresent) {
    const std::string dir = (fs::temp_directory_path() / "flexaidds_blacklist_only").string();
    fs::create_directories(dir);
    const std::string pdb = dir + "/only.pdb";
    const std::string sdf = dir + "/only_ligand.sdf";
    fs::remove(sdf);

    {
        std::ofstream ofs(pdb);
        ofs << "HETATM    1  PA  ATP A   1       0.000   0.000   0.000  1.00 20.00           P\n";
        ofs << "HETATM    2  PB  ATP A   1       1.500   0.000   0.000  1.00 20.00           P\n";
        ofs << "HETATM    3  C1  ATP A   1       3.000   0.000   0.000  1.00 20.00           C\n";
        ofs << "HETATM    4  N1  ATP A   1       4.500   0.000   0.000  1.00 20.00           N\n";
        ofs << "HETATM    5  O   HOH B   1      30.000  30.000  30.000  1.00  5.00           O\n";
        ofs << "END\n";
    }

    DatasetRunner runner(dir + "/cache");
    const bool ok = runner.extract_ligand(pdb, sdf);

    EXPECT_FALSE(ok)
        << "ATP is a blacklisted cofactor; with no real ligand extraction must fail";

    fs::remove_all(dir);
}
