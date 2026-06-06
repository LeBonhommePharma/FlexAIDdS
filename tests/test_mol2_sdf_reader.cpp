// tests/test_mol2_sdf_reader.cpp
// Unit tests for Mol2Reader and SdfReader — ligand file parsing
// Apache-2.0 © 2026 NRGlab, Université de Montréal

#include <gtest/gtest.h>
#include "../LIB/Mol2Reader.h"
#include "../LIB/SdfReader.h"
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <filesystem>
#include <string>

// ===========================================================================
// HELPER: Initialize FA_Global with minimum required fields
// ===========================================================================

static void init_fa_for_reader(FA_Global* FA, atom** atoms, resid** residue) {
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wnontrivial-memcall"
    std::memset(FA, 0, sizeof(FA_Global));
    #pragma clang diagnostic pop
    FA->MIN_NUM_ATOM     = 100;
    FA->MIN_NUM_RESIDUE  = 10;
    FA->MIN_FLEX_BONDS   = 5;
    FA->MIN_OPTRES       = 1;
    FA->atm_cnt          = 0;
    FA->atm_cnt_real     = 0;
    FA->res_cnt          = 0;
    FA->num_het          = 0;
    FA->num_het_atm      = 0;

    // PDB num → internal index mapping (same as read_pdb allocates)
    FA->num_atm = (int*)calloc(100000, sizeof(int));

    *atoms   = (atom*)calloc(FA->MIN_NUM_ATOM, sizeof(atom));
    *residue = (resid*)calloc(FA->MIN_NUM_RESIDUE, sizeof(resid));
}

static void cleanup_fa(FA_Global* FA, atom* atoms, resid* residue) {
    // Free residue sub-allocations created by the readers
    for (int r = 1; r <= FA->res_cnt; ++r) {
        free(residue[r].fatm);
        free(residue[r].latm);
        free(residue[r].bond);
    }
    free(FA->optres);
    free(FA->num_atm);
    free(atoms);
    free(residue);
}

// ===========================================================================
// MOL2 READER TESTS
// ===========================================================================

class Mol2ReaderTest : public ::testing::Test {
protected:
    std::string tmp_dir;

    void SetUp() override {
        tmp_dir = std::filesystem::temp_directory_path().string();
    }

    std::string write_mol2(const std::string& name, const std::string& content) {
        std::string path = tmp_dir + "/" + name;
        std::ofstream ofs(path);
        ofs << content;
        return path;
    }
};

TEST_F(Mol2ReaderTest, ReadsSimpleMolecule) {
    // Minimal water molecule in MOL2 format
    std::string mol2 = write_mol2("water.mol2",
        "@<TRIPOS>MOLECULE\n"
        "WAT\n"
        "3 2\n"
        "SMALL\n"
        "\n"
        "@<TRIPOS>ATOM\n"
        "1 O1   0.000  0.000  0.000 O.3   1 WAT  -0.834\n"
        "2 H1   0.957  0.000  0.000 H     1 WAT   0.417\n"
        "3 H2  -0.240  0.927  0.000 H     1 WAT   0.417\n"
        "@<TRIPOS>BOND\n"
        "1 1 2 1\n"
        "2 1 3 1\n"
    );

    FA_Global FA;
    atom* atoms = nullptr;
    resid* residue = nullptr;
    init_fa_for_reader(&FA, &atoms, &residue);

    int ok = read_mol2_ligand(&FA, &atoms, &residue, mol2.c_str());
    EXPECT_EQ(ok, 1);

    // Should have 3 atoms
    EXPECT_EQ(FA.num_het_atm, 3);
    EXPECT_EQ(FA.res_cnt, 1);

    // Atoms are stored 0-based (atoms[0] is the first atom; fatm/latm match).
    // Check coordinates of first atom (oxygen)
    EXPECT_NEAR(atoms[0].coor[0], 0.0f, 0.01f);
    EXPECT_NEAR(atoms[0].coor[1], 0.0f, 0.01f);
    EXPECT_NEAR(atoms[0].coor[2], 0.0f, 0.01f);

    // Check second atom coordinates (H1)
    EXPECT_NEAR(atoms[1].coor[0], 0.957f, 0.01f);

    // Check types: O.3 → type 10, H → type 22
    EXPECT_EQ(atoms[0].type, 10);
    EXPECT_EQ(atoms[1].type, 22);
    EXPECT_EQ(atoms[2].type, 22);

    // Check radii
    EXPECT_NEAR(atoms[0].radius, 1.52f, 0.01f);  // oxygen
    EXPECT_NEAR(atoms[1].radius, 1.20f, 0.01f);  // hydrogen

    // Check partial charges
    EXPECT_NEAR(atoms[0].charge, -0.834f, 0.01f);
    EXPECT_NEAR(atoms[1].charge,  0.417f, 0.01f);

    // Check bonds: O should have 2 bonds, each H should have 1
    EXPECT_EQ(atoms[0].bond[0], 2);
    EXPECT_EQ(atoms[1].bond[0], 1);
    EXPECT_EQ(atoms[2].bond[0], 1);

    // Residue should be set up as ligand
    EXPECT_EQ(residue[1].type, 1);
    EXPECT_EQ(FA.resligand, &residue[1]);

    cleanup_fa(&FA, atoms, residue);
    std::remove(mol2.c_str());
}

TEST_F(Mol2ReaderTest, ReadsDrugLikeMolecule) {
    // Aspirin-like structure (simplified) with aromatic and double bonds
    std::string mol2 = write_mol2("aspirin.mol2",
        "@<TRIPOS>MOLECULE\n"
        "ASP\n"
        "5 4\n"
        "SMALL\n"
        "\n"
        "@<TRIPOS>ATOM\n"
        "1 C1   0.000  0.000  0.000 C.ar  1 ASP  0.0\n"
        "2 C2   1.400  0.000  0.000 C.2   1 ASP  0.0\n"
        "3 O1   2.100  1.000  0.000 O.2   1 ASP -0.5\n"
        "4 O2   2.100 -1.000  0.000 O.3   1 ASP -0.3\n"
        "5 N1   0.000  1.400  0.000 N.am  1 ASP -0.2\n"
        "@<TRIPOS>BOND\n"
        "1 1 2 ar\n"
        "2 2 3 2\n"
        "3 2 4 1\n"
        "4 1 5 1\n"
    );

    FA_Global FA;
    atom* atoms = nullptr;
    resid* residue = nullptr;
    init_fa_for_reader(&FA, &atoms, &residue);

    int ok = read_mol2_ligand(&FA, &atoms, &residue, mol2.c_str());
    EXPECT_EQ(ok, 1);
    EXPECT_EQ(FA.num_het_atm, 5);

    // C.ar → type 3, C.2 → type 2, O.2 → type 11, O.3 → type 10, N.am → type 7
    EXPECT_EQ(atoms[0].type, 3);
    EXPECT_EQ(atoms[1].type, 2);
    EXPECT_EQ(atoms[2].type, 11);
    EXPECT_EQ(atoms[3].type, 10);
    EXPECT_EQ(atoms[4].type, 7);

    // C2 should have 3 bonds (to C1, O1, O2)
    EXPECT_EQ(atoms[1].bond[0], 3);

    cleanup_fa(&FA, atoms, residue);
    std::remove(mol2.c_str());
}

TEST_F(Mol2ReaderTest, FailsOnMissingFile) {
    FA_Global FA;
    atom* atoms = nullptr;
    resid* residue = nullptr;
    init_fa_for_reader(&FA, &atoms, &residue);

    int ok = read_mol2_ligand(&FA, &atoms, &residue, "/nonexistent/file.mol2");
    EXPECT_EQ(ok, 0);

    free(FA.num_atm);
    free(atoms);
    free(residue);
}

TEST_F(Mol2ReaderTest, FailsOnEmptyAtomBlock) {
    std::string mol2 = write_mol2("empty.mol2",
        "@<TRIPOS>MOLECULE\n"
        "EMPTY\n"
        "0 0\n"
        "SMALL\n"
        "\n"
        "@<TRIPOS>ATOM\n"
        "@<TRIPOS>BOND\n"
    );

    FA_Global FA;
    atom* atoms = nullptr;
    resid* residue = nullptr;
    init_fa_for_reader(&FA, &atoms, &residue);

    int ok = read_mol2_ligand(&FA, &atoms, &residue, mol2.c_str());
    EXPECT_EQ(ok, 0);

    free(FA.num_atm);
    free(atoms);
    free(residue);
    std::remove(mol2.c_str());
}

TEST_F(Mol2ReaderTest, HandlesUnknownAtomType) {
    std::string mol2 = write_mol2("unknown.mol2",
        "@<TRIPOS>MOLECULE\n"
        "UNK\n"
        "1 0\n"
        "SMALL\n"
        "\n"
        "@<TRIPOS>ATOM\n"
        "1 X1   1.0  2.0  3.0 Du    1 UNK  0.0\n"
    );

    FA_Global FA;
    atom* atoms = nullptr;
    resid* residue = nullptr;
    init_fa_for_reader(&FA, &atoms, &residue);

    int ok = read_mol2_ligand(&FA, &atoms, &residue, mol2.c_str());
    EXPECT_EQ(ok, 1);

    // Unknown type → dummy type 39
    EXPECT_EQ(atoms[0].type, 39);

    cleanup_fa(&FA, atoms, residue);
    std::remove(mol2.c_str());
}

TEST_F(Mol2ReaderTest, PDBNumbersStartAt90001) {
    std::string mol2 = write_mol2("numbering.mol2",
        "@<TRIPOS>MOLECULE\n"
        "NUM\n"
        "2 1\n"
        "SMALL\n"
        "\n"
        "@<TRIPOS>ATOM\n"
        "1 C1   0.0  0.0  0.0 C.3   1 NUM  0.0\n"
        "2 C2   1.5  0.0  0.0 C.3   1 NUM  0.0\n"
        "@<TRIPOS>BOND\n"
        "1 1 2 1\n"
    );

    FA_Global FA;
    atom* atoms = nullptr;
    resid* residue = nullptr;
    init_fa_for_reader(&FA, &atoms, &residue);

    int ok = read_mol2_ligand(&FA, &atoms, &residue, mol2.c_str());
    EXPECT_EQ(ok, 1);

    EXPECT_EQ(atoms[0].number, 90001);
    EXPECT_EQ(atoms[1].number, 90002);

    // Verify reverse mapping (0-based internal atom indices)
    EXPECT_EQ(FA.num_atm[90001], 0);
    EXPECT_EQ(FA.num_atm[90002], 1);

    cleanup_fa(&FA, atoms, residue);
    std::remove(mol2.c_str());
}

// ===========================================================================
// SDF READER TESTS
// ===========================================================================

class SdfReaderTest : public ::testing::Test {
protected:
    std::string tmp_dir;

    void SetUp() override {
        tmp_dir = std::filesystem::temp_directory_path().string();
    }

    std::string write_sdf(const std::string& name, const std::string& content) {
        std::string path = tmp_dir + "/" + name;
        std::ofstream ofs(path);
        ofs << content;
        return path;
    }
};

TEST_F(SdfReaderTest, ReadsSimpleMolecule) {
    // Methane: 1 carbon, 4 hydrogens
    std::string sdf = write_sdf("methane.sdf",
        "methane\n"
        "  test\n"
        "\n"
        "  5  4  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\n"
        "    0.6300    0.6300    0.6300 H   0  0  0  0  0  0\n"
        "   -0.6300   -0.6300    0.6300 H   0  0  0  0  0  0\n"
        "   -0.6300    0.6300   -0.6300 H   0  0  0  0  0  0\n"
        "    0.6300   -0.6300   -0.6300 H   0  0  0  0  0  0\n"
        "  1  2  1  0\n"
        "  1  3  1  0\n"
        "  1  4  1  0\n"
        "  1  5  1  0\n"
        "M  END\n"
        "$$$$\n"
    );

    FA_Global FA;
    atom* atoms = nullptr;
    resid* residue = nullptr;
    init_fa_for_reader(&FA, &atoms, &residue);

    int ok = read_sdf_ligand(&FA, &atoms, &residue, sdf.c_str());
    EXPECT_EQ(ok, 1);

    EXPECT_EQ(FA.num_het_atm, 5);
    EXPECT_EQ(FA.res_cnt, 1);

    // Carbon at origin (atoms stored 0-based)
    EXPECT_NEAR(atoms[0].coor[0], 0.0f, 0.01f);
    EXPECT_NEAR(atoms[0].coor[1], 0.0f, 0.01f);
    EXPECT_NEAR(atoms[0].coor[2], 0.0f, 0.01f);

    // Types: C → 1, H → 22
    EXPECT_EQ(atoms[0].type, 1);
    EXPECT_EQ(atoms[1].type, 22);

    // Carbon has 4 bonds
    EXPECT_EQ(atoms[0].bond[0], 4);
    // Each hydrogen has 1 bond
    EXPECT_EQ(atoms[1].bond[0], 1);
    EXPECT_EQ(atoms[2].bond[0], 1);

    // Radii
    EXPECT_NEAR(atoms[0].radius, 1.70f, 0.01f);  // carbon
    EXPECT_NEAR(atoms[1].radius, 1.20f, 0.01f);  // hydrogen

    // Residue setup
    EXPECT_EQ(residue[1].type, 1);
    EXPECT_EQ(FA.resligand, &residue[1]);

    cleanup_fa(&FA, atoms, residue);
    std::remove(sdf.c_str());
}

TEST_F(SdfReaderTest, ReadsHalogens) {
    // Test halogen type and radius mapping
    std::string sdf = write_sdf("halogens.sdf",
        "halogens\n"
        "\n"
        "\n"
        "  4  3  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\n"
        "    1.5000    0.0000    0.0000 F   0  0  0  0  0  0\n"
        "    0.0000    1.5000    0.0000 Cl  0  0  0  0  0  0\n"
        "    0.0000    0.0000    1.5000 Br  0  0  0  0  0  0\n"
        "  1  2  1  0\n"
        "  1  3  1  0\n"
        "  1  4  1  0\n"
        "M  END\n"
    );

    FA_Global FA;
    atom* atoms = nullptr;
    resid* residue = nullptr;
    init_fa_for_reader(&FA, &atoms, &residue);

    int ok = read_sdf_ligand(&FA, &atoms, &residue, sdf.c_str());
    EXPECT_EQ(ok, 1);
    EXPECT_EQ(FA.num_het_atm, 4);

    // F → type 13, Cl → type 14, Br → type 15 (atoms stored 0-based)
    EXPECT_EQ(atoms[1].type, 13);
    EXPECT_EQ(atoms[2].type, 14);
    EXPECT_EQ(atoms[3].type, 15);

    // Radii
    EXPECT_NEAR(atoms[1].radius, 1.47f, 0.01f);  // F
    EXPECT_NEAR(atoms[2].radius, 1.75f, 0.01f);  // Cl
    EXPECT_NEAR(atoms[3].radius, 1.85f, 0.01f);  // Br

    cleanup_fa(&FA, atoms, residue);
    std::remove(sdf.c_str());
}

TEST_F(SdfReaderTest, FailsOnMissingFile) {
    FA_Global FA;
    atom* atoms = nullptr;
    resid* residue = nullptr;
    init_fa_for_reader(&FA, &atoms, &residue);

    int ok = read_sdf_ligand(&FA, &atoms, &residue, "/nonexistent/file.sdf");
    EXPECT_EQ(ok, 0);

    free(FA.num_atm);
    free(atoms);
    free(residue);
}

TEST_F(SdfReaderTest, FailsOnInvalidAtomCount) {
    std::string sdf = write_sdf("bad_count.sdf",
        "bad\n"
        "\n"
        "\n"
        "  0  0  0  0  0  0  0  0  0  0999 V2000\n"
        "M  END\n"
    );

    FA_Global FA;
    atom* atoms = nullptr;
    resid* residue = nullptr;
    init_fa_for_reader(&FA, &atoms, &residue);

    int ok = read_sdf_ligand(&FA, &atoms, &residue, sdf.c_str());
    EXPECT_EQ(ok, 0);

    free(FA.num_atm);
    free(atoms);
    free(residue);
    std::remove(sdf.c_str());
}

TEST_F(SdfReaderTest, MoleculeNameExtracted) {
    std::string sdf = write_sdf("named.sdf",
        "Caffeine\n"
        "\n"
        "\n"
        "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\n"
        "M  END\n"
    );

    FA_Global FA;
    atom* atoms = nullptr;
    resid* residue = nullptr;
    init_fa_for_reader(&FA, &atoms, &residue);

    int ok = read_sdf_ligand(&FA, &atoms, &residue, sdf.c_str());
    EXPECT_EQ(ok, 1);

    // Residue name should be first 3 chars of molecule name
    EXPECT_STREQ(residue[1].name, "Caf");

    cleanup_fa(&FA, atoms, residue);
    std::remove(sdf.c_str());
}

TEST_F(SdfReaderTest, BondOutOfRangeIgnored) {
    // Bond referencing atom index > natoms should be silently skipped
    std::string sdf = write_sdf("bad_bond.sdf",
        "test\n"
        "\n"
        "\n"
        "  2  2  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\n"
        "    1.5000    0.0000    0.0000 C   0  0  0  0  0  0\n"
        "  1  2  1  0\n"
        "  1  9  1  0\n"
        "M  END\n"
    );

    FA_Global FA;
    atom* atoms = nullptr;
    resid* residue = nullptr;
    init_fa_for_reader(&FA, &atoms, &residue);

    int ok = read_sdf_ligand(&FA, &atoms, &residue, sdf.c_str());
    EXPECT_EQ(ok, 1);

    // Only the valid bond (1-2) should be recorded (atoms stored 0-based)
    EXPECT_EQ(atoms[0].bond[0], 1);
    EXPECT_EQ(atoms[1].bond[0], 1);

    cleanup_fa(&FA, atoms, residue);
    std::remove(sdf.c_str());
}

// ===========================================================================
// MAIN
// ===========================================================================

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
