// tests/test_mol2_sdf_reader.cpp
// Unit tests for Mol2Reader and SdfReader — ligand file parsing
// Apache-2.0 © 2026 NRGlab, Université de Montréal

#include <gtest/gtest.h>
#include "../LIB/Mol2Reader.h"
#include "../LIB/SdfReader.h"
#include "../LIB/fileio.h"
#include "../LIB/LigandRingFlex/LigandRingFlex.h"  // complete RingFlexGenes for delete
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <string>
#include <vector>

void assign_radii_types(FA_Global* FA, atom* atoms, resid* residue);

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
        free(residue[r].gpa);
    }
    // SdfReader may heap-allocate ring topology on FA->ring_flex_template
    // (even when FLEXAIDDS_RING_FLEX is off — detection always runs).
    delete FA->ring_flex_template;
    FA->ring_flex_template = nullptr;
    free(FA->optres);
    free(FA->num_atm);
    free(atoms);
    free(residue);
}

TEST(PdbReceptorTest, RetainedCofactorsHaveCanonicalTypesAndNonzeroRadii) {
    const std::string path =
        std::filesystem::temp_directory_path().string() + "/flexaids_cofactor_reader.pdb";
    {
        std::ofstream pdb(path);
        pdb << "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 10.00           C  \n"
            << "HETATM    2  CHA HEM A 500       2.000   0.000   0.000  1.00 10.00           C  \n"
            << "HETATM    3  FE  HEM A 500       3.000   0.000   0.000  1.00 10.00          FE  \n"
            << "HETATM    4  MG   MG B 401       4.000   0.000   0.000  1.00 10.00          MG  \n"
            << "END\n";
    }

    FA_Global FA{};
    FA.MIN_NUM_ATOM = 32;
    FA.MIN_NUM_RESIDUE = 8;
    FA.MIN_ROTAMER = 1;
    FA.MIN_FLEX_BONDS = 5;
    FA.ntypes = 40;
    atom* atoms = nullptr;
    resid* residue = nullptr;

    read_pdb(&FA, &atoms, &residue, const_cast<char*>(path.c_str()));
    ASSERT_EQ(FA.atm_cnt, 4);
    ASSERT_EQ(FA.res_cnt, 3);
    EXPECT_EQ(residue[FA.res_cnt].latm[0], 4);
    EXPECT_STREQ(atoms[2].element, "C");
    EXPECT_STREQ(atoms[3].element, "Fe");
    EXPECT_STREQ(atoms[4].element, "Mg");
    EXPECT_EQ(atoms[2].type, 3);
    EXPECT_EQ(atoms[3].type, 37);
    EXPECT_EQ(atoms[4].type, 28);

    assign_radii_types(&FA, atoms, residue);
    EXPECT_GT(atoms[2].radius, 1.0f);
    EXPECT_NEAR(atoms[3].radius, 0.61f, 1e-6f);
    EXPECT_NEAR(atoms[4].radius, 0.72f, 1e-6f);

    for (int r = 0; r <= FA.res_cnt; ++r) {
        free(residue[r].fatm);
        free(residue[r].latm);
        free(residue[r].bond);
    }
    free(FA.num_atm);
    free(atoms);
    free(residue);
    std::remove(path.c_str());
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

    // Atom arrays follow legacy FlexAID 1-based indexing; residues are also
    // 1-based.
    // Check coordinates of first atom (oxygen)
    EXPECT_NEAR(atoms[1].coor[0], 0.0f, 0.01f);
    EXPECT_NEAR(atoms[1].coor[1], 0.0f, 0.01f);
    EXPECT_NEAR(atoms[1].coor[2], 0.0f, 0.01f);

    // Check second atom coordinates (H1)
    EXPECT_NEAR(atoms[2].coor[0], 0.957f, 0.01f);

    // Check canonical VCT types: O.3 → 14, H → 39 (DUMMY, not scored)
    EXPECT_EQ(atoms[1].type, 14);
    EXPECT_EQ(atoms[2].type, 39);
    EXPECT_EQ(atoms[3].type, 39);

    // Check radii
    EXPECT_NEAR(atoms[1].radius, 1.52f, 0.01f);  // oxygen
    EXPECT_NEAR(atoms[2].radius, 1.20f, 0.01f);  // hydrogen

    // Check partial charges
    EXPECT_NEAR(atoms[1].charge, -0.834f, 0.01f);
    EXPECT_NEAR(atoms[2].charge,  0.417f, 0.01f);

    // Check bonds: O should have 2 bonds, each H should have 1
    EXPECT_EQ(atoms[1].bond[0], 2);
    EXPECT_EQ(atoms[2].bond[0], 1);
    EXPECT_EQ(atoms[3].bond[0], 1);

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

    // Canonical VCT table: C.ar → 4, C.2 → 2, O.2 → 13, O.3 → 14, N.am → 11
    EXPECT_EQ(atoms[1].type, 4);
    EXPECT_EQ(atoms[2].type, 2);
    EXPECT_EQ(atoms[3].type, 13);
    EXPECT_EQ(atoms[4].type, 14);
    EXPECT_EQ(atoms[5].type, 11);

    // C2 should have 3 bonds (to C1, O1, O2)
    EXPECT_EQ(atoms[2].bond[0], 3);

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
    EXPECT_EQ(atoms[1].type, 39);

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

    EXPECT_EQ(atoms[1].number, 90001);
    EXPECT_EQ(atoms[2].number, 90002);

    // Verify reverse mapping (legacy FlexAID atom arrays are 1-based)
    EXPECT_EQ(FA.num_atm[90001], 1);
    EXPECT_EQ(FA.num_atm[90002], 2);

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

static std::string make_single_atom_sdf(const char* molecule_name, const char* elem) {
    char atom_line[128];
    std::snprintf(atom_line, sizeof(atom_line),
                  "%10.4f%10.4f%10.4f %-3s 0  0  0  0  0  0\n",
                  0.0f, 0.0f, 0.0f, elem);

    return std::string(molecule_name) + "\n"
           "\n"
           "\n"
           "  1  0  0  0  0  0  0  0  0  0999 V2000\n" +
           std::string(atom_line) +
           "M  END\n"
           "$$$$\n";
}

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

    // Carbon at origin (legacy FlexAID 1-based atom storage)
    EXPECT_NEAR(atoms[1].coor[0], 0.0f, 0.01f);
    EXPECT_NEAR(atoms[1].coor[1], 0.0f, 0.01f);
    EXPECT_NEAR(atoms[1].coor[2], 0.0f, 0.01f);

    // Types: generic C → C.3 (3), H → DUMMY (39)
    EXPECT_EQ(atoms[1].type, 3);
    EXPECT_EQ(atoms[2].type, 39);

    // Carbon has 4 bonds
    EXPECT_EQ(atoms[1].bond[0], 4);
    // Each hydrogen has 1 bond
    EXPECT_EQ(atoms[2].bond[0], 1);
    EXPECT_EQ(atoms[3].bond[0], 1);

    // Radii
    EXPECT_NEAR(atoms[1].radius, 1.70f, 0.01f);  // carbon
    EXPECT_NEAR(atoms[2].radius, 1.20f, 0.01f);  // hydrogen

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

    // F → type 23, Cl → type 24, Br → type 25
    EXPECT_EQ(atoms[2].type, 23);
    EXPECT_EQ(atoms[3].type, 24);
    EXPECT_EQ(atoms[4].type, 25);

    // Radii
    EXPECT_NEAR(atoms[2].radius, 1.47f, 0.01f);  // F
    EXPECT_NEAR(atoms[3].radius, 1.75f, 0.01f);  // Cl
    EXPECT_NEAR(atoms[4].radius, 1.85f, 0.01f);  // Br

    cleanup_fa(&FA, atoms, residue);
    std::remove(sdf.c_str());
}

TEST_F(SdfReaderTest, MapsBareSdfElementsToCanonicalTypes) {
    struct Case {
        const char* elem;
        int expected_type;
    };

    const Case cases[] = {
        {"C", 3},
        {"N", 11},
        {"O", 14},
        {"S", 18},
        {"P", 22},
        {"F", 23},
        {"Cl", 24},
        {"Br", 25},
        {"I", 25},
        {"Se", 27},
        {"Mg", 28},
        {"Sr", 29},
        {"Cu", 30},
        {"Mn", 31},
        {"Hg", 32},
        {"Cd", 33},
        {"Ni", 34},
        {"Zn", 35},
        {"Ca", 36},
        {"Fe", 37},
        {"Co", 38},
        {"H", 39},
        {"Xx", 39},
    };

    for (const auto& tc : cases) {
        std::string sdf = write_sdf(std::string("bare_") + tc.elem + ".sdf",
                                    make_single_atom_sdf("bare", tc.elem));

        FA_Global FA;
        atom* atoms = nullptr;
        resid* residue = nullptr;
        init_fa_for_reader(&FA, &atoms, &residue);

        int ok = read_sdf_ligand(&FA, &atoms, &residue, sdf.c_str());
        EXPECT_EQ(ok, 1) << "failed to read element " << tc.elem;
        EXPECT_EQ(FA.num_het_atm, 1) << "unexpected atom count for " << tc.elem;
        EXPECT_EQ(atoms[1].type, tc.expected_type) << "element " << tc.elem;

        cleanup_fa(&FA, atoms, residue);
        std::remove(sdf.c_str());
    }
}

TEST_F(SdfReaderTest, BareNitrogenMapsToActiveScoringType) {
    // SDF/MOL atom blocks only carry an element symbol, not SYBYL hybridisation.
    // FlexAIDdS therefore maps bare N to N.am (11), not N.3 (8), because type 8
    // has no usable interactions in MC_st0r5.2_6.dat and would make the ligand
    // nitrogen effectively invisible to scoring.
    std::string sdf = write_sdf("nitrogen.sdf",
        "nitrogen\n"
        "\n"
        "\n"
        "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 N   0  0  0  0  0  0\n"
        "M  END\n"
    );

    FA_Global FA;
    atom* atoms = nullptr;
    resid* residue = nullptr;
    init_fa_for_reader(&FA, &atoms, &residue);

    int ok = read_sdf_ligand(&FA, &atoms, &residue, sdf.c_str());
    EXPECT_EQ(ok, 1);

    EXPECT_EQ(atoms[1].type, 11);
    EXPECT_NE(atoms[1].type, 8);

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

    // Only the valid bond (1-2) should be recorded.
    EXPECT_EQ(atoms[1].bond[0], 1);
    EXPECT_EQ(atoms[2].bond[0], 1);

    cleanup_fa(&FA, atoms, residue);
    std::remove(sdf.c_str());
}

TEST_F(SdfReaderTest, TopologyDerivedFrameAndTorsionPreserveLocalGeometry) {
    // The first three records are disconnected substituents. A direct reader
    // must choose its GPA from the bonded aromatic core, not record order.
    std::string sdf = write_sdf("topology_frame.sdf",
        "frame\n\n\n"
        " 12 13  0  0  0  0  0  0  0  0999 V2000\n"
        "    3.5000    0.0000    0.0000 C   0  0  0  0  0  0\n"
        "   -1.7500    3.0311    0.0000 C   0  0  0  0  0  0\n"
        "   -1.7500   -3.0311    0.0000 C   0  0  0  0  0  0\n"
        "    2.0000    0.0000    0.0000 C   0  0  0  0  0  0\n"
        "    1.0000    1.7321    0.0000 C   0  0  0  0  0  0\n"
        "   -1.0000    1.7321    0.0000 C   0  0  0  0  0  0\n"
        "   -2.0000    0.0000    0.0000 C   0  0  0  0  0  0\n"
        "   -1.0000   -1.7321    0.0000 C   0  0  0  0  0  0\n"
        "    1.0000   -1.7321    0.0000 C   0  0  0  0  0  0\n"
        "    1.0000    3.2321    0.0000 C   0  0  0  0  0  0\n"
        "    2.2000    4.0321    0.6000 C   0  0  0  0  0  0\n"
        "   -0.2000    4.0321   -0.6000 C   0  0  0  0  0  0\n"
        "  1  4  1  0\n"
        "  2  6  1  0\n"
        "  3  8  1  0\n"
        "  4  5  4  0\n"
        "  5  6  4  0\n"
        "  6  7  4  0\n"
        "  7  8  4  0\n"
        "  8  9  4  0\n"
        "  9  4  4  0\n"
        "  5 10  1  0\n"
        " 10 11  1  0\n"
        " 10 12  1  0\n"
        " 11 12  1  0\n"
        "M  END\n$$$$\n");

    FA_Global FA;
    atom* atoms = nullptr;
    resid* residue = nullptr;
    init_fa_for_reader(&FA, &atoms, &residue);
    ASSERT_EQ(read_sdf_ligand(&FA, &atoms, &residue, sdf.c_str()), 1);

    auto bonded = [&](int a, int b) {
        for (int i = 1; i <= atoms[a].bond[0]; ++i)
            if (atoms[a].bond[i] == b) return true;
        return false;
    };
    ASSERT_NE(residue[1].gpa, nullptr);
    EXPECT_TRUE(bonded(residue[1].gpa[0], residue[1].gpa[1]));
    EXPECT_TRUE(bonded(residue[1].gpa[1], residue[1].gpa[2]));
    EXPECT_GT(residue[1].gpa[0], 3);
    EXPECT_GT(residue[1].gpa[1], 3);
    EXPECT_GT(residue[1].gpa[2], 3);

    auto distance = [&](int a, int b) {
        const double dx = atoms[a].coor[0] - atoms[b].coor[0];
        const double dy = atoms[a].coor[1] - atoms[b].coor[1];
        const double dz = atoms[a].coor[2] - atoms[b].coor[2];
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    };
    auto angle = [&](int a, int b, int c) {
        double u[3], v[3];
        for (int i = 0; i < 3; ++i) {
            u[i] = atoms[a].coor[i] - atoms[b].coor[i];
            v[i] = atoms[c].coor[i] - atoms[b].coor[i];
        }
        const double un = std::sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
        const double vn = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        const double cosine = std::clamp(
            (u[0]*v[0] + u[1]*v[1] + u[2]*v[2]) / (un * vn), -1.0, 1.0);
        return std::acos(cosine) * 180.0 / 3.14159265358979323846;
    };

    struct BondMetric { int a, b; double value; };
    struct AngleMetric { int a, b, c; double value; };
    std::vector<BondMetric> bond_metrics;
    std::vector<AngleMetric> angle_metrics;
    for (int a = 1; a <= FA.num_het_atm; ++a) {
        for (int bi = 1; bi <= atoms[a].bond[0]; ++bi) {
            const int b = atoms[a].bond[bi];
            if (a < b) bond_metrics.push_back({a, b, distance(a, b)});
        }
        for (int i = 1; i <= atoms[a].bond[0]; ++i) {
            for (int j = i + 1; j <= atoms[a].bond[0]; ++j) {
                angle_metrics.push_back(
                    {atoms[a].bond[i], a, atoms[a].bond[j],
                     angle(atoms[a].bond[i], a, atoms[a].bond[j])});
            }
        }
    }

    auto expect_local_geometry = [&] {
        for (const auto& metric : bond_metrics)
            EXPECT_NEAR(distance(metric.a, metric.b), metric.value, 2e-3);
        for (const auto& metric : angle_metrics)
            EXPECT_NEAR(angle(metric.a, metric.b, metric.c), metric.value, 0.15);
    };

    int rebuild[MAX_ATM_HET] = {};
    int rebuilt = 0;
    buildlist(&FA, atoms, residue, 1, 0, &rebuilt, rebuild);
    atoms[residue[1].gpa[0]].dis += 4.0f;
    atoms[residue[1].gpa[0]].ang += 25.0f;
    atoms[residue[1].gpa[0]].dih += 35.0f;
    atoms[residue[1].gpa[1]].ang += 40.0f;
    atoms[residue[1].gpa[1]].dih -= 55.0f;
    atoms[residue[1].gpa[2]].dih += 70.0f;
    ASSERT_TRUE(buildcc(&FA, atoms, rebuilt, rebuild));
    expect_local_geometry();

    ASSERT_EQ(residue[1].fdih, 1);
    const int control = residue[1].bond[1];
    const int axis_child = atoms[control].rec[0];
    const int axis_parent = atoms[control].rec[1];
    EXPECT_TRUE(bonded(axis_child, axis_parent));
    atoms[control].dih += 60.0f;
    int current = control;
    int shifted = atoms[current].rec[3];
    while (shifted != 0 && shifted != control) {
        atoms[shifted].dih = atoms[current].dih + atoms[shifted].shift;
        current = shifted;
        shifted = atoms[current].rec[3];
    }
    rebuilt = 0;
    buildlist(&FA, atoms, residue, 1, 1, &rebuilt, rebuild);
    ASSERT_TRUE(buildcc(&FA, atoms, rebuilt, rebuild));
    expect_local_geometry();

    cleanup_fa(&FA, atoms, residue);
    std::remove(sdf.c_str());
}

// Nonterminal secondary/tertiary alkyl-amine C–N must remain a rotor.
// VCT type 11 is NOT amide proof (SDF generic N and MOL2 N.1/N.2/N.3 → type 11).
// Chain C–C–N–C–C so both ends of each C–N have heavy degree ≥ 2 (nonterminal).
TEST_F(SdfReaderTest, AlkylAmineCNRemainsRuntimeRotor) {
    std::string sdf = write_sdf("alkyl_amine_rotor.sdf",
        "amine\n\n\n"
        "  5  4  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\n"
        "    1.5400    0.0000    0.0000 C   0  0  0  0  0  0\n"
        "    2.2800    1.2600    0.0000 N   0  0  0  0  0  0\n"
        "    3.7200    1.2600    0.0000 C   0  0  0  0  0  0\n"
        "    4.4600    2.5200    0.0000 C   0  0  0  0  0  0\n"
        "  1  2  1  0\n"
        "  2  3  1  0\n"
        "  3  4  1  0\n"
        "  4  5  1  0\n"
        "M  END\n$$$$\n");

    FA_Global FA;
    atom* atoms = nullptr;
    resid* residue = nullptr;
    init_fa_for_reader(&FA, &atoms, &residue);
    ASSERT_EQ(read_sdf_ligand(&FA, &atoms, &residue, sdf.c_str()), 1);
    // Type 11 on N is expected (SDF generic N → 11) but must NOT freeze rotors.
    EXPECT_EQ(atoms[3].type, 11);

    auto el = [&](int idx) -> char {
        if (atoms[idx].element[0])
            return static_cast<char>(std::toupper(
                static_cast<unsigned char>(atoms[idx].element[0])));
        return '?';
    };
    int cn_rotors = 0;
    for (int d = 1; d <= residue[1].fdih; ++d) {
        const int control = residue[1].bond[d];
        ASSERT_GT(control, 0);
        const int child = atoms[control].rec[0];
        const int parent = atoms[control].rec[1];
        const bool cn =
            (el(child) == 'C' && el(parent) == 'N') ||
            (el(child) == 'N' && el(parent) == 'C');
        if (cn) ++cn_rotors;
    }
    EXPECT_GE(cn_rotors, 1)
        << "secondary alkyl-amine C–N must remain DirectLigandIC rotatable "
           "(type==11 must not freeze amine rotors); fdih="
        << residue[1].fdih;

    cleanup_fa(&FA, atoms, residue);
    std::remove(sdf.c_str());
}

// Resonance-locked C–N (amide) must not become a runtime DirectLigandIC rotor.
// Carbonyl C is identified as VCT type C.2 (SdfReader sets type=2 for C with
// double bond) — not "any C bonded to O", which would false-flag carbinolamines.
// Do NOT use atoms[n].type==11 as amide proof.
TEST_F(SdfReaderTest, AmideBondNotRuntimeRotor) {
    // CC(=O)NCC — amide C–N is a graph bridge but resonance-locked.
    std::string sdf = write_sdf("amide_rotor.sdf",
        "amide\n\n\n"
        "  5  4  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\n"
        "    1.5000    0.0000    0.0000 C   0  0  0  0  0  0\n"
        "    2.1000    1.1000    0.0000 O   0  0  0  0  0  0\n"
        "    2.1000   -1.2000    0.0000 N   0  0  0  0  0  0\n"
        "    3.5000   -1.2000    0.0000 C   0  0  0  0  0  0\n"
        "  1  2  1  0\n"
        "  2  3  2  0\n"
        "  2  4  1  0\n"
        "  4  5  1  0\n"
        "M  END\n$$$$\n");

    FA_Global FA;
    atom* atoms = nullptr;
    resid* residue = nullptr;
    init_fa_for_reader(&FA, &atoms, &residue);
    ASSERT_EQ(read_sdf_ligand(&FA, &atoms, &residue, sdf.c_str()), 1);

    auto el = [&](int idx) -> char {
        if (atoms[idx].element[0])
            return static_cast<char>(std::toupper(
                static_cast<unsigned char>(atoms[idx].element[0])));
        return '?';
    };
    for (int d = 1; d <= residue[1].fdih; ++d) {
        const int control = residue[1].bond[d];
        ASSERT_GT(control, 0);
        const int child = atoms[control].rec[0];
        const int parent = atoms[control].rec[1];
        const bool cn =
            (el(child) == 'C' && el(parent) == 'N') ||
            (el(child) == 'N' && el(parent) == 'C');
        if (!cn) continue;
        const int c_idx = el(child) == 'C' ? child : parent;
        // C.2 (carbonyl/sp2) is the typed carbonyl from double-bond perception.
        EXPECT_NE(atoms[c_idx].type, 2)
            << "amide C(=O)–N must not be a DirectLigandIC rotor (gene " << d
            << ", C type=" << atoms[c_idx].type << ")";
    }

    cleanup_fa(&FA, atoms, residue);
    std::remove(sdf.c_str());
}

// Real Astex 1M2Z: SDF read → buildlist → buildcc geometry round-trip.
// (Honest name: does not call ic2cf or pose writers — those need a full FA dock.)
// Historical rupture does NOT reproduce on modern DirectLigandIC + topology GPA.
TEST_F(SdfReaderTest, Real1M2Z_SdfReadBuildlistBuildccGeometryPreserved) {
    namespace fs = std::filesystem;
    const fs::path candidates[] = {
        fs::path("benchmarks/astex_diverse/astex_diverse/1M2Z/1M2Z_ligand.sdf"),
        fs::path("benchmarks/astex_diverse/data/astex_diverse/1M2Z/1M2Z_ligand.sdf"),
        fs::path(__FILE__).parent_path().parent_path() /
            "benchmarks/astex_diverse/astex_diverse/1M2Z/1M2Z_ligand.sdf",
        fs::path(__FILE__).parent_path().parent_path() /
            "benchmarks/astex_diverse/data/astex_diverse/1M2Z/1M2Z_ligand.sdf",
    };
    fs::path sdf;
    for (const auto& c : candidates) {
        if (fs::exists(c)) { sdf = c; break; }
    }
    if (sdf.empty()) {
        GTEST_SKIP() << "1M2Z_ligand.sdf not present in worktree";
    }

    FA_Global FA;
    atom* atoms = nullptr;
    resid* residue = nullptr;
    init_fa_for_reader(&FA, &atoms, &residue);
    ASSERT_EQ(read_sdf_ligand(&FA, &atoms, &residue, sdf.string().c_str()), 1)
        << "read_sdf_ligand failed for " << sdf;
    ASSERT_GE(FA.num_het_atm, 20);
    ASSERT_NE(residue[1].gpa, nullptr);

    auto distance = [&](int a, int b) {
        const double dx = atoms[a].coor[0] - atoms[b].coor[0];
        const double dy = atoms[a].coor[1] - atoms[b].coor[1];
        const double dz = atoms[a].coor[2] - atoms[b].coor[2];
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    };

    struct BondMetric { int a, b; double value; };
    std::vector<BondMetric> bonds;
    for (int a = 1; a <= FA.num_het_atm; ++a) {
        for (int bi = 1; bi <= atoms[a].bond[0]; ++bi) {
            const int b = atoms[a].bond[bi];
            if (a < b) bonds.push_back({a, b, distance(a, b)});
        }
    }
    ASSERT_FALSE(bonds.empty());

    // buildlist + buildcc (IC reconstruction from topology-derived tree).
    // buildcc is fail-closed: false means NaNs were written / reconstruction failed.
    int rebuild[MAX_ATM_HET] = {};
    int rebuilt = 0;
    buildlist(&FA, atoms, residue, 1, 0, &rebuilt, rebuild);
    ASSERT_GT(rebuilt, 0);
    ASSERT_TRUE(buildcc(&FA, atoms, rebuilt, rebuild))
        << "buildcc reconstruction failed (singular frame / non-finite)";

    double max_drift = 0.0;
    for (const auto& m : bonds) {
        const double d = distance(m.a, m.b);
        ASSERT_TRUE(std::isfinite(d)) << "non-finite bond after buildcc";
        max_drift = std::max(max_drift, std::abs(d - m.value));
        EXPECT_LT(std::abs(d - m.value), 0.05)
            << "bond " << m.a << "-" << m.b << " drifted " << (d - m.value);
    }
    EXPECT_LT(max_drift, 1e-3)
        << "1M2Z max bond drift " << max_drift
        << " (historical rupture does not reproduce; expect ~1e-5 A)";

    for (int a = 1; a <= FA.num_het_atm; ++a) {
        for (int k = 0; k < 3; ++k) {
            EXPECT_TRUE(std::isfinite(atoms[a].coor[k]))
                << "atom " << a << " coor non-finite";
        }
    }

    cleanup_fa(&FA, atoms, residue);
}

// ===========================================================================
// MAIN
// ===========================================================================

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
