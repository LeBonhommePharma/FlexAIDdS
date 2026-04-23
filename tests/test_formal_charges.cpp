// tests/test_formal_charges.cpp
// Unit tests for residue-based formal charge assignment (PDB input)
// Apache-2.0 (c) 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include "../LIB/assign_formal_charges.h"
#include "../LIB/flexaid.h"
#include <cstring>
#include <cstdlib>

using namespace formal_charges;

// ─── Test fixture: builds a minimal FA + atoms + residue structure ──────────
class FormalChargeTest : public ::testing::Test {
protected:
    FA_Global FA;
    atom atoms[50];
    resid residues[10];
    int fatm_storage[10];
    int latm_storage[10];

    void SetUp() override {
        std::memset(&FA, 0, sizeof(FA));
        std::memset(atoms, 0, sizeof(atoms));
        std::memset(residues, 0, sizeof(residues));
        FA.res_cnt = 0;
        FA.atm_cnt = 0;
        for (int i = 0; i < 10; i++) {
            residues[i].fatm = &fatm_storage[i];
            residues[i].latm = &latm_storage[i];
        }
    }

    // Helper: add a residue with given atoms
    void add_residue(const char* name, int type, const char* atom_names[],
                     int n_atoms, bool is_ter = false) {
        FA.res_cnt++;
        int r = FA.res_cnt;
        std::strncpy(residues[r].name, name, 3);
        residues[r].name[3] = '\0';
        residues[r].type = type;
        residues[r].ter = is_ter ? 1 : 0;
        residues[r].chn = 'A';
        residues[r].fatm[0] = FA.atm_cnt + 1;

        for (int i = 0; i < n_atoms; i++) {
            FA.atm_cnt++;
            std::strncpy(atoms[FA.atm_cnt].name, atom_names[i], 4);
            atoms[FA.atm_cnt].name[4] = '\0';
            atoms[FA.atm_cnt].charge = 0.0f;
            atoms[FA.atm_cnt].ofres = r;
            // Mark backbone atoms
            atoms[FA.atm_cnt].isbb = 0;
            if (std::strcmp(atom_names[i], " N  ") == 0 ||
                std::strcmp(atom_names[i], " CA ") == 0 ||
                std::strcmp(atom_names[i], " C  ") == 0 ||
                std::strcmp(atom_names[i], " O  ") == 0 ||
                std::strcmp(atom_names[i], " CB ") == 0 ||
                std::strcmp(atom_names[i], " OXT") == 0) {
                atoms[FA.atm_cnt].isbb = 1;
            }
        }
        residues[r].latm[0] = FA.atm_cnt;
    }
};

// ===========================================================================
// CHARGED AMINO ACIDS
// ===========================================================================

TEST_F(FormalChargeTest, AspartateCharges) {
    const char* asp_atoms[] = {" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " OD1", " OD2"};
    add_residue("ASP", 0, asp_atoms, 8);

    assign_formal_charges(&FA, atoms, residues);

    // Find OD1 and OD2
    for (int i = 1; i <= FA.atm_cnt; i++) {
        if (std::strncmp(atoms[i].name, " OD1", 4) == 0)
            EXPECT_NEAR(atoms[i].charge, -0.80f, 0.01f) << "ASP OD1 should be -0.80";
        if (std::strncmp(atoms[i].name, " OD2", 4) == 0)
            EXPECT_NEAR(atoms[i].charge, -0.80f, 0.01f) << "ASP OD2 should be -0.80";
        if (std::strncmp(atoms[i].name, " CG ", 4) == 0)
            EXPECT_NEAR(atoms[i].charge, +0.70f, 0.01f) << "ASP CG should be +0.70";
    }
}

TEST_F(FormalChargeTest, GlutamateCharges) {
    const char* glu_atoms[] = {" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD ", " OE1", " OE2"};
    add_residue("GLU", 0, glu_atoms, 9);

    assign_formal_charges(&FA, atoms, residues);

    for (int i = 1; i <= FA.atm_cnt; i++) {
        if (std::strncmp(atoms[i].name, " OE1", 4) == 0)
            EXPECT_NEAR(atoms[i].charge, -0.82f, 0.01f);
        if (std::strncmp(atoms[i].name, " OE2", 4) == 0)
            EXPECT_NEAR(atoms[i].charge, -0.82f, 0.01f);
        if (std::strncmp(atoms[i].name, " CD ", 4) == 0)
            EXPECT_NEAR(atoms[i].charge, +0.80f, 0.01f);
    }
}

TEST_F(FormalChargeTest, LysineCharge) {
    const char* lys_atoms[] = {" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD ", " CE ", " NZ "};
    add_residue("LYS", 0, lys_atoms, 9);

    assign_formal_charges(&FA, atoms, residues);

    for (int i = 1; i <= FA.atm_cnt; i++) {
        if (std::strncmp(atoms[i].name, " NZ ", 4) == 0)
            EXPECT_NEAR(atoms[i].charge, +1.00f, 0.01f) << "LYS NZ should be +1.0";
    }
}

TEST_F(FormalChargeTest, ArginineCharges) {
    const char* arg_atoms[] = {" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD ",
                               " NE ", " CZ ", " NH1", " NH2"};
    add_residue("ARG", 0, arg_atoms, 11);

    assign_formal_charges(&FA, atoms, residues);

    for (int i = 1; i <= FA.atm_cnt; i++) {
        if (std::strncmp(atoms[i].name, " NH1", 4) == 0)
            EXPECT_NEAR(atoms[i].charge, +0.45f, 0.01f);
        if (std::strncmp(atoms[i].name, " NH2", 4) == 0)
            EXPECT_NEAR(atoms[i].charge, +0.45f, 0.01f);
        if (std::strncmp(atoms[i].name, " CZ ", 4) == 0)
            EXPECT_NEAR(atoms[i].charge, +0.64f, 0.01f);
        if (std::strncmp(atoms[i].name, " NE ", 4) == 0)
            EXPECT_NEAR(atoms[i].charge, -0.54f, 0.01f);
    }
}

// ===========================================================================
// METAL IONS
// ===========================================================================

TEST_F(FormalChargeTest, CalciumIon) {
    const char* ca_atoms[] = {" CA "};
    add_residue("CA ", 0, ca_atoms, 1);

    assign_formal_charges(&FA, atoms, residues);

    EXPECT_NEAR(atoms[1].charge, +2.0f, 0.01f) << "Ca2+ should be +2.0";
}

TEST_F(FormalChargeTest, ZincIon) {
    const char* zn_atoms[] = {" ZN "};
    add_residue("ZN ", 0, zn_atoms, 1);

    assign_formal_charges(&FA, atoms, residues);

    EXPECT_NEAR(atoms[1].charge, +2.0f, 0.01f) << "Zn2+ should be +2.0";
}

TEST_F(FormalChargeTest, MagnesiumIon) {
    const char* mg_atoms[] = {" MG "};
    add_residue("MG ", 0, mg_atoms, 1);

    assign_formal_charges(&FA, atoms, residues);

    EXPECT_NEAR(atoms[1].charge, +2.0f, 0.01f) << "Mg2+ should be +2.0";
}

TEST_F(FormalChargeTest, ChlorideIon) {
    const char* cl_atoms[] = {" CL "};
    add_residue("CL ", 0, cl_atoms, 1);

    assign_formal_charges(&FA, atoms, residues);

    EXPECT_NEAR(atoms[1].charge, -1.0f, 0.01f) << "Cl- should be -1.0";
}

TEST_F(FormalChargeTest, IronIII) {
    const char* fe_atoms[] = {" FE "};
    add_residue("FE3", 0, fe_atoms, 1);

    assign_formal_charges(&FA, atoms, residues);

    EXPECT_NEAR(atoms[1].charge, +3.0f, 0.01f) << "Fe3+ should be +3.0";
}

// ===========================================================================
// BACKBONE
// ===========================================================================

TEST_F(FormalChargeTest, BackboneCarbonylOxygen) {
    const char* ala_atoms[] = {" N  ", " CA ", " C  ", " O  ", " CB "};
    add_residue("ALA", 0, ala_atoms, 5);

    assign_formal_charges(&FA, atoms, residues);

    for (int i = 1; i <= FA.atm_cnt; i++) {
        if (std::strncmp(atoms[i].name, " O  ", 4) == 0 && atoms[i].isbb)
            EXPECT_NEAR(atoms[i].charge, -0.57f, 0.01f) << "Backbone O should be -0.57";
    }
}

TEST_F(FormalChargeTest, CTerminusOXT) {
    const char* ala_atoms[] = {" N  ", " CA ", " C  ", " O  ", " CB ", " OXT"};
    add_residue("ALA", 0, ala_atoms, 6, /*is_ter=*/true);

    assign_formal_charges(&FA, atoms, residues);

    for (int i = 1; i <= FA.atm_cnt; i++) {
        if (std::strncmp(atoms[i].name, " OXT", 4) == 0)
            EXPECT_NEAR(atoms[i].charge, -0.83f, 0.01f) << "OXT should be -0.83";
        if (std::strncmp(atoms[i].name, " O  ", 4) == 0 && atoms[i].isbb)
            EXPECT_NEAR(atoms[i].charge, -0.83f, 0.01f) << "C-term O should match OXT at -0.83";
    }
}

// ===========================================================================
// PRESERVATION OF EXISTING CHARGES
// ===========================================================================

TEST_F(FormalChargeTest, DoesNotOverwriteExistingCharges) {
    const char* asp_atoms[] = {" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " OD1", " OD2"};
    add_residue("ASP", 0, asp_atoms, 8);

    // Pre-set a charge (simulating MOL2 input)
    for (int i = 1; i <= FA.atm_cnt; i++) {
        if (std::strncmp(atoms[i].name, " OD1", 4) == 0)
            atoms[i].charge = -0.99f;
    }

    assign_formal_charges(&FA, atoms, residues);

    for (int i = 1; i <= FA.atm_cnt; i++) {
        if (std::strncmp(atoms[i].name, " OD1", 4) == 0)
            EXPECT_NEAR(atoms[i].charge, -0.99f, 0.001f) << "Should not overwrite existing charge";
        if (std::strncmp(atoms[i].name, " OD2", 4) == 0)
            EXPECT_NEAR(atoms[i].charge, -0.80f, 0.01f) << "Should assign where charge was 0";
    }
}

TEST_F(FormalChargeTest, SkipsLigandResidues) {
    const char* lig_atoms[] = {" OD1", " OD2"};
    add_residue("ASP", 1, lig_atoms, 2);  // type=1 = ligand

    assign_formal_charges(&FA, atoms, residues);

    EXPECT_EQ(atoms[1].charge, 0.0f) << "Ligand atoms should not get charges";
    EXPECT_EQ(atoms[2].charge, 0.0f);
}

// ===========================================================================
// NEUTRAL AMINO ACIDS — POLAR ATOMS
// ===========================================================================

TEST_F(FormalChargeTest, CysteineThiol) {
    const char* cys_atoms[] = {" N  ", " CA ", " C  ", " O  ", " CB ", " SG "};
    add_residue("CYS", 0, cys_atoms, 6);

    assign_formal_charges(&FA, atoms, residues);

    for (int i = 1; i <= FA.atm_cnt; i++) {
        if (std::strncmp(atoms[i].name, " SG ", 4) == 0)
            EXPECT_NEAR(atoms[i].charge, -0.31f, 0.01f);
    }
}

TEST_F(FormalChargeTest, DeprotonatedCysteine) {
    const char* cym_atoms[] = {" N  ", " CA ", " C  ", " O  ", " CB ", " SG "};
    add_residue("CYM", 0, cym_atoms, 6);

    assign_formal_charges(&FA, atoms, residues);

    for (int i = 1; i <= FA.atm_cnt; i++) {
        if (std::strncmp(atoms[i].name, " SG ", 4) == 0)
            EXPECT_NEAR(atoms[i].charge, -0.80f, 0.01f) << "CYM SG should be strongly negative";
    }
}

TEST_F(FormalChargeTest, HistidineGeneric) {
    const char* his_atoms[] = {" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " ND1", " CE1", " NE2"};
    add_residue("HIS", 0, his_atoms, 9);

    assign_formal_charges(&FA, atoms, residues);

    for (int i = 1; i <= FA.atm_cnt; i++) {
        if (std::strncmp(atoms[i].name, " ND1", 4) == 0)
            EXPECT_NEAR(atoms[i].charge, -0.35f, 0.01f);
        if (std::strncmp(atoms[i].name, " NE2", 4) == 0)
            EXPECT_NEAR(atoms[i].charge, -0.35f, 0.01f);
    }
}

// ===========================================================================
// MULTIPLE RESIDUES
// ===========================================================================

TEST_F(FormalChargeTest, MixedReceptorWithMetal) {
    // ASP + Ca2+ + ALA: verify all get charges
    const char* asp_atoms[] = {" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " OD1", " OD2"};
    add_residue("ASP", 0, asp_atoms, 8);
    const char* ca_atoms[] = {" CA "};
    add_residue("CA ", 0, ca_atoms, 1);
    const char* ala_atoms[] = {" N  ", " CA ", " C  ", " O  ", " CB "};
    add_residue("ALA", 0, ala_atoms, 5);

    assign_formal_charges(&FA, atoms, residues);

    // Check ASP OD1
    bool found_od1 = false;
    for (int i = 1; i <= FA.atm_cnt; i++) {
        if (std::strncmp(atoms[i].name, " OD1", 4) == 0 &&
            std::strncmp(residues[atoms[i].ofres].name, "ASP", 3) == 0) {
            EXPECT_NEAR(atoms[i].charge, -0.80f, 0.01f);
            found_od1 = true;
        }
    }
    EXPECT_TRUE(found_od1);

    // Check Ca2+
    for (int i = 1; i <= FA.atm_cnt; i++) {
        if (std::strncmp(residues[atoms[i].ofres].name, "CA ", 3) == 0)
            EXPECT_NEAR(atoms[i].charge, +2.0f, 0.01f);
    }
}

// ===========================================================================
// NON-CHARGED RESIDUES — NO SIDE-CHAIN CHARGES
// ===========================================================================

TEST_F(FormalChargeTest, AlanineNoSideChainCharge) {
    const char* ala_atoms[] = {" N  ", " CA ", " C  ", " O  ", " CB "};
    add_residue("ALA", 0, ala_atoms, 5);

    assign_formal_charges(&FA, atoms, residues);

    // CB should remain 0.0 (no charge entry for ALA CB)
    for (int i = 1; i <= FA.atm_cnt; i++) {
        if (std::strncmp(atoms[i].name, " CB ", 4) == 0)
            EXPECT_EQ(atoms[i].charge, 0.0f) << "ALA CB should have no charge";
    }
}
