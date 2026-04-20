// tests/test_ptm_attachment.cpp
// Unit tests for PTMAttachment module (glycan/PTM on target chain)
// Apache-2.0 (c) 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include "../LIB/PTMAttachment/PTMAttachment.h"
#include <cstring>
#include <cmath>
#include <fstream>

using namespace ptm;

// ===========================================================================
// HELPERS
// ===========================================================================

/// Write a minimal glycan_conformers.json to a temp file for testing
static std::string write_test_json() {
    std::string path = "test_glycan_conformers.json";
    std::ofstream ofs(path);
    ofs << R"({
  "TestGlycan": {
    "description": "Test N-glycan",
    "attachment_residues": ["ASN"],
    "bond_atom": "ND2",
    "linkage": "N-glycosidic",
    "added_atoms": [
      {"name": "C1",  "element": "C", "radius": 1.70, "type_name": "C.3", "resp_charge": 0.312},
      {"name": "O5",  "element": "O", "radius": 1.52, "type_name": "O.3", "resp_charge": -0.458}
    ],
    "conformers": [
      {"phi": -75.0, "psi": 135.0, "omega": 60.0, "weight": 0.65},
      {"phi": -60.0, "psi": 180.0, "omega": -60.0, "weight": 0.35}
    ]
  },
  "TestPhos": {
    "description": "Test phosphorylation",
    "attachment_residues": ["SER", "THR", "TYR"],
    "bond_atom": "OG",
    "linkage": "phosphoester",
    "added_atoms": [
      {"name": "P",   "element": "P", "radius": 1.80, "type_name": "P.3", "resp_charge": 1.166},
      {"name": "O1P", "element": "O", "radius": 1.52, "type_name": "O.co2", "resp_charge": -0.776}
    ],
    "conformers": [
      {"phi": 180.0, "psi": 180.0, "weight": 1.0}
    ]
  }
})";
    ofs.close();
    return path;
}

/// Build a minimal mock atom array with an ASN residue containing ND2
static void setup_mock_receptor(FA_Global& FA, atom* atoms, resid* residue) {
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wnontrivial-memcall"
    std::memset(&FA, 0, sizeof(FA_Global));
    #pragma clang diagnostic pop
    std::memset(atoms, 0, sizeof(atom) * 100);
    std::memset(residue, 0, sizeof(resid) * 10);

    FA.atm_cnt = 5;
    FA.atm_cnt_real = 5;
    FA.res_cnt = 1;
    FA.MIN_NUM_ATOM = 100;

    // Residue 0: ASN 123
    std::strncpy(residue[0].name, "ASN", 3);
    residue[0].number = 123;
    residue[0].chn = 'A';
    residue[0].fatm = new int[1]{0};
    residue[0].latm = new int[1]{4};

    // Atom 0: N
    std::strncpy(atoms[0].name, "N", 4);
    atoms[0].coor[0] = 0.0f; atoms[0].coor[1] = 0.0f; atoms[0].coor[2] = 0.0f;
    atoms[0].radius = 1.55f; atoms[0].ofres = 0; atoms[0].number = 1;

    // Atom 1: CA
    std::strncpy(atoms[1].name, "CA", 4);
    atoms[1].coor[0] = 1.47f; atoms[1].coor[1] = 0.0f; atoms[1].coor[2] = 0.0f;
    atoms[1].radius = 1.70f; atoms[1].ofres = 0; atoms[1].number = 2;
    atoms[1].bond[0] = 1; atoms[1].bond[1] = 0;

    // Atom 2: CB
    std::strncpy(atoms[2].name, "CB", 4);
    atoms[2].coor[0] = 2.0f; atoms[2].coor[1] = 1.4f; atoms[2].coor[2] = 0.0f;
    atoms[2].radius = 1.70f; atoms[2].ofres = 0; atoms[2].number = 3;
    atoms[2].bond[0] = 1; atoms[2].bond[1] = 1;

    // Atom 3: CG
    std::strncpy(atoms[3].name, "CG", 4);
    atoms[3].coor[0] = 3.0f; atoms[3].coor[1] = 2.0f; atoms[3].coor[2] = 0.5f;
    atoms[3].radius = 1.70f; atoms[3].ofres = 0; atoms[3].number = 4;
    atoms[3].bond[0] = 1; atoms[3].bond[1] = 2;

    // Atom 4: ND2 (glycan attachment point)
    std::strncpy(atoms[4].name, "ND2", 4);
    atoms[4].coor[0] = 4.0f; atoms[4].coor[1] = 2.5f; atoms[4].coor[2] = 1.0f;
    atoms[4].radius = 1.55f; atoms[4].ofres = 0; atoms[4].number = 5;
    atoms[4].bond[0] = 1; atoms[4].bond[1] = 3;
}

static void cleanup_mock_receptor(resid* residue) {
    delete[] residue[0].fatm;
    delete[] residue[0].latm;
}

// ===========================================================================
// TEST: JSON LIBRARY LOADING
// ===========================================================================

class PTMLibraryTest : public ::testing::Test {
protected:
    std::string json_path;
    void SetUp() override { json_path = write_test_json(); }
    void TearDown() override { std::remove(json_path.c_str()); }
};

TEST_F(PTMLibraryTest, LoadLibraryReadsAllEntries) {
    auto lib = load_ptm_library(json_path.c_str());
    EXPECT_EQ(lib.size(), 2u);
    EXPECT_TRUE(lib.count("TestGlycan") > 0);
    EXPECT_TRUE(lib.count("TestPhos") > 0);
}

TEST_F(PTMLibraryTest, GlycanDefinitionFields) {
    auto lib = load_ptm_library(json_path.c_str());
    const auto& glycan = lib.at("TestGlycan");

    EXPECT_EQ(glycan.description, "Test N-glycan");
    EXPECT_EQ(glycan.bond_atom, "ND2");
    EXPECT_EQ(glycan.linkage, "N-glycosidic");
    EXPECT_EQ(glycan.attachment_residues.size(), 1u);
    EXPECT_EQ(glycan.attachment_residues[0], "ASN");
    EXPECT_EQ(glycan.added_atoms.size(), 2u);
    EXPECT_EQ(glycan.conformers.size(), 2u);
}

TEST_F(PTMLibraryTest, ConformerWeightsAreCorrect) {
    auto lib = load_ptm_library(json_path.c_str());
    const auto& conformers = lib.at("TestGlycan").conformers;

    EXPECT_NEAR(conformers[0].phi, -75.0, 1e-6);
    EXPECT_NEAR(conformers[0].psi, 135.0, 1e-6);
    EXPECT_NEAR(conformers[0].omega, 60.0, 1e-6);
    EXPECT_NEAR(conformers[0].weight, 0.65, 1e-6);
    EXPECT_NEAR(conformers[1].weight, 0.35, 1e-6);
}

TEST_F(PTMLibraryTest, AddedAtomsHaveRESPCharges) {
    auto lib = load_ptm_library(json_path.c_str());
    const auto& atoms = lib.at("TestGlycan").added_atoms;

    EXPECT_NEAR(atoms[0].resp_charge, 0.312f, 1e-3);
    EXPECT_NEAR(atoms[1].resp_charge, -0.458f, 1e-3);
}

TEST_F(PTMLibraryTest, NonexistentFileThrows) {
    EXPECT_THROW(load_ptm_library("/nonexistent/path.json"), std::runtime_error);
}

// ===========================================================================
// TEST: SPEC PARSING
// ===========================================================================

TEST(PTMSpecTest, ParseSingleSpec) {
    auto specs = parse_ptm_spec("ASN123:NMan9");
    ASSERT_EQ(specs.size(), 1u);
    EXPECT_EQ(specs[0].first, 123);
    EXPECT_EQ(specs[0].second, "NMan9");
}

TEST(PTMSpecTest, ParseMultipleSpecs) {
    auto specs = parse_ptm_spec("ASN123:NMan9,THR45:OGalNAc,SER10:Phosphate");
    ASSERT_EQ(specs.size(), 3u);
    EXPECT_EQ(specs[0].first, 123);
    EXPECT_EQ(specs[0].second, "NMan9");
    EXPECT_EQ(specs[1].first, 45);
    EXPECT_EQ(specs[1].second, "OGalNAc");
    EXPECT_EQ(specs[2].first, 10);
    EXPECT_EQ(specs[2].second, "Phosphate");
}

TEST(PTMSpecTest, ParseNumericOnlyResidueNumber) {
    auto specs = parse_ptm_spec("123:Phosphate");
    ASSERT_EQ(specs.size(), 1u);
    EXPECT_EQ(specs[0].first, 123);
}

TEST(PTMSpecTest, EmptySpecReturnsEmpty) {
    auto specs = parse_ptm_spec("");
    EXPECT_TRUE(specs.empty());

    specs = parse_ptm_spec(nullptr);
    EXPECT_TRUE(specs.empty());
}

TEST(PTMSpecTest, InvalidSpecThrows) {
    EXPECT_THROW(parse_ptm_spec("NOSEP"), std::runtime_error);
}

// ===========================================================================
// TEST: RESIDUE RESOLUTION
// ===========================================================================

TEST(PTMResolveTest, ResolvesKnownResidue) {
    resid residues[2];
    std::memset(residues, 0, sizeof(residues));
    residues[0].number = 123; residues[0].chn = 'A';
    residues[1].number = 45;  residues[1].chn = 'A';

    EXPECT_EQ(resolve_residue(residues, 2, 123), 0);
    EXPECT_EQ(resolve_residue(residues, 2, 45), 1);
}

TEST(PTMResolveTest, ReturnsMinusOneForUnknown) {
    resid residues[1];
    std::memset(residues, 0, sizeof(residues));
    residues[0].number = 123;

    EXPECT_EQ(resolve_residue(residues, 1, 999), -1);
}

// ===========================================================================
// TEST: VALIDATION & ATTACHMENT
// ===========================================================================

class PTMAttachmentTest : public ::testing::Test {
protected:
    FA_Global FA;
    atom atoms_arr[100];
    resid residues[10];
    std::string json_path;

    void SetUp() override {
        json_path = write_test_json();
        setup_mock_receptor(FA, atoms_arr, residues);
    }
    void TearDown() override {
        cleanup_mock_receptor(residues);
        std::remove(json_path.c_str());
    }
};

TEST_F(PTMAttachmentTest, FindAtomInResidue) {
    int idx = find_atom_in_residue(atoms_arr, residues, 0, "ND2");
    EXPECT_EQ(idx, 4);
}

TEST_F(PTMAttachmentTest, FindAtomReturnsMinusOneForMissing) {
    int idx = find_atom_in_residue(atoms_arr, residues, 0, "OG");
    EXPECT_EQ(idx, -1);
}

TEST_F(PTMAttachmentTest, ValidateAttachmentSucceeds) {
    auto lib = load_ptm_library(json_path.c_str());
    EXPECT_TRUE(validate_attachment(atoms_arr, residues, 0, lib.at("TestGlycan")));
}

TEST_F(PTMAttachmentTest, ValidateAttachmentFailsWrongResidue) {
    auto lib = load_ptm_library(json_path.c_str());
    // TestPhos attaches to SER/THR/TYR, not ASN
    EXPECT_FALSE(validate_attachment(atoms_arr, residues, 0, lib.at("TestPhos")));
}

TEST_F(PTMAttachmentTest, AttachModificationAddsAtoms) {
    auto lib = load_ptm_library(json_path.c_str());
    int atm_before = FA.atm_cnt_real;

    PTMSite site = attach_modification(&FA, atoms_arr, residues, 0, lib.at("TestGlycan"));

    EXPECT_EQ(site.mod_name, "TestGlycan");
    EXPECT_EQ(site.residue_index, 0);
    EXPECT_EQ(site.chain_atom, 4); // ND2
    EXPECT_EQ(site.n_added_atoms, 2);
    EXPECT_EQ(site.first_added_atom, atm_before);
    EXPECT_EQ(FA.atm_cnt_real, atm_before + 2);
}

TEST_F(PTMAttachmentTest, AddedAtomsHaveRESPFlags) {
    auto lib = load_ptm_library(json_path.c_str());
    PTMSite site = attach_modification(&FA, atoms_arr, residues, 0, lib.at("TestGlycan"));

    for (int i = 0; i < site.n_added_atoms; ++i) {
        int idx = site.first_added_atom + i;
        EXPECT_EQ(atoms_arr[idx].has_resp, 1);
        EXPECT_EQ(atoms_arr[idx].is_ptm, 1);
        EXPECT_EQ(atoms_arr[idx].ptm_parent, 4); // bonded to ND2
    }
}

TEST_F(PTMAttachmentTest, AddedAtomsHaveCorrectCharges) {
    auto lib = load_ptm_library(json_path.c_str());
    PTMSite site = attach_modification(&FA, atoms_arr, residues, 0, lib.at("TestGlycan"));

    EXPECT_NEAR(atoms_arr[site.first_added_atom].resp_charge, 0.312f, 1e-3);
    EXPECT_NEAR(atoms_arr[site.first_added_atom + 1].resp_charge, -0.458f, 1e-3);
}

TEST_F(PTMAttachmentTest, AddedAtomsArePositionedNearAttachmentPoint) {
    auto lib = load_ptm_library(json_path.c_str());
    PTMSite site = attach_modification(&FA, atoms_arr, residues, 0, lib.at("TestGlycan"));

    float* anchor = atoms_arr[site.chain_atom].coor;
    for (int i = 0; i < site.n_added_atoms; ++i) {
        int idx = site.first_added_atom + i;
        float dx = atoms_arr[idx].coor[0] - anchor[0];
        float dy = atoms_arr[idx].coor[1] - anchor[1];
        float dz = atoms_arr[idx].coor[2] - anchor[2];
        float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
        // Should be within a few bond lengths
        EXPECT_GT(dist, 0.5f);
        EXPECT_LT(dist, 15.0f);
    }
}

TEST_F(PTMAttachmentTest, ApplyConformerChangesPositions) {
    auto lib = load_ptm_library(json_path.c_str());
    PTMSite site = attach_modification(&FA, atoms_arr, residues, 0, lib.at("TestGlycan"));

    // Save conformer 0 positions
    float pos0_x = atoms_arr[site.first_added_atom].coor[0];
    float pos0_y = atoms_arr[site.first_added_atom].coor[1];

    // Apply conformer 1
    apply_conformer(atoms_arr, residues, site, lib.at("TestGlycan"), 1);

    float pos1_x = atoms_arr[site.first_added_atom].coor[0];
    float pos1_y = atoms_arr[site.first_added_atom].coor[1];

    // Different conformer should produce different positions
    EXPECT_TRUE(std::fabs(pos1_x - pos0_x) > 0.01f || std::fabs(pos1_y - pos0_y) > 0.01f);
}

TEST_F(PTMAttachmentTest, ApplyConformerOutOfRangeThrows) {
    auto lib = load_ptm_library(json_path.c_str());
    PTMSite site = attach_modification(&FA, atoms_arr, residues, 0, lib.at("TestGlycan"));

    EXPECT_THROW(
        apply_conformer(atoms_arr, residues, site, lib.at("TestGlycan"), 99),
        std::out_of_range);
}

// ===========================================================================
// TEST: HIGH-LEVEL apply_target_modifications
// ===========================================================================

TEST_F(PTMAttachmentTest, ApplyTargetModificationsEndToEnd) {
    auto lib = load_ptm_library(json_path.c_str());
    PTMState state = apply_target_modifications(
        &FA, atoms_arr, residues, lib, "ASN123:TestGlycan");

    EXPECT_EQ(state.sites.size(), 1u);
    EXPECT_EQ(state.total_added_atoms, 2);
    EXPECT_EQ(state.sites[0].mod_name, "TestGlycan");
}

TEST_F(PTMAttachmentTest, UnknownModificationThrows) {
    auto lib = load_ptm_library(json_path.c_str());
    EXPECT_THROW(
        apply_target_modifications(&FA, atoms_arr, residues, lib, "ASN123:FakeGlycan"),
        std::runtime_error);
}

TEST_F(PTMAttachmentTest, UnknownResidueThrows) {
    auto lib = load_ptm_library(json_path.c_str());
    EXPECT_THROW(
        apply_target_modifications(&FA, atoms_arr, residues, lib, "ASN999:TestGlycan"),
        std::runtime_error);
}

// ===========================================================================
// TEST: RESP CHARGE ON atom_struct
// ===========================================================================

TEST(RESPAtomTest, DefaultFieldsAreZero) {
    atom a;
    std::memset(&a, 0, sizeof(atom));

    EXPECT_FLOAT_EQ(a.resp_charge, 0.0f);
    EXPECT_EQ(a.has_resp, 0);
    EXPECT_EQ(a.is_ptm, 0);
    EXPECT_EQ(a.ptm_parent, 0);
}
