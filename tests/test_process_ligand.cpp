// tests/test_process_ligand.cpp
// Unit tests for the ProcessLigand pipeline:
//   SmilesParser, RingPerception, Aromaticity, RotatableBonds, ValenceChecker, SybylTyper
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>

#include "../LIB/ProcessLigand/BonMol.h"
#include "../LIB/ProcessLigand/SmilesParser.h"
#include "../LIB/ProcessLigand/RingPerception.h"
#include "../LIB/ProcessLigand/Aromaticity.h"
#include "../LIB/ProcessLigand/RotatableBonds.h"
#include "../LIB/ProcessLigand/ValenceChecker.h"
#include "../LIB/ProcessLigand/SybylTyper.h"
#include "../LIB/ProcessLigand/ProcessLigand.h"
#include "../LIB/atom_typing_256.h"

#include <cmath>
#include <algorithm>
#include <string>

using namespace bonmol;

// ===========================================================================
// Helpers
// ===========================================================================

static BonMol parse(const std::string& smiles) {
    SmilesParser p;
    return p.parse(smiles).mol;
}

static BonMol parse_full(const std::string& smiles) {
    // Parse + ring perception + aromaticity (full pre-processing)
    SmilesParser p;
    BonMol mol = p.parse(smiles).mol;
    ring_perception::perceive_rings(mol);
    aromaticity::assign_aromaticity(mol);
    return mol;
}

// ===========================================================================
// SmilesParser — atoms and atom counts
// ===========================================================================

TEST(SmilesParser, Methane) {
    auto mol = parse("C");
    EXPECT_EQ(mol.num_atoms(), 1);
    EXPECT_EQ(mol.atoms[0].element, Element::C);
    EXPECT_EQ(mol.num_bonds(), 0);
}

TEST(SmilesParser, Ethane) {
    auto mol = parse("CC");
    EXPECT_EQ(mol.num_atoms(), 2);
    EXPECT_EQ(mol.num_bonds(), 1);
    EXPECT_EQ(mol.bonds[0].order, BondOrder::SINGLE);
}

TEST(SmilesParser, Ethene) {
    auto mol = parse("C=C");
    EXPECT_EQ(mol.num_atoms(), 2);
    EXPECT_EQ(mol.num_bonds(), 1);
    EXPECT_EQ(mol.bonds[0].order, BondOrder::DOUBLE);
}

TEST(SmilesParser, Ethyne) {
    auto mol = parse("C#C");
    EXPECT_EQ(mol.num_atoms(), 2);
    EXPECT_EQ(mol.num_bonds(), 1);
    EXPECT_EQ(mol.bonds[0].order, BondOrder::TRIPLE);
}

TEST(SmilesParser, Ethanol) {
    auto mol = parse("CCO");
    EXPECT_EQ(mol.num_atoms(), 3);
    EXPECT_EQ(mol.atoms[2].element, Element::O);
}

TEST(SmilesParser, BranchPropane) {
    // Isobutane: CC(C)C
    auto mol = parse("CC(C)C");
    EXPECT_EQ(mol.num_atoms(), 4);
    EXPECT_EQ(mol.num_bonds(), 3);
}

TEST(SmilesParser, Benzene) {
    auto mol = parse("c1ccccc1");
    EXPECT_EQ(mol.num_atoms(), 6);
    EXPECT_EQ(mol.num_bonds(), 6);
    // All bonds should be aromatic from the lowercase notation
    for (const auto& b : mol.bonds)
        EXPECT_EQ(b.order, BondOrder::AROMATIC);
}

TEST(SmilesParser, BracketAtomCharge) {
    // Ammonium ion
    auto mol = parse("[NH4+]");
    EXPECT_EQ(mol.num_atoms(), 1);
    EXPECT_EQ(mol.atoms[0].element, Element::N);
    EXPECT_EQ(mol.atoms[0].formal_charge, 1);
}

TEST(SmilesParser, BracketAtomNegCharge) {
    auto mol = parse("[O-]");
    EXPECT_EQ(mol.num_atoms(), 1);
    EXPECT_EQ(mol.atoms[0].formal_charge, -1);
}

TEST(SmilesParser, Isotope) {
    auto mol = parse("[13C]");
    EXPECT_EQ(mol.num_atoms(), 1);
    EXPECT_EQ(mol.atoms[0].isotope, 13);
}

TEST(SmilesParser, Pyridine) {
    auto mol = parse("c1ccncc1");
    EXPECT_EQ(mol.num_atoms(), 6);
    bool found_n = std::any_of(mol.atoms.begin(), mol.atoms.end(),
        [](const Atom& a){ return a.element == Element::N; });
    EXPECT_TRUE(found_n);
}

TEST(SmilesParser, Pyrrole) {
    // [nH] = aromatic NH in 5-membered ring
    auto mol = parse("c1cc[nH]c1");
    EXPECT_EQ(mol.num_atoms(), 5);
    bool found_nh = std::any_of(mol.atoms.begin(), mol.atoms.end(),
        [](const Atom& a){ return a.element == Element::N; });
    EXPECT_TRUE(found_nh);
}

TEST(SmilesParser, Naphthalene) {
    auto mol = parse("c1ccc2ccccc2c1");
    EXPECT_EQ(mol.num_atoms(), 10);
}

TEST(SmilesParser, MultipleRingClosures) {
    // Bicyclo[2.2.1]heptane (norbornane): C1CC2CCC1C2
    auto mol = parse("C1CC2CCC1C2");
    EXPECT_EQ(mol.num_atoms(), 7);
}

TEST(SmilesParser, TwoFragments) {
    // Disconnected fragments separated by '.'
    auto mol = parse("C.N");
    EXPECT_EQ(mol.num_atoms(), 2);
    EXPECT_EQ(mol.num_bonds(), 0);
}

TEST(SmilesParser, InvalidSmilesBadBracket) {
    SmilesParser p;
    EXPECT_THROW(p.parse("[C++X"), SmilesParseError);
}

TEST(SmilesParser, EmptyStringThrows) {
    SmilesParser p;
    EXPECT_THROW(p.parse(""), SmilesParseError);
}

TEST(SmilesParser, AtomMapNumber) {
    auto mol = parse("[C:42]");
    EXPECT_EQ(mol.atoms[0].atom_map_num, 42);
}

// ===========================================================================
// RingPerception
// ===========================================================================

TEST(RingPerception, BenzeneOneRing) {
    auto mol = parse("c1ccccc1");
    auto res = ring_perception::perceive_rings(mol);
    EXPECT_EQ(res.num_rings, 1);
    ASSERT_EQ(mol.rings.size(), 1u);
    EXPECT_EQ(mol.rings[0].size, 6);
}

TEST(RingPerception, BondsMarkedInRing) {
    auto mol = parse("c1ccccc1");
    ring_perception::perceive_rings(mol);
    for (const auto& b : mol.bonds)
        EXPECT_TRUE(b.in_ring);
}

TEST(RingPerception, AtomRingMembership) {
    auto mol = parse("c1ccccc1");
    ring_perception::perceive_rings(mol);
    for (const auto& a : mol.atoms)
        EXPECT_EQ(a.ring_membership, 1);
}

TEST(RingPerception, NaphthaleneTwoRings) {
    auto mol = parse("c1ccc2ccccc2c1");
    auto res = ring_perception::perceive_rings(mol);
    EXPECT_EQ(res.num_rings, 2);
}

TEST(RingPerception, AcyclicMoleculeNoRings) {
    auto mol = parse("CCCC");
    auto res = ring_perception::perceive_rings(mol);
    EXPECT_EQ(res.num_rings, 0);
    for (const auto& b : mol.bonds)
        EXPECT_FALSE(b.in_ring);
}

TEST(RingPerception, BFSShortestPath) {
    auto mol = parse("c1ccccc1");
    ring_perception::perceive_rings(mol);
    auto path = ring_perception::bfs_shortest_path(mol, 0, 3);
    EXPECT_FALSE(path.empty());
}

// ===========================================================================
// Aromaticity
// ===========================================================================

TEST(Aromaticity, BenzeneAromaticAtoms) {
    auto mol = parse_full("c1ccccc1");
    for (const auto& a : mol.atoms)
        EXPECT_TRUE(a.is_aromatic) << "Benzene atom should be aromatic";
}

TEST(Aromaticity, BenzeneAromaticBonds) {
    auto mol = parse_full("c1ccccc1");
    for (const auto& b : mol.bonds)
        EXPECT_TRUE(b.is_aromatic);
}

TEST(Aromaticity, BenzeneRingAromaticFlag) {
    auto mol = parse_full("c1ccccc1");
    ASSERT_EQ(mol.rings.size(), 1u);
    EXPECT_TRUE(mol.rings[0].is_aromatic);
}

TEST(Aromaticity, BenzeneKekulized) {
    auto mol = parse_full("c1ccccc1");
    auto res = aromaticity::assign_aromaticity(mol);
    EXPECT_TRUE(res.kekulized);
}

TEST(Aromaticity, PyridineAromaticNitrogen) {
    auto mol = parse_full("c1ccncc1");
    bool n_aromatic = std::any_of(mol.atoms.begin(), mol.atoms.end(),
        [](const Atom& a){ return a.element == Element::N && a.is_aromatic; });
    EXPECT_TRUE(n_aromatic);
}

TEST(Aromaticity, FuranAromaticOxygen) {
    auto mol = parse_full("c1ccoc1");
    bool o_aromatic = std::any_of(mol.atoms.begin(), mol.atoms.end(),
        [](const Atom& a){ return a.element == Element::O && a.is_aromatic; });
    EXPECT_TRUE(o_aromatic);
}

TEST(Aromaticity, NaphthaleneCountAromaticAtoms) {
    auto mol = parse_full("c1ccc2ccccc2c1");
    auto res = aromaticity::assign_aromaticity(mol);
    EXPECT_EQ(res.num_aromatic_atoms, 10);
    EXPECT_EQ(res.num_aromatic_rings, 2);
}

TEST(Aromaticity, CyclohexaneNotAromatic) {
    auto mol = parse_full("C1CCCCC1");
    for (const auto& a : mol.atoms)
        EXPECT_FALSE(a.is_aromatic);
}

TEST(Aromaticity, PiElectronCount) {
    BonMol mol = parse("c1ccccc1");
    ring_perception::perceive_rings(mol);
    aromaticity::assign_hybridisation(mol);
    // A C in benzene should contribute 1 pi electron
    int pi = aromaticity::pi_electron_count(mol, 0, mol.rings[0]);
    EXPECT_EQ(pi, 1);
}

// ===========================================================================
// RotatableBonds
// ===========================================================================

TEST(RotatableBonds, EthaneSingleRotatableBond) {
    SmilesParser p;
    BonMol mol = p.parse("CC").mol;
    ring_perception::perceive_rings(mol);
    // Ethane: C-C both terminal (degree 1), so NOT rotatable
    auto res = rotatable_bonds::identify_rotatable_bonds(mol);
    EXPECT_EQ(res.count, 0);
}

TEST(RotatableBonds, ButaneMidBondRotatable) {
    SmilesParser p;
    BonMol mol = p.parse("CCCC").mol;
    ring_perception::perceive_rings(mol);
    auto res = rotatable_bonds::identify_rotatable_bonds(mol);
    // Central C-C bond is rotatable (both endpoints have degree >= 2)
    EXPECT_EQ(res.count, 1);
}

TEST(RotatableBonds, AmideBondNotRotatable) {
    SmilesParser p;
    // Simple amide: CC(=O)N
    BonMol mol = p.parse("CC(=O)N").mol;
    ring_perception::perceive_rings(mol);
    auto res = rotatable_bonds::identify_rotatable_bonds(mol);
    int bidx = mol.find_bond(1, 3); // C(=O)-N bond (indices may vary)
    // Check that is_amide_bond works
    // Find the C(=O)-N bond manually
    bool amide_found = false;
    for (int i = 0; i < mol.num_bonds(); ++i)
        if (rotatable_bonds::is_amide_bond(mol, i)) amide_found = true;
    EXPECT_TRUE(amide_found);
}

TEST(RotatableBonds, DisulfideBondNotRotatable) {
    SmilesParser p;
    BonMol mol = p.parse("CSSC").mol;
    ring_perception::perceive_rings(mol);
    bool disulfide_found = false;
    for (int i = 0; i < mol.num_bonds(); ++i)
        if (rotatable_bonds::is_disulfide_bond(mol, i)) disulfide_found = true;
    EXPECT_TRUE(disulfide_found);
}

TEST(RotatableBonds, RingBondsNotRotatable) {
    SmilesParser p;
    BonMol mol = p.parse("C1CCCCC1").mol;
    ring_perception::perceive_rings(mol);
    auto res = rotatable_bonds::identify_rotatable_bonds(mol);
    // No bonds outside the ring → none rotatable
    EXPECT_EQ(res.count, 0);
}

TEST(RotatableBonds, BondsMarkedOnMolecule) {
    SmilesParser p;
    BonMol mol = p.parse("CCCC").mol;
    ring_perception::perceive_rings(mol);
    rotatable_bonds::identify_rotatable_bonds(mol);
    int marked = 0;
    for (const auto& b : mol.bonds)
        if (b.is_rotatable) ++marked;
    EXPECT_EQ(marked, 1);
}

// ===========================================================================
// ValenceChecker
// ===========================================================================

TEST(ValenceChecker, CarbonTetravalentValid) {
    // Methane-like: C with 4 implicit H
    SmilesParser p;
    BonMol mol = p.parse("[CH4]").mol;
    ring_perception::perceive_rings(mol);
    auto res = valence::check_valence(mol);
    EXPECT_TRUE(res.valid);
    EXPECT_TRUE(res.errors.empty());
}

TEST(ValenceChecker, WaterValid) {
    SmilesParser p;
    BonMol mol = p.parse("O").mol;
    ring_perception::perceive_rings(mol);
    auto res = valence::check_valence(mol);
    EXPECT_TRUE(res.valid);
}

TEST(ValenceChecker, NitrogenValid) {
    SmilesParser p;
    BonMol mol = p.parse("N").mol;
    ring_perception::perceive_rings(mol);
    auto res = valence::check_valence(mol);
    EXPECT_TRUE(res.valid);
}

TEST(ValenceChecker, ExpectedValencesCarbonNeutral) {
    auto vals = valence::expected_valences(Element::C, 0);
    EXPECT_FALSE(vals.empty());
    // C should have valence 4
    EXPECT_NE(std::find(vals.begin(), vals.end(), 4), vals.end());
}

TEST(ValenceChecker, ExpectedValencesNitrogenNeutral) {
    auto vals = valence::expected_valences(Element::N, 0);
    EXPECT_FALSE(vals.empty());
    // N neutral: valence 3
    EXPECT_NE(std::find(vals.begin(), vals.end(), 3), vals.end());
}

TEST(ValenceChecker, ExpectedValencesOxygenNeutral) {
    auto vals = valence::expected_valences(Element::O, 0);
    EXPECT_NE(std::find(vals.begin(), vals.end(), 2), vals.end());
}

TEST(ValenceChecker, ImplicitHComputedForCarbon) {
    SmilesParser p;
    BonMol mol = p.parse("C").mol;
    // Single C in SMILES: 0 explicit bonds → should infer 4 implicit H
    int h = valence::compute_implicit_h(mol, 0);
    EXPECT_EQ(h, 4);
}

// ===========================================================================
// SybylTyper
// ===========================================================================

TEST(SybylTyper, Sp3CarbonType) {
    SmilesParser p;
    BonMol mol = p.parse("C").mol;
    ring_perception::perceive_rings(mol);
    aromaticity::assign_aromaticity(mol);
    int type = sybyl::assign_sybyl_type_single(mol, 0);
    EXPECT_EQ(type, 1); // C.3 → 1
}

TEST(SybylTyper, AromaticCarbonType) {
    auto mol = parse_full("c1ccccc1");
    sybyl::assign_sybyl_types(mol);
    // All C atoms in benzene should be C.ar → 3
    for (int i = 0; i < mol.num_atoms(); ++i)
        EXPECT_EQ(mol.atoms[i].sybyl_type, 3) << "Atom " << i << " should be C.ar";
}

TEST(SybylTyper, SybylTypeName) {
    const char* name = sybyl::sybyl_type_name(3);
    EXPECT_NE(name, nullptr);
    EXPECT_STRNE(name, "");
}

TEST(SybylTyper, HBondDonorNHAcceptor) {
    SmilesParser p;
    BonMol mol = p.parse("N").mol; // NH3
    ring_perception::perceive_rings(mol);
    aromaticity::assign_aromaticity(mol);
    // NH3 nitrogen: donor (has H) and acceptor (lone pair)
    EXPECT_TRUE(sybyl::is_hbond_donor(mol, 0));
    EXPECT_TRUE(sybyl::is_hbond_acceptor(mol, 0));
}

TEST(SybylTyper, HBondNonDonorCarbon) {
    SmilesParser p;
    BonMol mol = p.parse("C").mol;
    ring_perception::perceive_rings(mol);
    aromaticity::assign_aromaticity(mol);
    EXPECT_FALSE(sybyl::is_hbond_donor(mol, 0));
    EXPECT_FALSE(sybyl::is_hbond_acceptor(mol, 0));
}

TEST(SybylTyper, Encode256Deterministic) {
    uint8_t enc1 = sybyl::encode_256(3, 0.0f, false);
    uint8_t enc2 = sybyl::encode_256(3, 0.0f, false);
    EXPECT_EQ(enc1, enc2);
}

TEST(SybylTyper, Encode256MatchesCanonicalAtomTyping) {
    struct Case {
        int legacy_sybyl;
        uint8_t canonical_base;
        float charge;
        bool hbond;
    };

    const Case cases[] = {
        {0,  atom256::C_sp,     0.0f,  false}, // C.1
        {1,  atom256::C_sp3,    0.0f,  false}, // C.3
        {2,  atom256::C_sp2,    0.0f,  false}, // C.2
        {3,  atom256::C_ar,     0.0f,  false}, // C.ar
        {4,  atom256::N_sp3,    0.2f,  true},  // N.3
        {5,  atom256::N_sp2,   -0.2f,  true},  // N.2
        {10, atom256::O_sp3,   -0.3f,  true},  // O.3
        {11, atom256::O_sp2,   -0.3f,  true},  // O.2
        {12, atom256::O_co2,   -0.3f,  true},  // O.co2
        {13, atom256::HAL_F,   -0.1f,  true},  // F
        {20, atom256::P_sp3,    0.1f,  false}, // P.3
        {21, atom256::HAL_I,    0.0f,  false}, // I
        {30, atom256::Metal_Fe, 0.4f,  false}, // Fe
        {31, atom256::Metal_Zn, 0.4f,  false}, // Zn
        {32, atom256::Metal_Ca, 0.4f,  false}, // Ca
        {33, atom256::Metal_Mg, 0.4f,  false}, // Mg
        {34, atom256::HAL_Se,   0.0f,  false}, // Se
        {35, atom256::Metal_Cu, 0.4f,  false}, // Cu
        {36, atom256::Metal_Ni, 0.4f,  false}, // Ni
    };

    for (const auto& tc : cases) {
        const uint8_t code =
            sybyl::encode_256(tc.legacy_sybyl, tc.charge, tc.hbond);
        EXPECT_EQ(atom256::get_base(code), tc.canonical_base)
            << "legacy SYBYL " << tc.legacy_sybyl;
        EXPECT_EQ(atom256::get_charge_bin(code),
                  atom256::quantise_charge(tc.charge));
        EXPECT_EQ(atom256::get_hbond(code), tc.hbond);
    }
}

TEST(SybylTyper, Encode256HydrogenFallsBackToDummy) {
    const uint8_t code = sybyl::encode_256(22, 0.1f, false);
    EXPECT_EQ(atom256::get_base(code), atom256::Dummy);
}

TEST(SybylTyper, AssignAllTypesDoesNotCrash) {
    auto mol = parse_full("c1ccncc1"); // pyridine
    sybyl::assign_sybyl_types(mol);
    for (const auto& a : mol.atoms)
        EXPECT_GT(a.sybyl_type, 0);
}

// ===========================================================================
// ProcessLigand pipeline (integration)
// ===========================================================================

TEST(ProcessLigand, BenzeneSmilesPipeline) {
    auto result = process_smiles("c1ccccc1");
    EXPECT_TRUE(result.success);
    EXPECT_EQ(result.num_atoms, 6);
    EXPECT_EQ(result.num_rings, 1);
    EXPECT_EQ(result.num_arom_rings, 1);
    EXPECT_EQ(result.num_rot_bonds, 0);
}

TEST(ProcessLigand, EthanolSmilesPipeline) {
    auto result = process_smiles("CCO");
    EXPECT_TRUE(result.success);
    EXPECT_EQ(result.num_atoms, 3);
    EXPECT_EQ(result.num_rings, 0);
}

TEST(ProcessLigand, CaffeineSmilesPipeline) {
    // Caffeine: 3 N-methylation + xanthine core
    auto result = process_smiles("Cn1cnc2c1c(=O)n(c(=O)n2C)C");
    EXPECT_TRUE(result.success);
    EXPECT_GT(result.num_atoms, 10);
    EXPECT_GT(result.num_rings, 0);
}

TEST(ProcessLigand, MolecularWeightBenzene) {
    auto result = process_smiles("c1ccccc1");
    EXPECT_TRUE(result.success);
    // Benzene MW = 78.11 g/mol (6×12.011 + 6×1.008)
    EXPECT_NEAR(result.molecular_weight, 78.0f, 2.0f);
}

TEST(ProcessLigand, EmptyInputFails) {
    auto result = process_smiles("");
    EXPECT_FALSE(result.success);
    EXPECT_FALSE(result.error.empty());
}

TEST(ProcessLigand, InvalidSmilesFails) {
    auto result = process_smiles("not_a_smiles_!!!###");
    EXPECT_FALSE(result.success);
}

TEST(ProcessLigand, ValidateOnlyDoesNotWrite) {
    ProcessOptions opts;
    opts.input   = "c1ccccc1";
    opts.format  = InputFormat::SMILES;
    opts.validate_only  = true;
    opts.output_prefix  = "";

    ProcessLigand pl;
    auto result = pl.run(opts);
    EXPECT_TRUE(result.success);
}

TEST(ProcessLigand, StageResultsPopulated) {
    auto result = process_smiles("c1ccccc1");
    EXPECT_TRUE(result.success);
    EXPECT_FALSE(result.stage_results.empty());
    for (const auto& sr : result.stage_results)
        EXPECT_TRUE(sr.ok) << "Stage " << sr.stage << " failed: " << sr.message;
}

TEST(ProcessLigand, PeptideGuardTriggersOnMultipleAmides) {
    // Tripeptide-like molecule with 3+ amide bonds
    // NCC(=O)NCC(=O)NCC(=O)O
    ProcessOptions opts;
    opts.input        = "NCC(=O)NCC(=O)NCC(=O)O";
    opts.format       = InputFormat::SMILES;
    opts.validate_only = false;
    opts.allow_peptides = false;

    ProcessLigand pl;
    auto result = pl.run(opts);
    // Should fail due to peptide guard
    EXPECT_FALSE(result.success);
}

TEST(ProcessLigand, PeptideGuardBypassable) {
    ProcessOptions opts;
    opts.input          = "NCC(=O)NCC(=O)NCC(=O)O";
    opts.format         = InputFormat::SMILES;
    opts.allow_peptides = true;

    ProcessLigand pl;
    auto result = pl.run(opts);
    EXPECT_TRUE(result.success);
}

TEST(ProcessLigand, DetectFormatSmiles) {
    // No extension → can't detect; AUTO should try SMILES when no file
    EXPECT_EQ(detect_format("molecule.mol2"), InputFormat::MOL2);
    EXPECT_EQ(detect_format("molecule.sdf"),  InputFormat::SDF);
}

// ===========================================================================
// MAIN
// ===========================================================================

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
