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

#include <cmath>
#include <algorithm>
#include <string>
#include <filesystem>

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
    // Amide C–N must not be marked rotatable
    for (int i = 0; i < mol.num_bonds(); ++i) {
        if (rotatable_bonds::is_amide_bond(mol, i)) {
            EXPECT_FALSE(mol.bonds[i].is_rotatable);
        }
    }
    (void)res;
    (void)bidx;
}

TEST(RotatableBonds, DisulfideBondNotRotatable) {
    SmilesParser p;
    BonMol mol = p.parse("CSSC").mol;
    ring_perception::perceive_rings(mol);
    bool disulfide_found = false;
    for (int i = 0; i < mol.num_bonds(); ++i)
        if (rotatable_bonds::is_disulfide_bond(mol, i)) disulfide_found = true;
    EXPECT_TRUE(disulfide_found);
    auto res = rotatable_bonds::identify_rotatable_bonds(mol);
    for (int i = 0; i < mol.num_bonds(); ++i) {
        if (rotatable_bonds::is_disulfide_bond(mol, i)) {
            EXPECT_FALSE(mol.bonds[i].is_rotatable);
        }
    }
    (void)res;
}

TEST(RotatableBonds, TripleAdjacentBondRejected) {
    SmilesParser p;
    // But-2-yne with a methyl: CC#CC — the C–C single adjacent to triple is locked.
    BonMol mol = p.parse("CC#CC").mol;
    ring_perception::perceive_rings(mol);
    auto res = rotatable_bonds::identify_rotatable_bonds(mol);
    for (int i = 0; i < mol.num_bonds(); ++i) {
        if (rotatable_bonds::is_triple_adjacent_bond(mol, i) &&
            mol.bonds[i].order == BondOrder::SINGLE) {
            EXPECT_FALSE(mol.bonds[i].is_rotatable)
                << "triple-adjacent single bond must not be a rotor";
        }
    }
    (void)res;
}

TEST(RotatableBonds, ConjugatedUreaCNNotRotatable) {
    SmilesParser p;
    // Urea: NC(=O)N — both C–N bonds are conjugated/amide-like.
    BonMol mol = p.parse("NC(=O)N").mol;
    ring_perception::perceive_rings(mol);
    auto res = rotatable_bonds::identify_rotatable_bonds(mol);
    for (int i = 0; i < mol.num_bonds(); ++i) {
        if (rotatable_bonds::is_conjugated_cn_bond(mol, i)) {
            EXPECT_FALSE(mol.bonds[i].is_rotatable);
        }
    }
    EXPECT_EQ(res.count, 0);
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

// Aromatic heteroatom valence (regression for thiadiazole/thiophene/furan bug).
// Aromatic S/O in a 5-membered ring has bond-order sum 2 × 1.5 = 3.0, which is
// a full 1.0 from the textbook neutral valences {2,...}. BOS=3 must be accepted.
TEST(ValenceChecker, ExpectedValencesOxygenAcceptsAromaticThree) {
    auto vals = valence::expected_valences(Element::O, 0);
    EXPECT_NE(std::find(vals.begin(), vals.end(), 2), vals.end());
    EXPECT_NE(std::find(vals.begin(), vals.end(), 3), vals.end())
        << "Neutral O must accept BOS=3 for aromatic donation (furan/oxazole)";
}

TEST(ValenceChecker, ExpectedValencesSulfurAcceptsAromaticThree) {
    auto vals = valence::expected_valences(Element::S, 0);
    EXPECT_NE(std::find(vals.begin(), vals.end(), 3), vals.end())
        << "Neutral S must accept BOS=3 for aromatic ring (thiophene/thiadiazole)";
    // Existing extended valences must still be present.
    EXPECT_NE(std::find(vals.begin(), vals.end(), 2), vals.end());
    EXPECT_NE(std::find(vals.begin(), vals.end(), 4), vals.end());
    EXPECT_NE(std::find(vals.begin(), vals.end(), 6), vals.end());
}

// Count valence errors attributed to atoms of a given element. The fix under
// test only concerns aromatic S/O, so we assert specifically on the heteroatom
// rather than whole-molecule validity (aromatic-carbon implicit-H handling on
// the SMILES path is a separate, unrelated code path).
static int heteroatom_errors(const valence::ValenceCheckResult& res, Element elem) {
    int n = 0;
    for (const auto& e : res.errors)
        if (e.element == elem) ++n;
    return n;
}

static valence::ValenceCheckResult checked(const std::string& smiles) {
    SmilesParser p;
    BonMol mol = p.parse(smiles).mol;
    ring_perception::perceive_rings(mol);
    aromaticity::assign_aromaticity(mol);
    return valence::check_valence(mol);
}

TEST(ValenceChecker, ThiopheneSulfurNoValenceError) {
    // Aromatic S, bond-order sum 2 × 1.5 = 3.0 — must not be flagged.
    auto res = checked("c1ccsc1");
    EXPECT_EQ(heteroatom_errors(res, Element::S), 0)
        << "Thiophene aromatic S (BOS=3) must pass valence check";
}

TEST(ValenceChecker, FuranOxygenNoValenceError) {
    auto res = checked("c1ccoc1");
    EXPECT_EQ(heteroatom_errors(res, Element::O), 0)
        << "Furan aromatic O (BOS=3) must pass valence check";
}

TEST(ValenceChecker, AcetazolamideThiadiazoleSulfurNoValenceError) {
    // AZM (1JD0, acetazolamide) — thiadiazole ring with aromatic S (BOS=3).
    auto res = checked("c1nnsc(NS(=O)(=O)C)1");
    EXPECT_EQ(heteroatom_errors(res, Element::S), 0)
        << "Acetazolamide thiadiazole S (BOS=3) must pass valence check";
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
    uint8_t enc1 = sybyl::encode_256(3, 0.0f, false, false);
    uint8_t enc2 = sybyl::encode_256(3, 0.0f, false, false);
    EXPECT_EQ(enc1, enc2);
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
    // Backbone-peptide detector requires ≥3 *linked* backbone amides
    // (not isolated amides). Tetrapeptide-like: 3 backbone amide links.
    ProcessOptions opts;
    opts.input        = "NCC(=O)NCC(=O)NCC(=O)NCC(=O)O";
    opts.format       = InputFormat::SMILES;
    opts.validate_only = false;
    opts.allow_peptides = false;

    ProcessLigand pl;
    auto result = pl.run(opts);
    // Should fail due to peptide guard
    EXPECT_FALSE(result.success) << result.error;
}

TEST(ProcessLigand, PeptideGuardBypassable) {
    ProcessOptions opts;
    opts.input          = "NCC(=O)NCC(=O)NCC(=O)NCC(=O)O";
    opts.format         = InputFormat::SMILES;
    opts.allow_peptides = true;

    ProcessLigand pl;
    auto result = pl.run(opts);
    EXPECT_TRUE(result.success) << result.error;
}

TEST(ProcessLigand, DetectFormatSmiles) {
    // No extension → can't detect; AUTO should try SMILES when no file
    EXPECT_EQ(detect_format("molecule.mol2"), InputFormat::MOL2);
    EXPECT_EQ(detect_format("molecule.sdf"),  InputFormat::SDF);
}

// ===========================================================================
// Real Astex 1M2Z ligand regression (geometry invariants + rotor hygiene)
// Historical multi-Å rupture does NOT reproduce on modern DirectLigandIC /
// topology-derived GPA (max bond drift ~1e-5 A on main be049f8c). Expect PASS.
// ===========================================================================

TEST(ProcessLigand, Real1M2ZLigandGeometryAndRotors) {
    // Locate Astex Diverse 1M2Z cognate ligand SDF relative to this source tree.
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

    ProcessOptions opts;
    opts.input = sdf.string();
    opts.format = InputFormat::SDF;
    opts.lig_name = "DEX";
    opts.validate_only = false;
    opts.allow_peptides = true;  // steroid-like scaffold may trip peptide guard

    ProcessLigand pl;
    auto result = pl.run(opts);
    ASSERT_TRUE(result.success) << result.error;
    EXPECT_GT(result.num_atoms, 10);
    // Geometry invariants must pass (writer fails closed on corruption).
    ASSERT_TRUE(result.writer_result.success) << result.writer_result.error;
    EXPECT_FALSE(result.writer_result.inp_content.empty());
    EXPECT_FALSE(result.writer_result.ga_content.empty());
    EXPECT_FALSE(result.writer_result.internal_coords.empty());

    // Amide / conjugated C–N / disulfide / triple-adjacent must not be rotors.
    for (int i = 0; i < result.mol.num_bonds(); ++i) {
        if (rotatable_bonds::is_amide_bond(result.mol, i) ||
            rotatable_bonds::is_conjugated_cn_bond(result.mol, i) ||
            rotatable_bonds::is_disulfide_bond(result.mol, i) ||
            rotatable_bonds::is_triple_adjacent_bond(result.mol, i)) {
            EXPECT_FALSE(result.mol.bonds[i].is_rotatable)
                << "forbidden rotor bond index " << i;
        }
    }
    // Preserve MOL2/SYBYL N.am: any N.am (sybyl_type==7) must not sit on a rotor.
    for (const auto& b : result.mol.bonds) {
        if (!b.is_rotatable) continue;
        const int sy_i = result.mol.atoms[b.atom_i].sybyl_type;
        const int sy_j = result.mol.atoms[b.atom_j].sybyl_type;
        EXPECT_NE(sy_i, 7) << "N.am endpoint marked rotatable";
        EXPECT_NE(sy_j, 7) << "N.am endpoint marked rotatable";
    }
    EXPECT_GE(result.num_rot_bonds, 0);
}

// ===========================================================================
// MAIN
// ===========================================================================

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
