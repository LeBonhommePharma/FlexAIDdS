// tests/test_top_helpers.cpp
// First tests for logic that lived in LIB/top.cpp.
//
// top.cpp defines main() and exports nothing else, so none of this was
// reachable from a test target. The three functions under test moved to
// LIB/top_helpers.cpp unchanged; these are their first tests in the history of
// the file.
//
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>

#include "../LIB/top_helpers.h"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>

namespace fs = std::filesystem;

// ── sybyl_name_to_canonical_vct ────────────────────────────────────────────
//
// The mappings that are NOT the naive row for the type are the ones worth
// pinning. Each of these five redirects exists because the type's own row is
// all-zero in MC_st0r5.2_6.dat, so the naive mapping scored that chemistry as
// nothing at all against every partner. "Tidying" the table back to the
// obvious row numbers is a silent scoring regression -- sulfones and
// sulfonamides, among the most common motifs in drug-like ligands, go
// invisible to the complementarity function with no test failing and no
// crash. That is what this table exists to prevent.
struct SybylCase { const char* name; int expected; const char* why; };

TEST(SybylCanonicalVct, DeadRowSubstitutionsAreDeliberate) {
    const SybylCase cases[] = {
        {"O.ar",  14, "row 16 all-zero; ring O is ether-like -> O.3"},
        {"S.O2",  19, "row 20 all-zero; sulfone/sulfonamide -> S.O"},
        {"S.o2",  19, "lower-case spelling must map identically"},
        {"S.ar",  18, "row 21 all-zero; thiophene S is a thioether -> S.3"},
        {"Se",    18, "row 27 all-zero; Se arrives as selenomethionine -> S.3"},
        {"I",     25, "type-26 row has 3 live entries -> Br"},
        {"N.2",   10, "sp2 imine is an acceptor; N.am would reverse the H-bond sign"},
        {"N.3",   11, "row 8 all-zero -> N.am"},
        {"C.1",    2, "sp C rare in PDB sites; C.2 better sampled"},
    };
    for (const auto& c : cases) {
        EXPECT_EQ(sybyl_name_to_canonical_vct(c.name), c.expected)
            << c.name << ": " << c.why;
    }
}

TEST(SybylCanonicalVct, StraightforwardTypesMapToTheirOwnRow) {
    const SybylCase cases[] = {
        {"C.2", 2, ""},   {"C.3", 3, ""},   {"C.ar", 4, ""},  {"C.cat", 5, ""},
        {"N.1", 6, ""},   {"N.4", 9, ""},   {"N.ar", 10, ""}, {"N.am", 11, ""},
        {"N.pl3", 12, ""},{"O.2", 13, ""},  {"O.3", 14, ""},  {"O.co2", 15, ""},
        {"S.2", 17, ""},  {"S.3", 18, ""},  {"S.O", 19, ""},  {"P.3", 22, ""},
        {"F", 23, ""},    {"Cl", 24, ""},   {"Br", 25, ""},
        {"Mg", 28, ""},   {"Sr", 29, ""},   {"Cu", 30, ""},   {"Mn", 31, ""},
        {"Hg", 32, ""},   {"Cd", 33, ""},   {"Ni", 34, ""},   {"Zn", 35, ""},
        {"Ca", 36, ""},   {"Fe", 37, ""},   {"Co", 38, ""},   {"Co.oh", 38, ""},
    };
    for (const auto& c : cases) {
        EXPECT_EQ(sybyl_name_to_canonical_vct(c.name), c.expected) << c.name;
    }
}

TEST(SybylCanonicalVct, HydrogenAndUnknownAreDummy) {
    EXPECT_EQ(sybyl_name_to_canonical_vct("H"), FA_TYPE_DUMMY);
    EXPECT_EQ(sybyl_name_to_canonical_vct("X"), FA_TYPE_DUMMY);
    EXPECT_EQ(sybyl_name_to_canonical_vct(""), FA_TYPE_DUMMY);
    EXPECT_EQ(sybyl_name_to_canonical_vct("Du"), FA_TYPE_DUMMY);
    // Case matters for everything except the sulfoxide/sulfone pair, which is
    // explicitly spelled both ways. This pins that asymmetry rather than
    // leaving a future reader to guess whether it was intended.
    EXPECT_EQ(sybyl_name_to_canonical_vct("c.3"), FA_TYPE_DUMMY);
    EXPECT_EQ(sybyl_name_to_canonical_vct("S.o"), 19);
}

// ── is_valid_pdb_id ────────────────────────────────────────────────────────
//
// This validator runs on a user-supplied string that is then interpolated into
// an RCSB download URL and into a cache directory path. Alphanumeric-only is
// what keeps separators and traversal sequences out of both, so the rejection
// cases below are the load-bearing half of the test.
TEST(IsValidPdbId, AcceptsClassicAndExtendedCodes) {
    EXPECT_TRUE(is_valid_pdb_id("1P62"));
    EXPECT_TRUE(is_valid_pdb_id("1t40"));
    EXPECT_TRUE(is_valid_pdb_id("7XYZ9"));
    EXPECT_TRUE(is_valid_pdb_id("PDB12345"));  // 8 chars, the upper bound
}

TEST(IsValidPdbId, RejectsOutOfRangeLengths) {
    EXPECT_FALSE(is_valid_pdb_id(""));
    EXPECT_FALSE(is_valid_pdb_id("1P6"));        // 3, below the bound
    EXPECT_FALSE(is_valid_pdb_id("PDB123456"));  // 9, above the bound
}

TEST(IsValidPdbId, RejectsPathAndShellMetacharacters) {
    // Each of these would otherwise reach a URL or a filesystem path.
    const char* hostile[] = {
        "../..", "1P62/..", "..%2f", "1P62;rm", "1P62 x", "1P62\n",
        "$(id)", "`id`", "1P62|x", "1P62&x", "a/b/c", "a\\b",
    };
    for (const char* s : hostile) {
        EXPECT_FALSE(is_valid_pdb_id(s)) << "accepted hostile input: " << s;
    }
}

// ── detect_file_role ───────────────────────────────────────────────────────
namespace {
fs::path make_temp_dir() {
    const auto dir = fs::temp_directory_path() / "flexaids_top_helpers";
    fs::remove_all(dir);
    fs::create_directories(dir);
    return dir;
}

void write_file(const fs::path& p, const std::string& body) {
    std::ofstream o(p);
    o << body;
}
}  // namespace

TEST(DetectFileRole, ClassifiesByExtension) {
    const auto dir = make_temp_dir();
    struct { const char* name; const char* role; } cases[] = {
        {"lig.mol2", "ligand"},  {"lig.sdf", "ligand"},   {"lig.mol", "ligand"},
        {"lib.smi", "ligand"},   {"lib.smiles", "ligand"},
        {"cfg.json", "config"},  {"rec.cif", "receptor"}, {"rec.mmcif", "receptor"},
        {"old.inp", "legacy"},   {"old.dat", "legacy"},
    };
    for (const auto& c : cases) {
        const auto p = dir / c.name;
        write_file(p, "x\n");
        EXPECT_EQ(detect_file_role(p.string()), c.role) << c.name;
    }
    fs::remove_all(dir);
}

// The .pdb branch is the only one that reads content rather than the name, and
// the ATOM > 20 threshold is the whole of the receptor/ligand decision.
TEST(DetectFileRole, PdbIsDisambiguatedByAtomRecordCount) {
    const auto dir = make_temp_dir();

    std::string receptor;
    for (int i = 0; i < 25; ++i) receptor += "ATOM      1  CA  ALA A   1\n";
    const auto rec = dir / "rec.pdb";
    write_file(rec, receptor);
    EXPECT_EQ(detect_file_role(rec.string()), "receptor");

    std::string ligand;
    for (int i = 0; i < 12; ++i) ligand += "HETATM    1  C1  LIG A 900\n";
    const auto lig = dir / "lig.pdb";
    write_file(lig, ligand);
    EXPECT_EQ(detect_file_role(lig.string()), "ligand");

    // Exactly at the boundary: 20 ATOM records is NOT > 20, so a file with
    // only ATOM records and no HETATM falls through to the fallback branch.
    std::string boundary;
    for (int i = 0; i < 20; ++i) boundary += "ATOM      1  CA  ALA A   1\n";
    const auto edge = dir / "edge.pdb";
    write_file(edge, boundary);
    EXPECT_EQ(detect_file_role(edge.string()), "receptor");

    const auto empty = dir / "empty.pdb";
    write_file(empty, "REMARK nothing here\n");
    EXPECT_EQ(detect_file_role(empty.string()), "unknown");

    fs::remove_all(dir);
}

TEST(DetectFileRole, DirectoryIsALigandLibrary) {
    const auto dir = make_temp_dir();
    const auto sub = dir / "ligands";
    fs::create_directories(sub);
    EXPECT_EQ(detect_file_role(sub.string()), "ligand");
    fs::remove_all(dir);
}

// A nonexistent path may still be a SMILES string. The heuristic keys on
// chemistry characters and the absence of path separators, so these two cases
// bound it from both sides.
TEST(DetectFileRole, NonexistentPathMayBeSmiles) {
    EXPECT_EQ(detect_file_role("CC(=O)Oc1ccccc1C(=O)O"), "smiles");  // aspirin
    EXPECT_EQ(detect_file_role("c1ccccc1"), "smiles");
    EXPECT_EQ(detect_file_role(""), "unknown");
    EXPECT_EQ(detect_file_role("/no/such/file.xyz"), "unknown");
    // Contains a separator, so it is a path that does not exist -- not SMILES,
    // even though it carries chemistry characters.
    EXPECT_EQ(detect_file_role("dir/CC(=O)O"), "unknown");
}
