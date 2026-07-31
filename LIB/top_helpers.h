// Helpers lifted out of LIB/top.cpp so they can be tested.
//
// top.cpp defines main() and exports no other symbol, so a test target can
// neither link it (gtest_main supplies its own main) nor reach anything inside
// it (every helper was static). These three were pure functions of their
// arguments trapped in that TU; moving them here changes no behaviour and makes
// them addressable.

#ifndef FLEXAID_TOP_HELPERS_H
#define FLEXAID_TOP_HELPERS_H

#include <string>

// VCT type index for "not scored against the matrix" -- hydrogen and anything
// unrecognised.
constexpr int FA_TYPE_DUMMY = 39;

// SYBYL atom-type name -> canonical VCT row index.
//
// Several mappings are deliberately NOT the naive row for the type: O.ar, S.O2,
// S.ar, Se and I are redirected to live rows because their own rows are
// all-zero in MC_st0r5.2_6.dat, which made those chemistries invisible to the
// scorer. Mol2Reader.cpp carries a second copy of this table that must agree --
// see the parity test in tests/test_top_helpers.cpp.
int sybyl_name_to_canonical_vct(const char* s);

// RCSB PDB ID validation, applied to a user-supplied string before it is used
// to build a download URL and a cache directory path. Alphanumeric-only is what
// keeps path separators and traversal sequences out of both.
bool is_valid_pdb_id(const std::string& id);

// Classify a command-line argument as "receptor" / "ligand" / "config" /
// "legacy" / "smiles" / "unknown". Touches the filesystem: a path that does not
// exist may still be classified as a SMILES string, and .pdb is disambiguated
// by counting ATOM vs HETATM records.
std::string detect_file_role(const std::string& path);

#endif
