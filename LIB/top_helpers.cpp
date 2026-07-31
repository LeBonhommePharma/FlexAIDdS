#include "top_helpers.h"

#include <cctype>
#include <cstdio>
#include <cstring>
#include <filesystem>

// Bodies moved verbatim from LIB/top.cpp. No behaviour change -- the only edits
// are dropping `static` and the FA_TYPE_DUMMY constant moving to the header.

int sybyl_name_to_canonical_vct(const char* s) {
	if (!strcmp(s, "C.1"))   return 2;   // C.2 — sp C rare in PDB sites, C.2 better sampled
	if (!strcmp(s, "C.2"))   return 2;
	if (!strcmp(s, "C.3"))   return 3;
	if (!strcmp(s, "C.ar"))  return 4;
	if (!strcmp(s, "C.cat")) return 5;
	if (!strcmp(s, "N.1"))   return 6;
	if (!strcmp(s, "N.2"))   return 10;  // N.ar — sp2 imine is an acceptor; N.am (donor) reversed the H-bond sign
	if (!strcmp(s, "N.3"))   return 11;  // N.am — row 8 is all-zero in MC_st0r5.2_6.dat; matches Mol2Reader.cpp:44
	if (!strcmp(s, "N.4"))   return 9;
	if (!strcmp(s, "N.ar"))  return 10;
	if (!strcmp(s, "N.am"))  return 11;
	if (!strcmp(s, "N.pl3")) return 12;
	if (!strcmp(s, "O.2"))   return 13;
	if (!strcmp(s, "O.3"))   return 14;
	if (!strcmp(s, "O.co2")) return 15;
	// O.ar → O.3. Row 16 is all-zero in MC_st0r5.2_6.dat, so a furan / oxazole /
	// benzofuran oxygen scored exactly nothing against every partner. An
	// aromatic ring oxygen is divalent with no labile H — ether-like — so O.3
	// (row 14, 22 live entries) is the correct live surrogate, not the carbonyl
	// row O.2. Geometry is unaffected: both the vH recipe and the implicit-H
	// count for row 14 are gated on heavy_bonds<=1, and a ring O always has 2.
	if (!strcmp(s, "O.ar"))  return 14;
	if (!strcmp(s, "S.2"))   return 17;
	if (!strcmp(s, "S.3"))   return 18;
	if (!strcmp(s, "S.O") || !strcmp(s, "S.o"))   return 19;
	// S.O2 → S.O. Row 20 is all-zero, so every sulfone and sulfonamide sulfur —
	// one of the most common motifs in drug-like ligands — was invisible to the
	// scorer. S.O (row 19, 12 live entries) is the nearest live chemistry:
	// oxidised, tetrahedral, strongly polarised S=O.
	if (!strcmp(s, "S.O2") || !strcmp(s, "S.o2")) return 19;
	// S.ar → S.3. Row 21 is all-zero; a thiophene / thiazole sulfur is divalent
	// with no labile H, so the thioether row 18 (20 live entries) is the correct
	// live surrogate. Same heavy_bonds<=1 geometry gate as O.ar above.
	if (!strcmp(s, "S.ar"))  return 18;
	if (!strcmp(s, "P.3"))   return 22;
	if (!strcmp(s, "F"))     return 23;
	if (!strcmp(s, "Cl"))    return 24;
	if (!strcmp(s, "Br"))    return 25;
	if (!strcmp(s, "I"))     return 25;  // BR — iodo near-absent from PDB training; I/type-26 row has only 3 live entries
	// Se → S.3. Row 27 is all-zero. Selenium reaches the scorer almost entirely
	// as selenomethionine (MSE), a methionine surrogate used for phasing, so
	// scoring it on the thioether row reproduces the chemistry it stands in for
	// instead of making the side chain invisible.
	if (!strcmp(s, "Se"))    return 18;
	if (!strcmp(s, "Mg"))    return 28;
	if (!strcmp(s, "Sr"))    return 29;
	if (!strcmp(s, "Cu"))    return 30;
	if (!strcmp(s, "Mn"))    return 31;
	if (!strcmp(s, "Hg"))    return 32;
	if (!strcmp(s, "Cd"))    return 33;
	if (!strcmp(s, "Ni"))    return 34;
	if (!strcmp(s, "Zn"))    return 35;
	if (!strcmp(s, "Ca"))    return 36;
	if (!strcmp(s, "Fe"))    return 37;
	if (!strcmp(s, "Co.oh") || !strcmp(s, "Co")) return 38;
	return FA_TYPE_DUMMY; // H ("H") and anything unknown ("X")
}

// ── Idiotproof file role detection ──────────────────────────────────────────
// Returns: "receptor", "ligand", "config", "smiles", or "unknown"
std::string detect_file_role(const std::string& path) {
	// Not a file? Might be a SMILES string
	if (!std::filesystem::exists(path)) {
		// SMILES strings contain typical chemistry chars, no path separators
		if (!path.empty() &&
		    path.find('/') == std::string::npos &&
		    path.find('\\') == std::string::npos &&
		    (path.find('(') != std::string::npos ||
		     path.find('=') != std::string::npos ||
		     path.find('#') != std::string::npos ||
		     path.find('c') != std::string::npos ||
		     path.find('C') != std::string::npos ||
		     path.find('N') != std::string::npos ||
		     path.find('O') != std::string::npos)) {
			return "smiles";
		}
		return "unknown";
	}

	std::string ext;
	{
		auto dot = path.rfind('.');
		if (dot != std::string::npos) {
			ext = path.substr(dot);
			for (auto& c : ext) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
		}
	}

	// Ligand formats (single molecule or library)
	if (ext == ".mol2" || ext == ".sdf" || ext == ".mol") return "ligand";

	// SMILES file = ligand library
	if (ext == ".smi" || ext == ".smiles") return "ligand";

	// Config formats
	if (ext == ".json") return "config";

	// CIF/mmCIF — receptor (PDB archive format)
	if (ext == ".cif" || ext == ".mmcif") return "receptor";

	// Directory of ligand files = ligand library
	if (std::filesystem::is_directory(path)) return "ligand";

	// PDB could be receptor or ligand — peek at content
	if (ext == ".pdb" || ext == ".ent") {
		FILE* fp = fopen(path.c_str(), "r");
		if (!fp) return "unknown";
		int atom_count = 0;
		int hetatm_count = 0;
		char line[256];
		while (fgets(line, sizeof(line), fp) && (atom_count + hetatm_count) < 200) {
			if (strncmp(line, "ATOM  ", 6) == 0) atom_count++;
			else if (strncmp(line, "HETATM", 6) == 0) hetatm_count++;
		}
		fclose(fp);
		// Receptor: many ATOM records. Ligand PDB: mostly HETATM, few atoms.
		if (atom_count > 20) return "receptor";
		if (hetatm_count > 0 && atom_count <= 20) return "ligand";
		if (atom_count > 0) return "receptor"; // fallback
		return "unknown";
	}

	// Legacy input files
	if (ext == ".inp" || ext == ".dat") return "legacy";

	return "unknown";
}

/// Validate an RCSB PDB ID (classic 4-char alphanumeric; allow longer alphanumeric
/// accession-style codes used by some modern deposits, 4–8 chars).
bool is_valid_pdb_id(const std::string& id) {
	if (id.size() < 4 || id.size() > 8) return false;
	for (unsigned char c : id) {
		if (!std::isalnum(c)) return false;
	}
	return true;
}
