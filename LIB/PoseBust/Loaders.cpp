// Loaders.cpp — Clean-room SDF/PDB molecular I/O for PoseBust
//
// Apache-2.0. Implemented from first principles; no posebusters/RDKit code.
#include "Loaders.h"
#include "PdbCoords.h"

#include <array>
#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace flexaids::posebust {
namespace {

void set_err(std::string* err, const std::string& msg) {
    if (err != nullptr) {
        *err = msg;
    }
}

std::string trim_copy(std::string s) {
    auto not_space = [](unsigned char c) { return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
    return s;
}

std::string to_upper(std::string s) {
    for (char& c : s) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    return s;
}

// Safe fixed-width field extract (0-based start, length).
std::string field(const std::string& line, std::size_t start, std::size_t len) {
    if (start >= line.size()) {
        return {};
    }
    return trim_copy(line.substr(start, std::min(len, line.size() - start)));
}

bool parse_int_field(const std::string& s, int& out) {
    if (s.empty()) {
        return false;
    }
    try {
        std::size_t idx = 0;
        out = std::stoi(s, &idx);
        return idx > 0;
    } catch (...) {
        return false;
    }
}

bool parse_float_field(const std::string& s, float& out) {
    if (s.empty()) {
        return false;
    }
    try {
        std::size_t idx = 0;
        out = std::stof(s, &idx);
        // Require whole field consumption after trim so " -0.63-80.2" style
        // compact spills fail here and are handled by span regex parser.
        while (idx < s.size() && std::isspace(static_cast<unsigned char>(s[idx])))
            ++idx;
        return idx == s.size();
    } catch (...) {
        return false;
    }
}

// Standard 20 amino acids + common protein placeholders.
bool is_standard_aa(const std::string& res) {
    static const std::unordered_set<std::string> kAA = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "UNK", "MSE", "SEC", "PYL", "ASX", "GLX",
    };
    return kAA.count(to_upper(res)) > 0;
}

// Solvent / buffer / common counter-ions to skip for ligand loading.
bool is_solvent_or_ion(const std::string& res) {
    static const std::unordered_set<std::string> kSkip = {
        "WAT", "HOH", "H2O", "TIP", "TIP3", "SOL", "DOD",
        "SO4", "PO4", "NO3", "CO3",
        "NA",  "NA+", "CL",  "CL-", "K",   "K+",  "MG",  "MG2", "CA",  "CA2",
        "ZN",  "ZN2", "MN",  "FE",  "CU",  "NI",
        "GOL", "EDO", "PEG", "PGE", "ACT", "ACY", "DMS", "DMSO", "MPD",
        "TRS", "EPE", "HEZ", "BOG", "BME", "FMT", "ACE",
    };
    return kSkip.count(to_upper(res)) > 0;
}

// Common protein cofactors / crystallographic junk that are HETATM but are NOT
// the FlexAID docking ligand (1G9V: 128×HEM vs 25 docked ligand atoms).
bool is_common_cofactor(const std::string& res) {
    static const std::unordered_set<std::string> kCof = {
        "HEM", "HEC", "HEA", "HEB", "HEO", "HAS", "HDD", "HDE",
        "NAD", "NAP", "NDP", "NADP", "FAD", "FMN", "ADP", "ATP", "AMP",
        "GDP", "GTP", "GMP", "UDP", "UTP", "UMP", "CDP", "CTP", "CMP",
        "COA", "ACO", "SAM", "SAH", "PLP", "TPP", "THP", "B12", "CNC",
        "FES", "SF4", "F3S", "CLM", "CLA", "BCL", "CHL",
    };
    return kCof.count(to_upper(res)) > 0;
}

bool is_hetatm_junk(const std::string& res) {
    return is_solvent_or_ion(res) || is_common_cofactor(res);
}

// Element from PDB atom name / element column.
// Hardens against FlexAID misaligned writes: "Cl0" with shifted element "L",
// or hydrogen written as Du with name H*.
std::string element_from_pdb(const std::string& atom_name, const std::string& elem_col) {
    const std::string name = trim_copy(atom_name);
    const std::string name_up = to_upper(name);

    // Prefer two-letter elements from the atom name when unambiguous (Cl, Br, …).
    // This recovers Cl when a short name shifts the element column to "L".
    auto element_from_name = [&]() -> std::string {
        std::size_t i = 0;
        while (i < name.size() && std::isdigit(static_cast<unsigned char>(name[i]))) {
            ++i;
        }
        if (i >= name.size()) {
            return {};
        }
        std::string e;
        e.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(name[i]))));
        if (i + 1 < name.size() && std::isalpha(static_cast<unsigned char>(name[i + 1]))) {
            const char c2 =
                static_cast<char>(std::tolower(static_cast<unsigned char>(name[i + 1])));
            const std::string two = std::string(1, e[0]) + c2;
            if (atomic_number(two) > 0) {
                return two;
            }
        }
        return e;
    };

    const std::string from_name = element_from_name();

    if (!elem_col.empty()) {
        std::string e = trim_copy(elem_col);
        if (!e.empty()) {
            e[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(e[0])));
            for (std::size_t i = 1; i < e.size(); ++i) {
                e[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(e[i])));
            }
            if (e.size() == 2 && std::isupper(static_cast<unsigned char>(e[1]))) {
                e[1] = static_cast<char>(std::tolower(static_cast<unsigned char>(e[1])));
            }
            // Recover Cl/Br when element column is a single garbage letter (e.g. "L").
            if (e.size() == 1 && from_name.size() == 2 && atomic_number(from_name) > 0) {
                return from_name;
            }
            // Du/H dummy hydrogens written with element Du
            if ((to_upper(e) == "DU" || to_upper(e) == "D") && !name_up.empty() &&
                name_up[0] == 'H') {
                return "H";
            }
            if (atomic_number(e) > 0 || to_upper(e) == "DU") {
                return e;
            }
        }
    }

    if (!from_name.empty()) {
        return from_name;
    }
    return "C";
}

// Returns false if the atom should be dropped entirely (orphan dummy).
bool finalize_atom(Atom& a, const std::string& atom_name = {}) {
    const std::string name_up = to_upper(trim_copy(atom_name));
    const std::string el_up = to_upper(a.element);
    // Recover H from Du dummy labels left by some writers (H… Du).
    if ((el_up == "DU" || el_up == "D") && !name_up.empty() && name_up[0] == 'H') {
        a.element = "H";
    }
    a.atomic_num = atomic_number(a.element);
    const std::string el = to_upper(a.element);
    a.is_h = (a.atomic_num == 1) || (el == "H") || (el == "D") || (el == "T");
    // Drop pure dummy / unknown non-hydrogen placeholders (never count as heavy).
    if (a.atomic_num == 0 && !a.is_h) {
        return false;
    }
    return true;
}

// Infer single bonds from covalent radii (1.2 * (r1+r2)).
void infer_bonds(Molecule& mol) {
    mol.bonds.clear();
    const std::size_t n = mol.atoms.size();
    for (std::size_t i = 0; i < n; ++i) {
        const Atom& ai = mol.atoms[i];
        const float ri = covalent_radius(ai.atomic_num);
        const Vec3  pi{ai.x, ai.y, ai.z};
        for (std::size_t j = i + 1; j < n; ++j) {
            const Atom& aj = mol.atoms[j];
            const float rj = covalent_radius(aj.atomic_num);
            const float max_d = 1.2f * (ri + rj);
            const Vec3  pj{aj.x, aj.y, aj.z};
            if (dist2(pi, pj) <= max_d * max_d) {
                // Skip unphysical H–H bonds unless they are the only contacts
                // (still allow; order always 1 for inferred bonds).
                mol.bonds.push_back(Bond{static_cast<int>(i), static_cast<int>(j), 1});
            }
        }
    }
    mol.build_adjacency();
}

struct PdbAtomRec {
    bool        is_hetatm = false;
    int         serial    = 0;
    std::string name;
    std::string resname;
    std::string chain;
    int         resseq = 0;
    float       x = 0.f, y = 0.f, z = 0.f;
    std::string element;
};

bool parse_pdb_atom_line(const std::string& line, PdbAtomRec& rec) {
    if (line.size() < 54) {
        return false;
    }
    const std::string tag = field(line, 0, 6);
    if (tag != "ATOM" && tag != "HETATM") {
        return false;
    }
    rec.is_hetatm = (tag == "HETATM");
    if (!parse_int_field(field(line, 6, 5), rec.serial)) {
        rec.serial = 0;
    }
    rec.name    = field(line, 12, 4);
    rec.resname = field(line, 17, 3);
    rec.chain   = field(line, 21, 1);
    if (!parse_int_field(field(line, 22, 4), rec.resseq)) {
        rec.resseq = 0;
    }
    // Shared strict finite decoder (identical to DatasetRunner RMSD path).
    {
        std::array<float, 3> xyz{};
        if (!flexaids::pdb_coords::parse_xyz_span(line, xyz)) {
            return false;
        }
        rec.x = xyz[0];
        rec.y = xyz[1];
        rec.z = xyz[2];
    }
    // Element columns 77–78 (1-based) → 76–77 0-based; may be absent.
    rec.element = (line.size() >= 78) ? field(line, 76, 2) : std::string{};
    return true;
}

void molecule_from_pdb_recs(const std::vector<PdbAtomRec>& recs, Molecule& out,
                            const std::string& name) {
    out = Molecule{};
    out.name = name;
    out.atoms.reserve(recs.size());
    for (const PdbAtomRec& r : recs) {
        Atom a;
        a.id      = static_cast<int>(out.atoms.size()) + 1;
        a.element = element_from_pdb(r.name, r.element);
        a.x       = r.x;
        a.y       = r.y;
        a.z       = r.z;
        if (!finalize_atom(a, r.name)) {
            continue;
        }
        out.atoms.push_back(std::move(a));
    }
    // Re-number ids after filtering
    for (std::size_t i = 0; i < out.atoms.size(); ++i) {
        out.atoms[i].id = static_cast<int>(i) + 1;
    }
    infer_bonds(out);
}

}  // namespace

// ─── Public helpers ──────────────────────────────────────────────────────────

int atomic_number(const std::string& elem) {
    if (elem.empty()) {
        return 0;
    }
    // Normalize: first upper, rest lower.
    std::string e;
    e.reserve(elem.size());
    e.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(elem[0]))));
    for (std::size_t i = 1; i < elem.size(); ++i) {
        if (std::isalpha(static_cast<unsigned char>(elem[i]))) {
            e.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(elem[i]))));
        }
    }

    // Deuterium / tritium → H
    if (e == "D" || e == "T") {
        return 1;
    }

    static const std::pair<const char*, int> kTable[] = {
        {"H", 1},   {"B", 5},   {"C", 6},   {"N", 7},   {"O", 8},   {"F", 9},
        {"Si", 14}, {"P", 15},  {"S", 16},  {"Cl", 17}, {"Br", 35}, {"I", 53},
        // Extra common ligand metals / elements (harmless, radii provided below)
        {"Na", 11}, {"Mg", 12}, {"K", 19},  {"Ca", 20}, {"Mn", 25}, {"Fe", 26},
        {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30}, {"Se", 34}, {"As", 33},
    };
    for (const auto& [sym, z] : kTable) {
        if (e == sym) {
            return z;
        }
    }
    return 0;
}

float covalent_radius(int Z) {
    // Single-bond covalent radii (Å), Cordero et al. 2008-ish simplified set.
    switch (Z) {
        case 1:  return 0.31f;
        case 5:  return 0.84f;
        case 6:  return 0.76f;
        case 7:  return 0.71f;
        case 8:  return 0.66f;
        case 9:  return 0.57f;
        case 11: return 1.66f;
        case 12: return 1.41f;
        case 14: return 1.11f;
        case 15: return 1.07f;
        case 16: return 1.05f;
        case 17: return 1.02f;
        case 19: return 2.03f;
        case 20: return 1.76f;
        case 25: return 1.39f;
        case 26: return 1.32f;
        case 27: return 1.26f;
        case 28: return 1.24f;
        case 29: return 1.32f;
        case 30: return 1.22f;
        case 33: return 1.19f;
        case 34: return 1.20f;
        case 35: return 1.20f;
        case 53: return 1.39f;
        default: return 0.77f;  // carbon-like fallback
    }
}

// ─── SDF V2000 ───────────────────────────────────────────────────────────────

bool load_sdf(const std::string& path, Molecule& out, std::string* err) {
    out = Molecule{};
    std::ifstream in(path);
    if (!in) {
        set_err(err, "load_sdf: cannot open '" + path + "'");
        return false;
    }

    std::string line;
    // Header: line 1 name, line 2 program, line 3 comment
    if (!std::getline(in, line)) {
        set_err(err, "load_sdf: empty file");
        return false;
    }
    out.name = trim_copy(line);
    if (!std::getline(in, line)) {
        set_err(err, "load_sdf: truncated header (program line)");
        return false;
    }
    if (!std::getline(in, line)) {
        set_err(err, "load_sdf: truncated header (comment line)");
        return false;
    }

    // Counts line
    if (!std::getline(in, line)) {
        set_err(err, "load_sdf: missing counts line");
        return false;
    }
    // V2000: aaabbb… where aaa = nAtoms, bbb = nBonds (cols 1–3, 4–6)
    int n_atoms = 0;
    int n_bonds = 0;
    if (line.size() >= 6) {
        if (!parse_int_field(field(line, 0, 3), n_atoms) ||
            !parse_int_field(field(line, 3, 3), n_bonds)) {
            set_err(err, "load_sdf: bad counts line: '" + trim_copy(line) + "'");
            return false;
        }
    } else {
        // Whitespace-separated fallback
        std::istringstream iss(line);
        if (!(iss >> n_atoms >> n_bonds)) {
            set_err(err, "load_sdf: unparsable counts line");
            return false;
        }
    }
    if (n_atoms < 0 || n_bonds < 0) {
        set_err(err, "load_sdf: negative atom/bond count");
        return false;
    }

    out.atoms.reserve(static_cast<std::size_t>(n_atoms));
    for (int i = 0; i < n_atoms; ++i) {
        if (!std::getline(in, line)) {
            set_err(err, "load_sdf: truncated atom block at atom " + std::to_string(i + 1));
            return false;
        }
        Atom a;
        a.id = i + 1;
        // V2000 atom: xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aa…
        // cols 1–10 x, 11–20 y, 21–30 z, 32–34 symbol (1-based) → 0-based 0,10,20,31
        if (line.size() >= 33) {
            if (!parse_float_field(field(line, 0, 10), a.x) ||
                !parse_float_field(field(line, 10, 10), a.y) ||
                !parse_float_field(field(line, 20, 10), a.z)) {
                set_err(err, "load_sdf: bad coordinates on atom line " + std::to_string(i + 1));
                return false;
            }
            a.element = field(line, 31, 3);
        } else {
            std::istringstream iss(line);
            if (!(iss >> a.x >> a.y >> a.z >> a.element)) {
                set_err(err, "load_sdf: unparsable atom line " + std::to_string(i + 1));
                return false;
            }
        }
        if (a.element.empty()) {
            set_err(err, "load_sdf: empty element on atom " + std::to_string(i + 1));
            return false;
        }
        // Normalize element capitalization
        a.element[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(a.element[0])));
        for (std::size_t k = 1; k < a.element.size(); ++k) {
            a.element[k] =
                static_cast<char>(std::tolower(static_cast<unsigned char>(a.element[k])));
        }
        finalize_atom(a);
        out.atoms.push_back(std::move(a));
    }

    out.bonds.reserve(static_cast<std::size_t>(n_bonds));
    for (int i = 0; i < n_bonds; ++i) {
        if (!std::getline(in, line)) {
            set_err(err, "load_sdf: truncated bond block at bond " + std::to_string(i + 1));
            return false;
        }
        int a1 = 0, a2 = 0, order = 1;
        if (line.size() >= 9) {
            if (!parse_int_field(field(line, 0, 3), a1) ||
                !parse_int_field(field(line, 3, 3), a2) ||
                !parse_int_field(field(line, 6, 3), order)) {
                set_err(err, "load_sdf: bad bond line " + std::to_string(i + 1));
                return false;
            }
        } else {
            std::istringstream iss(line);
            if (!(iss >> a1 >> a2 >> order)) {
                set_err(err, "load_sdf: unparsable bond line " + std::to_string(i + 1));
                return false;
            }
        }
        if (a1 < 1 || a2 < 1 || a1 > n_atoms || a2 > n_atoms) {
            set_err(err, "load_sdf: bond " + std::to_string(i + 1) + " atom index out of range");
            return false;
        }
        if (order < 1) {
            order = 1;
        }
        // MDL: 1 single, 2 double, 3 triple, 4 aromatic; clamp unknowns to 1.
        if (order > 4) {
            order = 1;
        }
        out.bonds.push_back(Bond{a1 - 1, a2 - 1, order});
    }

    out.build_adjacency();
    return true;
}

bool write_sdf(const Molecule& mol, const std::string& path, std::string* err) {
    std::ofstream out(path);
    if (!out) {
        set_err(err, "write_sdf: cannot open '" + path + "' for write");
        return false;
    }

    // Header
    out << (mol.name.empty() ? "PoseBust" : mol.name) << "\n";
    out << "  FlexAIDdS PoseBust\n";
    out << "\n";

    const int n_atoms = static_cast<int>(mol.atoms.size());
    const int n_bonds = static_cast<int>(mol.bonds.size());
    // Counts line (V2000)
    char counts[64];
    std::snprintf(counts, sizeof(counts), "%3d%3d  0  0  0  0  0  0  0  0999 V2000\n", n_atoms,
                  n_bonds);
    out << counts;

    for (const Atom& a : mol.atoms) {
        char buf[128];
        std::snprintf(buf, sizeof(buf), "%10.4f%10.4f%10.4f %-3s 0  0  0  0  0  0  0  0  0  0  0  0\n",
                      static_cast<double>(a.x), static_cast<double>(a.y),
                      static_cast<double>(a.z), a.element.c_str());
        out << buf;
    }

    for (const Bond& b : mol.bonds) {
        char buf[64];
        const int order = (b.order >= 1 && b.order <= 4) ? b.order : 1;
        std::snprintf(buf, sizeof(buf), "%3d%3d%3d  0  0  0  0\n", b.a + 1, b.b + 1, order);
        out << buf;
    }

    out << "M  END\n";
    out << "$$$$\n";
    if (!out) {
        set_err(err, "write_sdf: write failed for '" + path + "'");
        return false;
    }
    return true;
}

// ─── PDB ligand / protein ────────────────────────────────────────────────────

bool load_pdb_ligand(const std::string& path, Molecule& out, std::string* err) {
    out = Molecule{};
    std::ifstream in(path);
    if (!in) {
        set_err(err, "load_pdb_ligand: cannot open '" + path + "'");
        return false;
    }

    std::vector<PdbAtomRec> hetatm;
    std::vector<PdbAtomRec> non_protein_atom;
    std::string line;
    while (std::getline(in, line)) {
        // Strip CR if present
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        PdbAtomRec rec;
        if (!parse_pdb_atom_line(line, rec)) {
            continue;
        }
        if (rec.is_hetatm) {
            if (!is_hetatm_junk(rec.resname)) {
                hetatm.push_back(std::move(rec));
            }
        } else {
            // ATOM: keep non-standard residues as ligand fallback
            if (!is_standard_aa(rec.resname) && !is_hetatm_junk(rec.resname)) {
                non_protein_atom.push_back(std::move(rec));
            }
        }
    }

    const std::vector<PdbAtomRec>* chosen = nullptr;
    if (!hetatm.empty()) {
        chosen = &hetatm;
    } else if (!non_protein_atom.empty()) {
        chosen = &non_protein_atom;
    } else {
        set_err(err, "load_pdb_ligand: no ligand-like atoms found in '" + path + "'");
        return false;
    }

    molecule_from_pdb_recs(*chosen, out, path);
    if (out.atoms.empty()) {
        set_err(err, "load_pdb_ligand: empty molecule after filtering");
        return false;
    }
    return true;
}

// ─── FlexAID pose ligand (CONECT / optimizable residue) ─────────────────────

bool load_pdb_flexaid_ligand(const std::string& path, Molecule& out, std::string* err) {
    out = Molecule{};
    std::ifstream in(path);
    if (!in) {
        set_err(err, "load_pdb_flexaid_ligand: cannot open '" + path + "'");
        return false;
    }

    std::vector<PdbAtomRec> all_atoms;
    std::vector<std::vector<int>> conect;  // each row: serials
    std::string opt_resname;
    int opt_resseq = -1;

    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();

        // REMARK optimizable residue RQ3   1
        if (line.size() >= 28 && line.compare(0, 6, "REMARK") == 0 &&
            line.find("optimizable residue") != std::string::npos) {
            // Parse after the phrase
            const auto pos = line.find("optimizable residue");
            std::string rest = line.substr(pos + std::string("optimizable residue").size());
            std::istringstream iss(rest);
            iss >> opt_resname >> opt_resseq;
            continue;
        }

        if (line.size() >= 11 && line.compare(0, 6, "CONECT") == 0) {
            std::vector<int> row;
            // Fixed-width serials of 5 chars starting at column 7 (1-based) = index 6
            for (std::size_t i = 6; i + 5 <= line.size(); i += 5) {
                std::string tok = trim_copy(line.substr(i, 5));
                if (tok.empty()) break;
                try {
                    row.push_back(std::stoi(tok));
                } catch (...) {
                    break;
                }
            }
            if (!row.empty()) conect.push_back(std::move(row));
            continue;
        }

        PdbAtomRec rec;
        if (parse_pdb_atom_line(line, rec)) {
            all_atoms.push_back(std::move(rec));
        }
    }

    std::unordered_map<int, PdbAtomRec*> by_serial;
    by_serial.reserve(all_atoms.size() * 2);
    for (auto& a : all_atoms) {
        by_serial[a.serial] = &a;
    }

    // Path 1: CONECT-referenced atoms (FlexAID ligand serials 90001+)
    if (!conect.empty()) {
        std::vector<int> ordered_serials;
        std::unordered_set<int> seen;
        for (const auto& row : conect) {
            for (int s : row) {
                if (seen.insert(s).second) ordered_serials.push_back(s);
            }
        }
        std::sort(ordered_serials.begin(), ordered_serials.end());

        std::vector<PdbAtomRec> lig;
        lig.reserve(ordered_serials.size());
        for (int s : ordered_serials) {
            auto it = by_serial.find(s);
            if (it == by_serial.end()) continue;
            // Never pull protein AA into ligand even if CONECT is weird
            if (!it->second->is_hetatm && is_standard_aa(it->second->resname)) continue;
            if (is_common_cofactor(it->second->resname)) continue;
            lig.push_back(*it->second);
        }
        if (lig.size() >= 3) {
            molecule_from_pdb_recs(lig, out, path);
            // Prefer CONECT bond graph over pure distance inference when possible
            out.bonds.clear();
            std::unordered_map<int, int> serial_to_idx;
            for (std::size_t i = 0; i < lig.size(); ++i) {
                serial_to_idx[lig[i].serial] = static_cast<int>(i);
            }
            std::set<std::pair<int, int>> bond_set;
            for (const auto& row : conect) {
                if (row.empty()) continue;
                const int a_ser = row[0];
                auto ia = serial_to_idx.find(a_ser);
                if (ia == serial_to_idx.end()) continue;
                for (std::size_t k = 1; k < row.size(); ++k) {
                    auto ib = serial_to_idx.find(row[k]);
                    if (ib == serial_to_idx.end()) continue;
                    int i = ia->second, j = ib->second;
                    if (i > j) std::swap(i, j);
                    if (i == j) continue;
                    if (bond_set.insert({i, j}).second) {
                        out.bonds.push_back(Bond{i, j, 1});  // order unknown from CONECT
                    }
                }
            }
            out.build_adjacency();
            if (!out.atoms.empty()) return true;
        }
    }

    // Path 2: REMARK optimizable residue
    if (!opt_resname.empty()) {
        std::vector<PdbAtomRec> lig;
        for (const auto& a : all_atoms) {
            if (to_upper(a.resname) == to_upper(opt_resname) &&
                (opt_resseq < 0 || a.resseq == opt_resseq) &&
                !is_common_cofactor(a.resname)) {
                lig.push_back(a);
            }
        }
        // FlexAID may split ligand across resnames (Q3 + RQ3) with same serial block
        if (lig.size() < 3 && !conect.empty()) {
            // already tried CONECT
        }
        if (lig.size() >= 3) {
            molecule_from_pdb_recs(lig, out, path);
            if (!out.atoms.empty()) return true;
        }
    }

    // Path 3: fallback without cofactors (still may be wrong multi-HETATM)
    if (!load_pdb_ligand(path, out, err)) {
        set_err(err, "load_pdb_flexaid_ligand: CONECT/optimizable/fallback all failed for '" +
                         path + "'");
        return false;
    }
    return true;
}

bool assign_topology_from_reference(Molecule& pred, const Molecule& reference,
                                    std::string* err) {
    // Graph-identity topology transfer only. Do NOT:
    //  - orphan explicit H/Du atoms on either side (counts must match including H)
    //  - accept repeated-element positional nearest-neighbour mapping without
    //    element sequence identity (audit P1).
    if (pred.atoms.size() != reference.atoms.size() || pred.atoms.empty()) {
        set_err(err, "assign_topology_from_reference: atom-count mismatch (pred=" +
                         std::to_string(pred.atoms.size()) + " ref=" +
                         std::to_string(reference.atoms.size()) +
                         ") — including explicit H/Du");
        return false;
    }

    auto elements_match = [&](int pi, int ri) -> bool {
        const Atom& pa = pred.atoms[static_cast<std::size_t>(pi)];
        const Atom& ra = reference.atoms[static_cast<std::size_t>(ri)];
        if (pa.atomic_num > 0 && ra.atomic_num > 0 && pa.atomic_num == ra.atomic_num) {
            return true;
        }
        return to_upper(pa.element) == to_upper(ra.element);
    };

    // Require full-atom element sequence identity (graph order).
    for (std::size_t i = 0; i < pred.atoms.size(); ++i) {
        if (!elements_match(static_cast<int>(i), static_cast<int>(i))) {
            set_err(err,
                    "assign_topology_from_reference: element sequence mismatch at "
                    "atom " +
                        std::to_string(i) + " (pred=" + pred.atoms[i].element +
                        " ref=" + reference.atoms[i].element +
                        ") — refuse positional remapping");
            return false;
        }
    }

    // Identity map: same atom order + element sequence.
    std::vector<Bond> new_bonds;
    for (const Bond& b : reference.bonds) {
        if (b.a < 0 || b.b < 0 ||
            static_cast<std::size_t>(b.a) >= pred.atoms.size() ||
            static_cast<std::size_t>(b.b) >= pred.atoms.size()) {
            continue;
        }
        new_bonds.push_back(Bond{b.a, b.b, b.order});
    }
    if (new_bonds.empty()) {
        set_err(err, "assign_topology_from_reference: no transferable bonds");
        return false;
    }
    pred.bonds = std::move(new_bonds);
    pred.build_adjacency();
    return true;
}

bool load_pdb_protein_heavy(const std::string& path, Molecule& out, std::string* err) {
    out = Molecule{};
    std::ifstream in(path);
    if (!in) {
        set_err(err, "load_pdb_protein_heavy: cannot open '" + path + "'");
        return false;
    }

    std::vector<PdbAtomRec> protein;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        PdbAtomRec rec;
        if (!parse_pdb_atom_line(line, rec)) {
            continue;
        }
        // Protein heavy atoms only: ATOM records with standard AA names, non-H.
        if (rec.is_hetatm) {
            continue;
        }
        if (!is_standard_aa(rec.resname)) {
            continue;
        }
        const std::string el = element_from_pdb(rec.name, rec.element);
        const int Z = atomic_number(el);
        if (Z == 1 || to_upper(el) == "H" || to_upper(el) == "D" || to_upper(el) == "T") {
            continue;
        }
        // Also skip hydrogens identified by atom name (e.g. " H  ", "HA ", "1HG1")
        std::string an = trim_copy(rec.name);
        if (!an.empty()) {
            // Strip leading digits
            std::size_t i = 0;
            while (i < an.size() && std::isdigit(static_cast<unsigned char>(an[i]))) {
                ++i;
            }
            if (i < an.size() &&
                (an[i] == 'H' || an[i] == 'h' || an[i] == 'D' || an[i] == 'd')) {
                // Avoid false positive for HE* metals — protein AA names never use He.
                // His HE1/HE2 are hydrogens; good. "HG" etc. start with H.
                continue;
            }
        }
        protein.push_back(std::move(rec));
    }

    if (protein.empty()) {
        set_err(err, "load_pdb_protein_heavy: no protein heavy atoms in '" + path + "'");
        return false;
    }

    molecule_from_pdb_recs(protein, out, path);
    // Protein connectivity for clash checks may use inferred bonds; keep them.
    // Drop pure-H if any slipped through.
    std::vector<Atom> heavy;
    heavy.reserve(out.atoms.size());
    for (Atom& a : out.atoms) {
        if (!a.is_h) {
            heavy.push_back(std::move(a));
        }
    }
    out.atoms = std::move(heavy);
    for (std::size_t i = 0; i < out.atoms.size(); ++i) {
        out.atoms[i].id = static_cast<int>(i) + 1;
    }
    infer_bonds(out);
    return true;
}

}  // namespace flexaids::posebust
