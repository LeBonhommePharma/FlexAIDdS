// =============================================================================
// receptor_prep.cpp — Binding-site rotamer pre-relaxation for apo receptors
//
// See receptor_prep.h for API documentation and design rationale.
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
// =============================================================================

#include "receptor_prep.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace receptor_prep {

// =============================================================================
// Math primitives — pure C++ stdlib, no Eigen/BLAS
// =============================================================================

struct Vec3 {
    float x{0}, y{0}, z{0};
    Vec3() = default;
    Vec3(float ax, float ay, float az) : x(ax), y(ay), z(az) {}
};

static inline Vec3 v_sub(Vec3 a, Vec3 b) { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
static inline Vec3 v_add(Vec3 a, Vec3 b) { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
static inline Vec3 v_scale(Vec3 a, float s) { return {a.x*s, a.y*s, a.z*s}; }
static inline float v_dot(Vec3 a, Vec3 b)   { return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline Vec3 v_cross(Vec3 a, Vec3 b) {
    return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}
static inline float v_norm(Vec3 a) {
    return std::sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}
static inline Vec3 v_normalize(Vec3 a) {
    float n = v_norm(a) + 1e-12f;
    return v_scale(a, 1.0f / n);
}
static inline float v_dist(Vec3 a, Vec3 b) { return v_norm(v_sub(a, b)); }

/// Dihedral angle (degrees) for four 3D positions p1-p2-p3-p4.
static float dihedral_deg(Vec3 p1, Vec3 p2, Vec3 p3, Vec3 p4) {
    Vec3 b1 = v_sub(p2, p1);
    Vec3 b2 = v_sub(p3, p2);
    Vec3 b3 = v_sub(p4, p3);
    Vec3 b2n = v_normalize(b2);
    Vec3 n1 = v_cross(b1, b2);
    Vec3 n2 = v_cross(b2, b3);
    Vec3 m1 = v_cross(n1, b2n);
    return std::atan2(v_dot(m1, n2), v_dot(n1, n2)) * (180.0f / 3.14159265358979f);
}

/// Rotate vector v about unit-axis k by theta_deg using Rodrigues' formula.
static Vec3 rodrigues(Vec3 v, Vec3 k, float theta_deg) {
    float t = theta_deg * (3.14159265358979f / 180.0f);
    float c = std::cos(t), s = std::sin(t);
    Vec3 kn = v_normalize(k);
    // R·v = v·cos(t) + (k×v)·sin(t) + k·(k·v)·(1-cos(t))
    return v_add(v_add(v_scale(v, c),
                       v_scale(v_cross(kn, v), s)),
                 v_scale(kn, v_dot(kn, v) * (1.0f - c)));
}

// =============================================================================
// PDB I/O — minimal, self-contained, no external headers
// =============================================================================

struct PrepAtom {
    int    serial{0};
    char   name[5]{};      // stripped (e.g. "CA", "OG1")
    char   resname[4]{};   // stripped residue name
    char   chain{' '};
    int    resseq{0};
    char   icode{' '};
    Vec3   xyz;
    bool   is_hetatm{false};
    char   raw[82]{};      // verbatim PDB line for output reconstruction
};

/// (chain, resseq, icode, resname) — uniquely identifies a residue.
using ResKey = std::tuple<char, int, char, std::string>;
static inline ResKey res_key(const PrepAtom& a) {
    return {a.chain, a.resseq, a.icode, std::string(a.resname)};
}

static std::string strip_ws(const std::string& s) {
    size_t a = s.find_first_not_of(' ');
    if (a == std::string::npos) return {};
    size_t b = s.find_last_not_of(' ');
    return s.substr(a, b - a + 1);
}

/// Parse ATOM records from a PDB file into PrepAtom vector.
/// Also returns raw lines (for verbatim-passthrough output).
static std::vector<PrepAtom> parse_pdb(const std::string& path,
                                        std::vector<std::string>& raw_lines) {
    std::vector<PrepAtom> atoms;
    std::ifstream ifs(path);
    if (!ifs) {
        std::cerr << "  [PREP] ERROR: cannot open " << path << "\n";
        return atoms;
    }
    std::string line;
    while (std::getline(ifs, line)) {
        raw_lines.push_back(line);
        // Pad to 80 chars for safe column extraction
        while (line.size() < 80) line += ' ';

        bool is_atom   = (line.compare(0, 6, "ATOM  ") == 0);
        bool is_hetatm = (line.compare(0, 6, "HETATM") == 0);
        if (!is_atom && !is_hetatm) continue;

        PrepAtom a;
        a.is_hetatm = is_hetatm;

        // Skip alternate conformers (B, C, …)
        char altloc = line[16];
        if (altloc != ' ' && altloc != 'A' && altloc != '\0') continue;

        try { a.serial = std::stoi(line.substr(6, 5)); }
        catch (...) { a.serial = 0; }

        // Atom name: cols 13-16 (0-indexed 12-15), strip spaces
        std::string aname = strip_ws(line.substr(12, 4));
        // Skip hydrogens
        if (!aname.empty() && (aname[0] == 'H' ||
            (aname.size() >= 2 && aname[0] >= '0' && aname[0] <= '9' &&
             aname[1] == 'H'))) continue;
        std::strncpy(a.name, aname.c_str(), 4);

        std::string rname = strip_ws(line.substr(17, 3));
        std::strncpy(a.resname, rname.c_str(), 3);

        a.chain = line[21];

        try { a.resseq = std::stoi(line.substr(22, 4)); }
        catch (...) { a.resseq = 0; }
        a.icode = line[26];

        try {
            a.xyz.x = std::stof(line.substr(30, 8));
            a.xyz.y = std::stof(line.substr(38, 8));
            a.xyz.z = std::stof(line.substr(46, 8));
        } catch (...) { continue; }

        // Store verbatim raw line (original length, not padded)
        std::strncpy(a.raw, raw_lines.back().c_str(), 80);
        a.raw[80] = '\0';

        atoms.push_back(a);
    }
    return atoms;
}

/// Compute centroid of all heavy atoms in a PDB file (used for oracle site).
static Vec3 site_centroid(const std::string& path) {
    std::vector<std::string> dummy;
    auto atoms = parse_pdb(path, dummy);
    Vec3 c{};
    int n = 0;
    for (const auto& a : atoms) {
        c.x += a.xyz.x; c.y += a.xyz.y; c.z += a.xyz.z;
        ++n;
    }
    if (n > 0) { c.x /= n; c.y /= n; c.z /= n; }
    return c;
}

/// Write a PDB file with updated atom coordinates.
/// For each atom whose serial is in coord_map, replaces columns 31-54 with
/// the new XYZ; all other records are written verbatim.
static bool write_pdb(const std::string& path,
                      const std::vector<std::string>& raw_lines,
                      const std::vector<PrepAtom>& atoms,
                      const std::map<int, Vec3>& updated) {
    // Build serial → raw-line-index map for fast lookup
    std::map<int, size_t> serial_to_idx;
    for (size_t i = 0; i < atoms.size(); ++i)
        serial_to_idx[atoms[i].serial] = i;

    // Build a serial → raw_lines index: we need to find the raw line
    // for each modified atom.  The raw lines vector includes ALL lines
    // (REMARK, HEADER, etc.), so we find them by string-scanning once.
    std::map<int, size_t> serial_to_rawline;
    for (size_t li = 0; li < raw_lines.size(); ++li) {
        const auto& rl = raw_lines[li];
        if (rl.size() < 11) continue;
        if (rl.compare(0, 4, "ATOM") != 0 && rl.compare(0, 6, "HETATM") != 0) continue;
        try {
            int s = std::stoi(rl.substr(6, 5));
            serial_to_rawline[s] = li;
        } catch (...) {}
    }

    std::ofstream ofs(path);
    if (!ofs) {
        std::cerr << "  [PREP] ERROR: cannot write " << path << "\n";
        return false;
    }

    for (const auto& line : raw_lines) {
        if (line.size() >= 11 &&
            (line.compare(0, 4, "ATOM") == 0 || line.compare(0, 6, "HETATM") == 0))
        {
            try {
                int s = std::stoi(line.substr(6, 5));
                auto it = updated.find(s);
                if (it != updated.end()) {
                    // Rebuild the coordinate columns in a copy
                    std::string out = line;
                    while (out.size() < 80) out += ' ';
                    char buf[25];
                    std::snprintf(buf, sizeof(buf), "%8.3f%8.3f%8.3f",
                                  (double)it->second.x,
                                  (double)it->second.y,
                                  (double)it->second.z);
                    out.replace(30, 24, buf);
                    ofs << out << "\n";
                    continue;
                }
            } catch (...) {}
        }
        ofs << line << "\n";
    }
    return true;
}

// =============================================================================
// Dunbrack 2010 backbone-independent mode rotamers (top-5 per residue)
// Source: Dunbrack RL 2002 Curr Opin Struct Biol + 2010 penultimate library.
// chi values in degrees; NaN = chi slot unused for this residue type.
// =============================================================================

static constexpr float kNaN = std::numeric_limits<float>::quiet_NaN();

struct Rotamer {
    float chi[4];   // chi1..chi4 in degrees; NaN if not applicable
    float prob;     // approximate backbone-independent probability
};

struct ResRotamers {
    const char*  resname;
    int          n_rot;
    Rotamer      rots[5];
};

// clang-format off
static const ResRotamers DUNBRACK[] = {
    // ── 1 chi (chi2-4 = NaN) ───────────────────────────────────────────────
    { "SER", 3, {{{ -65.4f,  kNaN,  kNaN,  kNaN}, 0.432f},
                 {{  64.1f,  kNaN,  kNaN,  kNaN}, 0.349f},
                 {{ 180.0f,  kNaN,  kNaN,  kNaN}, 0.219f}} },

    { "CYS", 3, {{{ -65.2f,  kNaN,  kNaN,  kNaN}, 0.505f},
                 {{  63.5f,  kNaN,  kNaN,  kNaN}, 0.290f},
                 {{ 180.0f,  kNaN,  kNaN,  kNaN}, 0.205f}} },

    { "THR", 3, {{{  62.1f,  kNaN,  kNaN,  kNaN}, 0.468f},
                 {{ -66.9f,  kNaN,  kNaN,  kNaN}, 0.300f},
                 {{ 180.0f,  kNaN,  kNaN,  kNaN}, 0.232f}} },

    { "VAL", 3, {{{  63.8f,  kNaN,  kNaN,  kNaN}, 0.401f},
                 {{ 179.2f,  kNaN,  kNaN,  kNaN}, 0.362f},
                 {{ -63.4f,  kNaN,  kNaN,  kNaN}, 0.237f}} },

    // ── 2 chi (chi3-4 = NaN) ───────────────────────────────────────────────
    { "ASP", 5, {{{ -70.1f,  -15.1f,  kNaN,  kNaN}, 0.310f},
                 {{ -70.1f,  165.0f,  kNaN,  kNaN}, 0.202f},
                 {{ 170.3f,   15.2f,  kNaN,  kNaN}, 0.132f},
                 {{  64.8f,   15.2f,  kNaN,  kNaN}, 0.128f},
                 {{ 170.3f,  165.0f,  kNaN,  kNaN}, 0.097f}} },

    { "ASN", 5, {{{ -67.5f,  -10.0f,  kNaN,  kNaN}, 0.285f},
                 {{ 179.8f,  -10.0f,  kNaN,  kNaN}, 0.215f},
                 {{  63.2f,  -10.0f,  kNaN,  kNaN}, 0.190f},
                 {{ -67.5f,  178.2f,  kNaN,  kNaN}, 0.109f},
                 {{ 179.8f,  178.2f,  kNaN,  kNaN}, 0.083f}} },

    { "HIS", 5, {{{ -65.0f,  -75.0f,  kNaN,  kNaN}, 0.217f},
                 {{ -65.0f,  -10.0f,  kNaN,  kNaN}, 0.198f},
                 {{ 179.8f, -165.0f,  kNaN,  kNaN}, 0.108f},
                 {{ -65.0f,  -60.0f,  kNaN,  kNaN}, 0.095f},
                 {{ 179.8f,  -90.0f,  kNaN,  kNaN}, 0.073f}} },

    { "PHE", 5, {{{ -65.1f,   90.0f,  kNaN,  kNaN}, 0.253f},
                 {{ 179.0f,   90.0f,  kNaN,  kNaN}, 0.196f},
                 {{ -65.1f,  -90.0f,  kNaN,  kNaN}, 0.126f},
                 {{  63.9f,   90.0f,  kNaN,  kNaN}, 0.107f},
                 {{ -65.1f,    0.0f,  kNaN,  kNaN}, 0.094f}} },

    { "TYR", 5, {{{ -65.1f,   90.0f,  kNaN,  kNaN}, 0.253f},
                 {{ 179.0f,   90.0f,  kNaN,  kNaN}, 0.196f},
                 {{ -65.1f,  -90.0f,  kNaN,  kNaN}, 0.126f},
                 {{  63.9f,   90.0f,  kNaN,  kNaN}, 0.107f},
                 {{ 179.0f,  -90.0f,  kNaN,  kNaN}, 0.094f}} },

    { "TRP", 5, {{{ -65.0f,  -90.0f,  kNaN,  kNaN}, 0.220f},
                 {{ -65.0f,   90.0f,  kNaN,  kNaN}, 0.169f},
                 {{ 179.8f,   90.0f,  kNaN,  kNaN}, 0.122f},
                 {{ 179.8f,  -90.0f,  kNaN,  kNaN}, 0.095f},
                 {{  63.5f,   90.0f,  kNaN,  kNaN}, 0.071f}} },

    { "LEU", 5, {{{ -65.2f,  174.9f,  kNaN,  kNaN}, 0.327f},
                 {{ -65.2f,  -63.6f,  kNaN,  kNaN}, 0.241f},
                 {{ 179.8f,  174.9f,  kNaN,  kNaN}, 0.168f},
                 {{  63.8f,  -63.6f,  kNaN,  kNaN}, 0.102f},
                 {{ 179.8f,  -63.6f,  kNaN,  kNaN}, 0.088f}} },

    { "ILE", 5, {{{ -65.2f,  169.8f,  kNaN,  kNaN}, 0.337f},
                 {{ -65.2f,  -63.6f,  kNaN,  kNaN}, 0.248f},
                 {{ 179.8f,  169.8f,  kNaN,  kNaN}, 0.190f},
                 {{  63.8f,  169.8f,  kNaN,  kNaN}, 0.095f},
                 {{  63.8f,  -63.6f,  kNaN,  kNaN}, 0.060f}} },

    // ── 3 chi (chi4 = NaN) ─────────────────────────────────────────────────
    { "GLU", 5, {{{ -67.5f,  180.0f,    0.0f,  kNaN}, 0.207f},
                 {{ -67.5f, -180.0f,  -20.0f,  kNaN}, 0.166f},
                 {{ 179.4f,  180.0f,    0.0f,  kNaN}, 0.136f},
                 {{  63.8f,  180.0f,    0.0f,  kNaN}, 0.090f},
                 {{ -67.5f,  -65.0f,    0.0f,  kNaN}, 0.079f}} },

    { "GLN", 5, {{{ -67.5f,  180.0f,    0.0f,  kNaN}, 0.161f},
                 {{ -67.5f, -180.0f,  -30.0f,  kNaN}, 0.134f},
                 {{ 179.4f,  180.0f,    0.0f,  kNaN}, 0.118f},
                 {{  63.8f,  180.0f,    0.0f,  kNaN}, 0.086f},
                 {{ -67.5f,  -65.0f,    0.0f,  kNaN}, 0.074f}} },

    { "MET", 5, {{{ -67.5f,  180.0f,  -65.0f,  kNaN}, 0.166f},
                 {{ -67.5f,  180.0f,  180.0f,  kNaN}, 0.142f},
                 {{ 179.4f,  180.0f,  -65.0f,  kNaN}, 0.108f},
                 {{ -67.5f,  -65.0f,  180.0f,  kNaN}, 0.085f},
                 {{ 179.4f,  180.0f,  180.0f,  kNaN}, 0.077f}} },

    // ── 4 chi ──────────────────────────────────────────────────────────────
    { "LYS", 5, {{{ -67.5f,  180.0f,  180.0f,  180.0f}, 0.098f},
                 {{ -67.5f,  180.0f,  180.0f,   65.0f}, 0.082f},
                 {{ -67.5f,  180.0f,   65.0f,  180.0f}, 0.071f},
                 {{ 179.4f,  180.0f,  180.0f,  180.0f}, 0.063f},
                 {{ -67.5f,  180.0f,  180.0f,  -65.0f}, 0.054f}} },

    { "ARG", 5, {{{ -67.6f,  180.0f,  180.0f,  180.0f}, 0.122f},
                 {{ -67.6f,  180.0f,  -80.0f,  180.0f}, 0.089f},
                 {{ -67.6f,  180.0f,  180.0f,  -80.0f}, 0.076f},
                 {{ 179.4f,  180.0f,  180.0f,  180.0f}, 0.068f},
                 {{ -67.6f,  -65.0f,  180.0f,  180.0f}, 0.059f}} },
};
// clang-format on

static constexpr int N_RES_TYPES =
    static_cast<int>(sizeof(DUNBRACK) / sizeof(DUNBRACK[0]));

static const ResRotamers* find_rotamers(const char* resname) {
    for (int i = 0; i < N_RES_TYPES; ++i)
        if (std::strcmp(DUNBRACK[i].resname, resname) == 0)
            return &DUNBRACK[i];
    return nullptr;
}

// =============================================================================
// Chi-angle geometry definitions
//
// For chi_i:  dihedral atoms = (d1, d2, d3, d4);  rotation axis = d2→d3.
// When setting chi_i, all atoms in mobile[] rotate about d2-d3.
// mobile[] is a nullptr-terminated list of atom names (stripped).
// =============================================================================

struct ChiDef {
    const char* d1;            // dihedral atom 1 (for reading current angle)
    const char* d2;
    const char* d3;
    const char* d4;
    const char* mobile[15];    // nullptr-terminated: atoms that MOVE for this chi
                               // (distal to d2-d3 bond; d3 itself is fixed)
};

struct ResChi {
    const char* resname;
    int         n_chi;
    ChiDef      chis[4];
};

// clang-format off
static const ResChi CHI_TABLE[] = {
    { "SER", 1, {
        // chi1: N-CA-CB-OG; rotate OG about CA-CB
        { "N","CA","CB","OG",  {"OG", nullptr} },
    }},
    { "CYS", 1, {
        { "N","CA","CB","SG",  {"SG", nullptr} },
    }},
    { "THR", 1, {
        // chi1: N-CA-CB-OG1; CG2 also rides the CB frame
        { "N","CA","CB","OG1", {"OG1","CG2", nullptr} },
    }},
    { "VAL", 1, {
        { "N","CA","CB","CG1", {"CG1","CG2", nullptr} },
    }},
    { "ASP", 2, {
        { "N",  "CA","CB","CG",  {"CG","OD1","OD2", nullptr} },
        { "CA", "CB","CG","OD1", {"OD1","OD2",       nullptr} },
    }},
    { "ASN", 2, {
        { "N",  "CA","CB","CG",  {"CG","OD1","ND2",       nullptr} },
        { "CA", "CB","CG","OD1", {"OD1","ND2",             nullptr} },
    }},
    { "HIS", 2, {
        { "N",  "CA","CB","CG",  {"CG","ND1","CD2","CE1","NE2", nullptr} },
        { "CA", "CB","CG","ND1", {"ND1","CD2","CE1","NE2",      nullptr} },
    }},
    { "PHE", 2, {
        { "N",  "CA","CB","CG",  {"CG","CD1","CD2","CE1","CE2","CZ",  nullptr} },
        { "CA", "CB","CG","CD1", {"CD1","CD2","CE1","CE2","CZ",       nullptr} },
    }},
    { "TYR", 2, {
        { "N",  "CA","CB","CG",  {"CG","CD1","CD2","CE1","CE2","CZ","OH", nullptr} },
        { "CA", "CB","CG","CD1", {"CD1","CD2","CE1","CE2","CZ","OH",      nullptr} },
    }},
    { "TRP", 2, {
        { "N",  "CA","CB","CG",  {"CG","CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2", nullptr} },
        { "CA", "CB","CG","CD1", {"CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2",      nullptr} },
    }},
    { "LEU", 2, {
        { "N",  "CA","CB","CG",  {"CG","CD1","CD2", nullptr} },
        { "CA", "CB","CG","CD1", {"CD1","CD2",       nullptr} },
    }},
    { "ILE", 2, {
        { "N",  "CA","CB","CG1", {"CG1","CG2","CD1", nullptr} },
        { "CA", "CB","CG1","CD1",{"CD1",              nullptr} },
    }},
    { "GLU", 3, {
        { "N",  "CA","CB","CG",  {"CG","CD","OE1","OE2",  nullptr} },
        { "CA", "CB","CG","CD",  {"CD","OE1","OE2",        nullptr} },
        { "CB", "CG","CD","OE1", {"OE1","OE2",              nullptr} },
    }},
    { "GLN", 3, {
        { "N",  "CA","CB","CG",  {"CG","CD","OE1","NE2",  nullptr} },
        { "CA", "CB","CG","CD",  {"CD","OE1","NE2",        nullptr} },
        { "CB", "CG","CD","OE1", {"OE1","NE2",              nullptr} },
    }},
    { "MET", 3, {
        { "N",  "CA","CB","CG",  {"CG","SD","CE",  nullptr} },
        { "CA", "CB","CG","SD",  {"SD","CE",        nullptr} },
        { "CB", "CG","SD","CE",  {"CE",              nullptr} },
    }},
    { "LYS", 4, {
        { "N",  "CA","CB","CG",  {"CG","CD","CE","NZ",  nullptr} },
        { "CA", "CB","CG","CD",  {"CD","CE","NZ",        nullptr} },
        { "CB", "CG","CD","CE",  {"CE","NZ",              nullptr} },
        { "CG", "CD","CE","NZ",  {"NZ",                   nullptr} },
    }},
    { "ARG", 4, {
        { "N",  "CA","CB","CG",  {"CG","CD","NE","CZ","NH1","NH2",  nullptr} },
        { "CA", "CB","CG","CD",  {"CD","NE","CZ","NH1","NH2",        nullptr} },
        { "CB", "CG","CD","NE",  {"NE","CZ","NH1","NH2",             nullptr} },
        { "CG", "CD","NE","CZ",  {"CZ","NH1","NH2",                  nullptr} },
    }},
};
// clang-format on

static constexpr int N_CHI_RES = static_cast<int>(sizeof(CHI_TABLE) / sizeof(CHI_TABLE[0]));

static const ResChi* find_chi(const char* resname) {
    for (int i = 0; i < N_CHI_RES; ++i)
        if (std::strcmp(CHI_TABLE[i].resname, resname) == 0)
            return &CHI_TABLE[i];
    return nullptr;
}

// =============================================================================
// VDW radii for clash scoring
// =============================================================================

static float vdw_radius(const char* atom_name) {
    // Identify element from stripped atom name: first alpha character
    const char* p = atom_name;
    while (*p && (*p < 'A' || *p > 'Z') && (*p < 'a' || *p > 'z')) ++p;
    if (!*p) return 1.70f;
    char elem = static_cast<char>(std::toupper((unsigned char)*p));
    switch (elem) {
        case 'C': return 1.70f;
        case 'N': return 1.55f;
        case 'O': return 1.52f;
        case 'S': return 1.80f;
        case 'P': return 1.80f;
        case 'F': return 1.47f;
        default:  return 1.70f;
    }
}

// =============================================================================
// Rotamer manipulation
// =============================================================================

using AtomCoords = std::map<std::string, Vec3>;   // {atom_name → xyz}

/// Apply a single Dunbrack rotamer to a coordinate map.
/// Returns false if required atoms are missing.
static bool apply_rotamer(AtomCoords& coords,
                           const ResChi& chi_res,
                           const Rotamer& rot) {
    int n_chi = chi_res.n_chi;
    for (int ci = 0; ci < n_chi; ++ci) {
        if (std::isnan(rot.chi[ci])) break;

        const ChiDef& cd = chi_res.chis[ci];
        float target = rot.chi[ci];

        // Check all required atoms are present
        for (const char* n : {cd.d1, cd.d2, cd.d3, cd.d4}) {
            if (coords.find(n) == coords.end()) return false;
        }

        // Compute current dihedral
        float current = dihedral_deg(coords.at(cd.d1), coords.at(cd.d2),
                                     coords.at(cd.d3), coords.at(cd.d4));
        float delta = target - current;
        // Normalise to (−180, 180]
        delta = std::fmod(delta + 180.0f, 360.0f) - 180.0f;
        if (std::fabs(delta) < 0.01f) continue;

        // Rotation axis = d2 → d3
        Vec3 origin = coords.at(cd.d2);
        Vec3 axis   = v_sub(coords.at(cd.d3), origin);

        // Rotate all mobile atoms about this axis
        for (int mi = 0; cd.mobile[mi] != nullptr; ++mi) {
            auto it = coords.find(cd.mobile[mi]);
            if (it == coords.end()) continue;
            Vec3 v = v_sub(it->second, origin);
            it->second = v_add(origin, rodrigues(v, axis, delta));
        }
    }
    return true;
}

// =============================================================================
// Clash scoring
// =============================================================================

struct EnvAtom {
    Vec3  xyz;
    float radius;
};

/// Summed VDW overlap: Σ max(0, rA + rB − tol − d) for mobile vs env.
static float clash_score(const AtomCoords& mobile_coords,
                          const std::vector<EnvAtom>& env,
                          float tol) {
    float score = 0.0f;
    for (const auto& [name, xyz] : mobile_coords) {
        float rm = vdw_radius(name.c_str());
        for (const auto& e : env) {
            float d = v_dist(xyz, e.xyz);
            float overlap = rm + e.radius - tol - d;
            if (overlap > 0.0f) score += overlap;
        }
    }
    return score;
}

/// Sidechain-only coords (atoms in any of the chi mobile sets for this residue).
static AtomCoords sidechain_coords(const AtomCoords& all_coords,
                                    const ResChi& chi_res) {
    std::set<std::string> sc_names;
    for (int ci = 0; ci < chi_res.n_chi; ++ci)
        for (int mi = 0; chi_res.chis[ci].mobile[mi] != nullptr; ++mi)
            sc_names.insert(chi_res.chis[ci].mobile[mi]);
    AtomCoords sc;
    for (const auto& nm : sc_names) {
        auto it = all_coords.find(nm);
        if (it != all_coords.end()) sc[nm] = it->second;
    }
    return sc;
}

// =============================================================================
// Main entry point
// =============================================================================

int prep_receptor_rotamers(const std::string& receptor_pdb,
                            const std::string& oracle_site_pdb,
                            const std::string& out_pdb,
                            float radius_ang,
                            int   top_n,
                            float vdw_tol)
{
    // ── 1. Load receptor ──────────────────────────────────────────────────
    std::vector<std::string> raw_lines;
    std::vector<PrepAtom> atoms = parse_pdb(receptor_pdb, raw_lines);
    if (atoms.empty()) {
        std::cerr << "  [PREP] ERROR: no atoms in " << receptor_pdb << "\n";
        return -1;
    }

    // ── 2. Oracle site centroid ───────────────────────────────────────────
    Vec3 centroid = site_centroid(oracle_site_pdb);
    std::cerr << "  [PREP] centroid = ("
              << centroid.x << ", " << centroid.y << ", " << centroid.z << ")\n";

    // ── 3. Identify pocket residues ───────────────────────────────────────
    // Build per-residue Cα coordinate
    std::map<ResKey, Vec3> ca_xyz;
    for (const auto& a : atoms) {
        if (!a.is_hetatm && std::strcmp(a.name, "CA") == 0)
            ca_xyz[res_key(a)] = a.xyz;
    }

    // Skip residues with no rotamer data (GLY, ALA, PRO, non-standard)
    static const std::set<std::string> NO_CHI_RES = {"GLY","ALA","PRO"};

    struct PocketRes { ResKey rk; float ca_dist; };
    std::vector<PocketRes> pocket;
    for (const auto& [rk, ca] : ca_xyz) {
        const std::string& resname = std::get<3>(rk);
        if (NO_CHI_RES.count(resname)) continue;
        if (!find_rotamers(resname.c_str())) continue;
        if (!find_chi(resname.c_str())) continue;
        float d = v_dist(ca, centroid);
        if (d <= radius_ang)
            pocket.push_back({rk, d});
    }

    std::cerr << "  [PREP] " << pocket.size()
              << " pocket residues within " << radius_ang << " Å\n";

    if (pocket.empty()) {
        // Nothing to do — write a copy and return
        std::ofstream ofs(out_pdb);
        for (const auto& l : raw_lines) ofs << l << "\n";
        return 0;
    }

    // ── 4. Build environment atom list (all residues) ─────────────────────
    // live_coords: serial → current xyz (updated as rotamers are accepted)
    std::map<int, Vec3> live_xyz;
    for (const auto& a : atoms) live_xyz[a.serial] = a.xyz;

    // Helper: build flat env list for a given residue (excluding that residue)
    auto build_env_for_res = [&](const ResKey& rk) -> std::vector<EnvAtom> {
        std::vector<EnvAtom> env;
        for (const auto& a : atoms) {
            if (res_key(a) == rk) continue;
            env.push_back({live_xyz[a.serial], vdw_radius(a.name)});
        }
        return env;
    };

    // Helper: build AtomCoords for a residue from live_xyz
    auto live_res_coords = [&](const ResKey& rk) -> AtomCoords {
        AtomCoords m;
        for (const auto& a : atoms)
            if (res_key(a) == rk)
                m[a.name] = live_xyz[a.serial];
        return m;
    };

    // ── 5. Compute initial clash scores, sort worst-first ─────────────────
    for (auto& pr : pocket) {
        const ResChi* chi_res = find_chi(std::get<3>(pr.rk).c_str());
        AtomCoords res_coords = live_res_coords(pr.rk);
        AtomCoords sc = sidechain_coords(res_coords, *chi_res);
        auto env = build_env_for_res(pr.rk);
        pr.ca_dist = clash_score(sc, env, vdw_tol);  // re-use ca_dist field as score
    }
    std::sort(pocket.begin(), pocket.end(),
              [](const PocketRes& a, const PocketRes& b) {
                  return a.ca_dist > b.ca_dist;  // worst (highest clash) first
              });

    // ── 6. Greedy single-residue rotamer optimisation ─────────────────────
    int n_modified = 0;
    top_n = std::min(top_n, 5);  // table has max 5 per residue

    for (auto& pr : pocket) {
        const ResKey& rk = pr.rk;
        const std::string& resname = std::get<3>(rk);
        char chain   = std::get<0>(rk);
        int  resseq  = std::get<1>(rk);

        const ResRotamers* rr  = find_rotamers(resname.c_str());
        const ResChi*      rc  = find_chi(resname.c_str());
        if (!rr || !rc) continue;

        auto env = build_env_for_res(rk);
        AtomCoords res_coords = live_res_coords(rk);
        AtomCoords sc_current = sidechain_coords(res_coords, *rc);
        float best_score = clash_score(sc_current, env, vdw_tol);
        float init_score = best_score;

        AtomCoords best_coords = res_coords;  // current as fallback

        int n_try = std::min(top_n, rr->n_rot);
        for (int ri = 0; ri < n_try; ++ri) {
            // Work on a copy of the residue coordinates
            AtomCoords trial = res_coords;
            if (!apply_rotamer(trial, *rc, rr->rots[ri])) continue;

            AtomCoords sc_trial = sidechain_coords(trial, *rc);
            float score = clash_score(sc_trial, env, vdw_tol);
            if (score < best_score - 1e-4f) {
                best_score = score;
                best_coords = trial;
            }
        }

        if (best_score < init_score - 1e-4f) {
            // Accept: update live_xyz for all atoms of this residue
            for (auto& a : atoms) {
                if (res_key(a) != rk) continue;
                auto it = best_coords.find(a.name);
                if (it != best_coords.end())
                    live_xyz[a.serial] = it->second;
            }
            ++n_modified;
            std::cerr << "  [PREP]   " << resname
                      << " " << chain << resseq
                      << "  clash " << init_score
                      << " -> " << best_score
                      << "  (Δ=" << (best_score - init_score) << ")\n";
        } else {
            std::cerr << "  [PREP]   " << resname
                      << " " << chain << resseq
                      << "  clash " << init_score
                      << " -> kept (no improvement)\n";
        }
    }

    std::cerr << "  [PREP] " << n_modified << "/" << pocket.size()
              << " pocket residues rotamer-optimised\n";

    // ── 7. Write output PDB ───────────────────────────────────────────────
    if (!write_pdb(out_pdb, raw_lines, atoms, live_xyz)) return -1;
    return n_modified;
}

} // namespace receptor_prep
