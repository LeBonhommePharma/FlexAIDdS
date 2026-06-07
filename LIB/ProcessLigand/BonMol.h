// BonMol.h — Eigen-based canonical molecule representation for FlexAIDdS
// ProcessLigand Phase 3 standalone ligand preprocessing engine
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
//
// This header defines the internal canonical molecule form used throughout
// ProcessLigand. All algorithms (ring perception, aromaticity, rotatable bonds,
// SYBYL typing) operate on BonMol. Conversion to FlexAID .inp/.ga format is
// handled by FlexAIDWriter via BonMol::to_flexaid().

#pragma once

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <cstdint>
#include <array>
#include <optional>
#include <span>
#include <stdexcept>
#include <algorithm>
#include <numeric>

namespace bonmol {

// ---------------------------------------------------------------------------
// Elemental types
// ---------------------------------------------------------------------------

enum class Element : uint8_t {
    Unknown = 0,
    H   = 1,  He  = 2,  Li  = 3,  Be  = 4,  B   = 5,
    C   = 6,  N   = 7,  O   = 8,  F   = 9,  Ne  = 10,
    Na  = 11, Mg  = 12, Al  = 13, Si  = 14, P   = 15,
    S   = 16, Cl  = 17, Ar  = 18, K   = 19, Ca  = 20,
    Fe  = 26, Ni  = 28, Cu  = 29, Zn  = 30, Se  = 34,
    Br  = 35, I   = 53
};

/// Convert element symbol string → Element enum.
/// Case-insensitive for the first character (capitalised symbol standard).
inline Element element_from_symbol(std::string_view sym) noexcept {
    if (sym.empty()) return Element::Unknown;
    // Two-character symbols first
    if (sym.size() >= 2) {
        char a = static_cast<char>(std::toupper(static_cast<unsigned char>(sym[0])));
        char b = static_cast<char>(std::tolower(static_cast<unsigned char>(sym[1])));
        if (a=='H' && b=='e') return Element::He;
        if (a=='L' && b=='i') return Element::Li;
        if (a=='B' && b=='e') return Element::Be;
        if (a=='N' && b=='e') return Element::Ne;
        if (a=='N' && b=='a') return Element::Na;
        if (a=='M' && b=='g') return Element::Mg;
        if (a=='A' && b=='l') return Element::Al;
        if (a=='S' && b=='i') return Element::Si;
        if (a=='A' && b=='r') return Element::Ar;
        if (a=='C' && b=='a') return Element::Ca;
        if (a=='F' && b=='e') return Element::Fe;
        if (a=='N' && b=='i') return Element::Ni;
        if (a=='C' && b=='u') return Element::Cu;
        if (a=='Z' && b=='n') return Element::Zn;
        if (a=='S' && b=='e') return Element::Se;
        if (a=='B' && b=='r') return Element::Br;
        if (a=='C' && b=='l') return Element::Cl;
    }
    char a = static_cast<char>(std::toupper(static_cast<unsigned char>(sym[0])));
    switch (a) {
        case 'H': return Element::H;
        case 'B': return Element::B;
        case 'C': return Element::C;
        case 'N': return Element::N;
        case 'O': return Element::O;
        case 'F': return Element::F;
        case 'P': return Element::P;
        case 'S': return Element::S;
        case 'K': return Element::K;
        case 'I': return Element::I;
        default:  return Element::Unknown;
    }
}

/// Return the standard atomic weight (approximate, for MW computation).
inline float atomic_mass(Element e) noexcept {
    switch (e) {
        case Element::H:  return 1.008f;
        case Element::He: return 4.003f;
        case Element::Li: return 6.941f;
        case Element::Be: return 9.012f;
        case Element::B:  return 10.811f;
        case Element::C:  return 12.011f;
        case Element::N:  return 14.007f;
        case Element::O:  return 15.999f;
        case Element::F:  return 18.998f;
        case Element::Na: return 22.990f;
        case Element::Mg: return 24.305f;
        case Element::Al: return 26.982f;
        case Element::Si: return 28.086f;
        case Element::P:  return 30.974f;
        case Element::S:  return 32.065f;
        case Element::Cl: return 35.453f;
        case Element::K:  return 39.098f;
        case Element::Ca: return 40.078f;
        case Element::Fe: return 55.845f;
        case Element::Ni: return 58.693f;
        case Element::Cu: return 63.546f;
        case Element::Zn: return 65.38f;
        case Element::Se: return 78.96f;
        case Element::Br: return 79.904f;
        case Element::I:  return 126.904f;
        default:          return 0.0f;
    }
}

/// Return the standard covalent radius in Angstroms.
inline float covalent_radius(Element e) noexcept {
    switch (e) {
        case Element::H:  return 0.31f;
        case Element::B:  return 0.84f;
        case Element::C:  return 0.76f;
        case Element::N:  return 0.71f;
        case Element::O:  return 0.66f;
        case Element::F:  return 0.57f;
        case Element::Si: return 1.11f;
        case Element::P:  return 1.07f;
        case Element::S:  return 1.05f;
        case Element::Cl: return 1.02f;
        case Element::Se: return 1.20f;
        case Element::Br: return 1.20f;
        case Element::I:  return 1.39f;
        case Element::Fe: return 1.32f;
        case Element::Zn: return 1.22f;
        default:          return 1.50f;
    }
}

// ---------------------------------------------------------------------------
// Bond order
// ---------------------------------------------------------------------------

enum class BondOrder : uint8_t {
    SINGLE   = 1,
    DOUBLE   = 2,
    TRIPLE   = 3,
    AROMATIC = 4
};

/// Effective π-bond weight: aromatic bonds count as 1.5 in order sums.
inline float bond_order_value(BondOrder bo) noexcept {
    switch (bo) {
        case BondOrder::SINGLE:   return 1.0f;
        case BondOrder::DOUBLE:   return 2.0f;
        case BondOrder::TRIPLE:   return 3.0f;
        case BondOrder::AROMATIC: return 1.5f;
    }
    return 1.0f;
}

// ---------------------------------------------------------------------------
// Hybridisation
// ---------------------------------------------------------------------------

enum class Hybridization : uint8_t {
    UNSET = 0,
    SP    = 1,
    SP2   = 2,
    SP3   = 3
};

// ---------------------------------------------------------------------------
// Per-atom data
// ---------------------------------------------------------------------------

struct Atom {
    Element       element          = Element::Unknown;
    Hybridization hybrid           = Hybridization::UNSET;
    int           formal_charge    = 0;    // integer formal charge
    float         partial_charge   = 0.0f; // AM1-BCC or RESP charge
    bool          is_aromatic      = false;
    bool          is_hbond_donor   = false;
    bool          is_hbond_acceptor = false;
    int           implicit_h_count = 0;    // computed H count
    int           ring_membership  = 0;    // number of SSSR rings containing this atom
    int           sybyl_type       = 0;    // FlexAID 40-type SYBYL code (1-30)
    uint8_t       type_256         = 0;    // 256-type encoding (atom_typing_256.h)
    int           atom_map_num     = 0;    // SMILES atom-map number (for tracking)
    int           isotope          = 0;    // isotope label (from SMILES bracket)
    char          name[5]          = {};   // PDB atom name (populated by FlexAIDWriter)
};

// ---------------------------------------------------------------------------
// Per-bond data
// ---------------------------------------------------------------------------

struct Bond {
    int       atom_i      = -1;
    int       atom_j      = -1;
    BondOrder order       = BondOrder::SINGLE;
    bool      is_aromatic = false;
    bool      is_rotatable = false;
    bool      in_ring     = false;
};

// ---------------------------------------------------------------------------
// Ring descriptor
// ---------------------------------------------------------------------------

struct Ring {
    std::vector<int> atom_indices; // ordered around the ring
    bool is_aromatic = false;
    int  size        = 0;          // == atom_indices.size()
};

// ---------------------------------------------------------------------------
// BonMol — core molecule class
// ---------------------------------------------------------------------------

class BonMol {
public:
    // -----------------------------------------------------------------------
    // Coordinate matrix — 3 × N column-major (SoA for SIMD-friendly access)
    // -----------------------------------------------------------------------
    Eigen::Matrix3Xf coords; // uninitialized until atoms are added

    // -----------------------------------------------------------------------
    // Atom, bond, ring lists
    // -----------------------------------------------------------------------
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    std::vector<Ring> rings;

    // -----------------------------------------------------------------------
    // Adjacency — adjacency[i] = sorted list of bonded atom indices
    // -----------------------------------------------------------------------
    std::vector<std::vector<int>> adjacency;

    // -----------------------------------------------------------------------
    // Bond lookup — bond_index[i][j] = index into bonds[], or -1
    // -----------------------------------------------------------------------
    // Stored as a flat upper-triangle encoded in a map for sparse graphs.
    // For dense lookups use find_bond().
    std::vector<std::vector<int>> bond_adj; // bond_adj[i][k] = bond index for adjacency[i][k]

    // -----------------------------------------------------------------------
    // Metadata
    // -----------------------------------------------------------------------
    std::string name;
    std::string smiles;        // input SMILES string (if any)
    int         total_charge   = 0;
    float       molecular_weight = 0.0f;

    // -----------------------------------------------------------------------
    // Construction API
    // -----------------------------------------------------------------------

    /// Add an atom with 3D coordinates. Returns the new atom index.
    int add_atom(Element elem, float x, float y, float z) {
        int idx = static_cast<int>(atoms.size());
        Atom a;
        a.element = elem;
        atoms.push_back(a);

        // Grow coordinate matrix
        coords.conservativeResize(3, idx + 1);
        coords.col(idx) = Eigen::Vector3f(x, y, z);

        adjacency.emplace_back();
        bond_adj.emplace_back();
        return idx;
    }

    /// Add an atom with no 3D coordinates (SMILES path). Coords are NaN.
    int add_atom(Element elem, int formal_charge = 0,
                 int isotope = 0, bool is_aromatic = false) {
        int idx = static_cast<int>(atoms.size());
        Atom a;
        a.element       = elem;
        a.formal_charge = formal_charge;
        a.isotope       = isotope;
        a.is_aromatic   = is_aromatic;
        atoms.push_back(a);

        coords.conservativeResize(3, idx + 1);
        coords.col(idx).fill(std::numeric_limits<float>::quiet_NaN());

        adjacency.emplace_back();
        bond_adj.emplace_back();
        return idx;
    }

    /// Add a bond between atom i and j with given order.
    void add_bond(int i, int j, BondOrder order, bool aromatic = false) {
        if (i == j) return;
        if (i < 0 || j < 0 ||
            i >= static_cast<int>(atoms.size()) ||
            j >= static_cast<int>(atoms.size())) {
            throw std::out_of_range("BonMol::add_bond: atom index out of range");
        }
        // Prevent duplicate bonds
        if (find_bond(i, j) >= 0) return;

        int bidx = static_cast<int>(bonds.size());
        Bond b;
        b.atom_i      = i;
        b.atom_j      = j;
        b.order       = order;
        b.is_aromatic = aromatic;
        bonds.push_back(b);

        adjacency[i].push_back(j);
        adjacency[j].push_back(i);
        bond_adj[i].push_back(bidx);
        bond_adj[j].push_back(bidx);
    }

    /// Find bond index between i and j. Returns -1 if not found.
    int find_bond(int i, int j) const noexcept {
        if (i < 0 || j < 0 ||
            i >= static_cast<int>(adjacency.size()) ||
            j >= static_cast<int>(adjacency.size())) return -1;
        for (size_t k = 0; k < adjacency[i].size(); ++k) {
            if (adjacency[i][k] == j) return bond_adj[i][k];
        }
        return -1;
    }

    /// Called after all atoms/bonds are added: computes MW and charge sum.
    /// Aromaticity, ring perception and SYBYL types are computed by the
    /// dedicated modules (RingPerception, Aromaticity, SybylTyper) called
    /// from the ProcessLigand pipeline.
    void finalize() {
        molecular_weight = 0.0f;
        total_charge     = 0;
        for (const auto& a : atoms) {
            molecular_weight += atomic_mass(a.element);
            // Add implicit H contribution
            molecular_weight += a.implicit_h_count * atomic_mass(Element::H);
            total_charge     += a.formal_charge;
        }
    }

    // -----------------------------------------------------------------------
    // Query API
    // -----------------------------------------------------------------------

    int num_atoms()  const noexcept { return static_cast<int>(atoms.size()); }
    int num_bonds()  const noexcept { return static_cast<int>(bonds.size()); }

    int num_heavy_atoms() const noexcept {
        return static_cast<int>(
            std::count_if(atoms.begin(), atoms.end(),
                [](const Atom& a){ return a.element != Element::H; }));
    }

    int num_rotatable_bonds() const noexcept {
        return static_cast<int>(
            std::count_if(bonds.begin(), bonds.end(),
                [](const Bond& b){ return b.is_rotatable; }));
    }

    /// Degree (heavy + H neighbours)
    int degree(int idx) const noexcept {
        if (idx < 0 || idx >= static_cast<int>(adjacency.size())) return 0;
        return static_cast<int>(adjacency[idx].size());
    }

    /// Heavy-atom degree (implicit H not counted, explicit H counted)
    int heavy_degree(int idx) const noexcept {
        if (idx < 0 || idx >= static_cast<int>(adjacency.size())) return 0;
        int cnt = 0;
        for (int nb : adjacency[idx])
            if (atoms[nb].element != Element::H) ++cnt;
        return cnt;
    }

    std::vector<int> get_neighbors(int atom_idx) const {
        if (atom_idx < 0 || atom_idx >= static_cast<int>(adjacency.size()))
            return {};
        return adjacency[atom_idx];
    }

    /// Sum of bond orders around atom idx (aromatic bonds count as 1.5).
    float bond_order_sum(int idx) const noexcept {
        float s = 0.0f;
        for (int bidx : bond_adj[idx])
            s += bond_order_value(bonds[bidx].order);
        return s;
    }

    /// Whether atom idx has 3D coordinates set (not NaN).
    bool has_coords(int idx) const noexcept {
        if (idx < 0 || idx >= coords.cols()) return false;
        return !std::isnan(coords(0, idx));
    }

    /// Geometric centre of mass (heavy atoms only).
    Eigen::Vector3f centroid() const {
        Eigen::Vector3f c = Eigen::Vector3f::Zero();
        int cnt = 0;
        for (int i = 0; i < num_atoms(); ++i) {
            if (atoms[i].element != Element::H && has_coords(i)) {
                c += coords.col(i);
                ++cnt;
            }
        }
        if (cnt > 0) c /= static_cast<float>(cnt);
        return c;
    }

    // -----------------------------------------------------------------------
    // FlexAID output structure
    // -----------------------------------------------------------------------
    struct FlexAIDOutput {
        std::string inp_content;   // complete .inp file text
        std::string ga_content;    // complete .ga  file text
        int         num_atoms     = 0;
        int         num_dihedrals = 0;
    };

    /// Delegate to FlexAIDWriter — defined in FlexAIDWriter.cpp.
    FlexAIDOutput to_flexaid(const std::string& lig_name = "LIG") const;

    // -----------------------------------------------------------------------
    // Validation
    // -----------------------------------------------------------------------
    struct ValidationResult {
        bool        valid               = false;
        std::string error;
        bool        has_peptide_backbone = false; // dual-use guard
    };

    /// Lightweight structural validation (full checks in ValenceChecker).
    ValidationResult validate() const {
        ValidationResult r;
        if (atoms.empty()) {
            r.error = "empty molecule";
            return r;
        }
        if (num_heavy_atoms() == 0) {
            r.error = "no heavy atoms";
            return r;
        }
        if (num_heavy_atoms() > 500) {
            r.error = "molecule too large (> 500 heavy atoms); possible macrocycle or peptide";
            return r;
        }
        // Peptide-backbone check.
        //
        // A genuine peptide is a *chain* of residues linked by backbone amide
        // bonds: ...Cα–C(=O)–N–Cα–C(=O)–N... Counting isolated amide bonds is a
        // poor discriminant: peptidomimetic drugs (cyclic ureas, carbamates,
        // anilides, bis-amide protease inhibitors) routinely carry 3+ isolated
        // amides yet are legitimate small-molecule docking targets — e.g. Astex
        // 1KZK/JE2 and 1R55/097, and HAP2 HIV-protease inhibitors. The old
        // "amide_count >= 3" test rejected exactly those (verified).
        //
        // Instead count *backbone* amides: an amide C(=O)–N whose nitrogen also
        // bonds an α-carbon that itself carries another carbonyl carbon. >=3
        // such linked amides indicates a true poly-peptide backbone, while
        // isolated/peptidomimetic amides score low and pass.
        auto is_carbonyl_c = [&](int c) -> bool {
            if (c < 0 || atoms[c].element != Element::C) return false;
            for (int nb : adjacency[c]) {
                int bi = find_bond(c, nb);
                if (bi >= 0 && bonds[bi].order == BondOrder::DOUBLE &&
                    atoms[nb].element == Element::O)
                    return true;
            }
            return false;
        };
        int backbone_amides = 0;
        for (const auto& b : bonds) {
            if (b.order != BondOrder::SINGLE) continue;
            int c_idx = -1, n_idx = -1;
            if (atoms[b.atom_i].element == Element::C &&
                atoms[b.atom_j].element == Element::N) { c_idx = b.atom_i; n_idx = b.atom_j; }
            else if (atoms[b.atom_j].element == Element::C &&
                     atoms[b.atom_i].element == Element::N) { c_idx = b.atom_j; n_idx = b.atom_i; }
            if (c_idx < 0 || !is_carbonyl_c(c_idx)) continue;
            // Backbone test: amide N must bond an α-carbon that bears ANOTHER
            // carbonyl carbon (the next residue's C=O).
            for (int alpha : adjacency[n_idx]) {
                if (alpha == c_idx || atoms[alpha].element != Element::C) continue;
                bool alpha_to_next_carbonyl = false;
                for (int nb : adjacency[alpha]) {
                    if (nb != c_idx && is_carbonyl_c(nb)) { alpha_to_next_carbonyl = true; break; }
                }
                if (alpha_to_next_carbonyl) { ++backbone_amides; break; }
            }
        }
        if (backbone_amides >= 3) {
            r.has_peptide_backbone = true;
            r.error = "molecule appears to be a peptide (>= 3 backbone amide bonds)";
            return r;
        }
        // Macrocycle check: ring size > 12
        for (const auto& ring : rings) {
            if (ring.size > 12) {
                r.error = "macrocycle detected (ring size " +
                          std::to_string(ring.size) + "); not supported";
                return r;
            }
        }
        r.valid = true;
        return r;
    }
};

// ---------------------------------------------------------------------------
// Factory functions (implemented in their respective .cpp files)
// ---------------------------------------------------------------------------

/// Build BonMol from an OpenSMILES string (no 3D coords).
BonMol from_smiles(const std::string& smiles);

/// Build BonMol from an SDF V2000 file (with 3D coords).
BonMol from_sdf(const std::string& filepath);

/// Build BonMol from a SYBYL MOL2 file (with 3D coords and SYBYL types).
BonMol from_mol2(const std::string& filepath);

} // namespace bonmol
