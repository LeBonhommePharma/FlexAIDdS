// FlexAIDWriter.cpp — FlexAID .inp and .ga file generation
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0

#include "FlexAIDWriter.h"
#include "SybylTyper.h"

#define _USE_MATH_DEFINES // MSVC compatibility for M_PI
#include <cmath>
#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif
#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <deque>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

namespace bonmol {
namespace writer {

// ---------------------------------------------------------------------------
// Geometry helpers
// ---------------------------------------------------------------------------

float FlexAIDWriter::dihedral_angle(const Eigen::Vector3f& a,
                                     const Eigen::Vector3f& b,
                                     const Eigen::Vector3f& c,
                                     const Eigen::Vector3f& d) const {
    Eigen::Vector3f b1 = b - a;
    Eigen::Vector3f b2 = c - b;
    Eigen::Vector3f b3 = d - c;

    Eigen::Vector3f n1 = b1.cross(b2);
    Eigen::Vector3f n2 = b2.cross(b3);

    float n1_norm = n1.norm();
    float n2_norm = n2.norm();
    if (n1_norm < 1e-6f || n2_norm < 1e-6f) return 0.0f;

    n1 /= n1_norm;
    n2 /= n2_norm;

    float cos_angle = n1.dot(n2);
    cos_angle = std::max(-1.0f, std::min(1.0f, cos_angle));
    float angle = std::acos(cos_angle) * 180.0f / static_cast<float>(M_PI);

    // Determine sign
    if (n1.cross(n2).dot(b2) < 0.0f) angle = -angle;
    return angle;
}

float FlexAIDWriter::bond_angle(const Eigen::Vector3f& a,
                                 const Eigen::Vector3f& b,
                                 const Eigen::Vector3f& c) const {
    Eigen::Vector3f v1 = (a - b).normalized();
    Eigen::Vector3f v2 = (c - b).normalized();
    float cos_a = v1.dot(v2);
    cos_a = std::max(-1.0f, std::min(1.0f, cos_a));
    return std::acos(cos_a) * 180.0f / static_cast<float>(M_PI);
}

// ---------------------------------------------------------------------------
// PDB ATOM record formatting
// ---------------------------------------------------------------------------

std::string FlexAIDWriter::format_pdb_atom(int serial, const std::string& atom_name,
                                            const std::string& res_name, int res_num,
                                            float x, float y, float z,
                                            float occ, float bfac) const {
    // PDB ATOM record: fixed-width columns
    // ATOM  SSSSS AAAA LRRR CNNNNI   XXXXXXXXYYYYYYYYZZZZZZZZ OOOOOO BBBBBB
    std::ostringstream ss;
    ss << std::left << std::setw(6) << "ATOM  "
       << std::right << std::setw(5) << serial
       << " "
       << std::left << std::setw(4) << atom_name.substr(0, 4)
       << " "
       << std::left << std::setw(3) << res_name.substr(0, 3)
       << " "
       << " "                                     // chain ID (blank for ligand)
       << std::right << std::setw(4) << res_num
       << "    "                                  // insertion code + padding
       << std::fixed << std::setprecision(3)
       << std::right << std::setw(8) << x
       << std::right << std::setw(8) << y
       << std::right << std::setw(8) << z
       << std::right << std::setw(6) << std::fixed << std::setprecision(2) << occ
       << std::right << std::setw(6) << std::fixed << std::setprecision(2) << bfac
       << "\n";
    return ss.str();
}

std::string FlexAIDWriter::format_atom_name(const BonMol& mol, int atom_idx,
                                             const std::string& prefix) const {
    // Generate unique atom name: element + serial, e.g. C1, N2, O3
    // If atom has a pre-assigned name, use it
    const Atom& a = mol.atoms[atom_idx];
    if (a.name[0] != '\0') {
        return std::string(a.name);
    }
    // Build from element
    std::string elem_str;
    switch (a.element) {
        case Element::H:  elem_str = "H";  break;
        case Element::C:  elem_str = "C";  break;
        case Element::N:  elem_str = "N";  break;
        case Element::O:  elem_str = "O";  break;
        case Element::F:  elem_str = "F";  break;
        case Element::P:  elem_str = "P";  break;
        case Element::S:  elem_str = "S";  break;
        case Element::Cl: elem_str = "CL"; break;
        case Element::Br: elem_str = "BR"; break;
        case Element::I:  elem_str = "I";  break;
        case Element::Fe: elem_str = "FE"; break;
        case Element::Zn: elem_str = "ZN"; break;
        default:          elem_str = "X";  break;
    }
    return elem_str + std::to_string(atom_idx + 1);
}

// ---------------------------------------------------------------------------
// Spanning tree (BFS) for Z-matrix ordering
// ---------------------------------------------------------------------------

FlexAIDWriter::SpanningTree FlexAIDWriter::build_spanning_tree(const BonMol& mol) const {
    const int N = mol.num_atoms();
    SpanningTree tree;
    tree.order.reserve(N);
    tree.parent.assign(N, -1);
    tree.parent2.assign(N, -1);
    tree.parent3.assign(N, -1);

    std::vector<bool> visited(N, false);

    // Start from atom 0 (or first heavy atom)
    int root = 0;
    for (int i = 0; i < N; ++i) {
        if (mol.atoms[i].element != Element::H) { root = i; break; }
    }

    std::deque<int> queue;
    queue.push_back(root);
    visited[root] = true;

    while (!queue.empty()) {
        int v = queue.front(); queue.pop_front();
        tree.order.push_back(v);

        // Set grandparent and great-grandparent
        int p1 = tree.parent[v];
        int p2 = (p1 >= 0) ? tree.parent[p1] : -1;
        int p3 = (p2 >= 0) ? tree.parent[p2] : -1;
        tree.parent2[v] = p2;
        tree.parent3[v] = p3;

        // Sort neighbours: heavy atoms first, then H
        std::vector<int> nbrs = mol.adjacency[v];
        std::sort(nbrs.begin(), nbrs.end(), [&](int a, int b){
            bool ah = mol.atoms[a].element == Element::H;
            bool bh = mol.atoms[b].element == Element::H;
            if (ah != bh) return !ah; // heavy first
            return a < b;
        });

        for (int nb : nbrs) {
            if (!visited[nb]) {
                visited[nb] = true;
                tree.parent[nb] = v;
                queue.push_back(nb);
            }
        }
    }

    // Handle disconnected components (shouldn't happen for valid ligands)
    for (int i = 0; i < N; ++i) {
        if (!visited[i]) {
            tree.order.push_back(i);
            tree.parent[i] = -1;
        }
    }

    return tree;
}

// ---------------------------------------------------------------------------
// Internal coordinate computation
// ---------------------------------------------------------------------------

std::vector<InternalCoord> FlexAIDWriter::compute_internal_coords(
        const BonMol& mol, const SpanningTree& tree) const {
    const int N = mol.num_atoms();
    std::vector<InternalCoord> ics(N);

    for (int k = 0; k < N; ++k) {
        int v = tree.order[k];
        InternalCoord& ic = ics[k];
        ic.atom_idx = v;
        ic.ref1 = tree.parent[v];
        ic.ref2 = tree.parent2[v];
        ic.ref3 = tree.parent3[v];

        ic.bond_length = 0.0f;
        ic.bond_angle  = 0.0f;
        ic.dihedral    = 0.0f;

        bool has_3d = mol.has_coords(v);

        if (ic.ref1 >= 0 && has_3d && mol.has_coords(ic.ref1)) {
            Eigen::Vector3f cv  = mol.coords.col(v);
            Eigen::Vector3f cr1 = mol.coords.col(ic.ref1);
            ic.bond_length = (cv - cr1).norm();

            if (ic.ref2 >= 0 && mol.has_coords(ic.ref2)) {
                Eigen::Vector3f cr2 = mol.coords.col(ic.ref2);
                ic.bond_angle = bond_angle(cv, cr1, cr2);

                if (ic.ref3 >= 0 && mol.has_coords(ic.ref3)) {
                    Eigen::Vector3f cr3 = mol.coords.col(ic.ref3);
                    ic.dihedral = dihedral_angle(cv, cr1, cr2, cr3);
                }
            }
        }
    }

    return ics;
}

// ---------------------------------------------------------------------------
// Geometry invariants (bond / angle / ring-closure) for production IC path
// ---------------------------------------------------------------------------

std::string FlexAIDWriter::check_geometry_invariants(
        const BonMol& mol,
        const std::vector<InternalCoord>& ics) const {
    // Covalent bond length bounds (Å) — generous to absorb crystal noise but
    // tight enough to catch IC corruption / broken reconstruction.
    constexpr float kBondMin = 0.70f;
    constexpr float kBondMax = 2.40f;
    // Valence angles (degrees)
    constexpr float kAngleMin = 50.0f;
    constexpr float kAngleMax = 170.0f;
    // Ring-closure: bonded ring atoms must still land near covalent distance.
    constexpr float kRingClosureMax = 2.60f;

    // 1) Bond lengths from topology + coordinates
    for (const auto& bond : mol.bonds) {
        const int i = bond.atom_i;
        const int j = bond.atom_j;
        if (!mol.has_coords(i) || !mol.has_coords(j)) continue;
        const float d = (mol.coords.col(i) - mol.coords.col(j)).norm();
        if (!std::isfinite(d) || d < kBondMin || d > kBondMax) {
            std::ostringstream oss;
            oss << "bond-length invariant failed atoms " << (i + 1) << "-"
                << (j + 1) << " d=" << d << " A (allowed [" << kBondMin
                << "," << kBondMax << "])";
            return oss.str();
        }
        if (bond.in_ring && d > kRingClosureMax) {
            std::ostringstream oss;
            oss << "ring-closure invariant failed atoms " << (i + 1) << "-"
                << (j + 1) << " d=" << d << " A";
            return oss.str();
        }
    }

    // 2) IC bond lengths / angles produced for reconstruction
    for (const auto& ic : ics) {
        if (ic.ref1 >= 0) {
            if (!std::isfinite(ic.bond_length) ||
                ic.bond_length < kBondMin || ic.bond_length > kBondMax) {
                // Root / partially defined atoms may legitimately have 0 length.
                if (ic.bond_length > 1e-6f) {
                    std::ostringstream oss;
                    oss << "IC bond-length invariant failed atom "
                        << (ic.atom_idx + 1) << " bl=" << ic.bond_length;
                    return oss.str();
                }
            }
        }
        if (ic.ref1 >= 0 && ic.ref2 >= 0 && ic.bond_angle > 1e-3f) {
            if (!std::isfinite(ic.bond_angle) ||
                ic.bond_angle < kAngleMin || ic.bond_angle > kAngleMax) {
                std::ostringstream oss;
                oss << "IC bond-angle invariant failed atom "
                    << (ic.atom_idx + 1) << " ang=" << ic.bond_angle;
                return oss.str();
            }
        }
        if (ic.ref1 >= 0 && ic.ref2 >= 0 && ic.ref3 >= 0) {
            if (!std::isfinite(ic.dihedral)) {
                std::ostringstream oss;
                oss << "IC dihedral non-finite atom " << (ic.atom_idx + 1);
                return oss.str();
            }
        }
    }

    // 3) Valence angles from Cartesian triples of bonded heavy atoms
    for (int c = 0; c < mol.num_atoms(); ++c) {
        if (mol.atoms[c].element == Element::H) continue;
        if (!mol.has_coords(c)) continue;
        const auto& nbrs = mol.adjacency[c];
        for (size_t a = 0; a < nbrs.size(); ++a) {
            for (size_t b = a + 1; b < nbrs.size(); ++b) {
                const int i = nbrs[a];
                const int j = nbrs[b];
                if (mol.atoms[i].element == Element::H ||
                    mol.atoms[j].element == Element::H) continue;
                if (!mol.has_coords(i) || !mol.has_coords(j)) continue;
                const float ang = bond_angle(mol.coords.col(i),
                                             mol.coords.col(c),
                                             mol.coords.col(j));
                if (!std::isfinite(ang) || ang < kAngleMin || ang > kAngleMax) {
                    // Linear sp centres can legitimately approach 180°; allow up to 180.
                    if (ang > 170.0f && ang <= 180.5f) continue;
                    std::ostringstream oss;
                    oss << "valence-angle invariant failed centre atom "
                        << (c + 1) << " ang=" << ang;
                    return oss.str();
                }
            }
        }
    }

    return {};
}

// ---------------------------------------------------------------------------
// Dihedral gene construction from rotatable bonds
// ---------------------------------------------------------------------------

std::vector<DihedralGene> FlexAIDWriter::build_dihedral_genes(const BonMol& mol) const {
    std::vector<DihedralGene> genes;

    for (int bidx = 0; bidx < mol.num_bonds(); ++bidx) {
        const Bond& bond = mol.bonds[bidx];
        if (!bond.is_rotatable) continue;

        int i = bond.atom_i;
        int j = bond.atom_j;

        // Find ref atoms: one neighbour of i (not j) and one neighbour of j (not i)
        int ref_a = -1, ref_b = -1;
        for (int nb : mol.adjacency[i]) {
            if (nb != j) { ref_a = nb; break; }
        }
        for (int nb : mol.adjacency[j]) {
            if (nb != i) { ref_b = nb; break; }
        }

        // Prefer heavy atoms for reference
        if (ref_a < 0 || ref_b < 0) continue; // terminal bond, skip

        DihedralGene gene;
        gene.bond_i   = i;
        gene.bond_j   = j;
        gene.ref_a    = ref_a;
        gene.ref_b    = ref_b;
        gene.min_deg  = -180.0f;
        gene.max_deg  =  180.0f;
        gene.step_deg =   1.0f;

        genes.push_back(gene);
    }

    return genes;
}

// ---------------------------------------------------------------------------
// Global positioning atom selection
// ---------------------------------------------------------------------------

std::array<int, 3> FlexAIDWriter::select_gpa(const BonMol& mol) const {
    std::array<int, 3> gpa = {0, 1, 2};
    const int N = mol.num_atoms();
    if (N < 3) return gpa;

    // Collect heavy atom indices with valid coordinates
    std::vector<int> heavy;
    heavy.reserve(N);
    for (int i = 0; i < N; ++i) {
        if (mol.atoms[i].element != Element::H && mol.has_coords(i))
            heavy.push_back(i);
    }
    if (heavy.size() < 3) {
        // Fall back to first 3 atoms
        for (int k = 0; k < std::min(3, N); ++k) gpa[k] = k;
        return gpa;
    }

    // Find the pair with maximum distance
    float max_dist = -1.0f;
    int ga = heavy[0], gb = heavy[1];
    for (size_t ii = 0; ii < heavy.size(); ++ii) {
        for (size_t jj = ii + 1; jj < heavy.size(); ++jj) {
            float d = (mol.coords.col(heavy[ii]) - mol.coords.col(heavy[jj])).norm();
            if (d > max_dist) { max_dist = d; ga = heavy[ii]; gb = heavy[jj]; }
        }
    }

    // Find atom maximally distant from midpoint of ga-gb (maximally non-collinear)
    Eigen::Vector3f mid = 0.5f * (mol.coords.col(ga) + mol.coords.col(gb));
    float max_d2 = -1.0f;
    int gc = -1;
    for (int hi : heavy) {
        if (hi == ga || hi == gb) continue;
        float d = (mol.coords.col(hi) - mid).norm();
        if (d > max_d2) { max_d2 = d; gc = hi; }
    }
    if (gc < 0) gc = heavy[2];

    gpa = {ga, gb, gc};
    return gpa;
}

// ---------------------------------------------------------------------------
// .inp file generator
// ---------------------------------------------------------------------------

std::string FlexAIDWriter::generate_inp(const BonMol& mol,
                                         const std::vector<InternalCoord>& ics,
                                         const std::array<int, 3>& gpa,
                                         const std::string& lig_name) const {
    std::ostringstream out;
    const int N = mol.num_atoms();

    // Header
    out << "REMARK  FlexAIDdS ProcessLigand v3.0  |  Apache-2.0  |  Le Bonhomme Pharma\n";
    out << "REMARK  Ligand: " << lig_name << "  Atoms: " << N
        << "  MW: " << std::fixed << std::setprecision(3) << mol.molecular_weight << "\n";
    out << "REMARK  GPA: " << gpa[0]+1 << " " << gpa[1]+1 << " " << gpa[2]+1 << "\n";

    // PDB ATOM records
    for (int k = 0; k < N; ++k) {
        int v = (k < static_cast<int>(ics.size())) ? ics[k].atom_idx : k;
        std::string aname = format_atom_name(mol, v, "");

        float x = 0.0f, y = 0.0f, z = 0.0f;
        if (mol.has_coords(v)) {
            x = mol.coords(0, v);
            y = mol.coords(1, v);
            z = mol.coords(2, v);
        }
        out << format_pdb_atom(v + 1, aname, lig_name, 1, x, y, z, 1.00f, 0.00f);
    }
    out << "END\n";

    // Internal coordinates section (FlexAID-specific extension)
    out << "REMARK  INTERNAL_COORDS\n";
    for (const auto& ic : ics) {
        out << "ICOORD "
            << std::setw(5) << ic.atom_idx + 1 << " "
            << std::setw(5) << (ic.ref1 >= 0 ? ic.ref1 + 1 : 0) << " "
            << std::setw(5) << (ic.ref2 >= 0 ? ic.ref2 + 1 : 0) << " "
            << std::setw(5) << (ic.ref3 >= 0 ? ic.ref3 + 1 : 0) << " "
            << std::fixed << std::setprecision(6)
            << std::setw(10) << ic.bond_length << " "
            << std::setw(10) << ic.bond_angle  << " "
            << std::setw(10) << ic.dihedral    << "\n";
    }

    // Atom types section
    out << "REMARK  ATOM_TYPES\n";
    for (int i = 0; i < N; ++i) {
        const Atom& a = mol.atoms[i];
        out << "ATYPE "
            << std::setw(5) << i + 1 << " "
            << std::setw(3) << a.sybyl_type << " "
            << std::setw(3) << static_cast<int>(a.type_256) << " "
            << std::fixed << std::setprecision(4)
            << std::setw(8) << a.partial_charge << " "
            << std::setw(8) << 0.0f              << "\n"; // resp_charge placeholder
    }

    // Connectivity section (mirrors FlexAID BOND records)
    out << "REMARK  CONNECTIVITY\n";
    for (const auto& bond : mol.bonds) {
        out << "BOND  "
            << std::setw(5) << bond.atom_i + 1 << " "
            << std::setw(5) << bond.atom_j + 1 << " "
            << std::setw(1) << static_cast<int>(bond.order) << "\n";
    }

    return out.str();
}

// ---------------------------------------------------------------------------
// .ga file generator
// ---------------------------------------------------------------------------

std::string FlexAIDWriter::generate_ga(const BonMol& mol,
                                        const std::vector<DihedralGene>& genes,
                                        const std::array<int, 3>& gpa) const {
    std::ostringstream out;

    out << "# FlexAIDdS GA gene descriptor\n";
    out << "# Generated by ProcessLigand v3.0 | Le Bonhomme Pharma | Apache-2.0\n";
    out << "#\n";

    // Global positioning atoms (3 atoms that anchor the pose)
    out << "GPA " << gpa[0]+1 << " " << gpa[1]+1 << " " << gpa[2]+1 << "\n";

    // Translation and rotation genes (6 fixed genes: tx,ty,tz,rx,ry,rz)
    out << "NUM_GENES " << (6 + static_cast<int>(genes.size())) << "\n";
    out << "#\n";
    out << "# Genes 1-3: translation (Angstroms)\n";
    out << "GENE 1 -20.0 20.0 0.1\n";
    out << "GENE 2 -20.0 20.0 0.1\n";
    out << "GENE 3 -20.0 20.0 0.1\n";
    out << "# Genes 4-6: rotation (degrees)\n";
    out << "GENE 4 -180.0 180.0 1.0\n";
    out << "GENE 5 -180.0 180.0 1.0\n";
    out << "GENE 6 -180.0 180.0 1.0\n";

    if (!genes.empty()) {
        out << "# Genes 7-" << (6 + static_cast<int>(genes.size()))
            << ": dihedral angles (degrees)\n";
    }

    for (int g = 0; g < static_cast<int>(genes.size()); ++g) {
        const DihedralGene& dg = genes[g];
        out << "GENE " << (7 + g) << " "
            << std::fixed << std::setprecision(1)
            << dg.min_deg << " " << dg.max_deg << " " << dg.step_deg << "\n";
        out << "DIHEDRAL " << (7 + g) << " "
            << dg.ref_a + 1 << " " << dg.bond_i + 1 << " "
            << dg.bond_j + 1 << " " << dg.ref_b + 1 << "\n";
    }

    out << "# END\n";
    return out.str();
}

// ---------------------------------------------------------------------------
// Main write entry point
// ---------------------------------------------------------------------------

FlexAIDWriterResult FlexAIDWriter::write(const BonMol& mol,
                                          const std::string& lig_name) const {
    FlexAIDWriterResult result;
    result.success = false;

    if (mol.num_atoms() == 0) {
        result.error = "empty molecule";
        return result;
    }

    // Check 3D coordinates are available
    bool any_3d = false;
    for (int i = 0; i < mol.num_atoms(); ++i) {
        if (mol.has_coords(i)) { any_3d = true; break; }
    }
    if (!any_3d) {
        result.error = "no 3D coordinates; cannot write .inp (SMILES-only input)";
        return result;
    }

    // Build spanning tree and internal coordinates
    SpanningTree tree = build_spanning_tree(mol);
    result.internal_coords = compute_internal_coords(mol, tree);

    // Fail-closed geometry invariants before emitting production .inp/.ga.
    {
        const std::string geom_err =
            check_geometry_invariants(mol, result.internal_coords);
        if (!geom_err.empty()) {
            result.error = geom_err;
            return result;
        }
    }

    // Dihedral genes
    result.dihedral_genes     = build_dihedral_genes(mol);
    result.num_dihedral_genes = static_cast<int>(result.dihedral_genes.size());

    // Global positioning atoms
    result.gpa = select_gpa(mol);

    // Generate file contents
    result.inp_content = generate_inp(mol, result.internal_coords, result.gpa, lig_name);
    result.ga_content  = generate_ga(mol, result.dihedral_genes, result.gpa);

    result.num_atoms = mol.num_atoms();
    result.success   = true;
    return result;
}

// ---------------------------------------------------------------------------
// BonMol::to_flexaid() delegation
// ---------------------------------------------------------------------------

} // namespace writer

// Implement the to_flexaid method declared in BonMol.h
BonMol::FlexAIDOutput BonMol::to_flexaid(const std::string& lig_name) const {
    writer::FlexAIDWriter fw;
    writer::FlexAIDWriterResult r = fw.write(*this, lig_name);
    FlexAIDOutput out;
    out.inp_content   = std::move(r.inp_content);
    out.ga_content    = std::move(r.ga_content);
    out.num_atoms     = r.num_atoms;
    out.num_dihedrals = r.num_dihedral_genes;
    return out;
}

} // namespace bonmol
