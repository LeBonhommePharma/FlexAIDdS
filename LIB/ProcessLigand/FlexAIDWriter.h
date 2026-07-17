// FlexAIDWriter.h — Generate FlexAID .inp and .ga input files from BonMol
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
//
// FlexAID expects two files per ligand:
//
//   <prefix>.inp — PDB-like ATOM records with:
//       columns 1-6:   record type ("ATOM  ")
//       columns 7-11:  atom serial number
//       columns 13-16: atom name (right-padded)
//       columns 17:    alternate location (blank)
//       columns 18-20: residue name (3 chars, right-padded)
//       columns 22:    chain ID (blank for ligand)
//       columns 23-26: residue sequence number
//       columns 27:    insertion code (blank)
//       columns 31-38: X coordinate (8.3f)
//       columns 39-46: Y coordinate (8.3f)
//       columns 47-54: Z coordinate (8.3f)
//       columns 55-60: occupancy (6.2f)
//       columns 61-66: B-factor (6.2f)
//       After PDB block: per-atom type and connectivity section
//
//   <prefix>.ga — Genetic algorithm gene descriptor:
//       NUM_GENES <n>
//       GENE <i> <min> <max> <step>  (for each dihedral gene)
//       GPA <i> <j> <k>              (global positioning atoms)
//
// Internal coordinates (bond lengths, angles, dihedrals) are computed from
// Cartesian coordinates using the Z-matrix convention:
//   atom 1: anchor (no internal coords)
//   atom 2: bond length to atom 1
//   atom 3: bond length to atom 2, angle to atom 1
//   atom ≥4: bond length, angle, dihedral
//
// Global positioning atoms (GPA) are 3 atoms that uniquely position the
// ligand in space (used by gaboom.cpp's ic2cf).
// Reconstruction atoms (rec[0-3]) define the local frame for each dihedral.

#pragma once

#include "BonMol.h"
#include <string>
#include <vector>
#include <array>
#include <optional>

namespace bonmol {
namespace writer {

struct InternalCoord {
    int   atom_idx;
    int   ref1, ref2, ref3; // reference atoms (-1 if undefined)
    float bond_length;      // distance to ref1 (Å)
    float bond_angle;       // angle ref2-ref1-atom (degrees)
    float dihedral;         // dihedral ref3-ref2-ref1-atom (degrees)
};

struct DihedralGene {
    int   bond_i, bond_j; // the rotatable bond atom indices
    int   ref_a, ref_b;   // additional reference atoms for dihedral definition
    float min_deg;        // gene lower bound (degrees)
    float max_deg;        // gene upper bound (degrees)
    float step_deg;       // GA step size (degrees)
};

struct FlexAIDWriterResult {
    std::string inp_content;      // complete .inp file text
    std::string ga_content;       // complete .ga file text
    int         num_atoms;
    int         num_dihedral_genes;
    std::vector<InternalCoord>  internal_coords;
    std::vector<DihedralGene>   dihedral_genes;
    std::array<int, 3>          gpa;           // global positioning atom indices
    bool                        success;
    std::string                 error;
};

class FlexAIDWriter {
public:
    FlexAIDWriter() = default;

    /// Generate both .inp and .ga file content from a BonMol.
    /// The molecule must have 3D coordinates (from SDF/MOL2).
    /// lig_name is the 3-character residue name used in the PDB records.
    FlexAIDWriterResult write(const BonMol& mol,
                              const std::string& lig_name = "LIG") const;

private:
    // -----------------------------------------------------------------------
    // Internal coordinate computation
    // -----------------------------------------------------------------------

    /// Build a spanning tree via BFS from atom 0 (or the heaviest connected atom).
    /// Returns the atom ordering and parent pointers for Z-matrix construction.
    struct SpanningTree {
        std::vector<int> order;   // BFS traversal order
        std::vector<int> parent;  // parent[i] = parent of order[i] in BFS (-1 for root)
        std::vector<int> parent2; // grandparent (for angle definition)
        std::vector<int> parent3; // great-grandparent (for dihedral definition)
    };
    SpanningTree build_spanning_tree(const BonMol& mol) const;

    /// Compute all internal coordinates from Cartesian coords and spanning tree.
    std::vector<InternalCoord> compute_internal_coords(
        const BonMol& mol, const SpanningTree& tree) const;

    /// Production geometry invariants around IC reconstruction:
    /// bond lengths, valence angles, and ring-closure distances.
    /// Returns empty string on success, else a human-readable failure reason.
    std::string check_geometry_invariants(
        const BonMol& mol,
        const std::vector<InternalCoord>& ics) const;

    // -----------------------------------------------------------------------
    // Dihedral gene construction
    // -----------------------------------------------------------------------

    std::vector<DihedralGene> build_dihedral_genes(const BonMol& mol) const;

    // -----------------------------------------------------------------------
    // Global positioning atoms selection
    // -----------------------------------------------------------------------

    /// Choose 3 atoms that maximally span the ligand (max pairwise distance).
    std::array<int, 3> select_gpa(const BonMol& mol) const;

    // -----------------------------------------------------------------------
    // File text generators
    // -----------------------------------------------------------------------

    std::string generate_inp(const BonMol& mol,
                             const std::vector<InternalCoord>& ic,
                             const std::array<int, 3>& gpa,
                             const std::string& lig_name) const;

    std::string generate_ga(const BonMol& mol,
                            const std::vector<DihedralGene>& genes,
                            const std::array<int, 3>& gpa) const;

    // -----------------------------------------------------------------------
    // Geometry helpers
    // -----------------------------------------------------------------------

    /// Compute dihedral angle (degrees) for atoms a-b-c-d.
    float dihedral_angle(const Eigen::Vector3f& a, const Eigen::Vector3f& b,
                         const Eigen::Vector3f& c, const Eigen::Vector3f& d) const;

    /// Compute bond angle (degrees) for atoms a-b-c (angle at b).
    float bond_angle(const Eigen::Vector3f& a, const Eigen::Vector3f& b,
                     const Eigen::Vector3f& c) const;

    /// Format PDB ATOM record line (80 chars).
    std::string format_pdb_atom(int serial, const std::string& name,
                                const std::string& res_name, int res_num,
                                float x, float y, float z,
                                float occ, float bfac) const;

    /// Format atom name for PDB (right-padded to 4 chars, element-aware).
    std::string format_atom_name(const BonMol& mol, int atom_idx,
                                 const std::string& prefix) const;
};

} // namespace writer
} // namespace bonmol
