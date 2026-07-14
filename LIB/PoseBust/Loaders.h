// Loaders.h — Clean-room SDF/PDB molecular I/O for PoseBust
//
// Apache-2.0. Implemented from first principles; no posebusters/RDKit code.
// Control flow uses bool returns + optional err strings (no exceptions required).
#pragma once

#include "Types.h"

#include <string>

namespace flexaids::posebust {

// Map element symbol to atomic number (0 if unknown).
// Supports H, C, N, O, F, P, S, Cl, Br, I, B, Si (and common aliases).
int atomic_number(const std::string& elem);

// Covalent radius (Å) for bond inference. Returns 0.77 for unknown Z.
float covalent_radius(int Z);

// Parse first molecule from a V2000 SDF/MOL file into `out`.
// On failure returns false and optionally writes a diagnostic to `err`.
bool load_sdf(const std::string& path, Molecule& out, std::string* err = nullptr);

// Load ligand-like atoms from a PDB (FALLBACK only — prefer FlexAID CONECT path).
//   1) HETATM excluding solvent/ions AND common cofactors (HEM, NAD, …)
//   2) if empty: non-protein ATOM residues (not standard 20 AA + UNK)
// Bonds are inferred from covalent radii unless assign_topology_from_reference is used.
// WARNING: On FlexAID complex PDBs, never use this alone — cofactors remain as HETATM.
bool load_pdb_ligand(const std::string& path, Molecule& out, std::string* err = nullptr);

// FlexAID elected-pose extraction (correct path for DatasetRunner):
//   1) Atoms referenced by CONECT records (serials, typically 90001+ for ligand)
//   2) Else REMARK "optimizable residue <name> <resseq>" atoms
//   3) Else load_pdb_ligand fallback with cofactor blacklist
// Coordinates come from the pose; bonds still need assign_topology_from_reference.
bool load_pdb_flexaid_ligand(const std::string& path, Molecule& out,
                             std::string* err = nullptr);

// Copy bond topology from reference (crystal SDF) onto pred when heavy-atom
// element sequences match (same count, same element order). Replaces inferred
// bonds. Returns false if topology cannot be assigned safely.
bool assign_topology_from_reference(Molecule& pred, const Molecule& reference,
                                    std::string* err = nullptr);

// Load protein ATOM heavy atoms only (standard amino-acid residue names).
bool load_pdb_protein_heavy(const std::string& path, Molecule& out,
                            std::string* err = nullptr);

// Write a minimal V2000 SDF for debugging.
bool write_sdf(const Molecule& mol, const std::string& path, std::string* err = nullptr);

}  // namespace flexaids::posebust
