// ProcessLigand.cpp — Unified ligand preprocessing pipeline
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0

#include "ProcessLigand.h"
#include "SmilesParser.h"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

// ---------------------------------------------------------------------------
// SDF and MOL2 readers parse directly into BonMol, including:
//   - SDF V2000: atom coords, bonds, and M CHG formal charge properties
//   - MOL2: atom coords, SYBYL types, bonds, and partial charges
// This avoids the heavy FA_Global dependency of the LIB/ readers
// (SdfReader.h / Mol2Reader.h) which target the legacy atom/resid arrays.
// ---------------------------------------------------------------------------

namespace bonmol {

// ---------------------------------------------------------------------------
// Format detection
// ---------------------------------------------------------------------------

InputFormat detect_format(const std::string& filepath) {
    // Check for bare SMILES: no path separators, no dot extension → SMILES
    // or explicit format override
    std::string lower = filepath;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c){ return std::tolower(c); });

    if (lower.size() >= 4 && lower.substr(lower.size() - 4) == ".sdf") return InputFormat::SDF;
    if (lower.size() >= 5 && lower.substr(lower.size() - 5) == ".mol2") return InputFormat::MOL2;
    if (lower.size() >= 4 && lower.substr(lower.size() - 4) == ".mol")  return InputFormat::SDF;
    if (lower.size() >= 7 && lower.substr(lower.size() - 7) == ".smiles") return InputFormat::SMILES;
    if (lower.size() >= 4 && lower.substr(lower.size() - 4) == ".smi")  return InputFormat::SMILES;

    // If it doesn't look like a file path, treat as SMILES string
    if (filepath.find('/') == std::string::npos &&
        filepath.find('\\') == std::string::npos &&
        filepath.find('.') == std::string::npos) {
        return InputFormat::SMILES;
    }

    return InputFormat::AUTO; // unknown — caller should handle
}

// ---------------------------------------------------------------------------
// Logging helper
// ---------------------------------------------------------------------------

void ProcessLigand::log(const std::string& msg) const {
    if (verbose_) std::cerr << "[ProcessLigand] " << msg << "\n";
}

// ---------------------------------------------------------------------------
// Minimal fallback SDF reader (standalone, no FlexAIDdS LIB dependency)
// Parses SDF V2000 into BonMol.
// ---------------------------------------------------------------------------

BonMol ProcessLigand::load_sdf(const std::string& filepath) {
    std::ifstream f(filepath);
    if (!f) throw std::runtime_error("cannot open SDF file: " + filepath);

    BonMol mol;
    std::string line;

    // Line 1: molecule name
    std::getline(f, line);
    mol.name = line;
    while (!mol.name.empty() && (mol.name.back() == '\r' || mol.name.back() == '\n'))
        mol.name.pop_back();

    // Lines 2-3: program/comment (skip)
    std::getline(f, line);
    std::getline(f, line);

    // Counts line (V2000)
    std::getline(f, line);
    if (line.size() < 6) throw std::runtime_error("SDF: invalid counts line");
    int num_atoms = std::stoi(line.substr(0, 3));
    int num_bonds = std::stoi(line.substr(3, 3));

    // Atom block
    std::vector<Element> atom_elements(num_atoms, Element::Unknown);
    for (int i = 0; i < num_atoms; ++i) {
        std::getline(f, line);
        if (line.size() < 39)
            throw std::runtime_error("SDF: atom line too short at index " + std::to_string(i));
        float x = std::stof(line.substr(0, 10));
        float y = std::stof(line.substr(10, 10));
        float z = std::stof(line.substr(20, 10));
        std::string sym = line.substr(31, 3);
        // Trim whitespace
        auto trim = [](std::string s) -> std::string {
            size_t start = s.find_first_not_of(" \t\r\n");
            size_t end   = s.find_last_not_of(" \t\r\n");
            if (start == std::string::npos) return "";
            return s.substr(start, end - start + 1);
        };
        sym = trim(sym);
        Element elem = element_from_symbol(sym);
        atom_elements[i] = elem;
        int idx = mol.add_atom(elem, x, y, z);
        (void)idx;
    }

    // Bond block
    for (int i = 0; i < num_bonds; ++i) {
        std::getline(f, line);
        if (line.size() < 9)
            throw std::runtime_error("SDF: bond line too short at index " + std::to_string(i));
        int a1    = std::stoi(line.substr(0, 3)) - 1;
        int a2    = std::stoi(line.substr(3, 3)) - 1;
        int btype = std::stoi(line.substr(6, 3));
        BondOrder order;
        bool arom = false;
        switch (btype) {
            case 1: order = BondOrder::SINGLE;   break;
            case 2: order = BondOrder::DOUBLE;   break;
            case 3: order = BondOrder::TRIPLE;   break;
            case 4: order = BondOrder::AROMATIC; arom = true; break;
            default: order = BondOrder::SINGLE;  break;
        }
        mol.add_bond(a1, a2, order, arom);
    }

    // Properties block: parse M  CHG lines for formal charges
    while (std::getline(f, line)) {
        // Strip trailing whitespace
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n' || line.back() == ' '))
            line.pop_back();
        if (line == "M  END" || line == "$$$$") break;

        // M  CHG  n  aaa vvv  aaa vvv ...
        // n = number of charge entries, each pair is atom_index(1-based) and charge
        if (line.size() >= 9 && line.substr(0, 6) == "M  CHG") {
            std::istringstream ss(line.substr(6));
            int n_entries = 0;
            ss >> n_entries;
            for (int i = 0; i < n_entries && ss; ++i) {
                int atom_idx = 0, charge = 0;
                ss >> atom_idx >> charge;
                if (atom_idx >= 1 && atom_idx <= num_atoms) {
                    mol.atoms[atom_idx - 1].formal_charge = charge;
                }
            }
        }
    }

    mol.finalize();
    return mol;
}

// ---------------------------------------------------------------------------
// Minimal fallback MOL2 reader
// ---------------------------------------------------------------------------

BonMol ProcessLigand::load_mol2(const std::string& filepath) {
    std::ifstream f(filepath);
    if (!f) throw std::runtime_error("cannot open MOL2 file: " + filepath);

    BonMol mol;
    std::string line;
    bool in_atom = false, in_bond = false;

    while (std::getline(f, line)) {
        // Trim trailing CR
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n'))
            line.pop_back();

        if (line.find("@<TRIPOS>MOLECULE") != std::string::npos) {
            in_atom = false; in_bond = false;
            std::getline(f, line); // molecule name
            while (!line.empty() && (line.back() == '\r' || line.back() == '\n'))
                line.pop_back();
            mol.name = line;
            continue;
        }
        if (line.find("@<TRIPOS>ATOM") != std::string::npos) {
            in_atom = true; in_bond = false; continue;
        }
        if (line.find("@<TRIPOS>BOND") != std::string::npos) {
            in_atom = false; in_bond = true; continue;
        }
        if (line.size() > 1 && line[0] == '@') {
            in_atom = false; in_bond = false; continue;
        }

        if (in_atom && !line.empty() && !std::isspace(line[0])) {
            std::istringstream ss(line);
            int atom_id; std::string aname, atype;
            float x, y, z;
            if (!(ss >> atom_id >> aname >> x >> y >> z >> atype)) continue;

            // Extract element from atype (e.g. "C.3" → "C", "N.ar" → "N")
            std::string elem_str = atype;
            auto dot = elem_str.find('.');
            if (dot != std::string::npos) elem_str = elem_str.substr(0, dot);

            Element elem = element_from_symbol(elem_str);
            int idx = mol.add_atom(elem, x, y, z);

            // Read optional substructure fields and partial charge
            // MOL2 format: atom_id name x y z type [subst_id subst_name charge]
            int res_id = 0; std::string res_name;
            float charge = 0.0f;
            if (ss >> res_id >> res_name >> charge) {
                mol.atoms[idx].partial_charge = charge;
            }

            // Decode SYBYL type from atype string
            // (mirrors Mol2Reader.cpp sybyl_to_flexaid_type mapping)
            auto sybyl_map = [](const std::string& t) -> int {
                if (t == "C.3")   return 1;  if (t == "C.2")   return 2;
                if (t == "C.ar")  return 3;  if (t == "C.1")   return 0;
                if (t == "N.3")   return 4;  if (t == "N.2")   return 5;
                if (t == "N.ar")  return 6;  if (t == "N.am")  return 7;
                if (t == "N.pl3") return 8;  if (t == "N.4")   return 9;
                if (t == "O.3")   return 10; if (t == "O.2")   return 11;
                if (t == "O.co2") return 12; if (t == "F")     return 13;
                if (t == "Cl")    return 14; if (t == "Br")    return 15;
                if (t == "S.3")   return 16; if (t == "S.2")   return 17;
                if (t == "S.O")   return 18; if (t == "S.O2")  return 19;
                if (t == "P.3")   return 20; if (t == "I")     return 21;
                if (t == "H")     return 22; if (t == "Fe")    return 30;
                return 1; // fallback
            };
            mol.atoms[idx].sybyl_type = sybyl_map(atype);
            // Mark aromatic from type
            if (atype == "C.ar" || atype == "N.ar")
                mol.atoms[idx].is_aromatic = true;

            // Atom name
            size_t copy_len = std::min(aname.size(), (size_t)4);
            std::memcpy(mol.atoms[idx].name, aname.c_str(), copy_len);
        }

        if (in_bond && !line.empty() && !std::isspace(line[0])) {
            std::istringstream ss(line);
            int bond_id, a1, a2; std::string btype;
            if (!(ss >> bond_id >> a1 >> a2 >> btype)) continue;
            a1--; a2--;
            BondOrder order = BondOrder::SINGLE;
            bool arom = false;
            if (btype == "2") order = BondOrder::DOUBLE;
            else if (btype == "3") order = BondOrder::TRIPLE;
            else if (btype == "ar") { order = BondOrder::AROMATIC; arom = true; }
            mol.add_bond(a1, a2, order, arom);
        }
    }

    mol.finalize();
    return mol;
}

// ---------------------------------------------------------------------------
// Stage 1: Parse
// ---------------------------------------------------------------------------

StageResult ProcessLigand::stage_parse(const ProcessOptions& opts, BonMol& mol) {
    StageResult r;
    r.stage = "parse";
    try {
        InputFormat fmt = opts.format;
        if (fmt == InputFormat::AUTO) fmt = detect_format(opts.input);

        switch (fmt) {
            case InputFormat::SMILES: {
                SmilesParser parser;
                auto parsed = parser.parse(opts.input);
                mol = std::move(parsed.mol);
                if (!parsed.warnings.empty()) {
                    std::ostringstream ws;
                    for (const auto& w : parsed.warnings) ws << w << "; ";
                    r.message = "SMILES warnings: " + ws.str();
                }
                break;
            }
            case InputFormat::SDF:
                mol = load_sdf(opts.input);
                break;
            case InputFormat::MOL2:
                mol = load_mol2(opts.input);
                break;
            default:
                // Try SDF then MOL2
                try { mol = load_sdf(opts.input); }
                catch (...) { mol = load_mol2(opts.input); }
                break;
        }

        if (!opts.lig_name.empty() && mol.name.empty())
            mol.name = opts.lig_name;

    } catch (const std::exception& e) {
        r.ok      = false;
        r.message = std::string("parse failed: ") + e.what();
    }
    return r;
}

// ---------------------------------------------------------------------------
// Stage 2: Validate
// ---------------------------------------------------------------------------

StageResult ProcessLigand::stage_validate(const ProcessOptions& opts, BonMol& mol,
                                           valence::ValenceCheckResult& valence_result) {
    StageResult r;
    r.stage = "validate";

    // Structural guards (peptide / macrocycle)
    auto v = mol.validate();
    if (!v.valid) {
        if (v.has_peptide_backbone && opts.allow_peptides) {
            r.message = "peptide warning suppressed (--allow-peptides)";
        } else if (!v.error.empty() && v.error.find("macrocycle") != std::string::npos
                   && opts.allow_macrocycles) {
            r.message = "macrocycle warning suppressed (--allow-macrocycles)";
        } else {
            r.ok      = false;
            r.message = v.error;
            return r;
        }
    }

    // Valence check
    valence_result = valence::check_valence(mol);
    if (!valence_result.valid) {
        std::ostringstream es;
        for (const auto& e : valence_result.errors) es << e.message << "; ";
        r.ok      = false;
        r.message = "valence errors: " + es.str();
        return r;
    }
    if (opts.strict_valence && !valence_result.warnings.empty()) {
        std::ostringstream es;
        for (const auto& e : valence_result.warnings) es << e.message << "; ";
        r.ok      = false;
        r.message = "valence warnings (strict mode): " + es.str();
    }
    return r;
}

// ---------------------------------------------------------------------------
// Stage 3: Ring perception
// ---------------------------------------------------------------------------

StageResult ProcessLigand::stage_rings(BonMol& mol,
                                        ring_perception::RingPerceptionResult& rr) {
    StageResult r;
    r.stage = "ring_perception";
    try {
        rr = ring_perception::perceive_rings(mol);
        r.message = std::to_string(rr.num_rings) + " rings found (circuit rank "
                  + std::to_string(rr.circuit_rank) + ")";
    } catch (const std::exception& e) {
        r.ok      = false;
        r.message = std::string("ring perception failed: ") + e.what();
    }
    return r;
}

// ---------------------------------------------------------------------------
// Stage 4: Aromaticity
// ---------------------------------------------------------------------------

StageResult ProcessLigand::stage_aromaticity(BonMol& mol,
                                              aromaticity::AromaticityResult& ar) {
    StageResult r;
    r.stage = "aromaticity";
    try {
        ar = aromaticity::assign_aromaticity(mol);
        r.message = std::to_string(ar.num_aromatic_rings) + " aromatic rings, "
                  + std::to_string(ar.num_aromatic_atoms) + " aromatic atoms";
        if (!ar.kekulized) r.message += " (Kekulization partial)";
    } catch (const std::exception& e) {
        r.ok      = false;
        r.message = std::string("aromaticity failed: ") + e.what();
    }
    return r;
}

// ---------------------------------------------------------------------------
// Stage 5: Rotatable bonds
// ---------------------------------------------------------------------------

StageResult ProcessLigand::stage_rotatable(BonMol& mol,
                                            rotatable_bonds::RotatableBondsResult& rr) {
    StageResult r;
    r.stage = "rotatable_bonds";
    try {
        rr = rotatable_bonds::identify_rotatable_bonds(mol);
        r.message = std::to_string(rr.count) + " rotatable bonds";
    } catch (const std::exception& e) {
        r.ok      = false;
        r.message = std::string("rotatable bond identification failed: ") + e.what();
    }
    return r;
}

// ---------------------------------------------------------------------------
// Stage 6: SYBYL typing
// ---------------------------------------------------------------------------

StageResult ProcessLigand::stage_typing(BonMol& mol) {
    StageResult r;
    r.stage = "sybyl_typing";
    try {
        sybyl::assign_sybyl_types(mol);
        r.message = "SYBYL types and 256-type encoding assigned";
    } catch (const std::exception& e) {
        r.ok      = false;
        r.message = std::string("SYBYL typing failed: ") + e.what();
    }
    return r;
}

// ---------------------------------------------------------------------------
// Stage 7: Write output
// ---------------------------------------------------------------------------

StageResult ProcessLigand::stage_write(const ProcessOptions& opts, const BonMol& mol,
                                        writer::FlexAIDWriterResult& wr) {
    StageResult r;
    r.stage = "write";

    if (opts.validate_only) {
        r.message = "validate-only mode; skipping output";
        return r;
    }

    writer::FlexAIDWriter fw;
    wr = fw.write(mol, opts.lig_name);

    if (!wr.success) {
        r.ok      = false;
        r.message = "writer error: " + wr.error;
        return r;
    }

    // Write to files if output prefix given
    if (!opts.output_prefix.empty()) {
        if (opts.write_inp) {
            std::string inp_path = opts.output_prefix + ".inp";
            std::ofstream of(inp_path);
            if (!of) {
                r.ok      = false;
                r.message = "cannot write " + inp_path;
                return r;
            }
            of << wr.inp_content;
        }
        if (opts.write_ga) {
            std::string ga_path = opts.output_prefix + ".ga";
            std::ofstream of(ga_path);
            if (!of) {
                r.ok      = false;
                r.message = "cannot write " + ga_path;
                return r;
            }
            of << wr.ga_content;
        }
    }

    r.message = "wrote " + std::to_string(wr.num_atoms) + " atoms, "
              + std::to_string(wr.num_dihedral_genes) + " dihedral genes";
    return r;
}

// ---------------------------------------------------------------------------
// Main pipeline runner
// ---------------------------------------------------------------------------

ProcessResult ProcessLigand::run(const ProcessOptions& opts) {
    verbose_ = opts.verbose;
    ProcessResult result;

    BonMol mol;

    // Stage 1: Parse
    log("Stage 1: parsing input");
    StageResult s1 = stage_parse(opts, mol);
    result.stage_results.push_back(s1);
    if (!s1.ok) { result.error = s1.message; return result; }

    // Stage 2: Validate
    log("Stage 2: validating structure");
    StageResult s2 = stage_validate(opts, mol, result.valence_result);
    result.stage_results.push_back(s2);
    if (!s2.ok) { result.error = s2.message; return result; }

    // Stage 3: Ring perception
    log("Stage 3: ring perception");
    StageResult s3 = stage_rings(mol, result.ring_result);
    result.stage_results.push_back(s3);
    if (!s3.ok) { result.error = s3.message; return result; }

    // Stage 4: Aromaticity
    log("Stage 4: aromaticity");
    StageResult s4 = stage_aromaticity(mol, result.arom_result);
    result.stage_results.push_back(s4);
    if (!s4.ok) { result.error = s4.message; return result; }

    // Stage 5: Rotatable bonds
    log("Stage 5: rotatable bonds");
    StageResult s5 = stage_rotatable(mol, result.rot_result);
    result.stage_results.push_back(s5);
    if (!s5.ok) { result.error = s5.message; return result; }

    // Stage 6: SYBYL typing
    log("Stage 6: SYBYL typing");
    StageResult s6 = stage_typing(mol);
    result.stage_results.push_back(s6);
    if (!s6.ok) { result.error = s6.message; return result; }

    // Stage 7: Write
    log("Stage 7: writing output");
    StageResult s7 = stage_write(opts, mol, result.writer_result);
    result.stage_results.push_back(s7);
    if (!s7.ok) { result.error = s7.message; return result; }

    // Populate summary
    result.mol              = std::move(mol);
    result.num_atoms        = result.mol.num_atoms();
    result.num_heavy_atoms  = result.mol.num_heavy_atoms();
    result.num_rings        = result.ring_result.num_rings;
    result.num_arom_rings   = result.arom_result.num_aromatic_rings;
    result.num_rot_bonds    = result.rot_result.count;
    result.molecular_weight = result.mol.molecular_weight;
    result.success          = true;

    return result;
}

ProcessResult ProcessLigand::process(const ProcessOptions& opts) {
    ProcessLigand pl;
    return pl.run(opts);
}

// ---------------------------------------------------------------------------
// Free function wrappers
// ---------------------------------------------------------------------------

ProcessResult process_smiles(const std::string& smiles,
                              const std::string& output_prefix,
                              const std::string& lig_name) {
    ProcessOptions opts;
    opts.input         = smiles;
    opts.format        = InputFormat::SMILES;
    opts.output_prefix = output_prefix;
    opts.lig_name      = lig_name;
    return ProcessLigand::process(opts);
}

ProcessResult process_sdf(const std::string& filepath,
                           const std::string& output_prefix,
                           const std::string& lig_name) {
    ProcessOptions opts;
    opts.input         = filepath;
    opts.format        = InputFormat::SDF;
    opts.output_prefix = output_prefix;
    opts.lig_name      = lig_name;
    return ProcessLigand::process(opts);
}

ProcessResult process_mol2(const std::string& filepath,
                            const std::string& output_prefix,
                            const std::string& lig_name) {
    ProcessOptions opts;
    opts.input         = filepath;
    opts.format        = InputFormat::MOL2;
    opts.output_prefix = output_prefix;
    opts.lig_name      = lig_name;
    return ProcessLigand::process(opts);
}

// ---------------------------------------------------------------------------
// BonMol factory free functions (declared in BonMol.h)
// These are thin wrappers around ProcessLigand's private loaders.
// ---------------------------------------------------------------------------

BonMol from_sdf(const std::string& filepath) {
    ProcessLigand pl;
    return pl.public_load_sdf(filepath);
}

BonMol from_mol2(const std::string& filepath) {
    ProcessLigand pl;
    return pl.public_load_mol2(filepath);
}

} // namespace bonmol
