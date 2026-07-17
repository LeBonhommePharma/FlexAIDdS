// ProcessLigand.h — Unified ligand preprocessing pipeline for FlexAIDdS
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
//
// Pipeline stages:
//   1. Parse input (SMILES → SmilesParser, SDF → SdfReader, MOL2 → Mol2Reader)
//   2. Validate structure (ValenceChecker + peptide/macrocycle guards)
//   3. Perceive rings (SSSR via RingPerception)
//   4. Assign aromaticity (Hückel via Aromaticity)
//   5. Identify rotatable bonds (RotatableBonds)
//   6. Assign SYBYL types + 256-type encoding (SybylTyper)
//   7. Generate FlexAID .inp and .ga files (FlexAIDWriter)
//
// The pipeline is designed to be run as:
//   ProcessLigand pl;
//   auto result = pl.run(options);
//
// or as a one-shot function:
//   auto result = process_ligand(options);

#pragma once

#include "BonMol.h"
#include "RingPerception.h"
#include "Aromaticity.h"
#include "RotatableBonds.h"
#include "ValenceChecker.h"
#include "SybylTyper.h"
#include "FlexAIDWriter.h"

#include <string>
#include <vector>
#include <optional>
#include <filesystem>

namespace bonmol {

// ---------------------------------------------------------------------------
// Input format
// ---------------------------------------------------------------------------

enum class InputFormat {
    AUTO,   // detect from file extension
    SMILES, // bare SMILES string (no file)
    SDF,    // MDL SDF V2000
    MOL2    // SYBYL MOL2
};

/// Detect input format from file extension.
InputFormat detect_format(const std::string& filepath);

// ---------------------------------------------------------------------------
// Pipeline options
// ---------------------------------------------------------------------------

struct ProcessOptions {
    // Input
    std::string  input;         // filepath or SMILES string
    InputFormat  format       = InputFormat::AUTO;
    std::string  lig_name     = "LIG"; // 3-char residue name for PDB records

    // Output
    std::string  output_prefix; // prefix for .inp and .ga files (empty = no file write)
    bool         write_inp    = true;
    bool         write_ga     = true;

    // Feature flags
    bool         validate_only  = false; // skip output, only validate
    bool         strict_valence = false; // fail on valence warnings (not just errors)
    bool         allow_macrocycles = false; // bypass macrocycle guard
    bool         allow_peptides    = false; // bypass peptide guard

    // Verbosity
    bool         verbose = false;
};

// ---------------------------------------------------------------------------
// Per-stage diagnostic data
// ---------------------------------------------------------------------------

struct StageResult {
    bool        ok     = true;
    std::string stage;
    std::string message;
};

// ---------------------------------------------------------------------------
// Full pipeline result
// ---------------------------------------------------------------------------

struct ProcessResult {
    bool success = false;
    std::string error;

    // Molecule after all stages
    BonMol mol;

    // Stage diagnostics
    std::vector<StageResult>                   stage_results;
    ring_perception::RingPerceptionResult      ring_result;
    aromaticity::AromaticityResult             arom_result;
    rotatable_bonds::RotatableBondsResult      rot_result;
    valence::ValenceCheckResult                valence_result;
    writer::FlexAIDWriterResult                writer_result;

    // Summary
    int num_atoms        = 0;
    int num_heavy_atoms  = 0;
    int num_rings        = 0;
    int num_arom_rings   = 0;
    int num_rot_bonds    = 0;
    float molecular_weight = 0.0f;
};

// ---------------------------------------------------------------------------
// Pipeline class
// ---------------------------------------------------------------------------

class ProcessLigand {
public:
    ProcessLigand() = default;

    /// Run the full pipeline. Returns a complete ProcessResult.
    ProcessResult run(const ProcessOptions& opts);

    /// One-shot static convenience wrapper.
    static ProcessResult process(const ProcessOptions& opts);

private:
    // -----------------------------------------------------------------------
    // Stage implementations
    // -----------------------------------------------------------------------

    /// Stage 1: Parse input into BonMol
    StageResult stage_parse(const ProcessOptions& opts, BonMol& mol);

    /// Stage 2: Structural validation
    StageResult stage_validate(const ProcessOptions& opts, BonMol& mol,
                               valence::ValenceCheckResult& valence_result);

    /// Stage 3: Ring perception
    StageResult stage_rings(BonMol& mol,
                            ring_perception::RingPerceptionResult& ring_result);

    /// Stage 4: Aromaticity assignment
    StageResult stage_aromaticity(BonMol& mol,
                                  aromaticity::AromaticityResult& arom_result);

    /// Stage 5: Rotatable bonds
    StageResult stage_rotatable(BonMol& mol,
                                rotatable_bonds::RotatableBondsResult& rot_result);

    /// Stage 6: SYBYL types + 256-type encoding
    StageResult stage_typing(BonMol& mol);

    /// Stage 7: Write output files
    StageResult stage_write(const ProcessOptions& opts, const BonMol& mol,
                            writer::FlexAIDWriterResult& writer_result,
                            bool enforce_geometry_invariants = true);

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------

    /// Load SDF file and return BonMol (delegates to SdfReader).
    BonMol load_sdf(const std::string& filepath);

    /// Load MOL2 file and return BonMol (delegates to Mol2Reader).
    BonMol load_mol2(const std::string& filepath);

public: // factory access — used by from_sdf/from_mol2 free functions
    BonMol public_load_sdf(const std::string& fp) { return load_sdf(fp); }
    BonMol public_load_mol2(const std::string& fp) { return load_mol2(fp); }

    bool verbose_ = false;
    void log(const std::string& msg) const;
};

// ---------------------------------------------------------------------------
// Free function convenience wrappers
// ---------------------------------------------------------------------------

/// Run the pipeline on a SMILES string.
ProcessResult process_smiles(const std::string& smiles,
                              const std::string& output_prefix = "",
                              const std::string& lig_name = "LIG");

/// Run the pipeline on an SDF file.
ProcessResult process_sdf(const std::string& filepath,
                           const std::string& output_prefix = "",
                           const std::string& lig_name = "LIG");

/// Run the pipeline on a MOL2 file.
ProcessResult process_mol2(const std::string& filepath,
                            const std::string& output_prefix = "",
                            const std::string& lig_name = "LIG");

} // namespace bonmol
