// Engine.h — Locked PoseBust evaluation API for DatasetRunner / tests
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "Types.h"

#include <string>

namespace flexaids::posebust {

/// Options for a single evaluation (defaults match dock-suite campaign use).
struct EvaluateOptions {
    Suite suite = Suite::Dock;
    /// Crop protein to heavy atoms within this distance (Å) of ligand COM.
    float protein_crop_A = 10.0f;
    /// Write JSON sidecar + extracted ligand SDF under this directory (empty = skip).
    std::string sidecar_dir;
    /// Optional PDB id for sidecar filenames.
    std::string pdb_id;
};

/// Evaluate a predicted ligand + protein (already loaded).
[[nodiscard]] PoseBustReport evaluate(const Molecule& ligand_pred,
                                      const Molecule& protein,
                                      const Molecule* ligand_true,
                                      const EvaluateOptions& opt = {});

/// Filesystem entry point used by DatasetRunner:
///   complex_pdb  — elected FlexAID pose (receptor+ligand)
///   receptor_pdb — apo/holo protein (may equal complex)
///   crystal_sdf  — optional topology/identity reference (empty ok)
PoseBustReport evaluate_paths(const std::string& complex_pdb,
                              const std::string& receptor_pdb,
                              const std::string& crystal_sdf,
                              const EvaluateOptions& opt = {});

/// Write full report as JSON (one object with checks[]).
bool write_report_json(const PoseBustReport& report, const std::string& path,
                       std::string* err = nullptr);

/// Resolve backend from FLEXAIDDS_POSEBUST / FLEXAIDDS_POSEBUST_BACKEND.
///   FLEXAIDDS_POSEBUST=0              → Off
///   FLEXAIDDS_POSEBUST_BACKEND=off    → Off
///   FLEXAIDDS_POSEBUST_BACKEND=bust   → upstream bust CLI cross-check
///   FLEXAIDDS_POSEBUST_BACKEND=native → NativePoseQC (default)
///   default                          → Native (clean-room PoseBusters path)
/// DatasetRunner.pb_pass is filled from the selected backend's full dock suite.
[[nodiscard]] Backend resolve_backend_from_env();

}  // namespace flexaids::posebust
