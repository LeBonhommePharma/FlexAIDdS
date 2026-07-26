// Engine.h — Locked PoseBust evaluation API for DatasetRunner / tests
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "Types.h"

#include <limits>
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
///   FLEXAIDDS_POSEBUST_BACKEND=bust   → upstream bust CLI (default)
///   FLEXAIDDS_POSEBUST_BACKEND=native → NativePoseQC diagnostic
///   default                           → BustCli (official PoseBusters)
/// DatasetRunner.pb_pass is filled from the selected backend's full dock suite.
[[nodiscard]] Backend resolve_backend_from_env();

// ─── Mandatory elected-pose validation (BindingMode rank-0 / elected) ────────

/// Durable outcome of PoseBust on the elected BindingMode pose.
/// Contract: pb_pass is never true unless pb_ran is true and validation completed
/// without hard error. success_pb is always success_rmsd && pb_pass after
/// finalize_success_pb().
struct ElectedPoseBustOutcome {
    bool        pb_ran  = false;
    bool        pb_pass = false;
    bool        success_pb = false;
    std::string pb_backend;  // bust_cli | native_pose_qc | native_pose_qc_fallback | error | skipped_*
    std::string pb_failed_keys;
    int         pb_n_pass   = 0;
    int         pb_n_fail   = 0;
    int         pb_n_checks = 0;
    bool        native_qc_ran  = false;
    bool        native_qc_pass = false;
    std::string native_qc_failed_keys;
    float       pb_min_lig_prot_dist =
        std::numeric_limits<float>::quiet_NaN();
    float       pb_volume_overlap =
        std::numeric_limits<float>::quiet_NaN();
    std::string elected_pose_path;
    std::string pose_sha256;
    std::string posebusters_pose_sha256;
    std::string posebusters_input_sha256;
    std::string error;

    /// success_pb := success_rmsd ∧ pb_pass; never implies pass without pb_ran.
    void finalize_success_pb(bool success_rmsd) noexcept {
        if (!pb_ran) {
            pb_pass = false;
        }
        success_pb = success_rmsd && pb_pass;
    }
};

/// Options for mandatory elected-pose validation.
struct ElectedPoseValidateOptions {
    Backend     backend = Backend::BustCli;
    std::string sidecar_dir;
    std::string pdb_id;
    /// When true (default), Backend::Off with a non-empty elected path still
    /// runs NativePoseQC (mandatory floor). Claim-ready still requires bust_cli.
    bool force_native_when_off = true;
    /// When BustCli is selected but `bust` is missing/fails to start, fall back
    /// to NativePoseQC so pb_ran can still be true (fail-closed on missing CLI
    /// only for STRICT claim_ready via pb_backend).
    bool native_fallback_if_bust_missing = true;
};

/// Validate the elected BindingMode pose (rank-0 / elected_pose.pdb).
///
/// This is the single post-election entry point used by DatasetRunner (and
/// callable from tests). Empty elected path → fail-closed (pb_pass=false).
/// Never returns pb_pass=true without a real backend run that passed.
[[nodiscard]] ElectedPoseBustOutcome validate_elected_pose(
    const std::string& elected_pose_path,
    const std::string& receptor_path,
    const std::string& crystal_sdf,
    const ElectedPoseValidateOptions& opt = {});

}  // namespace flexaids::posebust
