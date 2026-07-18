// BustCli.h — Authoritative upstream PoseBusters (`bust`) gate
//
// Native C++ QC (this library) is diagnostic only. Official pb_pass comes from
// the installed PoseBusters CLI (BSD-licensed, not vendored).
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "Types.h"

#include <string>
#include <vector>

namespace flexaids::posebust {

struct BustCliResult {
    bool        ran     = false;
    bool        pb_pass = false;   // all dock-suite booleans True (excl. RMSD column)
    std::string backend = "bust_cli";
    std::string error;
    std::string failed_keys;       // semicolon-separated failed columns
    int         n_pass  = 0;
    int         n_fail  = 0;
    int         n_checks = 0;
    std::string csv_path;          // validated CSV report if sidecar set
    std::string raw_csv;           // full raw stdout (preserved before schema fail)
    // Strict PoseBusters 0.6.5 redock contract + receipts.
    std::string schema_id = "posebusters-0.6.5-redock-csv-v1";
    std::string package_name = "posebusters";
    std::string package_version;       // distribution version (must be 0.6.5)
    std::string package_record_path;   // installed dist-info/RECORD manifest
    std::string package_record_sha256;
    std::string config_name = "redock";
    std::string config_path;           // exact installed redock.yml consumed
    std::string config_sha256;
    std::string bust_path;             // resolved console-script launcher
    std::string bust_sha256;           // launcher SHA-256 (not package identity)
    std::string bust_version;          // raw `bust --version` output
    std::vector<std::string> argv;     // exact argv passed to exec
    std::string argv_joined;           // human-readable argv (legacy receipt)
    int         exit_status = -1;
    std::string pred_sdf_path;
    std::string pred_sdf_sha256;
    std::string protein_pdb_path;
    std::string protein_pdb_sha256;
    std::string crystal_sdf_path;
    std::string crystal_sdf_sha256;
    std::string raw_csv_path;
    std::string raw_csv_sha256;        // SHA-256 of raw stdout sidecar
    std::string csv_sha256;            // SHA-256 of validated CSV report
};

/// Resolve bust binary: FLEXAIDDS_POSEBUSTERS_BIN, else PATH, else
/// $FLEXAIDDS_ROOT/.venv-posebusters/bin/bust, else repo-relative.
[[nodiscard]] std::string resolve_bust_binary();

/// Run upstream bust on predicted ligand SDF vs protein (+ optional crystal).
/// RMSD column is recorded but excluded from pb_pass (RMSD is success_rmsd).
/// Always preserves raw_csv before any schema-failure return.
BustCliResult run_upstream_bust(const std::string& pred_sdf,
                                const std::string& protein_pdb,
                                const std::string& crystal_sdf,
                                const std::string& sidecar_dir = {},
                                const std::string& stem = "pose");

/// Parse already-captured bust CSV body into a result (schema + pb_pass).
/// Used by unit tests and by run_upstream_bust after raw_csv is preserved.
/// raw_csv is set from `csv_body` when r.raw_csv is empty.
void apply_bust_csv_schema(const std::string& csv_body, BustCliResult& r);

/// Complete JSON receipt for package/version/config/command/input/output
/// identity. The receipt is deterministic and contains no inferred hashes.
[[nodiscard]] std::string format_bust_receipt_json(const BustCliResult& r);

/// SHA-256 hex of a file (empty on failure).
[[nodiscard]] std::string sha256_file(const std::string& path);

/// Copy file; returns false on failure.
bool copy_file_atomic(const std::string& src, const std::string& dst,
                      std::string* err = nullptr);

}  // namespace flexaids::posebust
