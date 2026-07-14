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
    std::string csv_path;          // written report if sidecar set
    std::string raw_csv;           // last line of bust csv
};

/// Resolve bust binary: FLEXAIDDS_POSEBUSTERS_BIN, else PATH, else
/// $FLEXAIDDS_ROOT/.venv-posebusters/bin/bust, else repo-relative.
[[nodiscard]] std::string resolve_bust_binary();

/// Run upstream bust on predicted ligand SDF vs protein (+ optional crystal).
/// RMSD column is recorded but excluded from pb_pass (RMSD is success_rmsd).
BustCliResult run_upstream_bust(const std::string& pred_sdf,
                                const std::string& protein_pdb,
                                const std::string& crystal_sdf,
                                const std::string& sidecar_dir = {},
                                const std::string& stem = "pose");

/// SHA-256 hex of a file (empty on failure).
[[nodiscard]] std::string sha256_file(const std::string& path);

/// Copy file; returns false on failure.
bool copy_file_atomic(const std::string& src, const std::string& dst,
                      std::string* err = nullptr);

}  // namespace flexaids::posebust
