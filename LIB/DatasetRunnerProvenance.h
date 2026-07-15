// =============================================================================
// DatasetRunnerProvenance.h — provenance.json writer for DatasetRunner
//
// Leaf extraction (P1 of the DatasetRunner split plan): hashes, matrix path
// resolution, and the per-run provenance.json document. No GA, ranking, or
// docking process control. Safe to unit-test with temp dirs.
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
// =============================================================================

#pragma once

#include <string>

namespace dataset {

/// Fields written to per-run provenance.json (key set is stable for audits).
struct RunProvenanceFields {
    std::string dataset;
    std::string matrix_path;
    std::string matrix_md5;
    std::string matrix_sha256;
    std::string binary_path;
    std::string binary_sha256;
    std::string git_commit;
    std::string oracle_site_dir;
    bool oracle_site_dir_set = false;
};

/// Escape a string for embedding inside a JSON string value.
/// Behaviour matches the former DatasetRunner::run() lambda (backslash, quote, newline).
std::string provenance_json_escape(const std::string& s);

/// First non-empty whitespace token of a shell command's stdout, or "".
std::string provenance_cmd_token(const std::string& cmd);

/// File MD5 via `md5 -q` then `md5sum`; empty if path missing or tools fail.
std::string provenance_file_md5(const std::string& path);

/// File SHA-256 via `shasum -a 256` then `sha256sum`; empty if path missing or tools fail.
std::string provenance_file_sha256(const std::string& path);

/// Resolve scoring matrix path with the same precedence as dock children:
/// 1) data_dir/MC_st0r5.2_6.dat when data_dir non-empty and file exists
/// 2) binary-adjacent MC_st0r5.2_6.dat
/// 3) binary/../WRK/MC_st0r5.2_6.dat development fallback
std::string resolve_scoring_matrix_path(const std::string& data_dir,
                                        const std::string& flexaidds_bin);

/// Build provenance fields (hashes + git commit) from run inputs.
/// oracle_site_dir_set is true when oracle_site_dir is non-empty.
RunProvenanceFields build_run_provenance(const std::string& dataset_name,
                                         const std::string& matrix_path,
                                         const std::string& binary_path,
                                         const std::string& oracle_site_dir);

/// Format provenance as the exact JSON document written by DatasetRunner::run.
std::string format_run_provenance_json(const RunProvenanceFields& p);

/// Write output_dir/provenance.json (creates directories). Best-effort.
/// Returns true on success. When log is true, emits the same cout/cerr lines
/// as the former inline block in DatasetRunner::run().
bool write_run_provenance_json(const std::string& output_dir,
                               const RunProvenanceFields& p,
                               bool log = true);

/// Convenience matching the former run() provenance block:
/// resolve matrix → hash → write provenance.json.
bool write_dataset_run_provenance(const std::string& output_dir,
                                  const std::string& dataset_name,
                                  const std::string& flexaidds_bin,
                                  const std::string& data_dir,
                                  const std::string& oracle_site_dir);

} // namespace dataset
