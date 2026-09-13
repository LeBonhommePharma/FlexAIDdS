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
///
/// eigen_version is recorded INLINE here, beside binary_sha256, rather than
/// being left to a join through the build directory.  That is a deliberate
/// consequence of a measurement: of 3044 stored receipts, 100% carried a
/// binary_sha256 but only 42 of 156 arms (27%) could be attributed to an
/// Eigen version retroactively, because attribution needed the engine's
/// build directory to still exist and CMakeCache.txt to still be readable.
/// A hash is a pointer into a mutable filesystem; once the build dir is
/// deleted the pointer dangles and the receipt can no longer be
/// interpreted.  A version STRING in the receipt survives deletion of
/// everything else, so it needs no second file to be readable.
///
/// This matters for entropy specifically: tENCoM classifies a zero mode by
/// the SIGN of a numerically-zero eigenvalue, so the eigensolver is part of
/// what produced a dS_vib mode count, and an arm that does not name its
/// solver cannot be compared with one from another machine.
struct RunProvenanceFields {
    std::string dataset;
    std::string matrix_path;
    std::string matrix_md5;
    std::string matrix_sha256;
    std::string binary_path;
    std::string binary_sha256;
    /// Short-form commit (git rev-parse --short), stamped at COMPILE time.
    /// Populated by build_run_provenance() from FLEXAIDS_GIT_COMMIT.
    ///
    /// It used to be measured at RUN time by shelling out to
    /// `git rev-parse HEAD`, which returned the empty string whenever the
    /// process CWD was outside a git work tree.  Benchmark runs execute in
    /// the results tree, not the repository, so the field was empty on
    /// 3017 of 3044 stored receipts (0.9% populated) while
    /// scripts/check_run_receipt.py lists it as REQUIRED.  The 27 that did
    /// populate carry 40-hex SHAs from runs that happened to start inside
    /// the repo; new receipts carry the 8-hex short form, which is a prefix
    /// of the full SHA and therefore still joinable against them.
    std::string git_commit;
    /// Eigen release string the engine was compiled against, e.g. the value
    /// of EIGEN_VERSION_STRING.  NOT the EIGEN_WORLD.MAJOR.MINOR macro
    /// triple, which is not the release number above Eigen 5 because WORLD
    /// is frozen at 3.  See cmake/FlexAIDEigen.cmake.
    std::string eigen_version;
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
