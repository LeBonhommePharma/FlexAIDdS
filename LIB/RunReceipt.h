// RunReceipt.h — Campaign / DatasetRunner RUN_RECEIPT.json provenance
//
// Aligns the C++ DatasetRunner path with C0 iCloud launch scripts
// (RUN_RECEIPT.json schema: matrix/binary hashes, GA knobs, ProtocolConfig).
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "ProtocolConfig.h"

#include <cstdint>
#include <string>

namespace flexaids {

/// Inputs for a campaign/run receipt. Hash fields may be precomputed or left
/// empty for the writer to fill when paths exist.
struct RunReceiptInput {
    std::string run_id;
    std::string started_utc;   ///< ISO-8601 UTC, e.g. 2026-07-15T02:54:52Z
    std::string output;
    std::string dataset;
    std::string mode;          ///< e.g. defined-cleft-redock
    double temperature_K{300.0};
    int pop{1000};
    int gen{2000};
    int restarts{5};
    std::uint64_t seed_base{0};
    bool seed_elitism{true};

    std::string matrix_path;
    std::string matrix_md5;
    std::string matrix_sha256;
    std::string binary_path;
    std::string binary_sha256;
    std::string runner_path;
    std::string runner_sha256;
    std::string git_commit;
    std::string oracle_site_dir;
    bool oracle_site_dir_set{false};

    ProtocolConfig protocol;
};

/// schema_version for RUN_RECEIPT.json (bump when keys change incompatibly).
inline constexpr int kRunReceiptSchemaVersion = 1;

/// Build a JSON object string (pretty-printed, trailing newline not included).
[[nodiscard]] std::string build_run_receipt_json(const RunReceiptInput& in);

/// Write RUN_RECEIPT.json (and optional legacy provenance.json) under output_dir.
/// Returns true if RUN_RECEIPT.json was written successfully.
bool write_run_receipt(const std::string& output_dir,
                       const RunReceiptInput& in,
                       bool also_write_provenance_json = true);

/// Current UTC time as YYYY-MM-DDTHH:MM:SSZ (best-effort; empty on failure).
[[nodiscard]] std::string utc_now_iso8601();

}  // namespace flexaids
