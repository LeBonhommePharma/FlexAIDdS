// FleetRunner.h - immutable Fleet result serialization for DatasetRunner.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "DatasetRunner.h"

#include <string>

namespace fleet {

struct ChunkMetadata {
    std::string campaign_id;
    std::string chunk_id;
    std::string attempt_id;
    std::string worker_id;
    std::string dataset;
    std::string command;
    std::string runner_path;
    std::string runner_sha256;
    std::string engine_path;
    std::string engine_sha256;
    std::string manifest_sha256;
};

class FleetRunner {
public:
    static std::string serialize_chunk_result(
        const ChunkMetadata& metadata,
        const dataset::BenchmarkReport& report,
        const dataset::DockingConfig& config,
        double duration_s);

    // Publishes with no-overwrite semantics. Returns false on any error and
    // leaves an existing destination untouched.
    static bool write_chunk_result_atomic(
        const std::string& destination,
        const std::string& contents,
        std::string* error = nullptr);
};

} // namespace fleet
