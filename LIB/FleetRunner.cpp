// FleetRunner.cpp - immutable Fleet result serialization for DatasetRunner.
// SPDX-License-Identifier: Apache-2.0

#include "FleetRunner.h"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <system_error>

#ifndef _WIN32
#include <fcntl.h>
#include <unistd.h>
#endif

namespace fs = std::filesystem;

namespace {

std::string json_escape(const std::string& value) {
    std::ostringstream out;
    for (const unsigned char c : value) {
        switch (c) {
            case '"': out << "\\\""; break;
            case '\\': out << "\\\\"; break;
            case '\b': out << "\\b"; break;
            case '\f': out << "\\f"; break;
            case '\n': out << "\\n"; break;
            case '\r': out << "\\r"; break;
            case '\t': out << "\\t"; break;
            default:
                if (c < 0x20) {
                    out << "\\u" << std::hex << std::setw(4)
                        << std::setfill('0') << static_cast<int>(c)
                        << std::dec << std::setfill(' ');
                } else {
                    out << static_cast<char>(c);
                }
        }
    }
    return out.str();
}

void json_string(std::ostringstream& out, const std::string& value) {
    out << '"' << json_escape(value) << '"';
}

template <typename T>
void json_number(std::ostringstream& out, T value) {
    if (std::isfinite(static_cast<double>(value))) {
        out << value;
    } else {
        out << "null";
    }
}

const char* json_bool(bool value) {
    return value ? "true" : "false";
}

const char* mode_name(dataset::BenchmarkMode mode) {
    switch (mode) {
        case dataset::BenchmarkMode::ORACLE_CEILING: return "oracle-ceiling";
        case dataset::BenchmarkMode::DEFINED_CLEFT_REDOCK: return "defined-cleft-redock";
        case dataset::BenchmarkMode::AUTONOMOUS: return "autonomous";
        case dataset::BenchmarkMode::UNSET: return "unset";
    }
    return "unset";
}

std::string utc_now() {
    const auto now = std::chrono::system_clock::now();
    const auto value = std::chrono::system_clock::to_time_t(now);
    std::tm tm{};
#ifdef _WIN32
    gmtime_s(&tm, &value);
#else
    gmtime_r(&value, &tm);
#endif
    std::ostringstream out;
    out << std::put_time(&tm, "%Y-%m-%dT%H:%M:%SZ");
    return out.str();
}

bool execution_completed(const dataset::DockingResult& result) {
    return result.docking_completed && result.docking_exit_code == 0 &&
           result.num_poses > 0 && !result.stuck;
}

bool validators_complete(const dataset::DockingResult& result) {
    if (!result.pb_ran || result.tencom_status != "ok" ||
        result.eigen_status != "ok" || result.pose_sha256.empty()) {
        return false;
    }
    return result.rmsd_pose_sha256 == result.pose_sha256 &&
           result.posebusters_pose_sha256 == result.pose_sha256 &&
           result.tencom_pose_sha256 == result.pose_sha256;
}

} // namespace

namespace fleet {

std::string FleetRunner::serialize_chunk_result(
    const ChunkMetadata& metadata,
    const dataset::BenchmarkReport& report,
    const dataset::DockingConfig& config,
    double duration_s) {
    int completed = 0;
    int successful_rmsd = 0;
    int successful_pb = 0;
    int claim_ready = 0;
    int validators_done = 0;

    for (const auto& result : report.results) {
        completed += execution_completed(result) ? 1 : 0;
        successful_rmsd += result.success_rmsd ? 1 : 0;
        successful_pb += result.success_pb ? 1 : 0;
        claim_ready += result.claim_ready ? 1 : 0;
        validators_done += validators_complete(result) ? 1 : 0;
    }

    std::ostringstream out;
    out << std::setprecision(9);
    out << "{\n  \"schema_version\": 1,\n";
    out << "  \"type\": \"flexaidds-fleet-chunk-result\",\n";
    out << "  \"created_at\": "; json_string(out, utc_now()); out << ",\n";
    out << "  \"campaign_id\": "; json_string(out, metadata.campaign_id); out << ",\n";
    out << "  \"chunk_id\": "; json_string(out, metadata.chunk_id); out << ",\n";
    out << "  \"attempt_id\": "; json_string(out, metadata.attempt_id); out << ",\n";
    out << "  \"worker_id\": "; json_string(out, metadata.worker_id); out << ",\n";
    out << "  \"dataset\": "; json_string(out, metadata.dataset); out << ",\n";
    out << "  \"duration_s\": "; json_number(out, duration_s); out << ",\n";

    out << "  \"provenance\": {\n";
    out << "    \"manifest_sha256\": "; json_string(out, metadata.manifest_sha256); out << ",\n";
    out << "    \"runner_path\": "; json_string(out, metadata.runner_path); out << ",\n";
    out << "    \"runner_sha256\": "; json_string(out, metadata.runner_sha256); out << ",\n";
    out << "    \"engine_path\": "; json_string(out, metadata.engine_path); out << ",\n";
    out << "    \"engine_sha256\": "; json_string(out, metadata.engine_sha256); out << ",\n";
    out << "    \"command\": "; json_string(out, metadata.command); out << "\n  },\n";

    out << "  \"protocol\": {\n";
    out << "    \"mode\": "; json_string(out, mode_name(config.mode)); out << ",\n";
    out << "    \"ga_population\": " << config.ga_population << ",\n";
    out << "    \"ga_generations\": " << config.ga_generations << ",\n";
    out << "    \"temperature_k\": "; json_number(out, config.temperature); out << ",\n";
    out << "    \"grid_spacing_a\": "; json_number(out, config.grid_spacing); out << ",\n";
    out << "    \"dataset_workers\": " << config.num_threads << ",\n";
    out << "    \"omp_threads_per_worker\": " << config.omp_threads_per_worker << ",\n";
    out << "    \"job_timeout_s\": " << config.per_job_timeout_s << ",\n";
    out << "    \"clustering\": "; json_string(out, config.clustering_algorithm); out << ",\n";
    out << "    \"gpu_enabled\": " << json_bool(config.use_gpu) << ",\n";
    out << "    \"gpu_backend\": "; json_string(out, config.gpu_backend); out << "\n  },\n";

    out << "  \"summary\": {\n";
    out << "    \"target_count\": " << report.results.size() << ",\n";
    out << "    \"execution_completed\": " << completed << ",\n";
    out << "    \"execution_failed\": " << (static_cast<int>(report.results.size()) - completed) << ",\n";
    out << "    \"success_rmsd\": " << successful_rmsd << ",\n";
    out << "    \"success_pb\": " << successful_pb << ",\n";
    out << "    \"validators_complete\": " << validators_done << ",\n";
    out << "    \"claim_ready\": " << claim_ready << "\n  },\n";

    out << "  \"targets\": [\n";
    for (std::size_t i = 0; i < report.results.size(); ++i) {
        const auto& result = report.results[i];
        out << "    {\n";
        out << "      \"pdb_id\": "; json_string(out, result.pdb_id); out << ",\n";
        out << "      \"execution_completed\": " << json_bool(execution_completed(result)) << ",\n";
        out << "      \"docking_exit_code\": " << result.docking_exit_code << ",\n";
        out << "      \"matrix_md5\": "; json_string(out, result.matrix_md5); out << ",\n";
        out << "      \"num_poses\": " << result.num_poses << ",\n";
        out << "      \"wall_time_s\": "; json_number(out, result.wall_time_s); out << ",\n";
        out << "      \"rmsd_hungarian_a\": "; json_number(out, result.rmsd_hungarian); out << ",\n";
        out << "      \"rmsd_serial_a\": "; json_number(out, result.rmsd_to_crystal); out << ",\n";
        out << "      \"success_rmsd\": " << json_bool(result.success_rmsd) << ",\n";
        out << "      \"posebusters_ran\": " << json_bool(result.pb_ran) << ",\n";
        out << "      \"posebusters_pass\": " << json_bool(result.pb_pass) << ",\n";
        out << "      \"posebusters_backend\": "; json_string(out, result.pb_backend); out << ",\n";
        out << "      \"posebusters_failed_keys\": "; json_string(out, result.pb_failed_keys); out << ",\n";
        out << "      \"success_pb\": " << json_bool(result.success_pb) << ",\n";
        out << "      \"tencom_status\": "; json_string(out, result.tencom_status); out << ",\n";
        out << "      \"eigen_status\": "; json_string(out, result.eigen_status); out << ",\n";
        out << "      \"eigen_n_modes\": " << result.eigen_n_modes << ",\n";
        out << "      \"validators_complete\": " << json_bool(validators_complete(result)) << ",\n";
        out << "      \"protocol_claim_eligible\": " << json_bool(result.protocol_claim_eligible) << ",\n";
        out << "      \"claim_ready\": " << json_bool(result.claim_ready) << ",\n";
        out << "      \"native_pose_seeded\": " << json_bool(result.native_pose_seeded) << ",\n";
        out << "      \"seed_echo\": " << json_bool(result.seed_echo) << ",\n";
        out << "      \"pose_source\": "; json_string(out, result.pose_source); out << ",\n";
        out << "      \"elected_pose_path\": "; json_string(out, result.elected_pose_path); out << ",\n";
        out << "      \"elected_pose_source\": "; json_string(out, result.elected_pose_source); out << ",\n";
        out << "      \"elected_restart\": " << result.elected_restart << ",\n";
        out << "      \"elected_cluster\": " << result.elected_cluster << ",\n";
        out << "      \"elected_cf\": "; json_number(out, result.elected_cf); out << ",\n";
        out << "      \"best_cluster_rmsd_a\": "; json_number(out, result.best_cluster_rmsd); out << ",\n";
        out << "      \"pose_sha256\": "; json_string(out, result.pose_sha256); out << ",\n";
        out << "      \"rmsd_pose_sha256\": "; json_string(out, result.rmsd_pose_sha256); out << ",\n";
        out << "      \"posebusters_pose_sha256\": "; json_string(out, result.posebusters_pose_sha256); out << ",\n";
        out << "      \"posebusters_input_sha256\": "; json_string(out, result.posebusters_input_sha256); out << ",\n";
        out << "      \"tencom_pose_sha256\": "; json_string(out, result.tencom_pose_sha256); out << ",\n";
        out << "      \"shannon_entropy_nats\": "; json_number(out, result.shannon_entropy); out << ",\n";
        out << "      \"elected_h_vib_nats\": "; json_number(out, result.elected_H_vib); out << ",\n";
        out << "      \"thermo_available\": " << json_bool(result.has_thermo) << "\n";
        out << "    }" << (i + 1 < report.results.size() ? "," : "") << "\n";
    }
    out << "  ]\n}\n";
    return out.str();
}

bool FleetRunner::write_chunk_result_atomic(
    const std::string& destination,
    const std::string& contents,
    std::string* error) {
    if (destination.empty() || destination == "-") {
        if (error) *error = "Fleet result requires a file destination";
        return false;
    }

    const fs::path target(destination);
    std::error_code ec;
    if (!target.parent_path().empty()) {
        fs::create_directories(target.parent_path(), ec);
        if (ec) {
            if (error) *error = "cannot create result directory: " + ec.message();
            return false;
        }
    }
    if (fs::exists(target, ec)) {
        if (error) *error = "refusing to overwrite immutable Fleet result";
        return false;
    }

    const fs::path temp = target.string() + ".tmp." +
        std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
    {
        std::ofstream stream(temp, std::ios::binary | std::ios::trunc);
        if (!stream) {
            if (error) *error = "cannot open temporary Fleet result";
            return false;
        }
        stream.write(contents.data(), static_cast<std::streamsize>(contents.size()));
        stream.flush();
        if (!stream) {
            if (error) *error = "cannot write temporary Fleet result";
            stream.close();
            fs::remove(temp, ec);
            return false;
        }
    }

#ifndef _WIN32
    const int fd = ::open(temp.c_str(), O_RDONLY);
    if (fd >= 0) {
        (void)::fsync(fd);
        (void)::close(fd);
    }
    if (::link(temp.c_str(), target.c_str()) != 0) {
        if (error) *error = fs::exists(target) ?
            "refusing to overwrite immutable Fleet result" :
            "cannot publish Fleet result";
        fs::remove(temp, ec);
        return false;
    }
    fs::remove(temp, ec);
    if (!target.parent_path().empty()) {
        const int dir_fd = ::open(target.parent_path().c_str(), O_RDONLY);
        if (dir_fd >= 0) {
            (void)::fsync(dir_fd);
            (void)::close(dir_fd);
        }
    }
#else
    fs::rename(temp, target, ec);
    if (ec) {
        if (error) *error = "cannot publish Fleet result: " + ec.message();
        fs::remove(temp, ec);
        return false;
    }
#endif
    return true;
}

} // namespace fleet
