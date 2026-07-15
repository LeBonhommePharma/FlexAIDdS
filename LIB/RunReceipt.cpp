// RunReceipt.cpp — RUN_RECEIPT.json builder/writer
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
// SPDX-License-Identifier: Apache-2.0

#include "RunReceipt.h"

#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace flexaids {
namespace {

namespace fs = std::filesystem;

std::string json_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 8);
    for (char c : s) {
        switch (c) {
            case '\\': out += "\\\\"; break;
            case '"':  out += "\\\""; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:   out += c; break;
        }
    }
    return out;
}

}  // namespace

std::string utc_now_iso8601() {
    using clock = std::chrono::system_clock;
    const auto now = clock::now();
    const std::time_t t = clock::to_time_t(now);
    std::tm tm{};
#if defined(_WIN32)
    gmtime_s(&tm, &t);
#else
    gmtime_r(&t, &tm);
#endif
    std::ostringstream o;
    o << std::put_time(&tm, "%Y-%m-%dT%H:%M:%SZ");
    return o.str();
}

std::string build_run_receipt_json(const RunReceiptInput& in) {
    std::ostringstream o;
    o.setf(std::ios::fixed);
    o.precision(6);
    o << "{\n";
    o << "  \"schema_version\": " << kRunReceiptSchemaVersion << ",\n";
    o << "  \"run_id\": \"" << json_escape(in.run_id) << "\",\n";
    o << "  \"started_utc\": \"" << json_escape(in.started_utc) << "\",\n";
    o << "  \"output\": \"" << json_escape(in.output) << "\",\n";
    o << "  \"dataset\": \"" << json_escape(in.dataset) << "\",\n";
    o << "  \"mode\": \"" << json_escape(in.mode) << "\",\n";
    o << "  \"temperature_K\": " << in.temperature_K << ",\n";
    o << "  \"pop\": " << in.pop << ",\n";
    o << "  \"gen\": " << in.gen << ",\n";
    o << "  \"restarts\": " << in.restarts << ",\n";
    o << "  \"seed_base\": " << in.seed_base << ",\n";
    o << "  \"seed_elitism\": " << (in.seed_elitism ? 1 : 0) << ",\n";
    o << "  \"matrix_path\": \"" << json_escape(in.matrix_path) << "\",\n";
    o << "  \"matrix_md5\": \"" << json_escape(in.matrix_md5) << "\",\n";
    o << "  \"matrix_sha256\": \"" << json_escape(in.matrix_sha256) << "\",\n";
    o << "  \"binary_path\": \"" << json_escape(in.binary_path) << "\",\n";
    o << "  \"binary_sha256\": \"" << json_escape(in.binary_sha256) << "\",\n";
    o << "  \"runner_path\": \"" << json_escape(in.runner_path) << "\",\n";
    o << "  \"runner_sha256\": \"" << json_escape(in.runner_sha256) << "\",\n";
    o << "  \"git_commit\": \"" << json_escape(in.git_commit) << "\",\n";
    o << "  \"oracle_site_dir\": \"" << json_escape(in.oracle_site_dir) << "\",\n";
    o << "  \"oracle_site_dir_set\": " << (in.oracle_site_dir_set ? "true" : "false") << ",\n";
    // Embed ProtocolConfig snapshot as a nested object (already JSON object text).
    o << "  \"protocol_config\": " << in.protocol.to_json() << "\n";
    o << "}";
    return o.str();
}

bool write_run_receipt(const std::string& output_dir,
                       const RunReceiptInput& in,
                       bool also_write_provenance_json) {
    std::error_code ec;
    fs::create_directories(output_dir, ec);

    const std::string body = build_run_receipt_json(in);
    const std::string receipt_path = output_dir + "/RUN_RECEIPT.json";
    {
        std::ofstream out(receipt_path);
        if (!out.is_open()) return false;
        out << body << "\n";
    }

    if (also_write_provenance_json) {
        // Legacy slim provenance.json — keep keys that older tools expect.
        const std::string prov_path = output_dir + "/provenance.json";
        std::ofstream pj(prov_path);
        if (pj.is_open()) {
            pj << "{\n"
               << "  \"dataset\": \"" << json_escape(in.dataset) << "\",\n"
               << "  \"matrix_path\": \"" << json_escape(in.matrix_path) << "\",\n"
               << "  \"matrix_md5\": \"" << json_escape(in.matrix_md5) << "\",\n"
               << "  \"matrix_sha256\": \"" << json_escape(in.matrix_sha256) << "\",\n"
               << "  \"binary_path\": \"" << json_escape(in.binary_path) << "\",\n"
               << "  \"binary_sha256\": \"" << json_escape(in.binary_sha256) << "\",\n"
               << "  \"git_commit\": \"" << json_escape(in.git_commit) << "\",\n"
               << "  \"oracle_site_dir\": \"" << json_escape(in.oracle_site_dir) << "\",\n"
               << "  \"oracle_site_dir_set\": "
               << (in.oracle_site_dir_set ? "true" : "false") << ",\n"
               << "  \"protocol_config\": " << in.protocol.to_json() << "\n"
               << "}\n";
        }
    }
    return true;
}

}  // namespace flexaids
