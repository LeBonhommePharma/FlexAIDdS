// =============================================================================
// DatasetRunnerProvenance.cpp — provenance.json writer for DatasetRunner
//
// Extracted from DatasetRunner.cpp (P1 leaf of the split plan). Behaviour is
// intentionally identical to the pre-split inline run() provenance block.
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
// =============================================================================

#include "DatasetRunnerProvenance.h"
#include "shell_exec.h"

#include <cctype>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

namespace dataset {

std::string provenance_json_escape(const std::string& s) {
    std::string r;
    for (char c : s) {
        if (c == '\\' || c == '"') {
            r += '\\';
            r += c;
        } else if (c == '\n') {
            r += "\\n";
        } else {
            r += c;
        }
    }
    return r;
}

std::string provenance_cmd_token(const std::string& cmd) {
    // Deprecated shell path: prefer provenance_file_md5/sha256 argv helpers.
    (void)cmd;
    return {};
}

std::string provenance_file_md5(const std::string& path) {
    using flexaids::shell_exec::is_safe_exec_path;
    using flexaids::shell_exec::run_argv_first_token;
    if (path.empty() || !fs::exists(path) || !is_safe_exec_path(path)) return "";
    std::string t = run_argv_first_token({"md5", "-q", path});
    if (t.empty()) t = run_argv_first_token({"md5sum", path});
    return t;
}

std::string provenance_file_sha256(const std::string& path) {
    using flexaids::shell_exec::is_safe_exec_path;
    using flexaids::shell_exec::run_argv_first_token;
    if (path.empty() || !fs::exists(path) || !is_safe_exec_path(path)) return "";
    std::string t = run_argv_first_token({"shasum", "-a", "256", path});
    if (t.empty()) t = run_argv_first_token({"sha256sum", path});
    return t;
}


std::string resolve_scoring_matrix_path(const std::string& data_dir,
                                        const std::string& flexaidds_bin) {
    // Same precedence as each dock child: explicit override, immutable data
    // staged beside the binary, then the mutable source-tree WRK fallback.
    if (!data_dir.empty()) {
        std::string cand = data_dir + "/MC_st0r5.2_6.dat";
        if (fs::exists(cand)) return cand;
    }
    if (!flexaidds_bin.empty()) {
        std::string bin_dir = flexaidds_bin;
        auto slash = bin_dir.rfind('/');
        if (slash != std::string::npos) bin_dir = bin_dir.substr(0, slash);
        const std::string staged = bin_dir + "/MC_st0r5.2_6.dat";
        const std::string source = bin_dir + "/../WRK/MC_st0r5.2_6.dat";
        if (fs::exists(staged)) return staged;
        if (fs::exists(source)) return source;
    }
    return "";
}

RunProvenanceFields build_run_provenance(const std::string& dataset_name,
                                         const std::string& matrix_path,
                                         const std::string& binary_path,
                                         const std::string& oracle_site_dir) {
    RunProvenanceFields p;
    p.dataset = dataset_name;
    p.matrix_path = matrix_path;
    p.matrix_md5 = provenance_file_md5(matrix_path);
    p.matrix_sha256 = provenance_file_sha256(matrix_path);
    p.binary_path = binary_path;
    p.binary_sha256 = provenance_file_sha256(binary_path);
    p.git_commit = provenance_cmd_token("git rev-parse HEAD 2>/dev/null");
    p.oracle_site_dir = oracle_site_dir;
    p.oracle_site_dir_set = !oracle_site_dir.empty();
    return p;
}

std::string format_run_provenance_json(const RunProvenanceFields& p) {
    // Exact key order and layout of the former DatasetRunner::run() writer.
    std::string j;
    j.reserve(512);
    j += "{\n";
    j += "  \"dataset\": \"";
    j += provenance_json_escape(p.dataset);
    j += "\",\n";
    j += "  \"matrix_path\": \"";
    j += provenance_json_escape(p.matrix_path);
    j += "\",\n";
    j += "  \"matrix_md5\": \"";
    j += provenance_json_escape(p.matrix_md5);
    j += "\",\n";
    j += "  \"matrix_sha256\": \"";
    j += provenance_json_escape(p.matrix_sha256);
    j += "\",\n";
    j += "  \"binary_path\": \"";
    j += provenance_json_escape(p.binary_path);
    j += "\",\n";
    j += "  \"binary_sha256\": \"";
    j += provenance_json_escape(p.binary_sha256);
    j += "\",\n";
    j += "  \"git_commit\": \"";
    j += provenance_json_escape(p.git_commit);
    j += "\",\n";
    j += "  \"oracle_site_dir\": \"";
    j += provenance_json_escape(p.oracle_site_dir);
    j += "\",\n";
    j += "  \"oracle_site_dir_set\": ";
    j += (p.oracle_site_dir_set ? "true" : "false");
    j += "\n";
    j += "}\n";
    return j;
}

bool write_run_provenance_json(const std::string& output_dir,
                               const RunProvenanceFields& p,
                               bool log) {
    std::error_code mk_ec;
    fs::create_directories(output_dir, mk_ec);
    const std::string prov_path = output_dir + "/provenance.json";
    std::ofstream pj(prov_path);
    if (!pj.is_open()) {
        if (log) {
            std::cerr << "[WARN] Could not write provenance.json to "
                      << prov_path << "\n";
        }
        return false;
    }
    pj << format_run_provenance_json(p);
    pj.close();
    if (log) {
        std::cout << "[DatasetRunner] Wrote provenance: " << prov_path << "\n";
    }
    return true;
}

bool write_dataset_run_provenance(const std::string& output_dir,
                                  const std::string& dataset_name,
                                  const std::string& flexaidds_bin,
                                  const std::string& data_dir,
                                  const std::string& oracle_site_dir) {
    const std::string matrix_path =
        resolve_scoring_matrix_path(data_dir, flexaidds_bin);
    const RunProvenanceFields fields =
        build_run_provenance(dataset_name, matrix_path, flexaidds_bin,
                             oracle_site_dir);
    return write_run_provenance_json(output_dir, fields, /*log=*/true);
}

} // namespace dataset
