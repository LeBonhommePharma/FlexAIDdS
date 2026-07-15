// BustCli.cpp — Authoritative upstream PoseBusters CLI driver
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0

#include "BustCli.h"

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace flexaids::posebust {
namespace {

namespace fs = std::filesystem;

[[nodiscard]] bool file_executable(const fs::path& p) {
    std::error_code ec;
    return fs::is_regular_file(p, ec) &&
           (fs::status(p, ec).permissions() & fs::perms::owner_exec) !=
               fs::perms::none;
}

std::string shell_quote(const std::string& s) {
    std::string o = "'";
    for (char c : s) {
        if (c == '\'') o += "'\\''";
        else o += c;
    }
    o += "'";
    return o;
}

std::string shell_quote(const std::string& s);

// Portable file SHA-256 via openssl (always present on macOS/Linux CI).
std::string sha256_via_openssl(const std::string& path) {
    const std::string cmd =
        "openssl dgst -sha256 -r " + shell_quote(path) + " 2>/dev/null";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) return {};
    char buf[128];
    std::string out;
    while (fgets(buf, sizeof(buf), pipe)) out += buf;
    pclose(pipe);
    // format: "<hex> *path" or "<hex> path"
    std::istringstream iss(out);
    std::string hex;
    iss >> hex;
    if (hex.size() != 64) return {};
    for (char& c : hex) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return hex;
}

// Columns that are metadata / optional, not part of official pb_pass dock gate.
bool is_excluded_from_pb_pass(std::string_view col) {
    if (col == "file" || col == "molecule" || col == "position") return true;
    // RMSD is success_rmsd, not pb_pass
    if (col.find("rmsd") != std::string_view::npos) return true;
    return false;
}

bool parse_truthy(std::string_view v) {
    return v == "True" || v == "true" || v == "1" || v == "TRUE";
}
bool parse_falsey(std::string_view v) {
    return v == "False" || v == "false" || v == "0" || v == "FALSE";
}

// RFC4180-enough CSV splitting for PoseBusters output: handles quoted fields
// and doubled quotes while keeping the implementation dependency-free.
std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> cols;
    std::string cur;
    bool in_quotes = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (c == '"') {
            if (in_quotes && i + 1 < line.size() && line[i + 1] == '"') {
                cur.push_back('"');
                ++i;
            } else {
                in_quotes = !in_quotes;
            }
        } else if (c == ',' && !in_quotes) {
            cols.push_back(cur);
            cur.clear();
        } else if (c != '\r') {
            cur.push_back(c);
        }
    }
    cols.push_back(cur);
    return cols;
}

}  // namespace

std::string resolve_bust_binary() {
    if (const char* e = std::getenv("FLEXAIDDS_POSEBUSTERS_BIN")) {
        if (e[0] && file_executable(e)) return e;
    }
    // PATH lookup
    if (const char* path = std::getenv("PATH")) {
        std::string p = path;
        std::size_t start = 0;
        while (start <= p.size()) {
            auto end = p.find(':', start);
            std::string dir = p.substr(start, end == std::string::npos ? std::string::npos : end - start);
            fs::path cand = fs::path(dir) / "bust";
            if (file_executable(cand)) return cand.string();
            if (end == std::string::npos) break;
            start = end + 1;
        }
    }
    if (const char* root = std::getenv("FLEXAIDDS_ROOT")) {
        fs::path cand = fs::path(root) / ".venv-posebusters" / "bin" / "bust";
        if (file_executable(cand)) return cand.string();
    }
    // Common relative to cwd / repo
    for (const char* rel : {".venv-posebusters/bin/bust",
                            "../.venv-posebusters/bin/bust"}) {
        if (file_executable(rel)) return fs::absolute(rel).string();
    }
    return {};
}

std::string sha256_file(const std::string& path) {
    return sha256_via_openssl(path);
}

bool copy_file_atomic(const std::string& src, const std::string& dst, std::string* err) {
    std::error_code ec;
    fs::path dpath(dst);
    if (dpath.has_parent_path()) fs::create_directories(dpath.parent_path(), ec);
    fs::path tmp = dpath;
    tmp += ".tmp";
    fs::copy_file(src, tmp, fs::copy_options::overwrite_existing, ec);
    if (ec) {
        if (err) *err = "copy_file: " + ec.message();
        return false;
    }
    fs::rename(tmp, dpath, ec);
    if (ec) {
        if (err) *err = "rename: " + ec.message();
        return false;
    }
    return true;
}

BustCliResult run_upstream_bust(const std::string& pred_sdf,
                                const std::string& protein_pdb,
                                const std::string& crystal_sdf,
                                const std::string& sidecar_dir,
                                const std::string& stem) {
    BustCliResult r;
    const std::string bust = resolve_bust_binary();
    if (bust.empty()) {
        r.error = "bust binary not found (set FLEXAIDDS_POSEBUSTERS_BIN)";
        r.backend = "bust_cli_missing";
        return r;
    }
    if (!fs::is_regular_file(pred_sdf) || !fs::is_regular_file(protein_pdb)) {
        r.error = "pred_sdf or protein_pdb missing";
        r.backend = "error";
        return r;
    }

    std::ostringstream cmd;
    cmd << shell_quote(bust) << " " << shell_quote(pred_sdf)
        << " -p " << shell_quote(protein_pdb);
    if (!crystal_sdf.empty() && fs::is_regular_file(crystal_sdf)) {
        cmd << " -l " << shell_quote(crystal_sdf);
    }
    cmd << " --outfmt csv 2>/dev/null";

    FILE* pipe = popen(cmd.str().c_str(), "r");
    if (!pipe) {
        r.error = "popen(bust) failed";
        r.backend = "error";
        return r;
    }
    char buf[4096];
    std::string out;
    while (fgets(buf, sizeof(buf), pipe)) out += buf;
    const int rc = pclose(pipe);
    r.raw_csv = out;
    r.ran = true;
    r.backend = "bust_cli";

    if (out.empty()) {
        r.error = "bust produced empty CSV (pclose_status=" + std::to_string(rc) + ")";
        r.pb_pass = false;
        return r;
    }

    std::istringstream iss(out);
    std::string header_line, data_line;
    if (!std::getline(iss, header_line) || !std::getline(iss, data_line)) {
        // single line?
        r.error = "bust CSV parse: need header+data";
        r.pb_pass = false;
        return r;
    }
    // strip trailing blank
    while (!data_line.empty() && (data_line.back() == '\n' || data_line.back() == '\r'))
        data_line.pop_back();

    auto headers = split_csv_line(header_line);
    auto values  = split_csv_line(data_line);
    if (values.size() < headers.size()) {
        // pad
        values.resize(headers.size());
    }

    bool all_ok = true;
    std::string failed;
    for (std::size_t i = 0; i < headers.size(); ++i) {
        const std::string& h = headers[i];
        const std::string& v = (i < values.size()) ? values[i] : "";
        if (is_excluded_from_pb_pass(h)) continue;
        if (!parse_truthy(v) && !parse_falsey(v)) {
            // non-boolean columns skipped
            continue;
        }
        ++r.n_checks;
        if (parse_truthy(v)) {
            ++r.n_pass;
        } else {
            ++r.n_fail;
            all_ok = false;
            if (!failed.empty()) failed += ';';
            failed += h;
        }
    }
    r.pb_pass = all_ok && r.n_checks > 0;
    r.failed_keys = failed;
    if (rc != 0) {
        if (!r.error.empty()) r.error += ";";
        r.error += "bust pclose_status=" + std::to_string(rc);
        r.pb_pass = false;
    }

    if (!sidecar_dir.empty()) {
        std::error_code ec;
        fs::create_directories(sidecar_dir, ec);
        r.csv_path = (fs::path(sidecar_dir) / (stem + "_bust.csv")).string();
        std::ofstream ofs(r.csv_path);
        if (ofs) ofs << out;
    }
    return r;
}

}  // namespace flexaids::posebust
