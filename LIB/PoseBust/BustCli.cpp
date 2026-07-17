// BustCli.cpp — Authoritative upstream PoseBusters CLI driver
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0

#include "BustCli.h"
#include "shell_exec.h"

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

using flexaids::shell_exec::is_safe_exec_path;
using flexaids::shell_exec::run_argv_capture;
using flexaids::shell_exec::run_argv_first_token;

[[nodiscard]] bool file_executable(const fs::path& p) {
    std::error_code ec;
    return fs::is_regular_file(p, ec) &&
           (fs::status(p, ec).permissions() & fs::perms::owner_exec) !=
               fs::perms::none;
}

// Portable file SHA-256 via openssl (always present on macOS/Linux CI).
// Argv exec — path never reaches a shell.
std::string sha256_via_openssl(const std::string& path) {
    if (!is_safe_exec_path(path)) return {};
    std::string hex = run_argv_first_token(
        {"openssl", "dgst", "-sha256", "-r", path});
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

    if (!is_safe_exec_path(bust) || !is_safe_exec_path(pred_sdf) ||
        !is_safe_exec_path(protein_pdb) ||
        (!crystal_sdf.empty() && !is_safe_exec_path(crystal_sdf))) {
        r.error = "unsafe path for bust argv (NUL/newline/control)";
        r.backend = "error";
        return r;
    }

    std::vector<std::string> argv = {bust, pred_sdf, "-p", protein_pdb};
    if (!crystal_sdf.empty() && fs::is_regular_file(crystal_sdf)) {
        argv.push_back("-l");
        argv.push_back(crystal_sdf);
    }
    argv.push_back("--outfmt");
    argv.push_back("csv");

    auto cap = run_argv_capture(argv, /*discard_stderr=*/true);
    if (!cap.ok) {
        r.error = "exec(bust) failed";
        r.backend = "error";
        return r;
    }
    const int rc = cap.exit_code;
    const std::string& out = cap.stdout_text;
    r.raw_csv = out;
    r.ran = true;
    r.backend = "bust_cli";

    if (out.empty()) {
        r.error = "bust produced empty CSV (exit_status=" + std::to_string(rc) + ")";
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
    // Fail-closed schema (audit P0 #17): no padding missing cells; column
    // count must match. Blank / NaN / non-boolean check columns fail pb_pass.
    if (values.size() != headers.size()) {
        r.error = "bust CSV schema: header/value column count mismatch (" +
                  std::to_string(headers.size()) + " vs " +
                  std::to_string(values.size()) + ")";
        r.pb_pass = false;
        r.n_checks = 0;
        r.n_pass = 0;
        r.n_fail = 0;
        r.failed_keys = "schema_column_count";
        return r;
    }

    auto is_nan_token = [](std::string_view v) -> bool {
        if (v.empty()) return false;
        std::string t(v);
        for (char& c : t) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        return t == "nan" || t == "na" || t == "none" || t == "null" ||
               t == "inf" || t == "+inf" || t == "-inf";
    };

    bool all_ok = true;
    std::string failed;
    for (std::size_t i = 0; i < headers.size(); ++i) {
        const std::string& h = headers[i];
        const std::string& v = values[i];
        if (is_excluded_from_pb_pass(h)) continue;

        // Mandatory boolean check columns: blank, NaN, or unparsed → fail.
        if (v.empty() || is_nan_token(v)) {
            ++r.n_checks;
            ++r.n_fail;
            all_ok = false;
            if (!failed.empty()) failed += ';';
            failed += h + (v.empty() ? ":blank" : ":uncomputed");
            continue;
        }
        if (!parse_truthy(v) && !parse_falsey(v)) {
            // Non-boolean in a non-excluded column is a schema failure.
            ++r.n_checks;
            ++r.n_fail;
            all_ok = false;
            if (!failed.empty()) failed += ';';
            failed += h + ":non_boolean";
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
    // Require a non-empty check set so a metadata-only row cannot pass.
    r.pb_pass = all_ok && r.n_checks > 0 && r.n_fail == 0;
    r.failed_keys = failed;
    if (rc != 0) {
        if (!r.error.empty()) r.error += ";";
        r.error += "bust pclose_status=" + std::to_string(rc);
        r.pb_pass = false;
    }
    if (r.n_checks == 0) {
        r.pb_pass = false;
        if (r.error.empty()) r.error = "bust CSV: no boolean check columns";
        if (r.failed_keys.empty()) r.failed_keys = "no_checks";
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
