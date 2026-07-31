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
#include <map>
#include <set>
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

// The four literal non-scored columns a PoseBusters 0.6.5 redock CSV emits.
// Deliberately an exact-name list, NOT substring matching: if membership were
// decided by `col.find("rmsd")`, a future column such as
// `rmsd_reference_source` would be silently subtracted from the canonical set
// and set equality below would be computed against a set that quietly shrank.
// Substring matching one layer down is still substring matching deciding the
// gate. 31 real header columns minus these four is exactly the canonical 27.
bool is_metadata_column(std::string_view col) {
    return col == "file" || col == "molecule" || col == "position" ||
           col == "rmsd_\u2264_2\u00e5";
}

// THE canonical PoseBusters 0.6.5 redock gate: the full set of 27 scored
// boolean checks, and the single authority for what pb_pass means.
//
// This list is authoritative in BOTH roles:
//   1. schema — the CSV's scored columns must equal this set exactly
//   2. scoring — pb_pass iterates this list, not the CSV's headers
// Consequently a schema change in either direction fails closed and demands a
// deliberate pin bump. A future column cannot silently strengthen the gate; an
// omitted column cannot silently weaken it. Gate membership is never decided
// by substring matching.
//
// Pin bump 2026-07-31: dropped "no_protein_clashes" — upstream 0.6.5 emits no
// such column, so every real bust run failed closed on schema and the campaign
// silently fell back to native_pose_qc. Filled out to the full 27 at the same
// time; a partial pin was itself an arbitrary subset.
//
// Water is deliberately IN. Upstream redock.yml selects both water checks, so
// removing them would produce a custom metric, not "PoseBusters pass" — even
// though water dominates observed failures (see the campaign note). Any
// non-water variant belongs beside pb_pass as a separately named diagnostic,
// never as a mutation of it.
//
// Verified against 34 real 0.6.5 bust CSVs (2 further files empty) under
// flexaidds_benchmark_results/: one single header layout, all 27 present in
// every file, no extras. Keep in CSV emission order so diffs against a real
// header are trivial to eyeball.
const std::vector<std::string>& mandatory_pb_check_columns() {
    static const std::vector<std::string> k = {
        "mol_pred_loaded",
        "mol_true_loaded",
        "mol_cond_loaded",
        "sanitization",
        "inchi_convertible",
        "all_atoms_connected",
        "no_radicals",
        "molecular_formula",
        "molecular_bonds",
        "double_bond_stereochemistry",
        "tetrahedral_chirality",
        "bond_lengths",
        "bond_angles",
        "internal_steric_clash",
        "aromatic_ring_flatness",
        "non-aromatic_ring_non-flatness",
        "double_bond_flatness",
        "internal_energy",
        "protein-ligand_maximum_distance",
        "minimum_distance_to_protein",
        "minimum_distance_to_organic_cofactors",
        "minimum_distance_to_inorganic_cofactors",
        "minimum_distance_to_waters",
        "volume_overlap_with_protein",
        "volume_overlap_with_organic_cofactors",
        "volume_overlap_with_inorganic_cofactors",
        "volume_overlap_with_waters",
    };
    return k;
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
#ifdef FLEXAIDDS_POSEBUSTERS_BIN_DEFAULT
    // Baked in by CMake's find_program(POSEBUSTERS_BIN ...) at configure time.
    if (file_executable(FLEXAIDDS_POSEBUSTERS_BIN_DEFAULT))
        return FLEXAIDDS_POSEBUSTERS_BIN_DEFAULT;
#endif
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

    // Receipts: path, binary hash, argv, exit, raw hash (always filled when possible).
    r.bust_path = bust;
    r.bust_sha256 = sha256_via_openssl(bust);
    {
        std::ostringstream aj;
        for (std::size_t i = 0; i < argv.size(); ++i) {
            if (i) aj << ' ';
            aj << argv[i];
        }
        r.argv_joined = aj.str();
    }
    // Best-effort version probe (non-fatal).
    {
        auto vcap = run_argv_capture({bust, "--version"}, /*discard_stderr=*/false);
        if (vcap.ok && !vcap.stdout_text.empty()) {
            r.bust_version = vcap.stdout_text;
            while (!r.bust_version.empty() &&
                   (r.bust_version.back() == '\n' || r.bust_version.back() == '\r'))
                r.bust_version.pop_back();
        }
    }

    auto cap = run_argv_capture(argv, /*discard_stderr=*/true);
    if (!cap.ok) {
        r.error = "exec(bust) failed";
        r.backend = "error";
        r.exit_status = -1;
        return r;
    }
    const int rc = cap.exit_code;
    const std::string& out = cap.stdout_text;
    // Always preserve raw CSV before any schema return.
    r.raw_csv = out;
    r.exit_status = rc;
    r.ran = true;
    r.backend = "bust_cli";
    if (!sidecar_dir.empty()) {
        std::error_code ec;
        fs::create_directories(sidecar_dir, ec);
        const std::string raw_path =
            (fs::path(sidecar_dir) / (stem + "_bust_raw.csv")).string();
        std::ofstream raw_ofs(raw_path);
        if (raw_ofs) raw_ofs << out;
        r.csv_path = raw_path;
        r.raw_csv_sha256 = sha256_via_openssl(raw_path);
        if (r.raw_csv_sha256.empty() && !out.empty()) {
            // Hash via temp if openssl path failed on empty write.
            r.raw_csv_sha256 = sha256_via_openssl(raw_path);
        }
    }

    // Schema validation (raw_csv already preserved above).
    apply_bust_csv_schema(out, r);
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

void apply_bust_csv_schema(const std::string& csv_body, BustCliResult& r) {
    if (r.raw_csv.empty()) r.raw_csv = csv_body;

    if (csv_body.empty()) {
        r.error = "bust produced empty CSV";
        r.pb_pass = false;
        return;
    }

    std::istringstream iss(csv_body);
    std::string header_line, data_line;
    if (!std::getline(iss, header_line) || !std::getline(iss, data_line)) {
        r.error = "bust CSV parse: need header+data";
        r.pb_pass = false;
        return;
    }
    while (!data_line.empty() &&
           (data_line.back() == '\n' || data_line.back() == '\r'))
        data_line.pop_back();

    auto headers = split_csv_line(header_line);
    auto values  = split_csv_line(data_line);
    // Fail-closed: no pad; column counts must match; no duplicate headers.
    if (values.size() != headers.size()) {
        r.error = "bust CSV schema: header/value column count mismatch (" +
                  std::to_string(headers.size()) + " vs " +
                  std::to_string(values.size()) + ")";
        r.pb_pass = false;
        r.n_checks = 0;
        r.n_pass = 0;
        r.n_fail = 0;
        r.failed_keys = "schema_column_count";
        return;
    }
    {
        std::set<std::string> seen;
        for (const auto& h : headers) {
            if (!seen.insert(h).second) {
                r.error = "bust CSV schema: duplicate header '" + h + "'";
                r.pb_pass = false;
                r.failed_keys = "duplicate_header:" + h;
                return;
            }
        }
    }
    // Version-pinned canonical check-set: enforce SET EQUALITY, not mere
    // presence. Drift in either direction fails closed and demands a
    // deliberate pin bump:
    //   missing canonical column  -> a future CSV cannot silently weaken pb_pass
    //   unexpected scored column  -> a future CSV cannot silently strengthen it
    {
        std::set<std::string> header_set(headers.begin(), headers.end());
        std::string missing;
        for (const auto& need : mandatory_pb_check_columns()) {
            if (header_set.count(need) == 0) {
                if (!missing.empty()) missing += ';';
                missing += need;
            }
        }
        if (!missing.empty()) {
            r.error = "bust CSV schema: missing mandatory check columns: " + missing;
            r.pb_pass = false;
            r.failed_keys = "mandatory_checks_missing:" + missing;
            return;
        }
        const auto& canon = mandatory_pb_check_columns();
        std::set<std::string> canon_set(canon.begin(), canon.end());
        std::string unexpected;
        for (const auto& h : headers) {
            if (is_metadata_column(h)) continue;
            if (canon_set.count(h) == 0) {
                if (!unexpected.empty()) unexpected += ';';
                unexpected += h;
            }
        }
        if (!unexpected.empty()) {
            r.error = "bust CSV schema: unexpected scored columns (pin bump "
                      "required): " + unexpected;
            r.pb_pass = false;
            r.failed_keys = "unexpected_scored_columns:" + unexpected;
            return;
        }
    }

    auto is_nan_token = [](std::string_view v) -> bool {
        if (v.empty()) return false;
        std::string t(v);
        for (char& c : t)
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        return t == "nan" || t == "na" || t == "none" || t == "null" ||
               t == "inf" || t == "+inf" || t == "-inf";
    };

    r.n_checks = 0;
    r.n_pass = 0;
    r.n_fail = 0;
    bool all_ok = true;
    std::string failed;
    // Iterate the CANONICAL list, not the CSV's headers. Set equality was
    // enforced above, so every canonical column is present exactly once and no
    // extra scored column exists — but driving the loop from the pin is what
    // makes that guarantee load-bearing rather than incidental.
    std::map<std::string, std::string> row;
    for (std::size_t i = 0; i < headers.size(); ++i) row[headers[i]] = values[i];
    for (const auto& h : mandatory_pb_check_columns()) {
        const std::string& v = row[h];

        if (v.empty() || is_nan_token(v)) {
            ++r.n_checks;
            ++r.n_fail;
            all_ok = false;
            if (!failed.empty()) failed += ';';
            failed += h + (v.empty() ? ":blank" : ":uncomputed");
            continue;
        }
        if (!parse_truthy(v) && !parse_falsey(v)) {
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
    r.pb_pass = all_ok && r.n_checks > 0 && r.n_fail == 0;
    r.failed_keys = failed;
    if (r.n_checks == 0) {
        r.pb_pass = false;
        if (r.error.empty()) r.error = "bust CSV: no boolean check columns";
        if (r.failed_keys.empty()) r.failed_keys = "no_checks";
    }
}

}  // namespace flexaids::posebust
