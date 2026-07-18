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

constexpr std::string_view kPoseBustersVersion = "0.6.5";
constexpr std::string_view kPoseBustersVersionOutput = "bust 0.6.5";
constexpr std::string_view kPoseBustersConfig = "redock";
constexpr std::string_view kPoseBustersRedockConfigSha256 =
    "4d551d898ff29a404f16e02ad5a7a2d4235e6b7b14e9a3e27f7c66b4d16b2da9";
constexpr std::string_view kPoseBustersSchema =
    "posebusters-0.6.5-redock-csv-v1";

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

// Exact boolean gate emitted by PoseBusters 0.6.5 config/redock.yml. The
// upstream config has no `no_protein_clashes` output: protein clashes are
// represented by minimum_distance_to_protein and volume_overlap_with_protein.
const std::vector<std::string>& required_pb_check_columns() {
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

// Header captured from the bundled upstream 0.6.5 1of6 redock example. Set
// equality is enforced, while column order remains parser-independent.
const std::vector<std::string>& expected_pb_csv_columns() {
    static const std::vector<std::string> k = {
        "file",
        "molecule",
        "position",
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
        "rmsd_≤_2å",
    };
    return k;
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

std::string identity_path(const std::string& path) {
    std::error_code ec;
    fs::path p = fs::absolute(fs::path(path), ec);
    if (ec) return fs::path(path).lexically_normal().string();
    return p.lexically_normal().string();
}

std::string launcher_python(const std::string& bust) {
    std::ifstream in(bust);
    std::string first;
    if (in && std::getline(in, first) && first.rfind("#!", 0) == 0) {
        std::string candidate = first.substr(2);
        while (!candidate.empty() &&
               (candidate.back() == '\r' || candidate.back() == ' ' ||
                candidate.back() == '\t')) {
            candidate.pop_back();
        }
        const auto first_non_space = candidate.find_first_not_of(" \t");
        if (first_non_space != std::string::npos)
            candidate.erase(0, first_non_space);
        if (candidate.find_first_of(" \t") == std::string::npos &&
            file_executable(candidate)) {
            return candidate;
        }
    }
    for (const char* name : {"python", "python3"}) {
        fs::path candidate = fs::path(bust).parent_path() / name;
        if (file_executable(candidate)) return candidate.string();
    }
    return {};
}

bool discover_posebusters_installation(const std::string& bust,
                                        BustCliResult& r) {
    const std::string python = launcher_python(bust);
    if (python.empty()) {
        r.error = "PoseBusters identity: cannot resolve launcher Python";
        r.failed_keys = "validator_identity:python";
        return false;
    }
    static const std::string kProbe =
        "import importlib.metadata as m,pathlib,posebusters;"
        "d=m.distribution('posebusters');"
        "f=next((x for x in (d.files or []) if str(x).endswith('.dist-info/RECORD')),None);"
        "print(d.version);"
        "print(d.locate_file(f) if f else '');"
        "print(pathlib.Path(posebusters.__file__).resolve().parent/'config'/'redock.yml')";
    auto cap = run_argv_capture({python, "-c", kProbe}, /*discard_stderr=*/true);
    if (!cap.ok || cap.exit_code != 0) {
        r.error = "PoseBusters identity: package metadata probe failed";
        r.failed_keys = "validator_identity:package_probe";
        return false;
    }
    std::istringstream lines(cap.stdout_text);
    std::string version, record_path, config_path;
    if (!std::getline(lines, version) || !std::getline(lines, record_path) ||
        !std::getline(lines, config_path)) {
        r.error = "PoseBusters identity: malformed package metadata probe";
        r.failed_keys = "validator_identity:package_probe_output";
        return false;
    }
    for (std::string* s : {&version, &record_path, &config_path}) {
        while (!s->empty() && s->back() == '\r') s->pop_back();
    }
    r.package_version = version;
    r.package_record_path = identity_path(record_path);
    r.package_record_sha256 = sha256_via_openssl(r.package_record_path);
    r.config_path = identity_path(config_path);
    r.config_sha256 = sha256_via_openssl(r.config_path);
    if (r.package_record_sha256.empty()) {
        r.error = "PoseBusters identity: package RECORD missing or unhashable";
        r.failed_keys = "validator_identity:package_record";
        return false;
    }
    if (r.config_sha256.empty()) {
        r.error = "PoseBusters identity: redock config missing or unhashable";
        r.failed_keys = "validator_identity:redock_config";
        return false;
    }
    return true;
}

std::string join_argv(const std::vector<std::string>& argv) {
    std::ostringstream out;
    for (std::size_t i = 0; i < argv.size(); ++i) {
        if (i) out << ' ';
        out << argv[i];
    }
    return out.str();
}

std::string json_escape(std::string_view value) {
    std::string out;
    out.reserve(value.size() + 8);
    static const char* hex = "0123456789abcdef";
    for (unsigned char c : value) {
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\b': out += "\\b"; break;
            case '\f': out += "\\f"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if (c < 0x20) {
                    out += "\\u00";
                    out.push_back(hex[(c >> 4) & 0x0f]);
                    out.push_back(hex[c & 0x0f]);
                } else {
                    out.push_back(static_cast<char>(c));
                }
        }
    }
    return out;
}

void append_error(BustCliResult& r, const std::string& message,
                  const std::string& failed_key) {
    if (!r.error.empty()) r.error += ';';
    r.error += message;
    if (!r.failed_keys.empty()) r.failed_keys += ';';
    r.failed_keys += failed_key;
    r.pb_pass = false;
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
    r.bust_path = identity_path(bust);
    r.bust_sha256 = sha256_via_openssl(r.bust_path);
    {
        auto vcap = run_argv_capture({bust, "--version"}, /*discard_stderr=*/false);
        if (vcap.ok && vcap.exit_code == 0) {
            r.bust_version = vcap.stdout_text;
            while (!r.bust_version.empty() &&
                   (r.bust_version.back() == '\n' || r.bust_version.back() == '\r'))
                r.bust_version.pop_back();
        }
    }
    if (!discover_posebusters_installation(bust, r)) {
        r.backend = "error";
        return r;
    }
    if (r.package_version != kPoseBustersVersion ||
        r.bust_version != kPoseBustersVersionOutput) {
        r.error = "PoseBusters identity: expected posebusters==0.6.5 and "
                  "'bust 0.6.5', got package='" + r.package_version +
                  "' launcher='" + r.bust_version + "'";
        r.failed_keys = "unsupported_posebusters_version";
        r.backend = "error";
        return r;
    }
    if (r.config_sha256 != kPoseBustersRedockConfigSha256) {
        r.error = "PoseBusters identity: redock.yml SHA-256 mismatch; expected " +
                  std::string(kPoseBustersRedockConfigSha256) + " got " +
                  r.config_sha256;
        r.failed_keys = "unsupported_posebusters_config";
        r.backend = "error";
        return r;
    }
    if (r.bust_sha256.empty()) {
        r.error = "PoseBusters identity: launcher unhashable";
        r.failed_keys = "validator_identity:launcher";
        r.backend = "error";
        return r;
    }

    if (!fs::is_regular_file(pred_sdf) || !fs::is_regular_file(protein_pdb) ||
        !fs::is_regular_file(crystal_sdf)) {
        r.error = "pred_sdf, protein_pdb, or crystal_sdf missing for redock";
        r.failed_keys = "validator_identity:input_missing";
        r.backend = "error";
        return r;
    }

    if (!is_safe_exec_path(bust) || !is_safe_exec_path(pred_sdf) ||
        !is_safe_exec_path(protein_pdb) ||
        !is_safe_exec_path(crystal_sdf) || !is_safe_exec_path(r.config_path)) {
        r.error = "unsafe path for bust argv (NUL/newline/control)";
        r.failed_keys = "validator_identity:unsafe_path";
        r.backend = "error";
        return r;
    }

    r.pred_sdf_path = identity_path(pred_sdf);
    r.pred_sdf_sha256 = sha256_via_openssl(r.pred_sdf_path);
    r.protein_pdb_path = identity_path(protein_pdb);
    r.protein_pdb_sha256 = sha256_via_openssl(r.protein_pdb_path);
    r.crystal_sdf_path = identity_path(crystal_sdf);
    r.crystal_sdf_sha256 = sha256_via_openssl(r.crystal_sdf_path);
    if (r.pred_sdf_sha256.empty() || r.protein_pdb_sha256.empty() ||
        r.crystal_sdf_sha256.empty()) {
        r.error = "PoseBusters identity: one or more redock inputs unhashable";
        r.failed_keys = "validator_identity:input_hash";
        r.backend = "error";
        return r;
    }

    r.argv = {bust, pred_sdf, "-p", protein_pdb, "-l", crystal_sdf,
              "--config", r.config_path, "--outfmt", "csv",
              "--max-workers", "0"};
    r.argv_joined = join_argv(r.argv);

    auto cap = run_argv_capture(r.argv, /*discard_stderr=*/true);
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
        {
            std::ofstream raw_ofs(raw_path);
            if (raw_ofs) raw_ofs << out;
        }
        r.raw_csv_path = identity_path(raw_path);
        r.raw_csv_sha256 = sha256_via_openssl(r.raw_csv_path);
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
        {
            std::ofstream ofs(r.csv_path);
            if (ofs) ofs << out;
        }
        r.csv_path = identity_path(r.csv_path);
        r.csv_sha256 = sha256_via_openssl(r.csv_path);
        if (r.raw_csv_sha256.empty() || r.csv_sha256.empty()) {
            append_error(r, "PoseBusters output sidecar missing or unhashable",
                         "validator_identity:output_hash");
        }
    }
    return r;
}

void apply_bust_csv_schema(const std::string& csv_body, BustCliResult& r) {
    if (r.raw_csv.empty()) r.raw_csv = csv_body;
    r.schema_id = std::string(kPoseBustersSchema);
    r.package_name = "posebusters";
    r.config_name = std::string(kPoseBustersConfig);

    if (r.package_version != kPoseBustersVersion) {
        r.error = "bust CSV schema requires PoseBusters 0.6.5, got '" +
                  r.package_version + "'";
        r.failed_keys = "unsupported_posebusters_version";
        r.pb_pass = false;
        return;
    }

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
    std::string extra_line;
    while (std::getline(iss, extra_line)) {
        while (!extra_line.empty() && extra_line.back() == '\r')
            extra_line.pop_back();
        if (!extra_line.empty()) {
            r.error = "bust CSV schema: expected exactly one data row";
            r.pb_pass = false;
            r.failed_keys = "schema_row_count";
            return;
        }
    }

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
    std::set<std::string> seen;
    std::map<std::string, std::size_t> column_index;
    {
        for (std::size_t i = 0; i < headers.size(); ++i) {
            const std::string& h = headers[i];
            if (h.empty()) {
                r.error = "bust CSV schema: blank header";
                r.pb_pass = false;
                r.failed_keys = "blank_header";
                return;
            }
            if (!seen.insert(h).second) {
                r.error = "bust CSV schema: duplicate header '" + h + "'";
                r.pb_pass = false;
                r.failed_keys = "duplicate_header:" + h;
                return;
            }
            column_index.emplace(h, i);
        }
    }

    // Version-pinned 0.6.5 redock schema. Required boolean checks receive the
    // explicit failure key; metadata/RMSD and unexpected columns are schema
    // failures but do not enter pb_pass check counts.
    {
        std::string missing_checks;
        for (const auto& need : required_pb_check_columns()) {
            if (seen.count(need) == 0) {
                if (!missing_checks.empty()) missing_checks += ';';
                missing_checks += need;
            }
        }
        if (!missing_checks.empty()) {
            r.error = "bust CSV schema: missing required 0.6.5 redock checks: " +
                      missing_checks;
            r.pb_pass = false;
            r.failed_keys = "required_checks_missing:" + missing_checks;
            return;
        }

        const std::set<std::string> expected(expected_pb_csv_columns().begin(),
                                             expected_pb_csv_columns().end());
        std::string missing_schema;
        for (const auto& need : expected) {
            if (seen.count(need) == 0) {
                if (!missing_schema.empty()) missing_schema += ';';
                missing_schema += need;
            }
        }
        if (!missing_schema.empty()) {
            r.error = "bust CSV schema: missing 0.6.5 columns: " + missing_schema;
            r.pb_pass = false;
            r.failed_keys = "schema_columns_missing:" + missing_schema;
            return;
        }
        std::string unexpected;
        for (const auto& got : seen) {
            if (expected.count(got) == 0) {
                if (!unexpected.empty()) unexpected += ';';
                unexpected += got;
            }
        }
        if (!unexpected.empty()) {
            r.error = "bust CSV schema: unexpected columns for 0.6.5: " + unexpected;
            r.pb_pass = false;
            r.failed_keys = "schema_columns_unexpected:" + unexpected;
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
    for (const std::string& h : required_pb_check_columns()) {
        const std::string& v = values[column_index.at(h)];

        if (v.empty() || is_nan_token(v)) {
            ++r.n_checks;
            ++r.n_fail;
            all_ok = false;
            if (!failed.empty()) failed += ';';
            failed += h + (v.empty() ? ":blank" : ":uncomputed");
            continue;
        }
        if (v != "True" && v != "False") {
            ++r.n_checks;
            ++r.n_fail;
            all_ok = false;
            if (!failed.empty()) failed += ';';
            failed += h + ":non_boolean";
            continue;
        }
        ++r.n_checks;
        if (v == "True") {
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

std::string format_bust_receipt_json(const BustCliResult& r) {
    auto q = [](std::string_view s) {
        return std::string("\"") + json_escape(s) + "\"";
    };
    std::ostringstream out;
    out << "{\n"
        << "  \"schema\": {\"id\": " << q(r.schema_id)
        << ", \"required_check_count\": " << required_pb_check_columns().size()
        << "},\n"
        << "  \"package\": {\"name\": " << q(r.package_name)
        << ", \"version\": " << q(r.package_version)
        << ", \"record_path\": " << q(r.package_record_path)
        << ", \"record_sha256\": " << q(r.package_record_sha256)
        << ", \"launcher_path\": " << q(r.bust_path)
        << ", \"launcher_sha256\": " << q(r.bust_sha256)
        << ", \"launcher_version_output\": " << q(r.bust_version) << "},\n"
        << "  \"config\": {\"name\": " << q(r.config_name)
        << ", \"path\": " << q(r.config_path)
        << ", \"sha256\": " << q(r.config_sha256) << "},\n"
        << "  \"command\": {\"argv\": [";
    for (std::size_t i = 0; i < r.argv.size(); ++i) {
        if (i) out << ", ";
        out << q(r.argv[i]);
    }
    out << "], \"exit_status\": " << r.exit_status << "},\n"
        << "  \"inputs\": {\n"
        << "    \"predicted_ligand\": {\"path\": " << q(r.pred_sdf_path)
        << ", \"sha256\": " << q(r.pred_sdf_sha256) << "},\n"
        << "    \"protein\": {\"path\": " << q(r.protein_pdb_path)
        << ", \"sha256\": " << q(r.protein_pdb_sha256) << "},\n"
        << "    \"crystal_ligand\": {\"path\": " << q(r.crystal_sdf_path)
        << ", \"sha256\": " << q(r.crystal_sdf_sha256) << "}\n"
        << "  },\n"
        << "  \"outputs\": {\n"
        << "    \"raw_csv\": {\"path\": " << q(r.raw_csv_path)
        << ", \"sha256\": " << q(r.raw_csv_sha256) << "},\n"
        << "    \"validated_csv\": {\"path\": " << q(r.csv_path)
        << ", \"sha256\": " << q(r.csv_sha256) << "}\n"
        << "  },\n"
        << "  \"result\": {\"backend\": " << q(r.backend)
        << ", \"ran\": " << (r.ran ? "true" : "false")
        << ", \"pb_pass\": " << (r.pb_pass ? "true" : "false")
        << ", \"failed_keys\": " << q(r.failed_keys)
        << ", \"error\": " << q(r.error) << "}\n"
        << "}\n";
    return out.str();
}

}  // namespace flexaids::posebust
