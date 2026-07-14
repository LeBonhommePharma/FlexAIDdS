// Engine.cpp — Native PoseBust evaluation orchestration (clean-room)
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
//
// Implements the locked API in Engine.h. PoseBusters-compatible check *keys*
// only; no posebusters/RDKit source.

#include "Engine.h"

#include "ChecksChemistry.h"
#include "ChecksGeometry.h"
#include "ChecksProtein.h"
#include "Loaders.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace flexaids::posebust {
namespace {

namespace fs = std::filesystem;

// ─── string helpers ──────────────────────────────────────────────────────────

[[nodiscard]] std::string json_escape(std::string_view s) {
    std::string out;
    out.reserve(s.size() + 8);
    for (unsigned char c : s) {
        switch (c) {
            case '"':  out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\b': out += "\\b"; break;
            case '\f': out += "\\f"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if (c < 0x20) {
                    char buf[8];
                    std::snprintf(buf, sizeof(buf), "\\u%04x",
                                  static_cast<unsigned>(c));
                    out += buf;
                } else {
                    out.push_back(static_cast<char>(c));
                }
                break;
        }
    }
    return out;
}

[[nodiscard]] bool iequals(std::string_view a, std::string_view b) noexcept {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const auto ca = static_cast<unsigned char>(a[i]);
        const auto cb = static_cast<unsigned char>(b[i]);
        if (std::tolower(ca) != std::tolower(cb)) return false;
    }
    return true;
}

[[nodiscard]] bool parse_tagged_float(const std::string& detail,
                                      std::string_view tag,
                                      float& out) {
    const auto pos = detail.find(tag);
    if (pos == std::string::npos) return false;
    const char* p = detail.c_str() + pos + tag.size();
    char* end = nullptr;
    const float v = std::strtof(p, &end);
    if (end == p || !std::isfinite(v)) return false;
    out = v;
    return true;
}

void fill_continuous_summaries(PoseBustReport& report) {
    for (const CheckItem& c : report.checks) {
        if (c.key == "no_clashes") {
            if (std::isfinite(c.metric)) {
                report.min_lig_prot_dist = c.metric;
            } else {
                float v = std::numeric_limits<float>::quiet_NaN();
                if (parse_tagged_float(c.detail, "min_dist=", v)) {
                    report.min_lig_prot_dist = v;
                }
            }
        } else if (c.key == "no_volume_clash") {
            if (std::isfinite(c.metric)) {
                report.volume_overlap = c.metric;
            } else {
                float v = std::numeric_limits<float>::quiet_NaN();
                if (parse_tagged_float(c.detail, "overlap_fraction=", v)) {
                    report.volume_overlap = v;
                }
            }
        }
    }
}

// ─── protein crop around ligand heavy COM ────────────────────────────────────

[[nodiscard]] bool is_heavy(const Atom& a) noexcept {
    if (a.is_h) return false;
    if (a.atomic_num == 1) return false;
    // Treat Z==0 with element H as hydrogen; otherwise keep (unknown metal etc.).
    if (a.atomic_num <= 0) {
        if (!a.element.empty()) {
            const char c0 = static_cast<char>(
                std::toupper(static_cast<unsigned char>(a.element[0])));
            if (c0 == 'H' &&
                (a.element.size() == 1 ||
                 !std::isalpha(static_cast<unsigned char>(a.element[1])))) {
                return false;
            }
        }
    }
    return true;
}

[[nodiscard]] bool ligand_heavy_com(const Molecule& ligand, Vec3& com_out) {
    double sx = 0.0, sy = 0.0, sz = 0.0;
    int n = 0;
    for (const Atom& a : ligand.atoms) {
        if (!is_heavy(a)) continue;
        sx += static_cast<double>(a.x);
        sy += static_cast<double>(a.y);
        sz += static_cast<double>(a.z);
        ++n;
    }
    if (n == 0) return false;
    const double inv = 1.0 / static_cast<double>(n);
    com_out = Vec3{static_cast<float>(sx * inv),
                   static_cast<float>(sy * inv),
                   static_cast<float>(sz * inv)};
    return true;
}

/// Keep protein heavy atoms within crop_A of ligand heavy COM; rebuild bonds.
[[nodiscard]] Molecule crop_protein_near_ligand(const Molecule& protein,
                                                const Molecule& ligand,
                                                float crop_A) {
    Molecule out;
    out.name = protein.name;

    if (protein.empty() || crop_A <= 0.f) {
        return out;
    }

    Vec3 com{};
    if (!ligand_heavy_com(ligand, com)) {
        // No ligand heavy atoms → empty crop (cannot define pocket centre).
        return out;
    }

    const float r2 = crop_A * crop_A;
    std::vector<int> old_to_new(protein.atoms.size(), -1);
    out.atoms.reserve(protein.atoms.size());

    for (std::size_t i = 0; i < protein.atoms.size(); ++i) {
        const Atom& a = protein.atoms[i];
        if (!is_heavy(a)) continue;
        const float d2 = dist2(a.pos(), com);
        if (d2 > r2) continue;
        old_to_new[i] = static_cast<int>(out.atoms.size());
        out.atoms.push_back(a);
    }

    out.bonds.reserve(protein.bonds.size());
    for (const Bond& b : protein.bonds) {
        if (b.a < 0 || b.b < 0) continue;
        if (static_cast<std::size_t>(b.a) >= old_to_new.size() ||
            static_cast<std::size_t>(b.b) >= old_to_new.size())
            continue;
        const int na = old_to_new[static_cast<std::size_t>(b.a)];
        const int nb = old_to_new[static_cast<std::size_t>(b.b)];
        if (na < 0 || nb < 0) continue;
        out.bonds.push_back(Bond{na, nb, b.order});
    }
    out.build_adjacency();
    return out;
}

void write_json_number_or_null(std::ostream& os, float v) {
    if (std::isfinite(v)) {
        os << v;
    } else {
        os << "null";
    }
}

[[nodiscard]] std::string sidecar_stem(const EvaluateOptions& opt) {
    if (!opt.pdb_id.empty()) return opt.pdb_id;
    return "pose";
}

bool write_sidecar(const Molecule& ligand,
                   const PoseBustReport& report,
                   const EvaluateOptions& opt,
                   std::string* err) {
    std::error_code ec;
    fs::create_directories(opt.sidecar_dir, ec);
    if (ec) {
        if (err) *err = "sidecar: cannot create directory '" + opt.sidecar_dir +
                        "': " + ec.message();
        return false;
    }

    const std::string stem = sidecar_stem(opt);
    const fs::path base(opt.sidecar_dir);
    const fs::path sdf_path  = base / (stem + "_ligand.sdf");
    const fs::path json_path = base / (stem + "_posebust.json");

    if (!write_sdf(ligand, sdf_path.string(), err)) return false;
    if (!write_report_json(report, json_path.string(), err)) return false;
    return true;
}

}  // namespace

// ─── public API ──────────────────────────────────────────────────────────────

PoseBustReport evaluate(const Molecule& ligand_pred,
                        const Molecule& protein,
                        const Molecule* ligand_true,
                        const EvaluateOptions& opt) {
    PoseBustReport report;
    report.ran     = true;
    report.backend = "native";
    report.n_ligand_atoms = static_cast<int>(ligand_pred.atoms.size());

    // Crop protein to pocket around ligand heavy COM (empty if no protein).
    const Molecule protein_cropped =
        protein.empty()
            ? Molecule{}
            : crop_protein_near_ligand(protein, ligand_pred, opt.protein_crop_A);
    report.n_protein_atoms_cropped = static_cast<int>(protein_cropped.atoms.size());

    // Loading reflects the *input* protein (pre-crop); pocket checks use crop.
    const Molecule* protein_for_loading =
        protein.empty() ? nullptr : &protein;

    // Always: loading + chemistry (full suite diagnostics)
    check_loading(&ligand_pred, protein_for_loading, report.checks);
    check_chemistry_sanity(ligand_pred, report.checks);

    // Geometry (ligand-only; soft fails do not block campaign gate)
    check_distance_geometry(ligand_pred, report.checks);
    check_flatness(ligand_pred, report.checks);

    // Protein-conditioned checks — high-impact campaign gate path
    if (!protein_cropped.empty()) {
        check_intermolecular_distance(ligand_pred, protein_cropped, report.checks);
        check_volume_overlap(ligand_pred, protein_cropped, report.checks);
        fill_continuous_summaries(report);
    } else if (!protein.empty()) {
        // Crop emptied — still emit fail-closed campaign keys via empty-protein path
        CheckItem miss;
        miss.key = "no_clashes";
        miss.label = "Minimum distance to protein";
        miss.passed = false;
        miss.detail = "protein crop empty (no heavy atoms within crop radius of ligand)";
        report.checks.push_back(miss);
        miss.key = "not_too_far_away";
        miss.label = "Protein-ligand maximum distance";
        miss.detail = "protein crop empty";
        report.checks.push_back(miss);
        miss.key = "no_volume_clash";
        miss.label = "Volume overlap with protein";
        miss.detail = "protein crop empty";
        report.checks.push_back(miss);
    }

    // Redock suite: identity vs crystal / true ligand
    if (opt.suite == Suite::Redock && ligand_true != nullptr) {
        check_identity_formula(ligand_pred, ligand_true, report.checks);
    }

    // Optional sidecar: extracted ligand SDF + JSON report
    if (!opt.sidecar_dir.empty()) {
        std::string side_err;
        if (!write_sidecar(ligand_pred, report, opt, &side_err)) {
            // Soft: sidecar I/O must not fail the campaign gate
            if (report.warning.empty())
                report.warning = side_err;
            else
                report.warning += "; " + side_err;
        }
    }

    return report;
}

PoseBustReport evaluate_paths(const std::string& complex_pdb,
                              const std::string& receptor_pdb,
                              const std::string& crystal_sdf,
                              const EvaluateOptions& opt) {
    PoseBustReport report;
    report.ran     = true;
    report.backend = "native_pose_qc";

    // 1) Coordinates from FlexAID pose via CONECT / optimizable residue
    //    (NOT all HETATM — that swallows HEM and cofactors).
    Molecule ligand;
    std::string err;
    if (!load_pdb_flexaid_ligand(complex_pdb, ligand, &err)) {
        report.backend = "error";
        report.error   = err.empty() ? "load_pdb_flexaid_ligand failed" : err;
        return report;
    }

    // 2) Topology from crystal SDF is MANDATORY (fail-closed).
    //    Never fall back to coordinate-inferred bonds for validation.
    Molecule crystal;
    if (crystal_sdf.empty()) {
        report.backend = "error";
        report.error =
            "crystal_sdf required for NativePoseQC (no inferred-bond fallback)";
        return report;
    }
    if (!load_sdf(crystal_sdf, crystal, &err)) {
        report.backend = "error";
        report.error   = err.empty() ? "load_sdf(crystal) failed" : err;
        return report;
    }
    std::string topo_err;
    if (!assign_topology_from_reference(ligand, crystal, &topo_err)) {
        report.backend = "error";
        report.error   = topo_err.empty()
                             ? "assign_topology_from_reference failed"
                             : topo_err;
        report.n_ligand_atoms = static_cast<int>(ligand.atoms.size());
        return report;
    }

    // 3) Protein from receptor apo preferred (no cofactors/ligand); complex fallback.
    Molecule protein;
    if (!receptor_pdb.empty()) {
        if (!load_pdb_protein_heavy(receptor_pdb, protein, &err)) {
            report.backend = "error";
            report.error   = err.empty() ? "load_pdb_protein_heavy(receptor) failed"
                                         : err;
            return report;
        }
    } else {
        std::string soft_err;
        if (!load_pdb_protein_heavy(complex_pdb, protein, &soft_err)) {
            protein = {};
        }
    }

    auto rep = evaluate(ligand, protein, &crystal, opt);
    rep.backend = "native_pose_qc";  // never claim "posebusters"
    return rep;
}

bool write_report_json(const PoseBustReport& report, const std::string& path,
                       std::string* err) {
    std::ofstream out(path);
    if (!out) {
        if (err) *err = "write_report_json: cannot open '" + path + "' for write";
        return false;
    }

    out << "{\n";
    out << "  \"ran\": " << (report.ran ? "true" : "false") << ",\n";
    out << "  \"backend\": \"" << json_escape(report.backend) << "\",\n";
    out << "  \"error\": \"" << json_escape(report.error) << "\",\n";
    out << "  \"warning\": \"" << json_escape(report.warning) << "\",\n";
    out << "  \"all_passed\": " << (report.all_passed() ? "true" : "false") << ",\n";
    // Diagnostic-only fields — NOT DatasetRunner success_pb (that is rmsd∧bust).
    out << "  \"native_qc_diagnostic_pass\": "
        << (report.native_qc_diagnostic_pass() ? "true" : "false") << ",\n";
    out << "  \"success_pb_campaign\": "
        << (report.success_pb_campaign() ? "true" : "false") << ",\n";
    out << "  \"success_pb_full\": " << (report.success_pb_full() ? "true" : "false")
        << ",\n";
    out << "  \"n_pass\": " << report.n_pass() << ",\n";
    out << "  \"n_fail\": " << report.n_fail() << ",\n";
    out << "  \"n_checks\": " << report.n_checks() << ",\n";
    out << "  \"n_ligand_atoms\": " << report.n_ligand_atoms << ",\n";
    out << "  \"n_protein_atoms_cropped\": " << report.n_protein_atoms_cropped << ",\n";
    out << "  \"failed_keys\": \"" << json_escape(report.failed_keys_csv())
        << "\",\n";
    out << "  \"failed_native_qc_keys\": \""
        << json_escape(report.failed_campaign_keys_csv()) << "\",\n";
    out << "  \"failed_campaign_keys\": \""
        << json_escape(report.failed_campaign_keys_csv()) << "\",\n";
    out << "  \"min_lig_prot_dist\": ";
    write_json_number_or_null(out, report.min_lig_prot_dist);
    out << ",\n";
    out << "  \"volume_overlap\": ";
    write_json_number_or_null(out, report.volume_overlap);
    out << ",\n";

    out << "  \"checks\": [\n";
    for (std::size_t i = 0; i < report.checks.size(); ++i) {
        const CheckItem& c = report.checks[i];
        out << "    {\n";
        out << "      \"key\": \"" << json_escape(c.key) << "\",\n";
        out << "      \"label\": \"" << json_escape(c.label) << "\",\n";
        out << "      \"passed\": " << (c.passed ? "true" : "false") << ",\n";
        out << "      \"detail\": \"" << json_escape(c.detail) << "\",\n";
        out << "      \"metric\": ";
        write_json_number_or_null(out, c.metric);
        out << ",\n";
        out << "      \"threshold\": ";
        write_json_number_or_null(out, c.threshold);
        out << ",\n";
        out << "      \"n_checked\": " << c.n_checked << ",\n";
        out << "      \"n_failed\": " << c.n_failed << "\n";
        out << "    }";
        if (i + 1 < report.checks.size()) out << ",";
        out << "\n";
    }
    out << "  ]\n";
    out << "}\n";

    if (!out) {
        if (err) *err = "write_report_json: write failed for '" + path + "'";
        return false;
    }
    return true;
}

Backend resolve_backend_from_env() {
    if (const char* v = std::getenv("FLEXAIDDS_POSEBUST")) {
        if (std::string_view(v) == "0") return Backend::Off;
    }
    if (const char* v = std::getenv("FLEXAIDDS_POSEBUST_BACKEND")) {
        if (iequals(v, "off")) return Backend::Off;
        if (iequals(v, "native") || iequals(v, "native_pose_qc"))
            return Backend::Native;
        if (iequals(v, "bust") || iequals(v, "bust_cli") || iequals(v, "posebusters"))
            return Backend::BustCli;
    }
    // Authoritative default: upstream bust — never claim native is PoseBusters.
    return Backend::BustCli;
}

}  // namespace flexaids::posebust
