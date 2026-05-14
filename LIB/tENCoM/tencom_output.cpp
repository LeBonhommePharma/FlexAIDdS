// tencom_output.cpp — BindingMode-like output for tENCoM vibrational analysis

#include "tencom_output.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstdio>
#include <cstring>

namespace tencom_output {

static std::string json_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        switch (c) {
        case '\\': out += "\\\\"; break;
        case '"': out += "\\\""; break;
        case '\n': out += "\\n"; break;
        case '\r': out += "\\r"; break;
        case '\t': out += "\\t"; break;
        default: out += c; break;
        }
    }
    return out;
}

// ─── FlexPopulation::sort_by_free_energy ────────────────────────────────────

void FlexPopulation::sort_by_free_energy() {
    // Keep mode 0 (reference) first, then sort remaining by delta_F_vib
    if (modes.size() <= 2) return;

    std::sort(modes.begin() + 1, modes.end(),
              [](const FlexMode& a, const FlexMode& b) {
                  return a.delta_F_vib < b.delta_F_vib;
              });
}

// ─── FlexPopulation::write_mode_pdb ─────────────────────────────────────────

void FlexPopulation::write_mode_pdb(const FlexMode& mode,
                                     const tencom_pdb::CalphaStructure& structure) const
{
    // Build output filename: prefix_mode_N.pdb
    std::string outfile = output_prefix + "_mode_" + std::to_string(mode.mode_id) + ".pdb";

    std::ofstream ofs(outfile);
    if (!ofs.is_open()) {
        std::cerr << "Error: cannot write to " << outfile << "\n";
        return;
    }

    // ── REMARK section — standardized KEY=VALUE format ─────────────────────
    // All REMARK lines use consistent KEY=VALUE pairs for machine parsing.
    ofs << "REMARK TENCOM_VERSION=1.0\n";
    ofs << "REMARK TOOL=FlexAIDdS_tENCoM\n";

    ofs << std::fixed;
    ofs << "REMARK MODE_ID=" << mode.mode_id << "\n";
    if (mode.mode_id == 0) {
        ofs << "REMARK MODE_TYPE=reference\n";
    } else {
        ofs << "REMARK MODE_TYPE=target\n";
    }
    ofs << "REMARK SOURCE=" << mode.pdb_path << "\n";

    ofs << std::setprecision(6);
    ofs << "REMARK S_VIB=" << mode.S_vib << "\n";
    ofs << "REMARK DELTA_S_VIB=" << mode.delta_S_vib << "\n";
    ofs << "REMARK S_VIB_STATUS="
        << (mode.frequency_calibrated ? "calibrated" : "model_scale_heuristic")
        << "\n";
    ofs << "REMARK EIGENVALUE_TO_OMEGA="
        << std::setprecision(8) << mode.eigenvalue_to_omega << "\n";
    ofs << "REMARK CALIBRATION_LABEL=" << mode.calibration_label << "\n";
    ofs << "REMARK CALIBRATION_PROVENANCE=" << mode.calibration_provenance << "\n";

    ofs << std::setprecision(4);
    ofs << "REMARK DELTA_F_VIB=" << mode.delta_F_vib << "\n";
    ofs << "REMARK N_MODES=" << mode.n_modes << "\n";
    ofs << "REMARK N_RESIDUES=" << mode.n_residues << "\n";
    ofs << "REMARK TEMPERATURE=" << temperature << "\n";
    ofs << "REMARK FULL_FLEXIBILITY=ON\n";

    // Eigenvalue summary (first few modes)
    if (!mode.mode_data.empty()) {
        int n_show = std::min(static_cast<int>(mode.mode_data.size()), 10);
        for (int i = 0; i < n_show; ++i) {
            const auto& mc = mode.mode_data[i];
            ofs << "REMARK EIGENVALUE_DIFF MODE=" << mc.mode_index
                << " DELTA_EIG=" << std::setprecision(6) << mc.delta_eigenvalue;
            if (!std::isnan(mc.overlap)) {
                ofs << " OVERLAP=" << std::setprecision(4) << mc.overlap;
            }
            ofs << "\n";
        }
        if (static_cast<int>(mode.mode_data.size()) > 10) {
            ofs << "REMARK EIGENVALUE_DIFF_REMAINING="
                << mode.mode_data.size() - 10 << "\n";
        }
    }

    // B-factor summary (per-residue)
    if (!mode.bfactors.empty()) {
        ofs << "REMARK BFACTORS";
        for (size_t i = 0; i < mode.bfactors.size(); ++i) {
            ofs << " " << std::setprecision(2) << mode.bfactors[i];
        }
        ofs << "\n";
    }

    if (!mode.delta_bfactors.empty()) {
        ofs << "REMARK DELTA_BFACTORS";
        for (size_t i = 0; i < mode.delta_bfactors.size(); ++i) {
            ofs << " " << std::setprecision(2) << mode.delta_bfactors[i];
        }
        ofs << "\n";
    }

    // Per-residue vibrational entropy
    if (!mode.per_residue_svib.empty()) {
        ofs << std::setprecision(6);
        ofs << "REMARK PER_RESIDUE_SVIB";
        for (double sv : mode.per_residue_svib) {
            ofs << " " << sv;
        }
        ofs << "\n";
    }

    if (!mode.per_residue_delta_svib.empty()) {
        ofs << std::setprecision(6);
        ofs << "REMARK PER_RESIDUE_DELTA_SVIB";
        for (double dsv : mode.per_residue_delta_svib) {
            ofs << " " << dsv;
        }
        ofs << "\n";
    }

    // Composition summary
    if (structure.n_protein > 0 || structure.n_dna > 0 || structure.n_rna > 0) {
        ofs << "REMARK COMPOSITION";
        if (structure.n_protein > 0) ofs << " PROTEIN=" << structure.n_protein;
        if (structure.n_dna > 0)     ofs << " DNA=" << structure.n_dna;
        if (structure.n_rna > 0)     ofs << " RNA=" << structure.n_rna;
        ofs << "\n";
    }

    // ── ATOM records (backbone representatives) ─────────────────────────────
    for (int ai = 1; ai <= structure.res_cnt; ++ai) {
        const atom& a = structure.atoms[ai];
        const resid& r = structure.residues[ai];

        // B-factor column: use computed B-factor if available, else 0
        float bf = 0.0f;
        if (ai - 1 < static_cast<int>(mode.bfactors.size())) {
            bf = mode.bfactors[ai - 1];
        }

        // Use correct atom name based on residue type
        const char* atom_label = " CA ";
        if (ai < static_cast<int>(structure.residue_types.size())) {
            auto rt = structure.residue_types[ai];
            if (rt == tencom_pdb::ResidueType::DNA ||
                rt == tencom_pdb::ResidueType::RNA) {
                atom_label = " C4'";
            }
        }

        char line[128];
        std::snprintf(line, sizeof(line),
            "ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           C",
            ai, atom_label, r.name, r.chn, r.number,
            a.coor[0], a.coor[1], a.coor[2],
            1.00f, bf);
        ofs << line << "\n";
    }

    ofs << "END\n";
    ofs.close();

    std::cout << "  Wrote: " << outfile << "\n";
}

// ─── FlexPopulation::print_summary ──────────────────────────────────────────

void FlexPopulation::print_summary(std::ostream& os) const {
    os << "\n=== FlexAIDdS tENCoM Vibrational Entropy Differential ===\n";
    os << "Temperature: " << temperature << " K\n";
    os << "Full flexibility: ON (all modes)\n\n";
    if (!modes.empty()) {
        const auto& calibration = modes.front();
        os << "Frequency calibration: "
           << (calibration.frequency_calibrated ? "calibrated" : "model-scale heuristic")
           << " (" << calibration.calibration_label
           << ", scale=" << calibration.eigenvalue_to_omega
           << " rad/s per sqrt(model unit))\n";
        if (!calibration.frequency_calibrated) {
            os << "S_vib values are relative/model-scale, not absolute thermodynamic claims.\n";
        }
        os << "\n";
    }

    os << std::left
       << std::setw(6)  << "Mode"
       << std::setw(40) << "Source"
       << std::setw(16) << "S_vib*"
       << std::setw(16) << "Delta_S_vib"
       << std::setw(16) << "Delta_F_vib"
       << std::setw(10) << "N_modes"
       << std::setw(10) << "N_res"
       << "\n";
    os << std::string(114, '-') << "\n";

    for (const auto& m : modes) {
        // Truncate source path for display
        std::string display_path = m.pdb_path;
        if (display_path.size() > 38) {
            display_path = "..." + display_path.substr(display_path.size() - 35);
        }

        os << std::left << std::fixed
           << std::setw(6)  << m.mode_id
           << std::setw(40) << display_path
           << std::setw(16) << std::setprecision(6) << m.S_vib
           << std::setw(16) << std::setprecision(6) << m.delta_S_vib
           << std::setw(16) << std::setprecision(4) << m.delta_F_vib
           << std::setw(10) << m.n_modes
           << std::setw(10) << m.n_residues
           << "\n";
    }

    // Mode overlap summary for targets
    for (const auto& m : modes) {
        if (m.mode_id == 0 || m.mode_data.empty()) continue;

        os << "\n  Mode " << m.mode_id << " (" << m.label << ") — top eigenvalue differentials:\n";
        int n_show = std::min(static_cast<int>(m.mode_data.size()), 5);
        for (int i = 0; i < n_show; ++i) {
            const auto& mc = m.mode_data[i];
            os << "    mode_" << mc.mode_index
               << "  eig_ref=" << std::setprecision(4) << mc.eigenvalue_ref
               << "  eig_tgt=" << mc.eigenvalue_tgt
               << "  delta=" << mc.delta_eigenvalue;
            if (!std::isnan(mc.overlap)) {
                os << "  overlap=" << mc.overlap;
            }
            os << "\n";
        }
    }

    os << "\n";
    if (!modes.empty() && !modes.front().frequency_calibrated) {
        os << "* heuristic model-scale entropy unless an eigenvalue-to-frequency"
              " calibration is supplied.\n\n";
    }
}

// ─── FlexPopulation::output_all ─────────────────────────────────────────────

void FlexPopulation::output_all(
    const std::vector<tencom_pdb::CalphaStructure>& structures) const
{
    print_summary();

    if (structures.size() != modes.size()) {
        std::cerr << "Warning: structure count (" << structures.size()
                  << ") != mode count (" << modes.size() << ")\n";
    }

    int n = std::min(structures.size(), modes.size());
    for (int i = 0; i < n; ++i) {
        write_mode_pdb(modes[i], structures[i]);
    }

    std::cout << "\ntENCoM analysis complete. "
              << modes.size() << " mode(s) written.\n";
}

// ─── FlexPopulation::write_json ──────────────────────────────────────────────

void FlexPopulation::write_json(
    const std::vector<tencom_pdb::CalphaStructure>& structures) const
{
    std::string outfile = output_prefix + "_results.json";
    std::ofstream ofs(outfile);
    if (!ofs.is_open()) {
        std::cerr << "Error: cannot write to " << outfile << "\n";
        return;
    }

    ofs << std::fixed;
    ofs << "{\n";
    ofs << "  \"tool\": \"FlexAIDdS_tENCoM\",\n";
    ofs << "  \"version\": \"1.0\",\n";
    ofs << "  \"temperature\": " << std::setprecision(1) << temperature << ",\n";
    ofs << "  \"full_flexibility\": true,\n";
    ofs << "  \"modes\": [\n";

    for (size_t mi = 0; mi < modes.size(); ++mi) {
        const auto& m = modes[mi];
        ofs << "    {\n";
        ofs << "      \"mode_id\": " << m.mode_id << ",\n";
        ofs << "      \"source\": \"" << json_escape(m.pdb_path) << "\",\n";
        ofs << "      \"label\": \"" << json_escape(m.label) << "\",\n";
        ofs << "      \"type\": \"" << (m.mode_id == 0 ? "reference" : "target") << "\",\n";
        ofs << "      \"S_vib\": " << std::setprecision(8) << m.S_vib << ",\n";
        ofs << "      \"delta_S_vib\": " << std::setprecision(8) << m.delta_S_vib << ",\n";
        ofs << "      \"delta_F_vib\": " << std::setprecision(6) << m.delta_F_vib << ",\n";
        ofs << "      \"S_vib_status\": \""
            << (m.frequency_calibrated ? "calibrated" : "model_scale_heuristic")
            << "\",\n";
        ofs << "      \"eigenvalue_to_omega\": " << std::setprecision(8)
            << m.eigenvalue_to_omega << ",\n";
        ofs << "      \"calibration_label\": \"" << json_escape(m.calibration_label) << "\",\n";
        ofs << "      \"calibration_provenance\": \""
            << json_escape(m.calibration_provenance) << "\",\n";
        ofs << "      \"n_modes\": " << m.n_modes << ",\n";
        ofs << "      \"n_residues\": " << m.n_residues << ",\n";

        // Composition
        if (mi < structures.size()) {
            const auto& s = structures[mi];
            ofs << "      \"composition\": {\"protein\": " << s.n_protein
                << ", \"dna\": " << s.n_dna
                << ", \"rna\": " << s.n_rna << "},\n";
        }

        // B-factors
        ofs << "      \"bfactors\": [";
        for (size_t i = 0; i < m.bfactors.size(); ++i) {
            if (i > 0) ofs << ", ";
            ofs << std::setprecision(4) << m.bfactors[i];
        }
        ofs << "],\n";

        // Delta B-factors
        ofs << "      \"delta_bfactors\": [";
        for (size_t i = 0; i < m.delta_bfactors.size(); ++i) {
            if (i > 0) ofs << ", ";
            ofs << std::setprecision(4) << m.delta_bfactors[i];
        }
        ofs << "],\n";

        // Per-residue S_vib
        ofs << "      \"per_residue_svib\": [";
        for (size_t i = 0; i < m.per_residue_svib.size(); ++i) {
            if (i > 0) ofs << ", ";
            ofs << std::setprecision(6) << m.per_residue_svib[i];
        }
        ofs << "],\n";

        ofs << "      \"per_residue_delta_svib\": [";
        for (size_t i = 0; i < m.per_residue_delta_svib.size(); ++i) {
            if (i > 0) ofs << ", ";
            ofs << std::setprecision(6) << m.per_residue_delta_svib[i];
        }
        ofs << "],\n";

        // Mode comparisons
        ofs << "      \"eigenvalue_diffs\": [";
        int n_show = std::min(static_cast<int>(m.mode_data.size()), 20);
        for (int i = 0; i < n_show; ++i) {
            const auto& mc = m.mode_data[i];
            if (i > 0) ofs << ", ";
            ofs << "{\"mode\": " << mc.mode_index
                << ", \"eig_ref\": " << std::setprecision(6) << mc.eigenvalue_ref
                << ", \"eig_tgt\": " << mc.eigenvalue_tgt
                << ", \"delta\": " << mc.delta_eigenvalue;
            if (!std::isnan(mc.overlap)) {
                ofs << ", \"overlap\": " << std::setprecision(4) << mc.overlap;
            }
            ofs << "}";
        }
        ofs << "]\n";

        ofs << "    }" << (mi + 1 < modes.size() ? "," : "") << "\n";
    }

    ofs << "  ]\n";
    ofs << "}\n";
    ofs.close();

    std::cout << "  Wrote: " << outfile << "\n";
}

// ─── FlexPopulation::write_csv ──────────────────────────────────────────────

void FlexPopulation::write_csv() const {
    std::string outfile = output_prefix + "_summary.csv";
    std::ofstream ofs(outfile);
    if (!ofs.is_open()) {
        std::cerr << "Error: cannot write to " << outfile << "\n";
        return;
    }

    // Header
    ofs << "mode_id,source,type,S_vib,S_vib_status,eigenvalue_to_omega,calibration_label,"
           "delta_S_vib,delta_F_vib,n_modes,n_residues\n";

    ofs << std::fixed;
    for (const auto& m : modes) {
        ofs << m.mode_id << ","
            << "\"" << m.pdb_path << "\","
            << (m.mode_id == 0 ? "reference" : "target") << ","
            << std::setprecision(8) << m.S_vib << ","
            << (m.frequency_calibrated ? "calibrated" : "model_scale_heuristic") << ","
            << std::setprecision(8) << m.eigenvalue_to_omega << ","
            << "\"" << m.calibration_label << "\","
            << std::setprecision(8) << m.delta_S_vib << ","
            << std::setprecision(6) << m.delta_F_vib << ","
            << m.n_modes << ","
            << m.n_residues << "\n";
    }

    ofs.close();
    std::cout << "  Wrote: " << outfile << "\n";
}

}  // namespace tencom_output
