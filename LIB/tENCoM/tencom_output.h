// tencom_output.h — BindingMode-like output for tENCoM vibrational analysis
//
// Provides FlexMode (analogous to BindingMode) and FlexPopulation (analogous
// to BindingPopulation) for outputting tENCoM results as PDB files with
// REMARK thermodynamic metadata — without any GA/chromosome dependencies.
//
// The reference structure is always mode 0.
// Full flexibility is always on by default (no option to disable).
#pragma once

#include "tencom_diff.h"
#include "pdb_calpha.h"

#include <string>
#include <vector>
#include <iostream>

namespace tencom_output {

// One flexibility mode — analogous to BindingMode but for global vibrational analysis
struct FlexMode {
    int         mode_id;        // 0 = reference, 1..N = targets
    std::string pdb_path;       // source PDB file path
    std::string label;          // human-readable label

    // Thermodynamic data
    double S_vib       = 0.0;   // vibrational entropy; heuristic unless calibrated
    double delta_S_vib = 0.0;   // relative to reference (0 for mode 0)
    double delta_F_vib = 0.0;   // -T × delta_S_vib (kcal/mol)
    double eigenvalue_to_omega = 1.0;
    bool   frequency_calibrated = false;
    std::string calibration_label = "model-scale";
    std::string calibration_provenance =
        "No mass/inertia or empirical eigenvalue-to-frequency calibration supplied";

    // Per-mode eigenvalue data
    std::vector<tencom_diff::ModeComparison> mode_data;

    // Per-residue B-factors and differentials
    std::vector<float> bfactors;
    std::vector<float> delta_bfactors;

    // Per-residue vibrational entropy decomposition
    std::vector<double> per_residue_svib;
    std::vector<double> per_residue_delta_svib;

    int n_modes   = 0;          // number of non-trivial normal modes
    int n_residues = 0;
};

// Population of FlexModes — analogous to BindingPopulation
struct FlexPopulation {
    double      temperature = 300.0;
    std::string output_prefix = "tencom";
    std::vector<FlexMode> modes;   // mode 0 = reference

    // Sort modes by delta_F_vib (ascending, most stabilizing first)
    void sort_by_free_energy();

    // Write PDB file for a given mode with REMARK thermodynamic metadata
    void write_mode_pdb(const FlexMode& mode,
                        const tencom_pdb::CalphaStructure& structure) const;

    // Write summary table to an output stream (default: stdout)
    void print_summary(std::ostream& os = std::cout) const;

    // Write all outputs: summary + PDB files
    void output_all(const std::vector<tencom_pdb::CalphaStructure>& structures) const;

    // Write JSON results file
    void write_json(const std::vector<tencom_pdb::CalphaStructure>& structures) const;

    // Write CSV summary table
    void write_csv() const;
};

}  // namespace tencom_output
