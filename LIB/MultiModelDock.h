// MultiModelDock.h — CCBM Conformer-Coupled Binding Mode orchestrator
//
// Reads multi-model PDB/CIF files, precomputes receptor strain energies,
// sets up FA_Global with model arrays, and provides output utilities for
// per-model top poses, global ensemble thermodynamics, conformer population
// weights, and entropy decomposition.
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "flexaid.h"
#include "gaboom.h"
#include "BindingMode.h"
#include "CifReader.h"
#include "statmech.h"

#include <string>
#include <vector>
#include <map>

namespace ccbm {

/// Receptor model with its own coordinate array and strain energy
struct ReceptorModel {
    int    model_num;      // PDB MODEL number (1-based)
    int    sequential_idx; // 0-based sequential index
    double strain_energy;  // E_strain(r) relative to reference conformer (kcal/mol)
    std::vector<float> atom_coords;  // flat [x0,y0,z0, x1,y1,z1, ...]
};

/// Per-model docking result summary
struct ModelResult {
    int    model_index;
    double best_score;     // best (lowest) total_energy
    double mean_score;     // Boltzmann-weighted mean energy
    int    n_poses;        // number of poses on this model
    double population;     // p(r) = marginal population weight
};

/// Ensemble thermodynamic summary
struct EnsembleThermo {
    double F;       // Helmholtz free energy (kcal/mol)
    double S;       // total configurational entropy (kcal/mol/K)
    double H;       // average enthalpy (kcal/mol)
    double Cv;      // heat capacity (kcal/mol/K)
    std::vector<double> conformer_populations;  // p(r) per model
    double S_receptor;   // Shannon entropy of conformer distribution
    double S_ligand;     // marginal ligand pose entropy
    double I_mutual;     // ligand-receptor mutual information
    double S_vibrational; // vibrational entropy (from ENCoM)
};

/// Full CCBM docking report
struct CCBMReport {
    std::string receptor_file;
    int n_models;
    int n_total_poses;
    EnsembleThermo thermo;
    std::vector<ModelResult> per_model;
    BindingMode::EntropyDecomposition entropy_decomp;
};


/// MultiModelDock: orchestrates multi-conformer ensemble docking
class MultiModelDock {
public:
    /// Construct with the receptor PDB/CIF file path
    explicit MultiModelDock(const std::string& receptor_path);

    /// Read multi-model file and populate receptor models.
    /// Returns number of models found (≥1).
    int load_models(FA_Global* FA, atom** atoms, resid** residue);

    /// Set strain energies from external source (e.g., ENCoM, AMBER).
    /// If not called, all strains default to 0.0 (uniform prior).
    /// energies[i] = E_strain for model i (kcal/mol, relative to reference).
    void set_strain_energies(const std::vector<double>& energies);

    /// Set strain energies from a file (one value per line).
    void set_strain_energies_from_file(const std::string& strain_file);

    /// Compute relative strain energies using a simple RMSD-based harmonic
    /// approximation: E_strain(r) = k * RMSD²(r, ref) / 2
    /// where ref is the first model and k is a spring constant (kcal/mol/Å²).
    void compute_harmonic_strain(double spring_constant = 1.0);

    /// Set up FA_Global with model arrays for GA integration
    void configure_fa(FA_Global* FA);

    /// Generate ensemble report from a completed BindingPopulation
    CCBMReport generate_report(const BindingPopulation& population) const;

    /// Write report to stdout
    static void print_report(const CCBMReport& report);

    /// Write report to file
    static void write_report(const CCBMReport& report, const std::string& outpath);

    /// Accessor
    int n_models() const { return static_cast<int>(models_.size()); }
    const std::vector<ReceptorModel>& models() const { return models_; }

private:
    std::string receptor_path_;
    std::vector<ReceptorModel> models_;
    bool is_cif_;  // true if .cif/.mmcif extension

    /// Detect file type by extension
    static bool is_cif_file(const std::string& path);
};

}  // namespace ccbm
