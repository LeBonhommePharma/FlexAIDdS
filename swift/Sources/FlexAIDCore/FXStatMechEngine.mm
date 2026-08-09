// FXStatMechEngine.mm — Objective-C++ implementation of the StatMechEngine C shim
//
// Bridges statmech::StatMechEngine (C++20) to plain C functions for Swift.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

#include "FXStatMechEngine.h"
#include "FXStatMechBridgeInternal.hpp"

#include <vector>
#include <cstring>

// ─── Memory helpers ─────────────────────────────────────────────────────────

extern "C" void fx_free_doubles(double* ptr) {
    delete[] ptr;
}

extern "C" void fx_free_wham_bins(FXWHAMBin* ptr) {
    delete[] ptr;
}

extern "C" void fx_free_pose_infos(FXPoseInfo* ptr) {
    delete[] ptr;
}

// ─── Lifecycle ──────────────────────────────────────────────────────────────

extern "C" FXStatMechEngineRef fx_statmech_create(double temperature_K) {
    return new FXStatMechEngineImpl(temperature_K);
}

extern "C" void fx_statmech_destroy(FXStatMechEngineRef engine) {
    delete engine;
}

extern "C" void fx_statmech_set_scientific_provenance(
    FXStatMechEngineRef engine,
    int32_t schema_version,
    int32_t energy_domain,
    int32_t ensemble_measure,
    int32_t reference_state,
    const char* energy_provenance,
    const char* measure_provenance,
    const char* reference_provenance
) {
    if (!engine) return;

    statmech::ScientificProvenance provenance;
    provenance.schema_version = schema_version;
    provenance.energy_domain = statmech_energy_domain_from_fx(energy_domain);
    provenance.ensemble_measure = statmech_ensemble_measure_from_fx(ensemble_measure);
    provenance.reference_state = statmech_reference_state_from_fx(reference_state);
    provenance.energy_provenance = energy_provenance ? energy_provenance : "";
    provenance.measure_provenance = measure_provenance ? measure_provenance : "";
    provenance.reference_provenance = reference_provenance ? reference_provenance : "";
    engine->engine.set_provenance(provenance);
}

// ─── Sample management ──────────────────────────────────────────────────────

extern "C" void fx_statmech_add_sample(FXStatMechEngineRef engine, double energy, int multiplicity) {
    if (engine) engine->engine.add_sample(energy, multiplicity);
}

extern "C" void fx_statmech_clear(FXStatMechEngineRef engine) {
    if (engine) engine->engine.clear();
}

extern "C" int fx_statmech_size(FXStatMechEngineRef engine) {
    return engine ? static_cast<int>(engine->engine.size()) : 0;
}

// ─── Thermodynamic computation ──────────────────────────────────────────────

extern "C" FXThermodynamics fx_statmech_compute(FXStatMechEngineRef engine) {
    if (!engine) {
        FXThermodynamics empty = {};
        return empty;
    }
    return fx_thermodynamics_from_cpp(engine->engine.compute());
}

extern "C" double* fx_statmech_boltzmann_weights(FXStatMechEngineRef engine, int* out_count) {
    if (!engine || !out_count) {
        if (out_count) *out_count = 0;
        return nullptr;
    }
    auto weights = engine->engine.boltzmann_weights();
    *out_count = static_cast<int>(weights.size());
    if (weights.empty()) return nullptr;

    double* result = new double[weights.size()];
    std::memcpy(result, weights.data(), weights.size() * sizeof(double));
    return result;
}

// ─── Comparative analysis ───────────────────────────────────────────────────

extern "C" double fx_statmech_delta_G(FXStatMechEngineRef engine, FXStatMechEngineRef reference) {
    if (!engine || !reference) return 0.0;
    return engine->engine.delta_G(reference->engine);
}

// ─── Static / pure functions ────────────────────────────────────────────────

extern "C" double fx_statmech_helmholtz(const double* energies, int count, double temperature) {
    if (!energies || count <= 0) return 0.0;
    std::span<const double> span(energies, static_cast<size_t>(count));
    return statmech::StatMechEngine::helmholtz(span, temperature);
}

extern "C" double fx_statmech_thermodynamic_integration(const FXTIPoint* points, int count) {
    if (!points || count <= 0) return 0.0;
    // Convert FXTIPoint array to statmech::TIPoint vector
    std::vector<statmech::TIPoint> ti_points(count);
    for (int i = 0; i < count; ++i) {
        ti_points[i].lambda = points[i].lambda;
        ti_points[i].dV_dlambda = points[i].dV_dlambda;
    }
    return statmech::StatMechEngine::thermodynamic_integration(ti_points);
}

extern "C" FXWHAMBin* fx_statmech_boltzmann_pmf(const double* energies, const double* coordinates,
                                                 int count, double temperature, int n_bins,
                                                 int max_iter, double tolerance, int* out_count) {
    if (!energies || !coordinates || count <= 0 || !out_count) {
        if (out_count) *out_count = 0;
        return nullptr;
    }

    std::span<const double> e_span(energies, static_cast<size_t>(count));
    std::span<const double> c_span(coordinates, static_cast<size_t>(count));

    auto bins = statmech::StatMechEngine::boltzmann_pmf(e_span, c_span, temperature,
                                                         n_bins, max_iter, tolerance);
    *out_count = static_cast<int>(bins.size());
    if (bins.empty()) return nullptr;

    FXWHAMBin* result = new FXWHAMBin[bins.size()];
    for (size_t i = 0; i < bins.size(); ++i) {
        result[i].coord_center = bins[i].coord_center;
        result[i].count        = bins[i].count;
        result[i].free_energy  = bins[i].free_energy;
    }
    return result;
}

// ─── Accessors ──────────────────────────────────────────────────────────────

extern "C" double fx_statmech_temperature(FXStatMechEngineRef engine) {
    return engine ? engine->engine.temperature() : 0.0;
}

extern "C" double fx_statmech_beta(FXStatMechEngineRef engine) {
    return engine ? engine->engine.beta() : 0.0;
}
