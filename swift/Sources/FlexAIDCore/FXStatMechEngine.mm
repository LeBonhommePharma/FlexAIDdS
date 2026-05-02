// FXStatMechEngine.mm — Objective-C++ implementation of the StatMechEngine C shim
//
// Bridges statmech::StatMechEngine (C++20) to plain C functions for Swift.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

#include "FXStatMechEngine.h"
#include "statmech.h"

#include <vector>
#include <cstring>

// ─── Opaque implementation ──────────────────────────────────────────────────

struct FXStatMechEngineImpl {
    statmech::StatMechEngine engine;
    explicit FXStatMechEngineImpl(double T) : engine(T) {}
};

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

// ─── Helper: C++ Thermodynamics → C FXThermodynamics ────────────────────────

static FXThermodynamics to_fx(const statmech::Thermodynamics& t) {
    FXThermodynamics fx;
    fx.temperature    = t.temperature;
    fx.log_Z          = t.log_Z;
    fx.free_energy    = t.free_energy;
    fx.mean_energy    = t.mean_energy;
    fx.mean_energy_sq = t.mean_energy_sq;
    fx.heat_capacity  = t.heat_capacity;
    fx.entropy        = t.entropy;
    fx.std_energy     = t.std_energy;
    return fx;
}

// ─── Lifecycle ──────────────────────────────────────────────────────────────

extern "C" FXStatMechEngineRef fx_statmech_create(double temperature_K) {
    return new FXStatMechEngineImpl(temperature_K);
}

extern "C" void fx_statmech_destroy(FXStatMechEngineRef engine) {
    delete engine;
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
    return to_fx(engine->engine.compute());
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
