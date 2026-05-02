// FXStatMechEngine.h — C interface to statmech::StatMechEngine
//
// Opaque-pointer pattern wrapping the C++20 StatMechEngine class.
// All functions are extern "C" safe for Swift consumption.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "FXTypes.h"

#ifdef __cplusplus
extern "C" {
#endif

// ─── Opaque handle ──────────────────────────────────────────────────────────

typedef struct FXStatMechEngineImpl* FXStatMechEngineRef;

// ─── Lifecycle ──────────────────────────────────────────────────────────────

FXStatMechEngineRef fx_statmech_create(double temperature_K);
void fx_statmech_destroy(FXStatMechEngineRef engine);

// ─── Sample management ──────────────────────────────────────────────────────

void fx_statmech_add_sample(FXStatMechEngineRef engine, double energy, int multiplicity);
void fx_statmech_clear(FXStatMechEngineRef engine);
int  fx_statmech_size(FXStatMechEngineRef engine);

// ─── Thermodynamic computation ──────────────────────────────────────────────

FXThermodynamics fx_statmech_compute(FXStatMechEngineRef engine);

// Returns heap-allocated array; caller must free via fx_free_doubles()
double* fx_statmech_boltzmann_weights(FXStatMechEngineRef engine, int* out_count);

// ─── Comparative analysis ───────────────────────────────────────────────────

// Delta-G relative to another engine's ensemble
double fx_statmech_delta_G(FXStatMechEngineRef engine, FXStatMechEngineRef reference);

// ─── Static / pure functions ────────────────────────────────────────────────

// Helmholtz free energy from raw energy array
double fx_statmech_helmholtz(const double* energies, int count, double temperature);

// Thermodynamic integration via trapezoidal rule
double fx_statmech_thermodynamic_integration(const FXTIPoint* points, int count);

// Boltzmann-reweighted PMF (free energy profile along a reaction coordinate)
// Returns heap-allocated array; caller must free via fx_free_wham_bins()
FXWHAMBin* fx_statmech_boltzmann_pmf(const double* energies, const double* coordinates,
                                      int count, double temperature, int n_bins,
                                      int max_iter, double tolerance, int* out_count);

// Backward-compatible alias
#define fx_statmech_wham fx_statmech_boltzmann_pmf

// ─── Accessors ──────────────────────────────────────────────────────────────

double fx_statmech_temperature(FXStatMechEngineRef engine);
double fx_statmech_beta(FXStatMechEngineRef engine);

#ifdef __cplusplus
}
#endif
