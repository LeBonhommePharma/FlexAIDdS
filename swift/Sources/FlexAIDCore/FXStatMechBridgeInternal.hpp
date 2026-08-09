// FXStatMechBridgeInternal.hpp — private C++/C thermodynamics bridge helpers
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "FXStatMechEngine.h"
#include "statmech.h"

#include <cstring>
#include <string>

struct FXStatMechEngineImpl {
    statmech::StatMechEngine engine;
    explicit FXStatMechEngineImpl(double temperature_K) : engine(temperature_K) {}
};

inline int32_t fx_energy_domain(statmech::EnergyDomain value) noexcept {
    switch (value) {
    case statmech::EnergyDomain::ContactFunctionArbitraryUnits:
        return FX_ENERGY_DOMAIN_CF_ARBITRARY_UNITS;
    case statmech::EnergyDomain::CalibratedKcalPerMol:
        return FX_ENERGY_DOMAIN_CALIBRATED_KCAL_PER_MOL;
    case statmech::EnergyDomain::ModelScale:
        return FX_ENERGY_DOMAIN_MODEL_SCALE;
    case statmech::EnergyDomain::Unclassified:
    default:
        return FX_ENERGY_DOMAIN_UNCLASSIFIED;
    }
}

inline int32_t fx_ensemble_measure(statmech::EnsembleMeasure value) noexcept {
    switch (value) {
    case statmech::EnsembleMeasure::OptimizerSamples:
        return FX_ENSEMBLE_MEASURE_OPTIMIZER_SAMPLES;
    case statmech::EnsembleMeasure::EnumeratedMicrostates:
        return FX_ENSEMBLE_MEASURE_ENUMERATED_MICROSTATES;
    case statmech::EnsembleMeasure::WeightedQuadrature:
        return FX_ENSEMBLE_MEASURE_WEIGHTED_QUADRATURE;
    case statmech::EnsembleMeasure::Unclassified:
    default:
        return FX_ENSEMBLE_MEASURE_UNCLASSIFIED;
    }
}

inline int32_t fx_reference_state(statmech::ReferenceState value) noexcept {
    switch (value) {
    case statmech::ReferenceState::BoundOnly:
        return FX_REFERENCE_STATE_BOUND_ONLY;
    case statmech::ReferenceState::MatchedAssociationCycle:
        return FX_REFERENCE_STATE_MATCHED_ASSOCIATION_CYCLE;
    case statmech::ReferenceState::None:
    default:
        return FX_REFERENCE_STATE_NONE;
    }
}

inline statmech::EnergyDomain statmech_energy_domain_from_fx(int32_t value) noexcept {
    switch (value) {
    case FX_ENERGY_DOMAIN_CF_ARBITRARY_UNITS:
        return statmech::EnergyDomain::ContactFunctionArbitraryUnits;
    case FX_ENERGY_DOMAIN_CALIBRATED_KCAL_PER_MOL:
        return statmech::EnergyDomain::CalibratedKcalPerMol;
    case FX_ENERGY_DOMAIN_MODEL_SCALE:
        return statmech::EnergyDomain::ModelScale;
    case FX_ENERGY_DOMAIN_UNCLASSIFIED:
    default:
        return statmech::EnergyDomain::Unclassified;
    }
}

inline statmech::EnsembleMeasure statmech_ensemble_measure_from_fx(int32_t value) noexcept {
    switch (value) {
    case FX_ENSEMBLE_MEASURE_OPTIMIZER_SAMPLES:
        return statmech::EnsembleMeasure::OptimizerSamples;
    case FX_ENSEMBLE_MEASURE_ENUMERATED_MICROSTATES:
        return statmech::EnsembleMeasure::EnumeratedMicrostates;
    case FX_ENSEMBLE_MEASURE_WEIGHTED_QUADRATURE:
        return statmech::EnsembleMeasure::WeightedQuadrature;
    case FX_ENSEMBLE_MEASURE_UNCLASSIFIED:
    default:
        return statmech::EnsembleMeasure::Unclassified;
    }
}

inline statmech::ReferenceState statmech_reference_state_from_fx(int32_t value) noexcept {
    switch (value) {
    case FX_REFERENCE_STATE_BOUND_ONLY:
        return statmech::ReferenceState::BoundOnly;
    case FX_REFERENCE_STATE_MATCHED_ASSOCIATION_CYCLE:
        return statmech::ReferenceState::MatchedAssociationCycle;
    case FX_REFERENCE_STATE_NONE:
    default:
        return statmech::ReferenceState::None;
    }
}

inline void fx_copy_provenance_text(
    char* destination,
    const std::string& source
) noexcept {
    // Evidence must cross the ABI losslessly. Truncation or an embedded NUL
    // could otherwise turn a hostile longer string into an apparently valid
    // `sha256:<digest>` receipt on the Swift side.
    if (source.size() >= FX_PROVENANCE_TEXT_CAPACITY ||
        source.find('\0') != std::string::npos) {
        destination[0] = '\0';
        return;
    }
    std::memcpy(destination, source.data(), source.size());
    destination[source.size()] = '\0';
}

inline FXScientificProvenance fx_scientific_provenance_from_cpp(
    const statmech::ScientificProvenance& provenance
) noexcept {
    FXScientificProvenance result{};
    result.schema_version = provenance.schema_version;
    result.energy_domain = fx_energy_domain(provenance.energy_domain);
    result.ensemble_measure = fx_ensemble_measure(provenance.ensemble_measure);
    result.reference_state = fx_reference_state(provenance.reference_state);
    fx_copy_provenance_text(result.energy_provenance, provenance.energy_provenance);
    fx_copy_provenance_text(result.measure_provenance, provenance.measure_provenance);
    fx_copy_provenance_text(result.reference_provenance, provenance.reference_provenance);
    return result;
}

inline FXThermodynamics fx_thermodynamics_from_cpp(
    const statmech::Thermodynamics& thermodynamics
) noexcept {
    FXThermodynamics result{};
    result.temperature = thermodynamics.temperature;
    result.log_Z = thermodynamics.log_Z;
    result.free_energy = thermodynamics.free_energy;
    result.mean_energy = thermodynamics.mean_energy;
    result.mean_energy_sq = thermodynamics.mean_energy_sq;
    result.heat_capacity = thermodynamics.heat_capacity;
    result.entropy = thermodynamics.entropy;
    result.std_energy = thermodynamics.std_energy;
    result.scientific_provenance =
        fx_scientific_provenance_from_cpp(thermodynamics.provenance);
    return result;
}
