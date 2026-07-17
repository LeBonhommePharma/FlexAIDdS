// PartitionFunctionBase.h — Abstract interface for NVT and μVT partition functions
//
// Thermodynamic justification
// ───────────────────────────
// Both the canonical partition function Z(N,V,T) and the grand partition
// function Ξ(μ,V,T) are exponential generating functions over microstates.
// They share:
//   • a temperature T (and β = 1/kT)
//   • a log-space value ln Z or ln Ξ (numerically via stabilized log-sum-exp)
//   • a free-energy map F = −kT ln(·)  (Helmholtz for Z; grand potential Ω for Ξ)
//
// Extracting a common abstract base lets GrandCanonicalEngine and adapters
// for StatMechEngine share a single polymorphic surface without duplicating
// the Shannon-stabilized log-sum-exp kernel (UnifiedHardwareDispatch).
//
// Ranking guardrail (AGENTS.md): PartitionFunctionBase is for the
// thermodynamic ledger only. GA pose ranking remains CF/contact-function
// proxy unless an experimental feature flag reweights search.
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "statmech.h"

#include <cmath>
#include <stdexcept>
#include <string_view>

namespace flexaids {

/// Ensemble family of a concrete partition-function engine.
enum class EnsembleKind {
    NVT,   ///< Canonical: fixed N, V, T → Z
    muVT,  ///< Grand-canonical: fixed μ, V, T → Ξ
};

/// Abstract partition-function surface shared by NVT and μVT engines.
///
/// Implementations must be thread-safe for concurrent *const* queries when
/// documented as such; mutations may require external synchronization.
class PartitionFunctionBase {
public:
    virtual ~PartitionFunctionBase() = default;

    PartitionFunctionBase(const PartitionFunctionBase&) = delete;
    PartitionFunctionBase& operator=(const PartitionFunctionBase&) = delete;
    PartitionFunctionBase(PartitionFunctionBase&&) = delete;
    PartitionFunctionBase& operator=(PartitionFunctionBase&&) = delete;

    /// Ensemble temperature in kelvin.
    [[nodiscard]] virtual double temperature() const noexcept = 0;

    /// β = 1/(k_B T) with k_B in kcal mol⁻¹ K⁻¹.
    [[nodiscard]] virtual double beta() const noexcept {
        return 1.0 / (statmech::kB_kcal * temperature());
    }

    /// ln of the partition function (ln Z for NVT, ln Ξ for μVT).
    [[nodiscard]] virtual double log_partition() const = 0;

    /// Free energy map: F = −kT ln Z (NVT) or Ω = −kT ln Ξ (μVT), kcal/mol.
    /// Thermodynamic: both are Legendre transforms of the appropriate potential.
    [[nodiscard]] virtual double free_energy() const {
        return -(1.0 / beta()) * log_partition();
    }

    /// Ensemble kind discriminator for diagnostic / serialization paths.
    [[nodiscard]] virtual EnsembleKind ensemble() const noexcept = 0;

    /// Human-readable engine name (for logs / audit trails).
    [[nodiscard]] virtual std::string_view engine_name() const noexcept = 0;

protected:
    PartitionFunctionBase() = default;

    static void require_positive_temperature(double T) {
        if (!(T > 0.0))
            throw std::invalid_argument("PartitionFunctionBase: temperature must be > 0 K");
    }
};

/// Lightweight NVT adapter around an already-computed (log_Z, F, S, …) result.
///
/// Does NOT own a StatMechEngine — use when you already called
/// StatMechEngine::compute() and only need the PartitionFunctionBase surface
/// (e.g. feeding Z(N) into GrandCanonicalEngine).
class CanonicalPartitionAdapter final : public PartitionFunctionBase {
public:
    CanonicalPartitionAdapter(double temperature_K, double log_Z)
        : T_(temperature_K)
        , log_Z_(log_Z)
    {
        require_positive_temperature(T_);
    }

    /// Build from a live StatMechEngine (calls compute() once).
    explicit CanonicalPartitionAdapter(const statmech::StatMechEngine& engine)
        : T_(engine.temperature())
    {
        require_positive_temperature(T_);
        if (engine.size() == 0)
            throw std::invalid_argument(
                "CanonicalPartitionAdapter: empty StatMechEngine ensemble");
        auto th = engine.compute();
        log_Z_ = th.log_Z;
        mean_energy_ = th.mean_energy;
        entropy_ = th.entropy;
        has_full_thermo_ = true;
    }

    [[nodiscard]] double temperature() const noexcept override { return T_; }
    [[nodiscard]] double log_partition() const override { return log_Z_; }
    [[nodiscard]] EnsembleKind ensemble() const noexcept override {
        return EnsembleKind::NVT;
    }
    [[nodiscard]] std::string_view engine_name() const noexcept override {
        return "CanonicalPartitionAdapter";
    }

    /// Optional entropy from full thermo (0 if constructed from log_Z only).
    [[nodiscard]] double entropy() const noexcept { return entropy_; }
    [[nodiscard]] double mean_energy() const noexcept { return mean_energy_; }
    [[nodiscard]] bool has_full_thermo() const noexcept { return has_full_thermo_; }

private:
    double T_ = 300.0;
    double log_Z_ = 0.0;
    double mean_energy_ = 0.0;
    double entropy_ = 0.0;
    bool has_full_thermo_ = false;
};

} // namespace flexaids
