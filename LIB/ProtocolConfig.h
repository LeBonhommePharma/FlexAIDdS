// ProtocolConfig.h — Typed protocol knobs with env-as-compat adapter
//
// Audit path: many FLEXAIDDS_* getenv() fragments currently configure
// benchmarks and the engine ad hoc. ProtocolConfig is the typed home for
// common science/benchmark knobs; getenv remains a compatibility adapter via
// from_env() until call sites are migrated (see docs/implementation/
// protocol-config.md).
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstdint>
#include <optional>
#include <string>

namespace flexaids {

/// Typed snapshot of high-traffic protocol / benchmark knobs.
/// Defaults match historical getenv fallbacks (byte-stable when env is unset).
struct ProtocolConfig {
    // ── Seed / restarts (DatasetRunner) ──────────────────────────────────
    std::uint64_t seed_base{0};       ///< FLEXAIDDS_SEED_BASE
    int restarts{5};                  ///< FLEXAIDDS_RESTARTS (min 1)
    /// When true and restarts > 1, launch restart workers in parallel.
    /// Env FLEXAIDDS_PARALLEL_RESTARTS: if unset, defaults to (restarts > 1).
    bool parallel_restarts{true};
    bool parallel_restarts_explicit{false}; ///< true if env overrode default

    // ── VCT / scoring knobs ──────────────────────────────────────────────
    double vct_r0{7.0};               ///< FLEXAIDDS_VCT_R0
    bool vct_normalize_contacts{false}; ///< FLEXAIDDS_VCT_NORM (presence)
    double vct_entropy_weight{0.0};   ///< FLEXAIDDS_VCT_ENTROPY_WEIGHT
    /// Niche-sharing exponent. nullopt → DatasetRunner uses pop-scaled 4.0.
    std::optional<double> sharing_alpha; ///< FLEXAIDDS_SHARING_ALPHA
    /// Boom inject fraction for legacy seed modes. nullopt → 1.0.
    std::optional<double> boom_frac;  ///< FLEXAIDDS_BOOM_FRAC
    int n_elite{1};                   ///< FLEXAIDDS_N_ELITE
    bool use_shannon{false};          ///< FLEXAIDDS_USE_SHANNON (presence)

    // ── Thermo engine (dock_config emission) ─────────────────────────────
    bool thermo_enabled{false};       ///< FLEXAIDDS_THERMO (presence)
    float t_eff{0.596f};              ///< FLEXAIDDS_T_EFF
    float tencom_scale{1.0f};         ///< FLEXAIDDS_TENCOM_SCALE

    // ── Paths ────────────────────────────────────────────────────────────
    std::string data_dir;             ///< FLEXAIDDS_DATA_DIR (empty = unset)

    // ── DatasetRunner pose-selector gates ────────────────────────────────
    bool cf_window_selector{false};   ///< FLEXAIDDS_CF_WINDOW_SELECTOR
    bool cluster_member_emit{false};  ///< FLEXAIDDS_CLUSTER_MEMBER_EMIT

    // ── Engine init (top.cpp) ────────────────────────────────────────────
    double hbond_weight{-2.5};        ///< FLEXAIDDS_HBOND_WEIGHT

    /// Built-in defaults (no env consultation).
    static ProtocolConfig defaults();

    /// Compatibility adapter: read FLEXAIDDS_* env vars into typed fields.
    /// Unset variables keep the same defaults as the pre-migration getenv sites.
    static ProtocolConfig from_env();

    /// Minimal JSON object of the typed fields (for audit / serialization).
    [[nodiscard]] std::string to_json() const;

    /// Parse a JSON object previously emitted by to_json() (or a superset).
    /// Missing keys keep defaults(). Throws std::runtime_error on bad JSON.
    static ProtocolConfig from_json(const std::string& json_text);

    /// Effective sharing alpha given population scaling (mirrors DatasetRunner).
    [[nodiscard]] double effective_sharing_alpha(int pop_base, int pop_scaled) const;

    /// Effective boom inject fraction: mode-aware caller still forces 0.0 for
    /// no-seed benchmark modes; this returns the env/default for seed modes.
    [[nodiscard]] double effective_boom_frac() const;
};

}  // namespace flexaids
