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
    /// True when FLEXAIDDS_N_ELITE was present (apply_config override gate).
    bool n_elite_set{false};
    /// True when FLEXAIDDS_VCT_ENTROPY_WEIGHT was present (apply_config gate).
    bool vct_entropy_weight_set{false};
    bool use_shannon{false};          ///< FLEXAIDDS_USE_SHANNON (presence)

    // ── Ranking / emission + GA ablation (config_parser / engine) ─────────
    /// nullopt = env unset (keep JSON/default). true/false = explicit override.
    std::optional<bool> force_cf_rank_emission;   ///< FLEXAIDDS_FORCE_CF_RANK_EMISSION
    /// nullopt = unset. false forces CF emission (historical classic=0 path).
    std::optional<bool> classic_entropy_ranking;  ///< FLEXAIDDS_CLASSIC_ENTROPY_RANKING
    /// nullopt = unset → keep JSON ga.entropy_weight; else override.
    std::optional<double> entropy_weight;         ///< FLEXAIDDS_ENTROPY_WEIGHT
    /// nullopt = unset → keep JSON ga.diversity_monitoring; else override.
    std::optional<bool> diversity_monitoring;     ///< FLEXAIDDS_DIVERSITY_MONITORING

    // ── Thermo engine (dock_config emission) ─────────────────────────────
    bool thermo_enabled{false};       ///< FLEXAIDDS_THERMO (presence)
    float t_eff{0.596f};              ///< FLEXAIDDS_T_EFF
    float tencom_scale{1.0f};         ///< FLEXAIDDS_TENCOM_SCALE

    // ── Paths ────────────────────────────────────────────────────────────
    std::string data_dir;             ///< FLEXAIDDS_DATA_DIR (empty = unset)
    std::string oracle_site_dir;      ///< FLEXAIDDS_ORACLE_SITE_DIR
    std::string oracle_site;          ///< FLEXAIDDS_ORACLE_SITE
    std::string cleft_sphere_file;    ///< FLEXAIDDS_CLEFT_SPHERE_FILE

    // ── DatasetRunner pose-selector gates ────────────────────────────────
    bool cf_window_selector{false};   ///< FLEXAIDDS_CF_WINDOW_SELECTOR
    bool cluster_member_emit{false};  ///< FLEXAIDDS_CLUSTER_MEMBER_EMIT
    bool seed_elitism{true};          ///< FLEXAIDDS_SEED_ELITISM (default ON)
    double seed_elitism_delta_cf{10.0}; ///< FLEXAIDDS_SEED_ELITISM_DELTA_CF
    bool freqsel{false};              ///< FLEXAIDDS_FREQSEL
    double freqsel_alpha{12.0};       ///< FLEXAIDDS_FREQSEL_ALPHA
    float freqsel_rmsd{1.5f};         ///< FLEXAIDDS_FREQSEL_RMSD
    bool consensus_scorer{false};     ///< FLEXAIDDS_CONSENSUS_SCORER
    /// tENCoM/H(ω) validator; default ON, disable with FLEXAIDDS_HVIB=0.
    bool hvib_enabled{true};

    // ── Budget / grid (DatasetRunner dock path) ──────────────────────────
    bool ring_flex{false};            ///< FLEXAIDDS_RING_FLEX
    /// 1=pop-scale (default), 0=legacy gen-scale, -1=fixed (off/none/fixed).
    int eval_scale_dihedral{1};       ///< FLEXAIDDS_EVAL_SCALE_DIHEDRAL
    bool budget_scale{true};          ///< FLEXAIDDS_BUDGET_SCALE (default ON)
    bool fine_grid{false};            ///< FLEXAIDDS_FINE_GRID (presence)
    int multi_cleft{0};               ///< FLEXAIDDS_MULTI_CLEFT (0 = off)
    bool cognate_site{false};         ///< FLEXAIDDS_COGNATE_SITE (presence)
    bool score_native{false};         ///< FLEXAIDDS_SCORE_NATIVE (truthy)
    bool native_only{false};          ///< FLEXAIDDS_NATIVE_ONLY (truthy)
    bool use_dp{false};               ///< FLEXAIDDS_USE_DP (=1)
    bool ignore_cache{false};         ///< FLEXAIDDS_IGNORE_CACHE
    bool thermo_csv{false};           ///< FLEXAIDDS_THERMO_CSV (presence)

    // ── Engine init (top.cpp) ────────────────────────────────────────────
    double hbond_weight{-2.5};        ///< FLEXAIDDS_HBOND_WEIGHT

    // ── GA engine (gaboom.cpp) ───────────────────────────────────────────
    bool no_sec{false};               ///< FLEXAIDDS_NO_SEC (presence)
    bool benchmark_mode{false};       ///< FLEXAIDDS_BENCHMARK (presence)
    double t_hot{0.0};                ///< FLEXAIDDS_T_HOT (0 = annealing off)
    /// 0 → keep compile-time GA_INSTREAM_INTERVAL; else override (>=1).
    int instream_interval{0};         ///< FLEXAIDDS_INSTREAM_INTERVAL
    bool chain_norm{false};           ///< FLEXAIDDS_CHAIN_NORM
    bool smfree_require_t{false};     ///< FLEXAIDDS_SMFREE_REQUIRE_T

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
