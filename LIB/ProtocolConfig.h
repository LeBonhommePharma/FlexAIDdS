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
    /// Reporting-only temperature for the whiteboard diagnostics (I_ES, ΔS_j,
    /// binding regime classifier, Boltzmann P_i) — ISMB 2017 calibration
    /// (kT_ISMB in ThermoWhiteboard.h). Per the whiteboard, T=21 is baked
    /// into the LEFT-hand quantity itself (ΔG₂₁, P_i(T=21)), not just
    /// substituted on the RHS — this field is that named constant.
    /// Independent of t_eff: never enters G_bind/CF scoring or GA selection.
    float report_T{21.0f};            ///< FLEXAIDDS_REPORT_T / --temperature

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
    // ── Two-gate spread guard (opt-in; 0 = disabled) ─────────────────────
    // Demotes a rank-0 cluster head only when it is BOTH spatially isolated
    // from its top-4 peers AND holds a minority of the merged population, AND
    // no quorum of restarts independently converges on it. The single-gate
    // version (d7ef67380, reverted in 024ba8068) demoted on isolation alone and
    // cost 64/85 Astex targets, so every field below defaults to a no-op.
    float cluster_spread_max{0.0f};        ///< FLEXAIDDS_CLUSTER_SPREAD_MAX (0 = off)
    /// Rank-0 must hold strictly less than this fraction of the merged
    /// population to be eligible for demotion.
    float cluster_pop_min_fraction{0.35f}; ///< FLEXAIDDS_CLUSTER_POP_MIN_FRACTION
    /// Å — a restart head within this radius of rank-0 counts as agreeing.
    float cluster_consensus_tau{2.0f};     ///< FLEXAIDDS_CLUSTER_CONSENSUS_TAU
    /// Minimum agreeing restarts that veto demotion.
    int   cluster_consensus_k{3};          ///< FLEXAIDDS_CLUSTER_CONSENSUS_K
    /// Pocket radius in Å used for the isolation threshold θ = 0.70·r.
    /// 0 = unknown → fall back to the population Q75 pairwise-RMSD spread.
    float cluster_pocket_radius{0.0f};     ///< FLEXAIDDS_CLUSTER_POCKET_RADIUS
    bool consensus_scorer{false};     ///< FLEXAIDDS_CONSENSUS_SCORER
    /// v135 crystal-blind basin recovery election (BCR-proxy). Default OFF —
    /// preserves claim ranking (AGENTS.md). Master switch enables:
    ///   • include freq=1 cluster heads (no freq>1 gate)
    ///   • score temperature τ in CF a.u. (not mixed kT=0.592 kcal-scale)
    bool election_v135{false};        ///< FLEXAIDDS_ELECTION_V135
    /// Score temperature τ in **CF arbitrary units** for Z-like composite.
    /// 0 = legacy hard-coded 0.592 (unit-mixed). When election_v135 and unset,
    /// from_env applies default τ=25 CF a.u.
    double election_score_tau{0.0};   ///< FLEXAIDDS_ELECTION_SCORE_TAU
    /// Include Frequency=1 heads in the election pool (helps BCR-like singletons).
    bool election_include_singletons{false}; ///< FLEXAIDDS_ELECTION_INCLUDE_SINGLETONS
    /// Rank already-clustered heads by Softβ Ĝ = H̃ − T S̃ (3Dsig / SoftBetaFreeEnergy).
    /// **Default OFF** for classic pilot and claim harness until explicitly opted in.
    /// Softβ reorders modes only — not a sampling method; cannot fix BCR=0.
    /// Opt in (either alias):
    ///   FLEXAIDDS_SOFTBETA_ELECTION=1  (preferred clear name)
    ///   FLEXAIDDS_ELECTION_SHANNON_F=1 (legacy alias, same bit)
    /// Force OFF path: FLEXAIDDS_ELECTION_LEGACY_ZH=1.
    /// When unset → Softβ S1 OFF (CF rank-0 / legacy ZH; no Softβ claim).
    bool election_shannon_free_energy{false}; ///< Softβ S1 (SOFTBETA_ELECTION / SHANNON_F)
    /// Soft-β temperature T for Ĝ. 0 → resolve at election:
    /// dock TEMPER / DockingConfig::temperature, else (legacy ZH) SCORE_TAU, else 298 K.
    /// Units: CF is scoring-proxy a.u.; FlexAID soft-β (β=1/T), **not** k_B·T kcal.
    /// TEMPER 21 on arm B is engine soft-T for FO/ACF, not physical kT.
    double election_soft_T{0.0};      ///< FLEXAIDDS_ELECTION_SOFT_T
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

/// FLEXAIDDS_THERMO_SCORE (truthy, default OFF): promote ΔG_eff = <CF> − T·H
/// from a reported diagnostic to the ranking criterion, in place of min(CF).
/// Read once on first use. Kept as a free function rather than a ProtocolConfig
/// field because it is consulted from scoring/reporting paths that do not build
/// a ProtocolConfig.
bool thermo_score_enabled();

/// True when the PoseBust pocket-presence penalty is configured to run.
///
/// The single source of truth is FA->pb_pocket_weight, which top.cpp populates
/// from FLEXAIDDS_PB_POCKET_WEIGHT and config_parser.cpp from the JSON key
/// scoring.pb_pocket_weight. Deliberately NOT a getenv() of its own: reading the
/// env independently would reproduce the FLEXAIDDS_PB_CLASH_ELECT_WEIGHT footgun
/// where the config path and the env path disagree and the term silently no-ops
/// for a JSON-only campaign.
///
/// Inline and branch-only so the hot path pays nothing when the term is off.
[[nodiscard]] inline bool pb_pocket_enabled(double pb_pocket_weight) noexcept {
    return pb_pocket_weight > 0.0;
}

}  // namespace flexaids
