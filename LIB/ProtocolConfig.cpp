// ProtocolConfig.cpp — env adapter + light JSON for typed protocol knobs
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
// SPDX-License-Identifier: Apache-2.0

#include "ProtocolConfig.h"
#include "json_value.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <stdexcept>

namespace flexaids {
namespace {

const char* env_raw(const char* name) {
    return std::getenv(name);
}

bool env_present(const char* name) {
    const char* e = env_raw(name);
    return e != nullptr && e[0] != '\0';
}

// Truthy presence: non-null and not an explicit "0"/false/off.
// Used for gates that historically treated any non-null getenv as ON, but
// also accept FLEXAIDDS_FOO=0 to disable when callers used atoi.
bool env_truthy_int(const char* name, bool default_value) {
    const char* e = env_raw(name);
    if (!e || e[0] == '\0') return default_value;
    return std::atoi(e) != 0;
}

// Presence-only historical gates: any non-null getenv is ON (even "0").
bool env_present_any(const char* name) {
    return env_raw(name) != nullptr;
}

std::optional<double> env_opt_double(const char* name) {
    const char* e = env_raw(name);
    if (!e || e[0] == '\0') return std::nullopt;
    try {
        return std::stod(e);
    } catch (...) {
        return std::nullopt;
    }
}

std::optional<int> env_opt_int(const char* name) {
    const char* e = env_raw(name);
    if (!e || e[0] == '\0') return std::nullopt;
    return std::atoi(e);
}

std::string json_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 8);
    for (char c : s) {
        switch (c) {
            case '\\': out += "\\\\"; break;
            case '"':  out += "\\\""; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:   out += c; break;
        }
    }
    return out;
}

void json_bool(std::ostringstream& o, const char* key, bool v, bool trailing_comma = true) {
    o << '"' << key << "\":" << (v ? "true" : "false");
    if (trailing_comma) o << ',';
}

}  // namespace

ProtocolConfig ProtocolConfig::defaults() {
    return ProtocolConfig{};
}

ProtocolConfig ProtocolConfig::from_env() {
    ProtocolConfig cfg = defaults();

    if (const char* e = env_raw("FLEXAIDDS_SEED_BASE")) {
        try {
            cfg.seed_base = std::stoull(e);
        } catch (...) {
            cfg.seed_base = 0;
        }
    }

    if (auto r = env_opt_int("FLEXAIDDS_RESTARTS")) {
        cfg.restarts = std::max(1, *r);
    }

    if (const char* e = env_raw("FLEXAIDDS_PARALLEL_RESTARTS")) {
        cfg.parallel_restarts_explicit = true;
        cfg.parallel_restarts = (std::atoi(e) != 0) && (cfg.restarts > 1);
    } else {
        cfg.parallel_restarts = (cfg.restarts > 1);
    }

    if (auto v = env_opt_double("FLEXAIDDS_VCT_R0")) {
        cfg.vct_r0 = *v;
    }
    cfg.vct_normalize_contacts = env_present("FLEXAIDDS_VCT_NORM");
    if (auto v = env_opt_double("FLEXAIDDS_VCT_ENTROPY_WEIGHT")) {
        cfg.vct_entropy_weight = *v;
        cfg.vct_entropy_weight_set = true;
    }

    cfg.sharing_alpha = env_opt_double("FLEXAIDDS_SHARING_ALPHA");
    cfg.boom_frac = env_opt_double("FLEXAIDDS_BOOM_FRAC");
    if (auto n = env_opt_int("FLEXAIDDS_N_ELITE")) {
        cfg.n_elite = *n;
        cfg.n_elite_set = true;
    }
    // Historical DatasetRunner: presence of FLEXAIDDS_USE_SHANNON enables.
    cfg.use_shannon = env_present("FLEXAIDDS_USE_SHANNON");

    // Ranking / emission + GA ablation (config_parser historical getenv sites).
    // Presence-only: unset keeps JSON/defaults; set applies override.
    {
        const char* e = env_raw("FLEXAIDDS_FORCE_CF_RANK_EMISSION");
        if (e && e[0] != '\0') {
            if (e[0] == '1' || e[0] == 't' || e[0] == 'T' || e[0] == 'y' || e[0] == 'Y')
                cfg.force_cf_rank_emission = true;
            else if (e[0] == '0' || e[0] == 'f' || e[0] == 'F' || e[0] == 'n' || e[0] == 'N')
                cfg.force_cf_rank_emission = false;
        }
    }
    {
        const char* e = env_raw("FLEXAIDDS_CLASSIC_ENTROPY_RANKING");
        if (e && e[0] != '\0') {
            if (e[0] == '0' || e[0] == 'f' || e[0] == 'F' || e[0] == 'n' || e[0] == 'N')
                cfg.classic_entropy_ranking = false;
            else if (e[0] == '1' || e[0] == 't' || e[0] == 'T' || e[0] == 'y' || e[0] == 'Y')
                cfg.classic_entropy_ranking = true;
        }
    }
    cfg.entropy_weight = env_opt_double("FLEXAIDDS_ENTROPY_WEIGHT");
    {
        const char* e = env_raw("FLEXAIDDS_DIVERSITY_MONITORING");
        if (e && e[0] != '\0') {
            cfg.diversity_monitoring = (std::atoi(e) != 0);
        }
    }

    cfg.thermo_enabled = env_present("FLEXAIDDS_THERMO");
    if (auto v = env_opt_double("FLEXAIDDS_T_EFF")) {
        cfg.t_eff = static_cast<float>(*v);
    }
    if (auto v = env_opt_double("FLEXAIDDS_TENCOM_SCALE")) {
        cfg.tencom_scale = static_cast<float>(*v);
    }

    if (const char* e = env_raw("FLEXAIDDS_DATA_DIR")) {
        if (e[0] != '\0') cfg.data_dir = e;
    }
    if (const char* e = env_raw("FLEXAIDDS_ORACLE_SITE_DIR")) {
        if (e[0] != '\0') cfg.oracle_site_dir = e;
    }
    if (const char* e = env_raw("FLEXAIDDS_ORACLE_SITE")) {
        if (e[0] != '\0') cfg.oracle_site = e;
    }
    if (const char* e = env_raw("FLEXAIDDS_CLEFT_SPHERE_FILE")) {
        if (e[0] != '\0') cfg.cleft_sphere_file = e;
    }

    cfg.cf_window_selector =
        env_truthy_int("FLEXAIDDS_CF_WINDOW_SELECTOR", /*default_value=*/false);
    cfg.cluster_member_emit =
        env_truthy_int("FLEXAIDDS_CLUSTER_MEMBER_EMIT", /*default_value=*/false);

    cfg.seed_elitism = env_truthy_int("FLEXAIDDS_SEED_ELITISM", /*default_value=*/true);
    if (auto v = env_opt_double("FLEXAIDDS_SEED_ELITISM_DELTA_CF")) {
        cfg.seed_elitism_delta_cf = *v;
    }
    cfg.freqsel = env_truthy_int("FLEXAIDDS_FREQSEL", /*default_value=*/false);
    if (auto v = env_opt_double("FLEXAIDDS_FREQSEL_ALPHA")) {
        cfg.freqsel_alpha = *v;
    }
    if (auto v = env_opt_double("FLEXAIDDS_FREQSEL_RMSD")) {
        cfg.freqsel_rmsd = static_cast<float>(*v);
    }
    cfg.consensus_scorer =
        env_truthy_int("FLEXAIDDS_CONSENSUS_SCORER", /*default_value=*/false);
    // v135 BCR-proxy election (default OFF — does not change claim ranking).
    cfg.election_v135 =
        env_truthy_int("FLEXAIDDS_ELECTION_V135", /*default_value=*/false);
    if (auto v = env_opt_double("FLEXAIDDS_ELECTION_SCORE_TAU")) {
        cfg.election_score_tau = *v;
    } else if (cfg.election_v135) {
        cfg.election_score_tau = 25.0;  // CF a.u. score temperature (not kcal)
    }
    cfg.election_include_singletons =
        env_truthy_int("FLEXAIDDS_ELECTION_INCLUDE_SINGLETONS",
                       /*default_value=*/cfg.election_v135);
    // Softβ S1 election (DatasetRunner): Ĝ = H̃ − T S̃ over already-clustered
    // heads. **Default OFF** — Softβ is reordering only, not sampling; pilot
    // harness must not claim Softβ S1 unless explicitly opted in.
    // Preferred: FLEXAIDDS_SOFTBETA_ELECTION=1
    // Legacy alias: FLEXAIDDS_ELECTION_SHANNON_F=1 (same bit)
    // Force OFF: FLEXAIDDS_ELECTION_LEGACY_ZH=1
    // Either ON alias wins if set truthy; LEGACY_ZH always forces false.
    {
        const bool legacy_zh =
            env_truthy_int("FLEXAIDDS_ELECTION_LEGACY_ZH", /*default_value=*/false);
        if (legacy_zh) {
            cfg.election_shannon_free_energy = false;
        } else {
            const bool softbeta =
                env_truthy_int("FLEXAIDDS_SOFTBETA_ELECTION",
                               /*default_value=*/false);
            const bool shannon_f =
                env_truthy_int("FLEXAIDDS_ELECTION_SHANNON_F",
                               /*default_value=*/false);
            cfg.election_shannon_free_energy = softbeta || shannon_f;
        }
    }
    if (auto v = env_opt_double("FLEXAIDDS_ELECTION_SOFT_T")) {
        cfg.election_soft_T = *v;
    }
    // HVIB default ON; only FLEXAIDDS_HVIB=0 disables (historical strcmp).
    {
        const char* e = env_raw("FLEXAIDDS_HVIB");
        cfg.hvib_enabled = !(e && std::strcmp(e, "0") == 0);
    }

    cfg.ring_flex = env_truthy_int("FLEXAIDDS_RING_FLEX", /*default_value=*/false);
    {
        const char* e = env_raw("FLEXAIDDS_EVAL_SCALE_DIHEDRAL");
        if (!e || e[0] == '\0') {
            cfg.eval_scale_dihedral = 1;
        } else if (std::strcmp(e, "off") == 0 ||
                   std::strcmp(e, "none") == 0 ||
                   std::strcmp(e, "fixed") == 0) {
            cfg.eval_scale_dihedral = -1;
        } else {
            cfg.eval_scale_dihedral = std::atoi(e);
        }
    }
    cfg.budget_scale = env_truthy_int("FLEXAIDDS_BUDGET_SCALE", /*default_value=*/true);
    cfg.fine_grid = env_present_any("FLEXAIDDS_FINE_GRID");
    if (auto n = env_opt_int("FLEXAIDDS_MULTI_CLEFT")) {
        cfg.multi_cleft = *n;
    }
    // COGNATE_SITE: historical presence gate (any non-null).
    cfg.cognate_site = env_present_any("FLEXAIDDS_COGNATE_SITE");
    {
        const char* e = env_raw("FLEXAIDDS_SCORE_NATIVE");
        cfg.score_native = e && e[0] != '\0' && std::strcmp(e, "0") != 0;
    }
    {
        const char* e = env_raw("FLEXAIDDS_NATIVE_ONLY");
        cfg.native_only = e && e[0] != '\0' && std::strcmp(e, "0") != 0;
    }
    {
        const char* e = env_raw("FLEXAIDDS_USE_DP");
        cfg.use_dp = e && std::strcmp(e, "1") == 0;
    }
    cfg.ignore_cache = env_truthy_int("FLEXAIDDS_IGNORE_CACHE", /*default_value=*/false);
    cfg.thermo_csv = env_present_any("FLEXAIDDS_THERMO_CSV");

    if (auto v = env_opt_double("FLEXAIDDS_HBOND_WEIGHT")) {
        cfg.hbond_weight = *v;
    }

    // gaboom / GA
    cfg.no_sec = env_present_any("FLEXAIDDS_NO_SEC");
    cfg.benchmark_mode = env_present_any("FLEXAIDDS_BENCHMARK");
    if (auto v = env_opt_double("FLEXAIDDS_T_HOT")) {
        cfg.t_hot = *v;
    }
    if (auto n = env_opt_int("FLEXAIDDS_INSTREAM_INTERVAL")) {
        if (*n >= 1) cfg.instream_interval = *n;
    }
    {
        const char* e = env_raw("FLEXAIDDS_CHAIN_NORM");
        cfg.chain_norm = e && e[0] != '\0' && e[0] != '0';
    }
    cfg.smfree_require_t = env_present("FLEXAIDDS_SMFREE_REQUIRE_T");

    return cfg;
}

double ProtocolConfig::effective_sharing_alpha(int pop_base, int pop_scaled) const {
    if (sharing_alpha.has_value()) return *sharing_alpha;
    if (pop_scaled <= 0) return 4.0;
    return 4.0 * static_cast<double>(pop_base) / static_cast<double>(pop_scaled);
}

double ProtocolConfig::effective_boom_frac() const {
    return boom_frac.value_or(1.0);
}

std::string ProtocolConfig::to_json() const {
    std::ostringstream o;
    o.setf(std::ios::fixed);
    o.precision(6);
    o << '{';
    o << "\"seed_base\":" << seed_base << ',';
    o << "\"restarts\":" << restarts << ',';
    json_bool(o, "parallel_restarts", parallel_restarts);
    o << "\"vct_r0\":" << vct_r0 << ',';
    json_bool(o, "vct_normalize_contacts", vct_normalize_contacts);
    o << "\"vct_entropy_weight\":" << vct_entropy_weight << ',';
    json_bool(o, "vct_entropy_weight_set", vct_entropy_weight_set);
    if (sharing_alpha) {
        o << "\"sharing_alpha\":" << *sharing_alpha << ',';
    } else {
        o << "\"sharing_alpha\":null,";
    }
    if (boom_frac) {
        o << "\"boom_frac\":" << *boom_frac << ',';
    } else {
        o << "\"boom_frac\":null,";
    }
    o << "\"n_elite\":" << n_elite << ',';
    json_bool(o, "n_elite_set", n_elite_set);
    if (force_cf_rank_emission) {
        json_bool(o, "force_cf_rank_emission", *force_cf_rank_emission);
    } else {
        o << "\"force_cf_rank_emission\":null,";
    }
    if (classic_entropy_ranking) {
        json_bool(o, "classic_entropy_ranking", *classic_entropy_ranking);
    } else {
        o << "\"classic_entropy_ranking\":null,";
    }
    if (entropy_weight) {
        o << "\"entropy_weight\":" << *entropy_weight << ',';
    } else {
        o << "\"entropy_weight\":null,";
    }
    if (diversity_monitoring) {
        json_bool(o, "diversity_monitoring", *diversity_monitoring);
    } else {
        o << "\"diversity_monitoring\":null,";
    }
    json_bool(o, "use_shannon", use_shannon);
    json_bool(o, "thermo_enabled", thermo_enabled);
    o << "\"t_eff\":" << t_eff << ',';
    o << "\"tencom_scale\":" << tencom_scale << ',';
    o << "\"data_dir\":\"" << json_escape(data_dir) << "\",";
    o << "\"oracle_site_dir\":\"" << json_escape(oracle_site_dir) << "\",";
    o << "\"oracle_site\":\"" << json_escape(oracle_site) << "\",";
    o << "\"cleft_sphere_file\":\"" << json_escape(cleft_sphere_file) << "\",";
    json_bool(o, "cf_window_selector", cf_window_selector);
    json_bool(o, "cluster_member_emit", cluster_member_emit);
    json_bool(o, "seed_elitism", seed_elitism);
    o << "\"seed_elitism_delta_cf\":" << seed_elitism_delta_cf << ',';
    json_bool(o, "freqsel", freqsel);
    o << "\"freqsel_alpha\":" << freqsel_alpha << ',';
    o << "\"freqsel_rmsd\":" << freqsel_rmsd << ',';
    json_bool(o, "consensus_scorer", consensus_scorer);
    json_bool(o, "election_v135", election_v135);
    o << "\"election_score_tau\":" << election_score_tau << ',';
    json_bool(o, "election_include_singletons", election_include_singletons);
    json_bool(o, "election_shannon_free_energy", election_shannon_free_energy);
    o << "\"election_soft_T\":" << election_soft_T << ',';
    json_bool(o, "hvib_enabled", hvib_enabled);
    json_bool(o, "ring_flex", ring_flex);
    o << "\"eval_scale_dihedral\":" << eval_scale_dihedral << ',';
    json_bool(o, "budget_scale", budget_scale);
    json_bool(o, "fine_grid", fine_grid);
    o << "\"multi_cleft\":" << multi_cleft << ',';
    json_bool(o, "cognate_site", cognate_site);
    json_bool(o, "score_native", score_native);
    json_bool(o, "native_only", native_only);
    json_bool(o, "use_dp", use_dp);
    json_bool(o, "ignore_cache", ignore_cache);
    json_bool(o, "thermo_csv", thermo_csv);
    o << "\"hbond_weight\":" << hbond_weight << ',';
    json_bool(o, "no_sec", no_sec);
    json_bool(o, "benchmark_mode", benchmark_mode);
    o << "\"t_hot\":" << t_hot << ',';
    o << "\"instream_interval\":" << instream_interval << ',';
    json_bool(o, "chain_norm", chain_norm);
    json_bool(o, "smfree_require_t", smfree_require_t, /*trailing_comma=*/false);
    o << '}';
    return o.str();
}

ProtocolConfig ProtocolConfig::from_json(const std::string& json_text) {
    ProtocolConfig cfg = defaults();
    const json::Value root = json::parse(json_text);
    if (!root.is_object()) {
        throw std::runtime_error("ProtocolConfig::from_json: root must be object");
    }

    if (!root["seed_base"].is_null())
        cfg.seed_base = static_cast<std::uint64_t>(root["seed_base"].as_int(0));
    if (!root["restarts"].is_null())
        cfg.restarts = std::max(1, root["restarts"].as_int(5));
    if (!root["parallel_restarts"].is_null())
        cfg.parallel_restarts = root["parallel_restarts"].as_bool(true);
    if (!root["vct_r0"].is_null())
        cfg.vct_r0 = root["vct_r0"].as_double(7.0);
    if (!root["vct_normalize_contacts"].is_null())
        cfg.vct_normalize_contacts = root["vct_normalize_contacts"].as_bool(false);
    if (!root["vct_entropy_weight"].is_null()) {
        cfg.vct_entropy_weight = root["vct_entropy_weight"].as_double(0.0);
        cfg.vct_entropy_weight_set = true;
    }
    if (!root["vct_entropy_weight_set"].is_null())
        cfg.vct_entropy_weight_set = root["vct_entropy_weight_set"].as_bool(false);
    if (!root["sharing_alpha"].is_null())
        cfg.sharing_alpha = root["sharing_alpha"].as_double(4.0);
    if (!root["boom_frac"].is_null())
        cfg.boom_frac = root["boom_frac"].as_double(1.0);
    if (!root["n_elite"].is_null()) {
        cfg.n_elite = root["n_elite"].as_int(1);
        cfg.n_elite_set = true;
    }
    if (!root["n_elite_set"].is_null())
        cfg.n_elite_set = root["n_elite_set"].as_bool(false);
    if (!root["force_cf_rank_emission"].is_null())
        cfg.force_cf_rank_emission = root["force_cf_rank_emission"].as_bool(false);
    if (!root["classic_entropy_ranking"].is_null())
        cfg.classic_entropy_ranking = root["classic_entropy_ranking"].as_bool(true);
    if (!root["entropy_weight"].is_null())
        cfg.entropy_weight = root["entropy_weight"].as_double(0.5);
    if (!root["diversity_monitoring"].is_null())
        cfg.diversity_monitoring = root["diversity_monitoring"].as_bool(true);
    if (!root["use_shannon"].is_null())
        cfg.use_shannon = root["use_shannon"].as_bool(false);
    if (!root["thermo_enabled"].is_null())
        cfg.thermo_enabled = root["thermo_enabled"].as_bool(false);
    if (!root["t_eff"].is_null())
        cfg.t_eff = root["t_eff"].as_float(0.596f);
    if (!root["tencom_scale"].is_null())
        cfg.tencom_scale = root["tencom_scale"].as_float(1.0f);
    if (!root["data_dir"].is_null())
        cfg.data_dir = root["data_dir"].as_string("");
    if (!root["oracle_site_dir"].is_null())
        cfg.oracle_site_dir = root["oracle_site_dir"].as_string("");
    if (!root["oracle_site"].is_null())
        cfg.oracle_site = root["oracle_site"].as_string("");
    if (!root["cleft_sphere_file"].is_null())
        cfg.cleft_sphere_file = root["cleft_sphere_file"].as_string("");
    if (!root["cf_window_selector"].is_null())
        cfg.cf_window_selector = root["cf_window_selector"].as_bool(false);
    if (!root["cluster_member_emit"].is_null())
        cfg.cluster_member_emit = root["cluster_member_emit"].as_bool(false);
    if (!root["seed_elitism"].is_null())
        cfg.seed_elitism = root["seed_elitism"].as_bool(true);
    if (!root["seed_elitism_delta_cf"].is_null())
        cfg.seed_elitism_delta_cf = root["seed_elitism_delta_cf"].as_double(10.0);
    if (!root["freqsel"].is_null())
        cfg.freqsel = root["freqsel"].as_bool(false);
    if (!root["freqsel_alpha"].is_null())
        cfg.freqsel_alpha = root["freqsel_alpha"].as_double(12.0);
    if (!root["freqsel_rmsd"].is_null())
        cfg.freqsel_rmsd = root["freqsel_rmsd"].as_float(1.5f);
    if (!root["consensus_scorer"].is_null())
        cfg.consensus_scorer = root["consensus_scorer"].as_bool(false);
    if (!root["election_v135"].is_null())
        cfg.election_v135 = root["election_v135"].as_bool(false);
    if (!root["election_score_tau"].is_null())
        cfg.election_score_tau = root["election_score_tau"].as_double(0.0);
    if (!root["election_include_singletons"].is_null())
        cfg.election_include_singletons =
            root["election_include_singletons"].as_bool(false);
    if (!root["election_shannon_free_energy"].is_null())
        cfg.election_shannon_free_energy =
            root["election_shannon_free_energy"].as_bool(false);
    if (!root["election_soft_T"].is_null())
        cfg.election_soft_T = root["election_soft_T"].as_double(0.0);
    if (!root["hvib_enabled"].is_null())
        cfg.hvib_enabled = root["hvib_enabled"].as_bool(true);
    if (!root["ring_flex"].is_null())
        cfg.ring_flex = root["ring_flex"].as_bool(false);
    if (!root["eval_scale_dihedral"].is_null())
        cfg.eval_scale_dihedral = root["eval_scale_dihedral"].as_int(1);
    if (!root["budget_scale"].is_null())
        cfg.budget_scale = root["budget_scale"].as_bool(true);
    if (!root["fine_grid"].is_null())
        cfg.fine_grid = root["fine_grid"].as_bool(false);
    if (!root["multi_cleft"].is_null())
        cfg.multi_cleft = root["multi_cleft"].as_int(0);
    if (!root["cognate_site"].is_null())
        cfg.cognate_site = root["cognate_site"].as_bool(false);
    if (!root["score_native"].is_null())
        cfg.score_native = root["score_native"].as_bool(false);
    if (!root["native_only"].is_null())
        cfg.native_only = root["native_only"].as_bool(false);
    if (!root["use_dp"].is_null())
        cfg.use_dp = root["use_dp"].as_bool(false);
    if (!root["ignore_cache"].is_null())
        cfg.ignore_cache = root["ignore_cache"].as_bool(false);
    if (!root["thermo_csv"].is_null())
        cfg.thermo_csv = root["thermo_csv"].as_bool(false);
    if (!root["hbond_weight"].is_null())
        cfg.hbond_weight = root["hbond_weight"].as_double(-2.5);
    if (!root["no_sec"].is_null())
        cfg.no_sec = root["no_sec"].as_bool(false);
    if (!root["benchmark_mode"].is_null())
        cfg.benchmark_mode = root["benchmark_mode"].as_bool(false);
    if (!root["t_hot"].is_null())
        cfg.t_hot = root["t_hot"].as_double(0.0);
    if (!root["instream_interval"].is_null())
        cfg.instream_interval = root["instream_interval"].as_int(0);
    if (!root["chain_norm"].is_null())
        cfg.chain_norm = root["chain_norm"].as_bool(false);
    if (!root["smfree_require_t"].is_null())
        cfg.smfree_require_t = root["smfree_require_t"].as_bool(false);

    return cfg;
}

}  // namespace flexaids
