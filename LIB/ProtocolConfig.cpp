// ProtocolConfig.cpp — env adapter + light JSON for typed protocol knobs
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
// SPDX-License-Identifier: Apache-2.0

#include "ProtocolConfig.h"
#include "json_value.h"

#include <algorithm>
#include <cstdlib>
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
    }

    cfg.sharing_alpha = env_opt_double("FLEXAIDDS_SHARING_ALPHA");
    cfg.boom_frac = env_opt_double("FLEXAIDDS_BOOM_FRAC");
    if (auto n = env_opt_int("FLEXAIDDS_N_ELITE")) {
        cfg.n_elite = *n;
    }
    // Historical DatasetRunner: presence of FLEXAIDDS_USE_SHANNON enables.
    cfg.use_shannon = env_present("FLEXAIDDS_USE_SHANNON");

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

    cfg.cf_window_selector =
        env_truthy_int("FLEXAIDDS_CF_WINDOW_SELECTOR", /*default_value=*/false);
    cfg.cluster_member_emit =
        env_truthy_int("FLEXAIDDS_CLUSTER_MEMBER_EMIT", /*default_value=*/false);

    if (auto v = env_opt_double("FLEXAIDDS_HBOND_WEIGHT")) {
        cfg.hbond_weight = *v;
    }

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
    o << "\"parallel_restarts\":" << (parallel_restarts ? "true" : "false") << ',';
    o << "\"vct_r0\":" << vct_r0 << ',';
    o << "\"vct_normalize_contacts\":"
      << (vct_normalize_contacts ? "true" : "false") << ',';
    o << "\"vct_entropy_weight\":" << vct_entropy_weight << ',';
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
    o << "\"use_shannon\":" << (use_shannon ? "true" : "false") << ',';
    o << "\"thermo_enabled\":" << (thermo_enabled ? "true" : "false") << ',';
    o << "\"t_eff\":" << t_eff << ',';
    o << "\"tencom_scale\":" << tencom_scale << ',';
    o << "\"data_dir\":\"" << json_escape(data_dir) << "\",";
    o << "\"cf_window_selector\":" << (cf_window_selector ? "true" : "false") << ',';
    o << "\"cluster_member_emit\":" << (cluster_member_emit ? "true" : "false") << ',';
    o << "\"hbond_weight\":" << hbond_weight;
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
    if (!root["vct_entropy_weight"].is_null())
        cfg.vct_entropy_weight = root["vct_entropy_weight"].as_double(0.0);
    if (!root["sharing_alpha"].is_null())
        cfg.sharing_alpha = root["sharing_alpha"].as_double(4.0);
    if (!root["boom_frac"].is_null())
        cfg.boom_frac = root["boom_frac"].as_double(1.0);
    if (!root["n_elite"].is_null())
        cfg.n_elite = root["n_elite"].as_int(1);
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
    if (!root["cf_window_selector"].is_null())
        cfg.cf_window_selector = root["cf_window_selector"].as_bool(false);
    if (!root["cluster_member_emit"].is_null())
        cfg.cluster_member_emit = root["cluster_member_emit"].as_bool(false);
    if (!root["hbond_weight"].is_null())
        cfg.hbond_weight = root["hbond_weight"].as_double(-2.5);

    return cfg;
}

}  // namespace flexaids
