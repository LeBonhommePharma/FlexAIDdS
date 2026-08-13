// flexaidds_flags.cpp — resolve + dump for the unified gate registry
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0

#include "flexaidds_flags.h"

#include <cctype>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace flexaidds {
namespace flags {
namespace {

enum class Kind : unsigned char { Bool, Enum, Value, Compile };

struct Gate {
    const char* name = nullptr;
    Kind kind = Kind::Bool;
    bool requested = false;
    bool active = false;
    std::string value;
    std::string reason;
};

struct Registry {
    std::vector<Gate> gates;
    std::unordered_map<std::string, std::size_t> index;  // normalized name → idx
    std::unordered_map<std::string, std::size_t> alias;  // normalized alias → idx
};

Registry g_reg;
std::mutex g_mu;
bool g_resolved = false;

std::string norm_key(const char* s) {
    std::string o;
    if (!s) return o;
    o.reserve(std::strlen(s));
    for (const char* p = s; *p; ++p) {
        unsigned char c = static_cast<unsigned char>(*p);
        if (c == '-') c = '_';
        o.push_back(static_cast<char>(std::toupper(c)));
    }
    return o;
}

// Truthiness: unset / empty / leading '0' → off; any other non-empty → on.
// Matches existing gates that test `e && e[0] != '\0' && e[0] != '0'`.
bool env_truthy(const char* e) {
    return e && e[0] != '\0' && e[0] != '0';
}

const char* env_raw(const char* name) {
    return std::getenv(name);
}

void add_gate(const char* name, Kind kind) {
    Gate g;
    g.name = name;
    g.kind = kind;
    const std::size_t idx = g_reg.gates.size();
    g_reg.gates.push_back(std::move(g));
    g_reg.index.emplace(norm_key(name), idx);
}

void add_alias(const char* alias, const char* canonical) {
    const auto it = g_reg.index.find(norm_key(canonical));
    if (it == g_reg.index.end()) return;
    g_reg.alias.emplace(norm_key(alias), it->second);
}

Gate* find_gate(const char* name) {
    if (!name || !*name) return nullptr;
    const std::string k = norm_key(name);
    if (auto it = g_reg.index.find(k); it != g_reg.index.end())
        return &g_reg.gates[it->second];
    if (auto it = g_reg.alias.find(k); it != g_reg.alias.end())
        return &g_reg.gates[it->second];
    static const char* kPrefix[] = {"FLEXAIDDS_", "FLEXAID_", "FLEXAIDS_"};
    for (const char* pre : kPrefix) {
        if (auto it = g_reg.index.find(std::string(pre) + k); it != g_reg.index.end())
            return &g_reg.gates[it->second];
    }
    return nullptr;
}

void seed_runtime_gates() {
    // Required / high-traffic gates first so dump() leads with them.
    static const char* kBool[] = {
        "FLEXAIDDS_RIGID_FASTPATH",
        "FLEXAIDDS_HOIST_RECEPTOR_INDEX",
        "FLEXAIDDS_CONTACTS_EPOCH",
        "FLEXAIDDS_GET_YVAL_LUT",
        "FLEXAIDDS_FIXED_ORDER_LSE",
        "FLEXAIDDS_PARALLEL_REPRODUCE",
        "FLEXAIDDS_RNG_STREAM_FIX",
        "FLEXAID_DETERMINISTIC",
        "FLEXAIDDS_FORCE_CPU",
        "FLEXAIDDS_FORCE_RIGID",
        "FLEXAIDDS_MEDOID_REFINE",
        "FLEXAIDDS_CLEFT_SORT",
        "FLEXAIDDS_WAL_COERCIVE",
        "FLEXAIDDS_SOFTCORE_WAL",
        "FLEXAIDDS_NO_SAS",
        "FLEXAIDDS_POSEBUST",
        "FLEXAIDS_SOA_ASSERT",
        "FLEXAIDDS_FLAGS_DUMP",
        // remaining LIB/ getenv + ProtocolConfig switches
        "FLEXAIDS_VH_DEBUG",
        "FLEXAIDDS_ADAPTIVE_GENERATIONS",
        "FLEXAIDDS_BASIN_REINJECT",
        "FLEXAIDDS_BENCHMARK",
        "FLEXAIDDS_BUDGET_SCALE",
        "FLEXAIDDS_CF_WINDOW_SELECTOR",
        "FLEXAIDDS_CHAIN_NORM",
        "FLEXAIDDS_CLASSIC_ENTROPY_RANKING",
        "FLEXAIDDS_CLUSTER_MEMBER_EMIT",
        "FLEXAIDDS_COGNATE_SITE",
        "FLEXAIDDS_CONSENSUS_SCORER",
        "FLEXAIDDS_DEBUG_TYPES",
        "FLEXAIDDS_DIST_WEIGHT_CON",
        "FLEXAIDDS_DIVERSITY_MONITORING",
        "FLEXAIDDS_DUMP_POP",
        "FLEXAIDDS_ELECT_LEGACY_ACF",
        "FLEXAIDDS_ELECTION_INCLUDE_SINGLETONS",
        "FLEXAIDDS_ELECTION_LEGACY_ZH",
        "FLEXAIDDS_ELECTION_SHANNON_F",
        "FLEXAIDDS_ELECTION_V135",
        "FLEXAIDDS_FINE_GRID",
        "FLEXAIDDS_FORCE_CF_RANK_EMISSION",
        "FLEXAIDDS_FRAME_CHART_STRICT",
        "FLEXAIDDS_FREQSEL",
        "FLEXAIDDS_GENTRACE",
        "FLEXAIDDS_HBOND_RANK",
        "FLEXAIDDS_HVIB",
        "FLEXAIDDS_IGNORE_CACHE",
        "FLEXAIDDS_MEMETIC",
        "FLEXAIDDS_MUTATION_GRANULAR",
        "FLEXAIDDS_NAN_RANK_GUARD",
        "FLEXAIDDS_NATIVE_ONLY",
        "FLEXAIDDS_NICHE_CARTESIAN",
        "FLEXAIDDS_NO_SEC",
        "FLEXAIDDS_NO_TENCOM",
        "FLEXAIDDS_PARALLEL_RESTARTS",
        "FLEXAIDDS_PB_AWARE_PROMOTION",
        "FLEXAIDDS_PB_CLASH_DEBUG",
        "FLEXAIDDS_PB_CLASH_PHASE2_PASS",
        "FLEXAIDDS_PB_METAL_CARVEOUT",
        "FLEXAIDDS_PB_POCKET_DEBUG",
        "FLEXAIDDS_PB_VDW_CACHED",
        "FLEXAIDDS_PHENOTYPE_UNIQUE",
        "FLEXAIDDS_POSEBUSTERS_REQUIRE_CLI",
        "FLEXAIDDS_RECEPTOR_ROTAMER_PREP",
        "FLEXAIDDS_RING_FLEX",
        "FLEXAIDDS_SCORE_NATIVE",
        "FLEXAIDDS_SEED_ELITISM",
        "FLEXAIDDS_SHANNON_ROBUST",
        "FLEXAIDDS_SMFREE_REQUIRE_T",
        "FLEXAIDDS_SOFTBETA_ELECTION",
        "FLEXAIDDS_THERMO",
        "FLEXAIDDS_THERMO_CSV",
        "FLEXAIDDS_THERMO_SCORE",
        "FLEXAIDDS_USE_DP",
        "FLEXAIDDS_USE_ELEC",
        "FLEXAIDDS_USE_SHANNON",
        "FLEXAIDDS_VCT_NORM",
        "FLEXAIDDS_WALL_PILOT_PASS",
    };
    static const char* kEnum[] = {
        "FLEXAIDDS_CLUSTER_REP",
        "FLEXAIDDS_SEARCH",
        "FLEXAIDDS_POSEBUST_BACKEND",
        "FLEXAIDDS_NEW_SEARCH_ARCH",
        "FLEXAIDDS_FLAGS",
    };
    static const char* kValue[] = {
        "FLEXAID_SEED",
        "FLEXAIDDS_WAL_STIFF",
        "FLEXAIDDS_SOFTCORE_FLOOR",
        "FLEXAIDDS_CON_R0",
        "FLEXAIDDS_COM_FLOOR",
        "FLEXAIDDS_POLAR_DESOLV_WEIGHT",
        "FLEXAIDDS_PB_CLASH_WEIGHT",
        "FLEXAIDDS_PB_CLASH_EXP",
        "FLEXAIDDS_PB_CLASH_RATIO",
        "FLEXAIDDS_PB_POCKET_WEIGHT",
        "FLEXAIDDS_PB_POCKET_RADIUS",
        "FLEXAIDDS_PB_CLASH_ELECT_WEIGHT",
        "FLEXAIDDS_CMAES_MAX_EVALS",
        "FLEXAIDDS_ADAPTIVE_EPS",
        "FLEXAIDDS_BASIN_SIGMA_ANG",
        "FLEXAIDDS_BINARY",
        "FLEXAIDDS_BOOM_FRAC",
        "FLEXAIDDS_BOOM_INTERVAL",
        "FLEXAIDDS_BUILD",
        "FLEXAIDDS_CLEFT_SPHERE_FILE",
        "FLEXAIDDS_CLUSTER_CONSENSUS_K",
        "FLEXAIDDS_CLUSTER_CONSENSUS_TAU",
        "FLEXAIDDS_CLUSTER_POCKET_RADIUS",
        "FLEXAIDDS_CLUSTER_POP_MIN_FRACTION",
        "FLEXAIDDS_CLUSTER_SPREAD_MAX",
        "FLEXAIDDS_COARSE_GRID_STEP",
        "FLEXAIDDS_COARSE_ORIENTATIONS",
        "FLEXAIDDS_DATA_DIR",
        "FLEXAIDDS_ELECTION_SCORE_TAU",
        "FLEXAIDDS_ELECTION_SOFT_T",
        "FLEXAIDDS_ENTROPY_WEIGHT",
        "FLEXAIDDS_EVAL_SCALE_DIHEDRAL",
        "FLEXAIDDS_FREQSEL_ALPHA",
        "FLEXAIDDS_FREQSEL_RMSD",
        "FLEXAIDDS_GENTRACE_EVERY",
        "FLEXAIDDS_GRID_CACHE_DIR",
        "FLEXAIDDS_HBOND_WEIGHT",
        "FLEXAIDDS_HTTP_RETRIES",
        "FLEXAIDDS_INCHI_BIN",
        "FLEXAIDDS_INSTREAM_INTERVAL",
        "FLEXAIDDS_LIGAND_BATCH",
        "FLEXAIDDS_MAX_CONCURRENT_RESTARTS",
        "FLEXAIDDS_MULTI_CLEFT",
        "FLEXAIDDS_N_ELITE",
        "FLEXAIDDS_NICHE_SIGMA_ANG",
        "FLEXAIDDS_ORACLE_SITE",
        "FLEXAIDDS_ORACLE_SITE_DIR",
        "FLEXAIDDS_POSEBUSTERS_BIN",
        "FLEXAIDDS_PRIORITY_TARGETS",
        "FLEXAIDDS_REDOCK_CACHE",
        "FLEXAIDDS_REPO",
        "FLEXAIDDS_REPORT_T",
        "FLEXAIDDS_RESTARTS",
        "FLEXAIDDS_RMSDST",
        "FLEXAIDDS_ROOT",
        "FLEXAIDDS_SAS_WEIGHT",
        "FLEXAIDDS_SEED_BASE",
        "FLEXAIDDS_SEED_ELITISM_DELTA_CF",
        "FLEXAIDDS_SHARING_ALPHA",
        "FLEXAIDDS_SIGMA_SCALE",
        "FLEXAIDDS_T_EFF",
        "FLEXAIDDS_T_HOT",
        "FLEXAIDDS_TENCOM_SCALE",
        "FLEXAIDDS_VCT_ENTROPY_WEIGHT",
        "FLEXAIDDS_VCT_R0",
        "FLEXAIDDS_VIB_ENTROPY_BINS",
    };

    g_reg.gates.clear();
    g_reg.index.clear();
    g_reg.alias.clear();
    g_reg.gates.reserve(200);

    for (const char* n : kBool) add_gate(n, Kind::Bool);
    for (const char* n : kEnum) add_gate(n, Kind::Enum);
    for (const char* n : kValue) add_gate(n, Kind::Value);

    add_alias("hoist", "FLEXAIDDS_HOIST_RECEPTOR_INDEX");
    add_alias("epoch", "FLEXAIDDS_CONTACTS_EPOCH");
    add_alias("fastpath", "FLEXAIDDS_RIGID_FASTPATH");
    add_alias("rng-stream-fix", "FLEXAIDDS_RNG_STREAM_FIX");
    add_alias("rng_stream_fix", "FLEXAIDDS_RNG_STREAM_FIX");
}

void add_compile(const char* name, bool on, const char* note) {
    add_gate(name, Kind::Compile);
    Gate* g = find_gate(name);
    if (!g) return;
    g->requested = on;
    g->active = on;
    g->value = on ? "1" : "0";
    if (note && *note) g->reason = note;  // documentation; cleared if active
    if (on) g->reason.clear();
    else if (note && *note) g->reason = note;
}

void seed_compile_gates() {
#ifdef FLEXAIDS_USE_AVX2
    add_compile("FLEXAIDS_USE_AVX2", true, "cmake -DFLEXAIDS_USE_AVX2");
#else
    add_compile("FLEXAIDS_USE_AVX2", false, "cmake -DFLEXAIDS_USE_AVX2 (read-only)");
#endif
#ifdef FLEXAIDS_USE_AVX512
    add_compile("FLEXAIDS_USE_AVX512", true, "cmake -DFLEXAIDS_USE_AVX512");
#else
    add_compile("FLEXAIDS_USE_AVX512", false, "cmake -DFLEXAIDS_USE_AVX512 (read-only)");
#endif
#ifdef FLEXAIDS_USE_CUDA
    add_compile("FLEXAIDS_USE_CUDA", true, "cmake -DFLEXAIDS_USE_CUDA");
#else
    add_compile("FLEXAIDS_USE_CUDA", false, "cmake -DFLEXAIDS_USE_CUDA (read-only)");
#endif
#ifdef FLEXAIDS_USE_METAL
    add_compile("FLEXAIDS_USE_METAL", true, "cmake -DFLEXAIDS_USE_METAL");
#else
    add_compile("FLEXAIDS_USE_METAL", false, "cmake -DFLEXAIDS_USE_METAL (read-only)");
#endif
#ifdef FLEXAIDS_USE_ROCM
    add_compile("FLEXAIDS_USE_ROCM", true, "cmake -DFLEXAIDS_USE_ROCM");
#else
    add_compile("FLEXAIDS_USE_ROCM", false, "cmake -DFLEXAIDS_USE_ROCM (read-only)");
#endif
#ifdef FLEXAIDS_USE_WEBGPU
    add_compile("FLEXAIDS_USE_WEBGPU", true, "cmake -DFLEXAIDS_USE_WEBGPU");
#else
    add_compile("FLEXAIDS_USE_WEBGPU", false, "cmake -DFLEXAIDS_USE_WEBGPU (read-only)");
#endif
#ifdef FLEXAIDS_USE_OPENMP
    add_compile("FLEXAIDS_USE_OPENMP", true, "cmake -DFLEXAIDS_USE_OPENMP");
#else
    add_compile("FLEXAIDS_USE_OPENMP", false, "cmake -DFLEXAIDS_USE_OPENMP (read-only)");
#endif
#ifdef FLEXAIDS_USE_EIGEN
    add_compile("FLEXAIDS_USE_EIGEN", true, "cmake -DFLEXAIDS_USE_EIGEN");
#else
    add_compile("FLEXAIDS_USE_EIGEN", false, "cmake -DFLEXAIDS_USE_EIGEN (read-only)");
#endif
#ifdef FLEXAIDS_USE_MPI
    add_compile("FLEXAIDS_USE_MPI", true, "cmake -DFLEXAIDS_USE_MPI");
#else
    add_compile("FLEXAIDS_USE_MPI", false, "cmake -DFLEXAIDS_USE_MPI (read-only)");
#endif
#ifdef FLEXAIDS_USE_NEON
    add_compile("FLEXAIDS_USE_NEON", true, "cmake -DFLEXAIDS_USE_NEON");
#else
    add_compile("FLEXAIDS_USE_NEON", false, "cmake -DFLEXAIDS_USE_NEON (read-only)");
#endif
#ifdef FLEXAIDS_USE_256_MATRIX
    add_compile("FLEXAIDS_USE_256_MATRIX", true, "cmake -DFLEXAIDS_USE_256_MATRIX");
#else
    add_compile("FLEXAIDS_USE_256_MATRIX", false, "cmake -DFLEXAIDS_USE_256_MATRIX (read-only)");
#endif
#ifdef FLEXAIDS_USE_SOA_DISTANCES
    add_compile("FLEXAIDS_USE_SOA_DISTANCES", true, "cmake -DFLEXAIDS_USE_SOA_DISTANCES");
#else
    add_compile("FLEXAIDS_USE_SOA_DISTANCES", false,
                "cmake -DFLEXAIDS_USE_SOA_DISTANCES (read-only)");
#endif
#ifdef BUILD_FLEXAIDDS_FAST
    add_compile("BUILD_FLEXAIDDS_FAST", true, "cmake -DBUILD_FLEXAIDDS_FAST (read-only)");
#else
    add_compile("BUILD_FLEXAIDDS_FAST", false, "cmake -DBUILD_FLEXAIDDS_FAST (read-only)");
#endif
#ifdef FLEXAID_PGO
    add_compile("FLEXAID_PGO", true, "cmake -DFLEXAID_PGO (read-only)");
#else
    add_compile("FLEXAID_PGO", false, "cmake -DFLEXAID_PGO=off|generate|use (read-only)");
#endif
}

void warn_once_line(const char* msg) {
    std::fprintf(stderr, "WARNING: %s\n", msg);
}

void apply_overlay_token(const char* token, std::size_t n) {
    if (!token || n == 0) return;
    std::string t(token, n);
    Gate* g = find_gate(t.c_str());
    if (!g) {
        std::fprintf(stderr,
                     "WARNING: FLEXAIDDS_FLAGS token '%s' is not a known gate; ignored.\n",
                     t.c_str());
        return;
    }
    if (g->kind == Kind::Compile) return;  // overlay does not flip compile defs
    g->requested = true;
    if (g->value.empty()) g->value = "1";
}

void apply_flags_overlay(const char* list) {
    if (!list || !*list) return;
    // A lone "0" is off (no tokens). Any other non-empty list is a token list.
    if (list[0] == '0' && list[1] == '\0') return;
    const char* p = list;
    while (*p) {
        while (*p == ' ' || *p == ',' || *p == ';' || *p == '|' || *p == '\t') ++p;
        if (!*p) break;
        const char* start = p;
        while (*p && *p != ' ' && *p != ',' && *p != ';' && *p != '|' && *p != '\t') ++p;
        apply_overlay_token(start, static_cast<std::size_t>(p - start));
    }
}

void read_env_into_gates() {
    for (Gate& g : g_reg.gates) {
        if (g.kind == Kind::Compile) continue;
        const char* e = env_raw(g.name);
        if (!e) continue;
        g.value = e;
        switch (g.kind) {
            case Kind::Bool:
                g.requested = env_truthy(e);
                break;
            case Kind::Enum:
            case Kind::Value:
                g.requested = (e[0] != '\0');
                break;
            case Kind::Compile:
                break;
        }
    }

#ifdef FLEXAID_DETERMINISTIC
    if (Gate* d = find_gate("FLEXAID_DETERMINISTIC")) {
        d->requested = true;
        if (d->value.empty()) d->value = "1";
    }
#endif
}

void disable_loser(Gate* loser, const char* why) {
    if (!loser) return;
    loser->active = false;
    loser->reason = why;
    if (loser->requested) warn_once_line(why);
}

void apply_implications_and_exclusions() {
    for (Gate& g : g_reg.gates) {
        if (g.kind == Kind::Compile) continue;
        g.active = g.requested;
        if (g.active) g.reason.clear();
    }

    // FLEXAIDDS_RIGID_FASTPATH implies hoist; both stay active (superset).
    Gate* fast = find_gate("FLEXAIDDS_RIGID_FASTPATH");
    Gate* hoist = find_gate("FLEXAIDDS_HOIST_RECEPTOR_INDEX");
    if (fast && hoist && fast->requested) {
        fast->active = true;
        hoist->active = true;
        if (!hoist->requested && hoist->reason.empty())
            hoist->reason = "implied by FLEXAIDDS_RIGID_FASTPATH";
    }

    // CLUSTER_REP (any explicit value) wins over MEDOID_REFINE.
    Gate* crep = find_gate("FLEXAIDDS_CLUSTER_REP");
    Gate* med = find_gate("FLEXAIDDS_MEDOID_REFINE");
    if (crep && med && crep->requested && med->requested) {
        disable_loser(med,
                      "FLEXAIDDS_MEDOID_REFINE disabled by FLEXAIDDS_CLUSTER_REP "
                      "(ClusterRepMode.h policy)");
    }

    // CLEFT_SORT is superseded: stay in the registry, become inactive.
    Gate* cleft = find_gate("FLEXAIDDS_CLEFT_SORT");
    if (cleft && cleft->requested) {
        disable_loser(cleft, "FLEXAIDDS_CLEFT_SORT superseded, ignored");
    }

    // Two WAL models: coercive is more specific and wins if both are set.
    Gate* wal_c = find_gate("FLEXAIDDS_WAL_COERCIVE");
    Gate* wal_s = find_gate("FLEXAIDDS_SOFTCORE_WAL");
    if (wal_c && wal_s && wal_c->requested && wal_s->requested) {
        wal_c->active = true;
        disable_loser(wal_s,
                      "FLEXAIDDS_SOFTCORE_WAL disabled by FLEXAIDDS_WAL_COERCIVE "
                      "(two WAL models; coercive is more specific)");
    }

    // FLEXAIDDS_SEARCH=cmaes is one backend; other SEARCH values stay requested
    // but are not the cmaes backend. The SEARCH gate itself remains active so
    // the chosen backend can be read via value().
    Gate* search = find_gate("FLEXAIDDS_SEARCH");
    if (search && search->requested) {
        search->active = true;
        // Document non-cmaes / non-ga leftovers as "not a recognised backend"
        // without dropping the flag from the registry.
        const std::string& v = search->value;
        const bool cmaes = (v == "cmaes" || v == "CMAES" || v == "Cmaes");
        const bool ga = (v.empty() || v == "ga" || v == "GA" || v == "Ga");
        if (!cmaes && !ga) {
            search->reason = "FLEXAIDDS_SEARCH: unrecognised backend (engine keeps GA); "
                             "cmaes is the only opt-in backend";
        }
    }

    // FORCE_CPU: runtime ablation of GPU dispatch (compile flags stay requested).
    Gate* cpu = find_gate("FLEXAIDDS_FORCE_CPU");
    if (cpu && cpu->requested) {
        static const char* kGpu[] = {
            "FLEXAIDS_USE_CUDA",
            "FLEXAIDS_USE_METAL",
            "FLEXAIDS_USE_ROCM",
            "FLEXAIDS_USE_WEBGPU",
        };
        for (const char* n : kGpu) {
            Gate* gpu = find_gate(n);
            if (gpu && gpu->requested) {
                disable_loser(gpu,
                              "GPU dispatch disabled by FLEXAIDDS_FORCE_CPU "
                              "(runtime force; compile -D flag retained)");
            }
        }
    }

    // FLEXAID_DETERMINISTIC does NOT disable PARALLEL_REPRODUCE (different axes).
}

void dump_locked(FILE* out) {
    if (!out) out = stderr;
    std::fprintf(out, "# FlexAIDdS flag registry\n");
    std::fprintf(out, "# name  requested  active  value  reason\n");
    for (const Gate& g : g_reg.gates) {
        std::fprintf(out, "%-36s  requested=%d  active=%d  value=%s  %s\n",
                     g.name, g.requested ? 1 : 0, g.active ? 1 : 0,
                     g.value.empty() ? "-" : g.value.c_str(),
                     g.reason.empty() ? "-" : g.reason.c_str());
    }
    std::fprintf(out,
                 "# note: compile AVX512 vs AVX2 is resolved by cmake "
                 "(AVX512 wins when both requested); both stay in the registry.\n");
}

void do_resolve_locked() {
    seed_runtime_gates();
    seed_compile_gates();
    read_env_into_gates();
    apply_flags_overlay(env_raw("FLEXAIDDS_FLAGS"));
    apply_implications_and_exclusions();

    Gate* dump_flag = find_gate("FLEXAIDDS_FLAGS_DUMP");
    if (dump_flag && dump_flag->requested) dump_locked(stderr);
}

}  // namespace

void resolve_once() {
    std::lock_guard<std::mutex> lock(g_mu);
    if (g_resolved) return;
    do_resolve_locked();
    g_resolved = true;
}

void reset_for_tests() {
    std::lock_guard<std::mutex> lock(g_mu);
    g_resolved = false;
    g_reg.gates.clear();
    g_reg.index.clear();
    g_reg.alias.clear();
}

bool requested(const char* name) {
    resolve_once();
    std::lock_guard<std::mutex> lock(g_mu);
    if (const Gate* g = find_gate(name)) return g->requested;
    return env_truthy(env_raw(name));
}

bool active(const char* name) {
    resolve_once();
    std::lock_guard<std::mutex> lock(g_mu);
    if (const Gate* g = find_gate(name)) return g->active;
    return env_truthy(env_raw(name));
}

const char* value(const char* name) {
    resolve_once();
    std::lock_guard<std::mutex> lock(g_mu);
    if (const Gate* g = find_gate(name)) return g->value.c_str();
    return "";
}

const char* reason(const char* name) {
    resolve_once();
    std::lock_guard<std::mutex> lock(g_mu);
    if (const Gate* g = find_gate(name)) return g->reason.c_str();
    return "";
}

void dump(FILE* out) {
    resolve_once();
    std::lock_guard<std::mutex> lock(g_mu);
    dump_locked(out);
}

void apply_to_environ() {
    resolve_once();
    std::lock_guard<std::mutex> lock(g_mu);
    for (const Gate& g : g_reg.gates) {
        if (g.kind == Kind::Compile) continue;
        if (std::strcmp(g.name, "FLEXAIDDS_FLAGS") == 0) continue;
        if (std::strcmp(g.name, "FLEXAIDDS_FLAGS_DUMP") == 0) continue;
        if (g.active) {
            if (!env_raw(g.name) || g.kind == Kind::Bool) {
                // Overlay / implication: publish a value existing getenv() sites see.
                // Do not overwrite a caller-supplied non-empty value.
                if (!env_raw(g.name) || env_raw(g.name)[0] == '\0') {
                    const char* v = g.value.empty() ? "1" : g.value.c_str();
                    setenv(g.name, v, 0);
                }
            }
        } else if (g.requested && (g.kind == Kind::Bool || g.kind == Kind::Enum)) {
            // Mutual-exclusion loser or superseded gate: hide from getenv()
            // so legacy call sites follow the winner. The flag stays in the
            // registry (requested() remains true).
            unsetenv(g.name);
        }
    }
}

}  // namespace flags
}  // namespace flexaidds
