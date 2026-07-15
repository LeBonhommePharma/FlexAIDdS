// tests/test_protocol_config.cpp — ProtocolConfig from_env + JSON round-trip
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>

#include "ProtocolConfig.h"
#include "RunReceipt.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <unistd.h>

namespace {

// RAII env override that restores the previous value on scope exit.
class ScopedEnv {
public:
    ScopedEnv(const char* name, const char* value) : name_(name) {
        const char* prev = std::getenv(name);
        if (prev) {
            had_prev_ = true;
            prev_ = prev;
        }
        if (value) {
#if defined(_WIN32)
            _putenv_s(name, value);
#else
            setenv(name, value, 1);
#endif
        } else {
#if defined(_WIN32)
            _putenv_s(name, "");
#else
            unsetenv(name);
#endif
        }
    }

    ~ScopedEnv() {
#if defined(_WIN32)
        if (had_prev_) _putenv_s(name_.c_str(), prev_.c_str());
        else _putenv_s(name_.c_str(), "");
#else
        if (had_prev_) setenv(name_.c_str(), prev_.c_str(), 1);
        else unsetenv(name_.c_str());
#endif
    }

    ScopedEnv(const ScopedEnv&) = delete;
    ScopedEnv& operator=(const ScopedEnv&) = delete;

private:
    std::string name_;
    std::string prev_;
    bool had_prev_{false};
};

// Clear the knobs ProtocolConfig::from_env() reads so tests are hermetic.
struct ClearProtocolEnv {
    ClearProtocolEnv()
        : seed_base("FLEXAIDDS_SEED_BASE", nullptr)
        , restarts("FLEXAIDDS_RESTARTS", nullptr)
        , parallel("FLEXAIDDS_PARALLEL_RESTARTS", nullptr)
        , vct_r0("FLEXAIDDS_VCT_R0", nullptr)
        , vct_norm("FLEXAIDDS_VCT_NORM", nullptr)
        , vct_ew("FLEXAIDDS_VCT_ENTROPY_WEIGHT", nullptr)
        , sharing("FLEXAIDDS_SHARING_ALPHA", nullptr)
        , boom("FLEXAIDDS_BOOM_FRAC", nullptr)
        , n_elite("FLEXAIDDS_N_ELITE", nullptr)
        , shannon("FLEXAIDDS_USE_SHANNON", nullptr)
        , thermo("FLEXAIDDS_THERMO", nullptr)
        , t_eff("FLEXAIDDS_T_EFF", nullptr)
        , tencom("FLEXAIDDS_TENCOM_SCALE", nullptr)
        , data_dir("FLEXAIDDS_DATA_DIR", nullptr)
        , oracle_dir("FLEXAIDDS_ORACLE_SITE_DIR", nullptr)
        , oracle_site("FLEXAIDDS_ORACLE_SITE", nullptr)
        , cleft("FLEXAIDDS_CLEFT_SPHERE_FILE", nullptr)
        , cf_win("FLEXAIDDS_CF_WINDOW_SELECTOR", nullptr)
        , cluster("FLEXAIDDS_CLUSTER_MEMBER_EMIT", nullptr)
        , seed_elit("FLEXAIDDS_SEED_ELITISM", nullptr)
        , seed_delta("FLEXAIDDS_SEED_ELITISM_DELTA_CF", nullptr)
        , freqsel("FLEXAIDDS_FREQSEL", nullptr)
        , freq_a("FLEXAIDDS_FREQSEL_ALPHA", nullptr)
        , freq_r("FLEXAIDDS_FREQSEL_RMSD", nullptr)
        , consensus("FLEXAIDDS_CONSENSUS_SCORER", nullptr)
        , hvib("FLEXAIDDS_HVIB", nullptr)
        , ring("FLEXAIDDS_RING_FLEX", nullptr)
        , eval_scale("FLEXAIDDS_EVAL_SCALE_DIHEDRAL", nullptr)
        , budget("FLEXAIDDS_BUDGET_SCALE", nullptr)
        , fine("FLEXAIDDS_FINE_GRID", nullptr)
        , multi("FLEXAIDDS_MULTI_CLEFT", nullptr)
        , cognate("FLEXAIDDS_COGNATE_SITE", nullptr)
        , score_n("FLEXAIDDS_SCORE_NATIVE", nullptr)
        , native_o("FLEXAIDDS_NATIVE_ONLY", nullptr)
        , use_dp("FLEXAIDDS_USE_DP", nullptr)
        , ignore("FLEXAIDDS_IGNORE_CACHE", nullptr)
        , thermo_csv("FLEXAIDDS_THERMO_CSV", nullptr)
        , hbond("FLEXAIDDS_HBOND_WEIGHT", nullptr)
        , no_sec("FLEXAIDDS_NO_SEC", nullptr)
        , bench("FLEXAIDDS_BENCHMARK", nullptr)
        , t_hot("FLEXAIDDS_T_HOT", nullptr)
        , instream("FLEXAIDDS_INSTREAM_INTERVAL", nullptr)
        , chain("FLEXAIDDS_CHAIN_NORM", nullptr)
        , smfree("FLEXAIDDS_SMFREE_REQUIRE_T", nullptr)
        , force_cf("FLEXAIDDS_FORCE_CF_RANK_EMISSION", nullptr)
        , classic("FLEXAIDDS_CLASSIC_ENTROPY_RANKING", nullptr)
        , ent_w("FLEXAIDDS_ENTROPY_WEIGHT", nullptr)
        , div_m("FLEXAIDDS_DIVERSITY_MONITORING", nullptr) {}

    ScopedEnv seed_base, restarts, parallel, vct_r0, vct_norm, vct_ew,
              sharing, boom, n_elite, shannon, thermo, t_eff, tencom,
              data_dir, oracle_dir, oracle_site, cleft, cf_win, cluster,
              seed_elit, seed_delta, freqsel, freq_a, freq_r, consensus,
              hvib, ring, eval_scale, budget, fine, multi, cognate, score_n,
              native_o, use_dp, ignore, thermo_csv, hbond, no_sec, bench,
              t_hot, instream, chain, smfree, force_cf, classic, ent_w, div_m;
};

}  // namespace

TEST(ProtocolConfig, DefaultsMatchHistoricalFallbacks) {
    ClearProtocolEnv clear;
    const auto d = flexaids::ProtocolConfig::defaults();
    const auto e = flexaids::ProtocolConfig::from_env();

    EXPECT_EQ(e.seed_base, d.seed_base);
    EXPECT_EQ(e.restarts, 5);
    EXPECT_TRUE(e.parallel_restarts);  // restarts default 5 > 1
    EXPECT_DOUBLE_EQ(e.vct_r0, 7.0);
    EXPECT_FALSE(e.vct_normalize_contacts);
    EXPECT_DOUBLE_EQ(e.vct_entropy_weight, 0.0);
    EXPECT_FALSE(e.sharing_alpha.has_value());
    EXPECT_FALSE(e.boom_frac.has_value());
    EXPECT_EQ(e.n_elite, 1);
    EXPECT_FALSE(e.use_shannon);
    EXPECT_FALSE(e.thermo_enabled);
    EXPECT_FLOAT_EQ(e.t_eff, 0.596f);
    EXPECT_FLOAT_EQ(e.tencom_scale, 1.0f);
    EXPECT_TRUE(e.data_dir.empty());
    EXPECT_FALSE(e.cf_window_selector);
    EXPECT_FALSE(e.cluster_member_emit);
    EXPECT_DOUBLE_EQ(e.hbond_weight, -2.5);

    // Chunk 2 defaults
    EXPECT_TRUE(e.seed_elitism);
    EXPECT_DOUBLE_EQ(e.seed_elitism_delta_cf, 10.0);
    EXPECT_FALSE(e.freqsel);
    EXPECT_DOUBLE_EQ(e.freqsel_alpha, 12.0);
    EXPECT_FLOAT_EQ(e.freqsel_rmsd, 1.5f);
    EXPECT_FALSE(e.consensus_scorer);
    EXPECT_FALSE(e.election_v135);
    EXPECT_DOUBLE_EQ(e.election_score_tau, 0.0);
    EXPECT_FALSE(e.election_include_singletons);
    EXPECT_TRUE(e.election_shannon_free_energy);  // 3Dsig default ON
    EXPECT_DOUBLE_EQ(e.election_soft_T, 0.0);
    EXPECT_TRUE(e.hvib_enabled);
    EXPECT_FALSE(e.ring_flex);
    EXPECT_EQ(e.eval_scale_dihedral, 1);
    EXPECT_TRUE(e.budget_scale);
    EXPECT_FALSE(e.fine_grid);
    EXPECT_EQ(e.multi_cleft, 0);
    EXPECT_FALSE(e.cognate_site);
    EXPECT_FALSE(e.score_native);
    EXPECT_FALSE(e.native_only);
    EXPECT_FALSE(e.use_dp);
    EXPECT_FALSE(e.ignore_cache);
    EXPECT_FALSE(e.thermo_csv);
    EXPECT_FALSE(e.no_sec);
    EXPECT_FALSE(e.benchmark_mode);
    EXPECT_DOUBLE_EQ(e.t_hot, 0.0);
    EXPECT_EQ(e.instream_interval, 0);
    EXPECT_FALSE(e.chain_norm);
    EXPECT_FALSE(e.smfree_require_t);

    // Chunk 3 defaults (optional ranking/ablation knobs)
    EXPECT_FALSE(e.n_elite_set);
    EXPECT_FALSE(e.vct_entropy_weight_set);
    EXPECT_FALSE(e.force_cf_rank_emission.has_value());
    EXPECT_FALSE(e.classic_entropy_ranking.has_value());
    EXPECT_FALSE(e.entropy_weight.has_value());
    EXPECT_FALSE(e.diversity_monitoring.has_value());

    EXPECT_DOUBLE_EQ(e.effective_sharing_alpha(1000, 1000), 4.0);
    EXPECT_DOUBLE_EQ(e.effective_sharing_alpha(1000, 2000), 2.0);
    EXPECT_DOUBLE_EQ(e.effective_boom_frac(), 1.0);
}

TEST(ProtocolConfig, FromEnvOverridesKeyVars) {
    ClearProtocolEnv clear;
    ScopedEnv seed("FLEXAIDDS_SEED_BASE", "42");
    ScopedEnv restarts("FLEXAIDDS_RESTARTS", "3");
    ScopedEnv vct("FLEXAIDDS_VCT_R0", "5.5");
    ScopedEnv alpha("FLEXAIDDS_SHARING_ALPHA", "2.5");
    ScopedEnv elite("FLEXAIDDS_N_ELITE", "4");
    ScopedEnv data("FLEXAIDDS_DATA_DIR", "/tmp/flexaidds-data");
    ScopedEnv hbond("FLEXAIDDS_HBOND_WEIGHT", "-3.1");

    const auto cfg = flexaids::ProtocolConfig::from_env();
    EXPECT_EQ(cfg.seed_base, 42u);
    EXPECT_EQ(cfg.restarts, 3);
    EXPECT_DOUBLE_EQ(cfg.vct_r0, 5.5);
    ASSERT_TRUE(cfg.sharing_alpha.has_value());
    EXPECT_DOUBLE_EQ(*cfg.sharing_alpha, 2.5);
    EXPECT_DOUBLE_EQ(cfg.effective_sharing_alpha(1000, 2000), 2.5);
    EXPECT_EQ(cfg.n_elite, 4);
    EXPECT_TRUE(cfg.n_elite_set);
    EXPECT_EQ(cfg.data_dir, "/tmp/flexaidds-data");
    EXPECT_DOUBLE_EQ(cfg.hbond_weight, -3.1);
}

TEST(ProtocolConfig, PresenceFlagsAndParallelRestarts) {
    ClearProtocolEnv clear;
    ScopedEnv norm("FLEXAIDDS_VCT_NORM", "1");
    ScopedEnv shannon("FLEXAIDDS_USE_SHANNON", "1");
    ScopedEnv thermo("FLEXAIDDS_THERMO", "1");
    ScopedEnv teff("FLEXAIDDS_T_EFF", "0.8");
    ScopedEnv tscale("FLEXAIDDS_TENCOM_SCALE", "1.5");
    ScopedEnv cf("FLEXAIDDS_CF_WINDOW_SELECTOR", "1");
    ScopedEnv cl("FLEXAIDDS_CLUSTER_MEMBER_EMIT", "0");
    ScopedEnv restarts("FLEXAIDDS_RESTARTS", "1");
    ScopedEnv par("FLEXAIDDS_PARALLEL_RESTARTS", "1");

    const auto cfg = flexaids::ProtocolConfig::from_env();
    EXPECT_TRUE(cfg.vct_normalize_contacts);
    EXPECT_TRUE(cfg.use_shannon);
    EXPECT_TRUE(cfg.thermo_enabled);
    EXPECT_FLOAT_EQ(cfg.t_eff, 0.8f);
    EXPECT_FLOAT_EQ(cfg.tencom_scale, 1.5f);
    EXPECT_TRUE(cfg.cf_window_selector);
    EXPECT_FALSE(cfg.cluster_member_emit);
    EXPECT_EQ(cfg.restarts, 1);
    // Explicit parallel=1 but restarts==1 → parallel stays false.
    EXPECT_FALSE(cfg.parallel_restarts);
    EXPECT_TRUE(cfg.parallel_restarts_explicit);
}

TEST(ProtocolConfig, Chunk2PoseAndBudgetKnobs) {
    ClearProtocolEnv clear;
    ScopedEnv elit("FLEXAIDDS_SEED_ELITISM", "0");
    ScopedEnv delta("FLEXAIDDS_SEED_ELITISM_DELTA_CF", "7.5");
    ScopedEnv freq("FLEXAIDDS_FREQSEL", "1");
    ScopedEnv fa("FLEXAIDDS_FREQSEL_ALPHA", "9.0");
    ScopedEnv fr("FLEXAIDDS_FREQSEL_RMSD", "2.0");
    ScopedEnv ring("FLEXAIDDS_RING_FLEX", "1");
    ScopedEnv eval("FLEXAIDDS_EVAL_SCALE_DIHEDRAL", "off");
    ScopedEnv budget("FLEXAIDDS_BUDGET_SCALE", "0");
    ScopedEnv fine("FLEXAIDDS_FINE_GRID", "1");
    ScopedEnv multi("FLEXAIDDS_MULTI_CLEFT", "8");
    ScopedEnv hvib("FLEXAIDDS_HVIB", "0");
    ScopedEnv cons("FLEXAIDDS_CONSENSUS_SCORER", "1");
    ScopedEnv v135("FLEXAIDDS_ELECTION_V135", "1");
    ScopedEnv hot("FLEXAIDDS_T_HOT", "1500");
    ScopedEnv nosec("FLEXAIDDS_NO_SEC", "1");
    ScopedEnv bench("FLEXAIDDS_BENCHMARK", "1");

    const auto cfg = flexaids::ProtocolConfig::from_env();
    EXPECT_FALSE(cfg.seed_elitism);
    EXPECT_DOUBLE_EQ(cfg.seed_elitism_delta_cf, 7.5);
    EXPECT_TRUE(cfg.freqsel);
    EXPECT_DOUBLE_EQ(cfg.freqsel_alpha, 9.0);
    EXPECT_FLOAT_EQ(cfg.freqsel_rmsd, 2.0f);
    EXPECT_TRUE(cfg.ring_flex);
    EXPECT_EQ(cfg.eval_scale_dihedral, -1);
    EXPECT_FALSE(cfg.budget_scale);
    EXPECT_TRUE(cfg.fine_grid);
    EXPECT_EQ(cfg.multi_cleft, 8);
    EXPECT_FALSE(cfg.hvib_enabled);
    EXPECT_TRUE(cfg.consensus_scorer);
    EXPECT_TRUE(cfg.election_v135);
    EXPECT_DOUBLE_EQ(cfg.election_score_tau, 25.0);  // CF a.u. default under v135
    EXPECT_TRUE(cfg.election_include_singletons);
    EXPECT_DOUBLE_EQ(cfg.t_hot, 1500.0);
    EXPECT_TRUE(cfg.no_sec);
    EXPECT_TRUE(cfg.benchmark_mode);
}

TEST(ProtocolConfig, ElectionV135TauOverride) {
    ClearProtocolEnv clear;
    ScopedEnv v135("FLEXAIDDS_ELECTION_V135", "1");
    ScopedEnv tau("FLEXAIDDS_ELECTION_SCORE_TAU", "40.0");
    ScopedEnv sing("FLEXAIDDS_ELECTION_INCLUDE_SINGLETONS", "0");
    const auto cfg = flexaids::ProtocolConfig::from_env();
    EXPECT_TRUE(cfg.election_v135);
    EXPECT_DOUBLE_EQ(cfg.election_score_tau, 40.0);
    EXPECT_FALSE(cfg.election_include_singletons);
}

// FLEXAIDDS_ELECTION_SOFT_T=0 (unset) means "resolve at election": dock TEMPER
// first, then 298 fallback. The env knob itself stays 0 in ProtocolConfig.
// Positive values override dock T in select_pose_freq_gated_pooled (source=env).
TEST(ProtocolConfig, ElectionSoftTFromEnv) {
    ClearProtocolEnv clear;
    {
        const auto unset = flexaids::ProtocolConfig::from_env();
        EXPECT_DOUBLE_EQ(unset.election_soft_T, 0.0);
        EXPECT_TRUE(unset.election_shannon_free_energy);  // 3Dsig default ON
    }
    {
        ScopedEnv soft("FLEXAIDDS_ELECTION_SOFT_T", "21.0");
        const auto cfg = flexaids::ProtocolConfig::from_env();
        EXPECT_DOUBLE_EQ(cfg.election_soft_T, 21.0);
        EXPECT_TRUE(cfg.election_shannon_free_energy);
    }
    {
        ScopedEnv soft("FLEXAIDDS_ELECTION_SOFT_T", "298.15");
        const auto cfg = flexaids::ProtocolConfig::from_env();
        EXPECT_DOUBLE_EQ(cfg.election_soft_T, 298.15);
    }
}

TEST(ProtocolConfig, ElectionSoftTJsonRoundTrip) {
    ClearProtocolEnv clear;
    ScopedEnv soft("FLEXAIDDS_ELECTION_SOFT_T", "21.0");
    const auto a = flexaids::ProtocolConfig::from_env();
    EXPECT_DOUBLE_EQ(a.election_soft_T, 21.0);
    const std::string json = a.to_json();
    EXPECT_NE(json.find("\"election_soft_T\":"), std::string::npos);
    const auto b = flexaids::ProtocolConfig::from_json(json);
    EXPECT_DOUBLE_EQ(b.election_soft_T, 21.0);
}

TEST(ProtocolConfig, JsonRoundTrip) {
    ClearProtocolEnv clear;
    ScopedEnv seed("FLEXAIDDS_SEED_BASE", "99");
    ScopedEnv restarts("FLEXAIDDS_RESTARTS", "7");
    ScopedEnv alpha("FLEXAIDDS_SHARING_ALPHA", "3.25");
    ScopedEnv data("FLEXAIDDS_DATA_DIR", "/opt/share/flexaidds");
    ScopedEnv hot("FLEXAIDDS_T_HOT", "900");
    ScopedEnv multi("FLEXAIDDS_MULTI_CLEFT", "4");

    const auto a = flexaids::ProtocolConfig::from_env();
    const std::string json = a.to_json();
    EXPECT_NE(json.find("\"seed_base\":99"), std::string::npos);
    EXPECT_NE(json.find("\"restarts\":7"), std::string::npos);
    EXPECT_NE(json.find("\"t_hot\":"), std::string::npos);
    EXPECT_NE(json.find("\"multi_cleft\":4"), std::string::npos);

    const auto b = flexaids::ProtocolConfig::from_json(json);
    EXPECT_EQ(b.seed_base, a.seed_base);
    EXPECT_EQ(b.restarts, a.restarts);
    ASSERT_TRUE(b.sharing_alpha.has_value());
    EXPECT_DOUBLE_EQ(*b.sharing_alpha, 3.25);
    EXPECT_EQ(b.data_dir, "/opt/share/flexaidds");
    EXPECT_DOUBLE_EQ(b.t_hot, 900.0);
    EXPECT_EQ(b.multi_cleft, 4);
}

TEST(ProtocolConfig, Chunk3RankingAndAblationKnobs) {
    ClearProtocolEnv clear;
    ScopedEnv force("FLEXAIDDS_FORCE_CF_RANK_EMISSION", "1");
    ScopedEnv classic("FLEXAIDDS_CLASSIC_ENTROPY_RANKING", "0");
    ScopedEnv ew("FLEXAIDDS_ENTROPY_WEIGHT", "0.25");
    ScopedEnv div("FLEXAIDDS_DIVERSITY_MONITORING", "0");
    ScopedEnv vctew("FLEXAIDDS_VCT_ENTROPY_WEIGHT", "0.7");
    ScopedEnv elite("FLEXAIDDS_N_ELITE", "3");

    const auto cfg = flexaids::ProtocolConfig::from_env();
    ASSERT_TRUE(cfg.force_cf_rank_emission.has_value());
    EXPECT_TRUE(*cfg.force_cf_rank_emission);
    ASSERT_TRUE(cfg.classic_entropy_ranking.has_value());
    EXPECT_FALSE(*cfg.classic_entropy_ranking);
    ASSERT_TRUE(cfg.entropy_weight.has_value());
    EXPECT_DOUBLE_EQ(*cfg.entropy_weight, 0.25);
    ASSERT_TRUE(cfg.diversity_monitoring.has_value());
    EXPECT_FALSE(*cfg.diversity_monitoring);
    EXPECT_TRUE(cfg.vct_entropy_weight_set);
    EXPECT_DOUBLE_EQ(cfg.vct_entropy_weight, 0.7);
    EXPECT_TRUE(cfg.n_elite_set);
    EXPECT_EQ(cfg.n_elite, 3);
}

TEST(RunReceipt, BuildJsonHasRequiredKeys) {
    ClearProtocolEnv clear;
    flexaids::RunReceiptInput in;
    in.run_id = "unit_test_run";
    in.started_utc = "2026-07-15T00:00:00Z";
    in.output = "/tmp/flexaidds_receipt_test";
    in.dataset = "astex_diverse";
    in.mode = "defined-cleft-redock";
    in.temperature_K = 298.0;
    in.pop = 1000;
    in.gen = 6000;
    in.restarts = 5;
    in.seed_base = 42;
    in.seed_elitism = false;
    in.matrix_path = "/data/MC_st0r5.2_6.dat";
    in.matrix_md5 = "deadbeef";
    in.matrix_sha256 = "cafebabe";
    in.binary_path = "/bin/FlexAIDdS";
    in.binary_sha256 = "0123456789abcdef";
    in.git_commit = "abc123";
    in.protocol = flexaids::ProtocolConfig::from_env();
    in.protocol.seed_base = 42;
    in.protocol.restarts = 5;

    const std::string json = flexaids::build_run_receipt_json(in);
    EXPECT_NE(json.find("\"schema_version\": 1"), std::string::npos);
    EXPECT_NE(json.find("\"run_id\": \"unit_test_run\""), std::string::npos);
    EXPECT_NE(json.find("\"mode\": \"defined-cleft-redock\""), std::string::npos);
    EXPECT_NE(json.find("\"temperature_K\":"), std::string::npos);
    EXPECT_NE(json.find("\"pop\": 1000"), std::string::npos);
    EXPECT_NE(json.find("\"gen\": 6000"), std::string::npos);
    EXPECT_NE(json.find("\"restarts\": 5"), std::string::npos);
    EXPECT_NE(json.find("\"seed_base\": 42"), std::string::npos);
    EXPECT_NE(json.find("\"seed_elitism\": 0"), std::string::npos);
    EXPECT_NE(json.find("\"matrix_md5\": \"deadbeef\""), std::string::npos);
    EXPECT_NE(json.find("\"matrix_sha256\": \"cafebabe\""), std::string::npos);
    EXPECT_NE(json.find("\"binary_sha256\": \"0123456789abcdef\""), std::string::npos);
    EXPECT_NE(json.find("\"git_commit\": \"abc123\""), std::string::npos);
    EXPECT_NE(json.find("\"protocol_config\":"), std::string::npos);
    EXPECT_NE(json.find("\"seed_base\":42"), std::string::npos);
}

TEST(RunReceipt, WriteReceiptFiles) {
    ClearProtocolEnv clear;
    namespace fs = std::filesystem;
    const auto tmp = fs::temp_directory_path() /
        ("flexaidds_receipt_" + std::to_string(::getpid()));
    fs::create_directories(tmp);

    flexaids::RunReceiptInput in;
    in.run_id = "write_test";
    in.started_utc = flexaids::utc_now_iso8601();
    in.output = tmp.string();
    in.dataset = "smoke";
    in.mode = "autonomous";
    in.temperature_K = 300.0;
    in.pop = 500;
    in.gen = 1000;
    in.restarts = 2;
    in.seed_elitism = true;
    in.protocol = flexaids::ProtocolConfig::defaults();

    ASSERT_TRUE(flexaids::write_run_receipt(tmp.string(), in, true));
    EXPECT_TRUE(fs::exists(tmp / "RUN_RECEIPT.json"));
    EXPECT_TRUE(fs::exists(tmp / "provenance.json"));

    std::ifstream rf(tmp / "RUN_RECEIPT.json");
    std::string body((std::istreambuf_iterator<char>(rf)),
                     std::istreambuf_iterator<char>());
    EXPECT_NE(body.find("\"run_id\": \"write_test\""), std::string::npos);
    EXPECT_NE(body.find("\"protocol_config\":"), std::string::npos);

    fs::remove_all(tmp);
}
