// tests/test_protocol_config.cpp — ProtocolConfig from_env + JSON round-trip
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>

#include "ProtocolConfig.h"

#include <cstdlib>
#include <string>

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
        , cf_win("FLEXAIDDS_CF_WINDOW_SELECTOR", nullptr)
        , cluster("FLEXAIDDS_CLUSTER_MEMBER_EMIT", nullptr)
        , hbond("FLEXAIDDS_HBOND_WEIGHT", nullptr) {}

    ScopedEnv seed_base, restarts, parallel, vct_r0, vct_norm, vct_ew,
              sharing, boom, n_elite, shannon, thermo, t_eff, tencom,
              data_dir, cf_win, cluster, hbond;
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

TEST(ProtocolConfig, JsonRoundTrip) {
    ClearProtocolEnv clear;
    ScopedEnv seed("FLEXAIDDS_SEED_BASE", "99");
    ScopedEnv restarts("FLEXAIDDS_RESTARTS", "7");
    ScopedEnv alpha("FLEXAIDDS_SHARING_ALPHA", "3.25");
    ScopedEnv data("FLEXAIDDS_DATA_DIR", "/opt/share/flexaidds");

    const auto a = flexaids::ProtocolConfig::from_env();
    const std::string json = a.to_json();
    EXPECT_NE(json.find("\"seed_base\":99"), std::string::npos);
    EXPECT_NE(json.find("\"restarts\":7"), std::string::npos);

    const auto b = flexaids::ProtocolConfig::from_json(json);
    EXPECT_EQ(b.seed_base, a.seed_base);
    EXPECT_EQ(b.restarts, a.restarts);
    ASSERT_TRUE(b.sharing_alpha.has_value());
    EXPECT_DOUBLE_EQ(*b.sharing_alpha, 3.25);
    EXPECT_EQ(b.data_dir, "/opt/share/flexaidds");
}
