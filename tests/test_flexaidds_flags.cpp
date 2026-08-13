// tests/test_flexaidds_flags.cpp
// Unified FLEXAIDDS_* / FLEXAID_* / FLEXAIDS_* gate overlay.
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0

#include "flexaidds_flags.h"
#include "EnvFlags.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <gtest/gtest.h>
#include <string>
#include <utility>
#include <vector>

namespace {

void clear_env(const char* key) {
#if defined(_WIN32)
    _putenv_s(key, "");
#else
    unsetenv(key);
#endif
}

void set_env(const char* key, const char* val) {
#if defined(_WIN32)
    _putenv_s(key, val);
#else
    setenv(key, val, 1);
#endif
}

void resolve_env() {
    flexaidds::flags::reset_for_tests();
    flexaidds::flags::resolve_once();
}

class FlexaiddsFlagsEnv : public ::testing::Test {
protected:
    void SetUp() override { clear_all(); }
    void TearDown() override { clear_all(); }

    static void clear_all() {
        static const char* kKeys[] = {
            "FLEXAIDDS_RIGID_FASTPATH",
            "FLEXAIDDS_HOIST_RECEPTOR_INDEX",
            "FLEXAIDDS_CONTACTS_EPOCH",
            "FLEXAIDDS_PARALLEL_REPRODUCE",
            "FLEXAIDDS_RNG_STREAM_FIX",
            "FLEXAIDDS_VORONOI_KEYED_JITTER",
            "FLEXAIDDS_GET_YVAL_LUT",
            "FLEXAIDDS_CA_REC_FLAT",
            "FLEXAIDDS_NICHE_HASH",
            "FLEXAIDDS_FIXED_ORDER_LSE",
            "FLEXAID_DETERMINISTIC",
            "FLEXAID_SEED",
            "FLEXAIDDS_FORCE_CPU",
            "FLEXAIDDS_FORCE_RIGID",
            "FLEXAIDDS_SEARCH",
            "FLEXAIDDS_CLUSTER_REP",
            "FLEXAIDDS_MEDOID_REFINE",
            "FLEXAIDDS_CLEFT_SORT",
            "FLEXAIDDS_WAL_COERCIVE",
            "FLEXAIDDS_SOFTCORE_WAL",
            "FLEXAIDDS_WAL_STIFF",
            "FLEXAIDDS_NO_SAS",
            "FLEXAIDDS_POSEBUST",
            "FLEXAIDDS_POSEBUST_BACKEND",
            "FLEXAIDS_SOA_ASSERT",
            "FLEXAIDDS_FLAGS",
            "FLEXAIDDS_FLAGS_DUMP",
        };
        for (const char* k : kKeys) clear_env(k);
        flexaidds::flags::reset_for_tests();
    }
};

std::string dump_to_string() {
    FILE* f = std::tmpfile();
    EXPECT_NE(f, nullptr);
    if (!f) return {};
    flexaidds::flags::dump(f);
    std::rewind(f);
    std::string out;
    char buf[512];
    while (std::fgets(buf, sizeof(buf), f)) out += buf;
    std::fclose(f);
    return out;
}

}  // namespace

TEST_F(FlexaiddsFlagsEnv, TruthTableOffFalseNoNeverEnable) {
    // B2: one parser through env_bool, flags::active (env_truthy), and the
    // engine reader flexaids::rigid_fastpath_requested().
    const std::vector<std::pair<const char*, bool>> kCases = {
        {"off", false}, {"false", false}, {"no", false}, {"0", false},
        {"", false},
        {"1", true}, {"true", true}, {"on", true}, {"yes", true},
        {"OFF", false}, {"False", false}, {"NO", false},
        {"TRUE", true}, {"On", true}, {"Yes", true},
    };
    for (const auto& [val, expect_on] : kCases) {
        if (val[0] == '\0') clear_env("FLEXAIDDS_RIGID_FASTPATH");
        else set_env("FLEXAIDDS_RIGID_FASTPATH", val);
        resolve_env();
        EXPECT_EQ(flexaids::env_bool_str(val, false), expect_on)
            << "env_bool_str(" << val << ")";
        EXPECT_EQ(flexaids::env_bool("FLEXAIDDS_RIGID_FASTPATH", false), expect_on)
            << "env_bool(" << val << ")";
        EXPECT_EQ(flexaidds::flags::active("FLEXAIDDS_RIGID_FASTPATH"), expect_on)
            << "flags::active(" << val << ")";
        EXPECT_EQ(flexaids::rigid_fastpath_requested(), expect_on)
            << "rigid_fastpath_requested(" << val << ")";
        // Overlay parser is shared: another acceleration gate must agree.
        if (val[0] == '\0') clear_env("FLEXAIDDS_GET_YVAL_LUT");
        else set_env("FLEXAIDDS_GET_YVAL_LUT", val);
        resolve_env();
        EXPECT_EQ(flexaidds::flags::active("FLEXAIDDS_GET_YVAL_LUT"), expect_on)
            << "GET_YVAL_LUT flags::active(" << val << ")";
        clear_env("FLEXAIDDS_GET_YVAL_LUT");
    }
}

TEST_F(FlexaiddsFlagsEnv, OverlayApplyReachesEngineReaders) {
    // B1: no individual env; overlay implies fastpath+hoist; after
    // apply_to_environ the live engine readers must see ON.
    set_env("FLEXAIDDS_FLAGS", "fastpath");
    resolve_env();
    flexaidds::flags::apply_to_environ();
    EXPECT_TRUE(flexaidds::flags::active("FLEXAIDDS_RIGID_FASTPATH"));
    EXPECT_TRUE(flexaidds::flags::active("FLEXAIDDS_HOIST_RECEPTOR_INDEX"));
    ASSERT_NE(std::getenv("FLEXAIDDS_RIGID_FASTPATH"), nullptr);
    ASSERT_NE(std::getenv("FLEXAIDDS_HOIST_RECEPTOR_INDEX"), nullptr);
    EXPECT_TRUE(flexaids::rigid_fastpath_requested());
    EXPECT_TRUE(flexaids::hoist_receptor_index_env());
    EXPECT_TRUE(flexaids::env_bool("FLEXAIDDS_RIGID_FASTPATH", false));
}

TEST_F(FlexaiddsFlagsEnv, VoronoiKeyedJitterDefaultOffIndependentOfRngFix) {
    resolve_env();
    EXPECT_FALSE(flexaidds::flags::active("FLEXAIDDS_VORONOI_KEYED_JITTER"));
    set_env("FLEXAIDDS_RNG_STREAM_FIX", "1");
    resolve_env();
    EXPECT_TRUE(flexaidds::flags::active("FLEXAIDDS_RNG_STREAM_FIX"));
    EXPECT_FALSE(flexaidds::flags::active("FLEXAIDDS_VORONOI_KEYED_JITTER"));
}

TEST_F(FlexaiddsFlagsEnv, DefaultRigidFastpathInactive) {
    resolve_env();
    EXPECT_FALSE(flexaidds::flags::requested("FLEXAIDDS_RIGID_FASTPATH"));
    EXPECT_FALSE(flexaidds::flags::active("FLEXAIDDS_RIGID_FASTPATH"));
    EXPECT_FALSE(flexaidds::flags::rigid_fastpath());
    EXPECT_FALSE(flexaidds::flags::active("FLEXAIDDS_HOIST_RECEPTOR_INDEX"));
}

TEST_F(FlexaiddsFlagsEnv, RigidFastpathImpliesHoist) {
    set_env("FLEXAIDDS_RIGID_FASTPATH", "1");
    resolve_env();
    EXPECT_TRUE(flexaidds::flags::requested("FLEXAIDDS_RIGID_FASTPATH"));
    EXPECT_TRUE(flexaidds::flags::active("FLEXAIDDS_RIGID_FASTPATH"));
    EXPECT_TRUE(flexaidds::flags::rigid_fastpath());
    EXPECT_FALSE(flexaidds::flags::requested("FLEXAIDDS_HOIST_RECEPTOR_INDEX"));
    EXPECT_TRUE(flexaidds::flags::active("FLEXAIDDS_HOIST_RECEPTOR_INDEX"));
    EXPECT_TRUE(flexaidds::flags::hoist_receptor_index());
}

TEST_F(FlexaiddsFlagsEnv, WalCoerciveWinsOverSoftcore) {
    set_env("FLEXAIDDS_WAL_COERCIVE", "1");
    set_env("FLEXAIDDS_SOFTCORE_WAL", "1");
    resolve_env();
    EXPECT_TRUE(flexaidds::flags::requested("FLEXAIDDS_WAL_COERCIVE"));
    EXPECT_TRUE(flexaidds::flags::requested("FLEXAIDDS_SOFTCORE_WAL"));
    EXPECT_TRUE(flexaidds::flags::active("FLEXAIDDS_WAL_COERCIVE"));
    EXPECT_FALSE(flexaidds::flags::active("FLEXAIDDS_SOFTCORE_WAL"));
    EXPECT_NE(std::strstr(flexaidds::flags::reason("FLEXAIDDS_SOFTCORE_WAL"),
                          "WAL_COERCIVE"),
              nullptr);
}

TEST_F(FlexaiddsFlagsEnv, ClusterRepWinsOverMedoidRefine) {
    set_env("FLEXAIDDS_CLUSTER_REP", "lowcf");
    set_env("FLEXAIDDS_MEDOID_REFINE", "1");
    resolve_env();
    EXPECT_TRUE(flexaidds::flags::requested("FLEXAIDDS_CLUSTER_REP"));
    EXPECT_TRUE(flexaidds::flags::active("FLEXAIDDS_CLUSTER_REP"));
    EXPECT_TRUE(flexaidds::flags::requested("FLEXAIDDS_MEDOID_REFINE"));
    EXPECT_FALSE(flexaidds::flags::active("FLEXAIDDS_MEDOID_REFINE"));
}

TEST_F(FlexaiddsFlagsEnv, CleftSortSuperseded) {
    set_env("FLEXAIDDS_CLEFT_SORT", "1");
    resolve_env();
    EXPECT_TRUE(flexaidds::flags::requested("FLEXAIDDS_CLEFT_SORT"));
    EXPECT_FALSE(flexaidds::flags::active("FLEXAIDDS_CLEFT_SORT"));
    EXPECT_NE(std::strstr(flexaidds::flags::reason("FLEXAIDDS_CLEFT_SORT"),
                          "superseded"),
              nullptr);
}

TEST_F(FlexaiddsFlagsEnv, FlagsOverlayEnablesEpochAndFastpathWithoutClearingOthers) {
    set_env("FLEXAIDDS_RNG_STREAM_FIX", "1");
    set_env("FLEXAIDDS_CONTACTS_EPOCH", "1");
    set_env("FLEXAIDDS_FLAGS", "epoch,fastpath");
    resolve_env();
    EXPECT_TRUE(flexaidds::flags::requested("FLEXAIDDS_CONTACTS_EPOCH"));
    EXPECT_TRUE(flexaidds::flags::active("FLEXAIDDS_CONTACTS_EPOCH"));
    EXPECT_TRUE(flexaidds::flags::active("epoch"));
    EXPECT_TRUE(flexaidds::flags::requested("FLEXAIDDS_RIGID_FASTPATH"));
    EXPECT_TRUE(flexaidds::flags::active("FLEXAIDDS_RIGID_FASTPATH"));
    EXPECT_TRUE(flexaidds::flags::active("fastpath"));
    EXPECT_TRUE(flexaidds::flags::active("FLEXAIDDS_HOIST_RECEPTOR_INDEX"));
    EXPECT_TRUE(flexaidds::flags::requested("FLEXAIDDS_RNG_STREAM_FIX"));
    EXPECT_TRUE(flexaidds::flags::active("FLEXAIDDS_RNG_STREAM_FIX"));
}

TEST_F(FlexaiddsFlagsEnv, DumpIncludesRigidFastpath) {
    resolve_env();
    const std::string text = dump_to_string();
    EXPECT_NE(text.find("RIGID_FASTPATH"), std::string::npos);
    EXPECT_NE(text.find("FLEXAIDDS_RIGID_FASTPATH"), std::string::npos);
}

TEST_F(FlexaiddsFlagsEnv, OverlayAliasesAreCaseInsensitive) {
    set_env("FLEXAIDDS_FLAGS", "Hoist, EPOCH, FastPath, rng-stream-fix");
    resolve_env();
    EXPECT_TRUE(flexaidds::flags::requested("hoist"));
    EXPECT_TRUE(flexaidds::flags::active("FLEXAIDDS_HOIST_RECEPTOR_INDEX"));
    EXPECT_TRUE(flexaidds::flags::active("FLEXAIDDS_CONTACTS_EPOCH"));
    EXPECT_TRUE(flexaidds::flags::active("FLEXAIDDS_RIGID_FASTPATH"));
    EXPECT_TRUE(flexaidds::flags::active("FLEXAIDDS_RNG_STREAM_FIX"));
}

TEST_F(FlexaiddsFlagsEnv, ApplyToEnvironPublishesOverlayAndHidesLoser) {
    set_env("FLEXAIDDS_FLAGS", "fastpath");
    set_env("FLEXAIDDS_WAL_COERCIVE", "1");
    set_env("FLEXAIDDS_SOFTCORE_WAL", "1");
    set_env("FLEXAIDDS_CLEFT_SORT", "1");
    resolve_env();
    flexaidds::flags::apply_to_environ();

    const char* fp = std::getenv("FLEXAIDDS_RIGID_FASTPATH");
    ASSERT_NE(fp, nullptr);
    EXPECT_STREQ(fp, "1");
    const char* hoist = std::getenv("FLEXAIDDS_HOIST_RECEPTOR_INDEX");
    ASSERT_NE(hoist, nullptr);
    EXPECT_STREQ(hoist, "1");
    EXPECT_EQ(std::getenv("FLEXAIDDS_SOFTCORE_WAL"), nullptr);
    EXPECT_EQ(std::getenv("FLEXAIDDS_CLEFT_SORT"), nullptr);
    EXPECT_STREQ(std::getenv("FLEXAIDDS_WAL_COERCIVE"), "1");
    EXPECT_TRUE(flexaidds::flags::requested("FLEXAIDDS_SOFTCORE_WAL"));
    EXPECT_FALSE(flexaidds::flags::active("FLEXAIDDS_SOFTCORE_WAL"));
}
