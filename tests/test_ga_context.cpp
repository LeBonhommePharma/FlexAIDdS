// test_ga_context.cpp — Unit tests for GAContext re-entrant GA state
//
// Verifies that multiple GAContext instances maintain independent state,
// enabling parallel GA execution in ParallelDock and ParallelCampaign.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>
#include "GAContext.h"
#include "GPUContextPool.h"
#include "RngSeed.h"

#include <cstdlib>

static void set_test_env(const char* key, const char* value) {
#ifdef _WIN32
    _putenv_s(key, value);
#else
    setenv(key, value, 1);
#endif
}

static void unset_test_env(const char* key) {
#ifdef _WIN32
    _putenv_s(key, "");
#else
    unsetenv(key);
#endif
}

TEST(GAContextTest, DefaultConstruction) {
    GAContext ctx;
    EXPECT_EQ(ctx.gen_id, 0);
    EXPECT_EQ(ctx.nrejected, 0);
    EXPECT_FALSE(ctx.dispatch_logged);
    EXPECT_EQ(ctx.tqcm, nullptr);
    EXPECT_EQ(ctx.tqcm_ntypes, 0);
}

TEST(GAContextTest, IndependentCounters) {
    GAContext ctx1;
    GAContext ctx2;

    ctx1.gen_id = 42;
    ctx1.nrejected = 10;
    ctx1.dispatch_logged = true;
    ctx1.tqcm_ntypes = 256;

    // ctx2 must be unaffected
    EXPECT_EQ(ctx2.gen_id, 0);
    EXPECT_EQ(ctx2.nrejected, 0);
    EXPECT_FALSE(ctx2.dispatch_logged);
    EXPECT_EQ(ctx2.tqcm_ntypes, 0);
}

TEST(GAContextTest, MoveSemantics) {
    GAContext ctx1;
    ctx1.gen_id = 100;
    ctx1.nrejected = 5;
    ctx1.tqcm_ntypes = 64;

    GAContext ctx2 = std::move(ctx1);
    EXPECT_EQ(ctx2.gen_id, 100);
    EXPECT_EQ(ctx2.nrejected, 5);
    EXPECT_EQ(ctx2.tqcm_ntypes, 64);
}

TEST(GAContextTest, NotCopyable) {
    // GAContext should not be copyable (has unique_ptr member)
    EXPECT_FALSE(std::is_copy_constructible_v<GAContext>);
    EXPECT_FALSE(std::is_copy_assignable_v<GAContext>);
}

TEST(GAContextTest, IsMovable) {
    EXPECT_TRUE(std::is_move_constructible_v<GAContext>);
    EXPECT_TRUE(std::is_move_assignable_v<GAContext>);
}

TEST(RngSeedTest, FlexaidSeedMakesStreamSeedsRepeatable) {
    set_test_env("FLEXAID_SEED", "42");
    auto s1 = flexaids_rng::seed_from_env_or_random(123);
    auto s2 = flexaids_rng::seed_from_env_or_random(123);
    auto s3 = flexaids_rng::seed_from_env_or_random(124);
    unset_test_env("FLEXAID_SEED");

    EXPECT_EQ(s1, s2);
    EXPECT_NE(s1, s3);
}

#ifdef FLEXAIDS_USE_CUDA
TEST(GPUContextPoolTest, SingletonInstance) {
    auto& pool1 = GPUContextPool::instance();
    auto& pool2 = GPUContextPool::instance();
    EXPECT_EQ(&pool1, &pool2);
}
#endif

#ifndef FLEXAIDS_USE_CUDA
#ifndef FLEXAIDS_USE_METAL
TEST(GPUContextPoolTest, SingletonInstanceNoGPU) {
    // Pool should still be constructible without GPU backends
    auto& pool1 = GPUContextPool::instance();
    auto& pool2 = GPUContextPool::instance();
    EXPECT_EQ(&pool1, &pool2);
}
#endif
#endif
