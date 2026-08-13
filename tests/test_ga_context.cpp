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
#include "EnvFlags.h"

#include <cstdlib>
#include <cmath>

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

// Sets an env var for one test and always clears it, including when an
// assertion fails part-way. Tests share a process, so a leaked behaviour flag
// would silently flip whichever test ran next.
namespace {
struct ScopedEnv {
    const char* key;
    ScopedEnv(const char* k, const char* v) : key(k) { set_test_env(k, v); }
    ~ScopedEnv() { unset_test_env(key); }
    ScopedEnv(const ScopedEnv&) = delete;
    ScopedEnv& operator=(const ScopedEnv&) = delete;
};
}  // namespace

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

TEST(RngSeedTest, MasterSeedOverridesRandomDevice) {
    unset_test_env("FLEXAID_SEED");
    flexaids_rng::set_master_seed(11);
    auto s1 = flexaids_rng::seed_from_env_or_random(0x0C0A11ULL);
    auto s2 = flexaids_rng::seed_from_env_or_random(0x0C0A11ULL);
    EXPECT_EQ(s1, s2);
    EXPECT_TRUE(flexaids_rng::has_master_seed());
    EXPECT_EQ(flexaids_rng::master_seed(), 11u);
}

TEST(RngSeedTest, LazyThreadRngRespectsMasterSeedEpoch) {
    flexaids_rng::set_master_seed(11);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    const double a = dist(flexaids_rng::lazy_thread_rng(0x9A800DULL));
    flexaids_rng::set_master_seed(11);
    const double b = dist(flexaids_rng::lazy_thread_rng(0x9A800DULL));
    const double c = dist(flexaids_rng::lazy_thread_rng(0x9A800DULL));
    EXPECT_DOUBLE_EQ(a, b);
    EXPECT_NE(a, c);
}

TEST(RngSeedTest, LazyThreadRngDefaultOffReproducesLegacySequence) {
    // DEFAULT OFF must be the pre-fix behaviour bit-for-bit, because every
    // frozen reference number in this repo was produced under it.
    //
    // Legacy semantics: ONE generator per thread, re-seeded whenever the
    // requested stream differs from the cached one. So interleaving streams
    // collapses each to its FIRST draw forever — that collapse is the thing
    // being asserted, not a nicety.
    unset_test_env("FLEXAIDDS_RNG_STREAM_FIX");
    constexpr std::uint64_t kGa = 0x9A800DULL;
    constexpr std::uint64_t kVcontacts = 0x0C0A11ULL;
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    set_test_env("FLEXAID_SEED", "12345");
    flexaids_rng::set_master_seed(12345);  // bumps the epoch: flag is re-read

    const double a0 = dist(flexaids_rng::lazy_thread_rng(kGa));
    const double b0 = dist(flexaids_rng::lazy_thread_rng(kVcontacts));
    const double a1 = dist(flexaids_rng::lazy_thread_rng(kGa));

    // Absolute check: each draw equals the FIRST draw of a freshly seeded
    // generator for that stream.
    std::mt19937 ref_ga = flexaids_rng::make_thread_rng(kGa);
    std::mt19937 ref_vc = flexaids_rng::make_thread_rng(kVcontacts);
    const double first_ga = dist(ref_ga);
    const double first_vc = dist(ref_vc);

    EXPECT_DOUBLE_EQ(a0, first_ga);
    EXPECT_DOUBLE_EQ(b0, first_vc);
    EXPECT_DOUBLE_EQ(a1, first_ga) << "OFF must re-seed on stream switch";
    EXPECT_DOUBLE_EQ(a1, a0) << "legacy collapse to first draw not reproduced";

    unset_test_env("FLEXAID_SEED");
}

TEST(RngSeedTest, LazyThreadRngKeepsIndependentStreams) {
    // Opt-in path: FLEXAIDDS_RNG_STREAM_FIX=1.
    // Interleaving two stream ids on one thread must not re-seed a shared
    // generator. Stream A after a B draw must continue A's sequence.
    ScopedEnv fix("FLEXAIDDS_RNG_STREAM_FIX", "1");
    constexpr std::uint64_t kGa = 0x9A800DULL;
    constexpr std::uint64_t kPucker = 0x5A6A9ULL;
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    flexaids_rng::set_master_seed(11);
    const double a0 = dist(flexaids_rng::lazy_thread_rng(kGa));
    const double b0 = dist(flexaids_rng::lazy_thread_rng(kPucker));
    const double a1 = dist(flexaids_rng::lazy_thread_rng(kGa));

    flexaids_rng::set_master_seed(11);
    const double a0_only = dist(flexaids_rng::lazy_thread_rng(kGa));
    const double a1_only = dist(flexaids_rng::lazy_thread_rng(kGa));

    EXPECT_DOUBLE_EQ(a0, a0_only);
    EXPECT_DOUBLE_EQ(a1, a1_only);
    EXPECT_NE(a1, a0);
    EXPECT_NE(b0, a0);
}

TEST(RngSeedTest, LazyThreadRngHeldReferenceSurvivesNewStream) {
    // Opt-in path: FLEXAIDDS_RNG_STREAM_FIX=1.
    // Callers keep auto& rng across other stream lookups (SugarPucker,
    // Vcontacts). Insert must not invalidate that reference.
    ScopedEnv fix("FLEXAIDDS_RNG_STREAM_FIX", "1");
    constexpr std::uint64_t kGa = 0x9A800DULL;
    constexpr std::uint64_t kPucker = 0x5A6A9ULL;
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    flexaids_rng::set_master_seed(11);
    auto& ga = flexaids_rng::lazy_thread_rng(kGa);
    const double a0 = dist(ga);
    (void)dist(flexaids_rng::lazy_thread_rng(kPucker));
    const double a1 = dist(ga);

    flexaids_rng::set_master_seed(11);
    auto& ga_only = flexaids_rng::lazy_thread_rng(kGa);
    const double b0 = dist(ga_only);
    const double b1 = dist(ga_only);
    EXPECT_DOUBLE_EQ(a0, b0);
    EXPECT_DOUBLE_EQ(a1, b1);
}

TEST(RngSeedTest, LazyThreadRngThreeStreamInterleaveKeepsEachSequence) {
    // F9: GA / Vcontacts / FOPTICS streams on one thread, FIX=1.
    ScopedEnv fix("FLEXAIDDS_RNG_STREAM_FIX", "1");
    constexpr std::uint64_t kGa = 0x9A800DULL;
    constexpr std::uint64_t kVc = 0x0C0A11ULL;
    constexpr std::uint64_t kFo = 0xF0701C5ULL;
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    flexaids_rng::set_master_seed(12345);

    const double g0 = dist(flexaids_rng::lazy_thread_rng(kGa));
    const double v0 = dist(flexaids_rng::lazy_thread_rng(kVc));
    const double f0 = dist(flexaids_rng::lazy_thread_rng(kFo));
    const double g1 = dist(flexaids_rng::lazy_thread_rng(kGa));
    const double v1 = dist(flexaids_rng::lazy_thread_rng(kVc));

    flexaids_rng::set_master_seed(12345);
    const double g0s = dist(flexaids_rng::lazy_thread_rng(kGa));
    const double g1s = dist(flexaids_rng::lazy_thread_rng(kGa));
    const double v0s = dist(flexaids_rng::lazy_thread_rng(kVc));
    const double v1s = dist(flexaids_rng::lazy_thread_rng(kVc));
    const double f0s = dist(flexaids_rng::lazy_thread_rng(kFo));

    EXPECT_DOUBLE_EQ(g0, g0s);
    EXPECT_DOUBLE_EQ(g1, g1s);
    EXPECT_DOUBLE_EQ(v0, v0s);
    EXPECT_DOUBLE_EQ(v1, v1s);
    EXPECT_DOUBLE_EQ(f0, f0s);
    EXPECT_NE(g0, v0);
    EXPECT_NE(g0, f0);
}

TEST(RngSeedTest, VoronoiKeyedJitterIndependentOfRngStreamFix) {
    unset_test_env("FLEXAIDDS_RNG_STREAM_FIX");
    unset_test_env("FLEXAIDDS_VORONOI_KEYED_JITTER");
    flexaids_rng::set_master_seed(12345);
    EXPECT_FALSE(flexaids_rng::voronoi_keyed_jitter_enabled());
    {
        ScopedEnv fix("FLEXAIDDS_RNG_STREAM_FIX", "1");
        flexaids_rng::set_master_seed(12345);
        EXPECT_TRUE(flexaids_rng::rng_stream_fix_enabled());
        EXPECT_FALSE(flexaids_rng::voronoi_keyed_jitter_enabled());
    }
    flexaids_rng::set_master_seed(12345);
    {
        ScopedEnv jit("FLEXAIDDS_VORONOI_KEYED_JITTER", "1");
        flexaids_rng::set_master_seed(12345);
        EXPECT_TRUE(flexaids_rng::voronoi_keyed_jitter_enabled());
        EXPECT_FALSE(flexaids_rng::rng_stream_fix_enabled());
    }
    {
        ScopedEnv off("FLEXAIDDS_VORONOI_KEYED_JITTER", "off");
        EXPECT_FALSE(flexaids_rng::voronoi_keyed_jitter_enabled());
    }
}

TEST(RngSeedTest, KeyedJitterDependsOnPoseAtomNotDrawOrder) {
    ScopedEnv seed("FLEXAID_SEED", "12345");
    flexaids_rng::set_master_seed(12345);
    const auto id_a = flexaids_rng::pose_atom_identity(90001, 1.0f, 2.0f, 3.0f);
    const auto id_b = flexaids_rng::pose_atom_identity(90002, 1.0f, 2.0f, 3.0f);
    const auto id_a2 = flexaids_rng::pose_atom_identity(90001, 1.0f, 2.0f, 3.0f);
    EXPECT_EQ(id_a, id_a2);
    EXPECT_NE(id_a, id_b);

    const float ax0 = flexaids_rng::keyed_jitter(id_a, 0);
    const float ay0 = flexaids_rng::keyed_jitter(id_a, 1);
    const float bx0 = flexaids_rng::keyed_jitter(id_b, 0);
    // Consume the shared thread RNG so a thread-keyed implementation would move.
    (void)flexaids_rng::lazy_thread_rng(0x0C0A11ULL)();
    const float ax1 = flexaids_rng::keyed_jitter(id_a, 0);
    const float ay1 = flexaids_rng::keyed_jitter(id_a, 1);

    EXPECT_FLOAT_EQ(ax0, ax1);
    EXPECT_FLOAT_EQ(ay0, ay1);
    EXPECT_NE(ax0, bx0);
    EXPECT_LE(std::fabs(ax0), 0.005f);
    EXPECT_LE(std::fabs(ay0), 0.005f);
}

TEST(EnvFlagsTest, ParallelReproduceDefaultOff) {
    unset_test_env("FLEXAIDDS_PARALLEL_REPRODUCE");
    EXPECT_FALSE(flexaids::parallel_reproduce_enabled());
    {
        ScopedEnv on("FLEXAIDDS_PARALLEL_REPRODUCE", "1");
        EXPECT_TRUE(flexaids::parallel_reproduce_enabled());
    }
    {
        ScopedEnv off("FLEXAIDDS_PARALLEL_REPRODUCE", "0");
        EXPECT_FALSE(flexaids::parallel_reproduce_enabled());
    }
    EXPECT_FALSE(flexaids::parallel_reproduce_enabled());
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
