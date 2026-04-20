/// @file test_knowledge_pool.cpp
/// @brief Comprehensive tests for TargetKnowledgeBase and SharedPosePool.
///
/// Tests thread safety, serialization round-trips, edge cases, invariants,
/// and the MPIMergeable concept satisfaction.

#include <gtest/gtest.h>
#include <thread>
#include <vector>
#include <cmath>
#include <numeric>
#include <limits>
#include <random>

#include "../LIB/TargetKnowledgeBase.h"
#include "../LIB/SharedPosePool.h"
#include "../LIB/MPIMergeable.h"

using namespace target;

// ============================================================================
// Compile-time: MPIMergeable concept satisfaction
// ============================================================================

static_assert(transport::MPIMergeable<SharedPosePool>,
    "SharedPosePool must satisfy MPIMergeable");
static_assert(transport::MPIMergeable<TargetKnowledgeBase>,
    "TargetKnowledgeBase must satisfy MPIMergeable");

// ============================================================================
// TargetKnowledgeBase — Construction
// ============================================================================

TEST(TargetKnowledgeBase, DefaultConstruction) {
    TargetKnowledgeBase kb;
    EXPECT_EQ(kb.n_models(), 0);
    EXPECT_EQ(kb.n_grid_points(), 0);
    EXPECT_EQ(kb.conformer_observation_count(), 0);
    EXPECT_EQ(kb.grid_observation_count(), 0);
    EXPECT_EQ(kb.hit_count(), 0);
    EXPECT_TRUE(kb.check_invariants());
}

TEST(TargetKnowledgeBase, ParameterizedConstruction) {
    TargetKnowledgeBase kb(5, 100);
    EXPECT_EQ(kb.n_models(), 5);
    EXPECT_EQ(kb.n_grid_points(), 100);
    EXPECT_TRUE(kb.check_invariants());
}

TEST(TargetKnowledgeBase, ConstructionWithZeroDimensions) {
    TargetKnowledgeBase kb(0, 0);
    EXPECT_EQ(kb.n_models(), 0);
    EXPECT_EQ(kb.n_grid_points(), 0);
    EXPECT_TRUE(kb.check_invariants());
}

TEST(TargetKnowledgeBase, NegativeDimensionsThrow) {
    EXPECT_THROW(TargetKnowledgeBase(-1, 10), std::invalid_argument);
    EXPECT_THROW(TargetKnowledgeBase(10, -1), std::invalid_argument);
    EXPECT_THROW(TargetKnowledgeBase(-5, -5), std::invalid_argument);
}

// ============================================================================
// TargetKnowledgeBase — Conformer accumulation
// ============================================================================

TEST(TargetKnowledgeBase, ConformerAccumulation) {
    TargetKnowledgeBase kb(3, 0);

    EXPECT_TRUE(kb.accumulate_conformer_weights({0.5, 0.3, 0.2}));
    EXPECT_EQ(kb.conformer_observation_count(), 1);

    EXPECT_TRUE(kb.accumulate_conformer_weights({0.6, 0.2, 0.2}));
    EXPECT_EQ(kb.conformer_observation_count(), 2);
    EXPECT_TRUE(kb.check_invariants());
}

TEST(TargetKnowledgeBase, ConformerDimensionMismatch) {
    TargetKnowledgeBase kb(3, 0);

    // Wrong size
    EXPECT_FALSE(kb.accumulate_conformer_weights({0.5, 0.5}));
    EXPECT_EQ(kb.conformer_observation_count(), 0);

    // Empty
    EXPECT_FALSE(kb.accumulate_conformer_weights({}));
    EXPECT_EQ(kb.conformer_observation_count(), 0);
}

TEST(TargetKnowledgeBase, ConformerAutoInitDimension) {
    TargetKnowledgeBase kb;  // n_models = 0 initially
    EXPECT_EQ(kb.n_models(), 0);

    EXPECT_TRUE(kb.accumulate_conformer_weights({0.4, 0.3, 0.2, 0.1}));
    EXPECT_EQ(kb.n_models(), 4);
    EXPECT_EQ(kb.conformer_observation_count(), 1);

    // Now different size should fail
    EXPECT_FALSE(kb.accumulate_conformer_weights({0.5, 0.5}));
}

TEST(TargetKnowledgeBase, ConformerPosteriorWithPrior) {
    TargetKnowledgeBase kb(2, 0);

    // No observations: posterior should be uniform (from prior only)
    auto prior = kb.conformer_posterior();
    ASSERT_EQ(prior.size(), 2u);
    EXPECT_NEAR(prior[0], 0.5, 1e-10);
    EXPECT_NEAR(prior[1], 0.5, 1e-10);

    // After one observation favoring model 0
    kb.accumulate_conformer_weights({0.8, 0.2});
    auto post = kb.conformer_posterior();
    ASSERT_EQ(post.size(), 2u);
    // posterior = (sum + 1) / total
    // model 0: (0.8 + 1) / (0.8 + 0.2 + 2) = 1.8 / 3.0 = 0.6
    // model 1: (0.2 + 1) / 3.0 = 1.2 / 3.0 = 0.4
    EXPECT_NEAR(post[0], 0.6, 1e-10);
    EXPECT_NEAR(post[1], 0.4, 1e-10);

    // Sum should be 1.0
    double sum = post[0] + post[1];
    EXPECT_NEAR(sum, 1.0, 1e-10);
}

TEST(TargetKnowledgeBase, ConformerPosteriorEmpty) {
    TargetKnowledgeBase kb;
    auto post = kb.conformer_posterior();
    EXPECT_TRUE(post.empty());
}

// ============================================================================
// TargetKnowledgeBase — Grid energy accumulation
// ============================================================================

TEST(TargetKnowledgeBase, GridEnergyAccumulation) {
    TargetKnowledgeBase kb(0, 3);

    EXPECT_TRUE(kb.accumulate_grid_energies({-2.0f, 0.5f, -0.5f}));
    EXPECT_EQ(kb.grid_observation_count(), 1);

    auto means = kb.grid_mean_energies();
    ASSERT_EQ(means.size(), 3u);
    EXPECT_FLOAT_EQ(means[0], -2.0f);
    EXPECT_FLOAT_EQ(means[1], 0.5f);
    EXPECT_FLOAT_EQ(means[2], -0.5f);

    // Second ligand
    EXPECT_TRUE(kb.accumulate_grid_energies({-4.0f, 1.5f, -1.5f}));
    EXPECT_EQ(kb.grid_observation_count(), 2);

    means = kb.grid_mean_energies();
    EXPECT_FLOAT_EQ(means[0], -3.0f);   // (-2 + -4) / 2
    EXPECT_FLOAT_EQ(means[1], 1.0f);    // (0.5 + 1.5) / 2
    EXPECT_FLOAT_EQ(means[2], -1.0f);   // (-0.5 + -1.5) / 2
}

TEST(TargetKnowledgeBase, GridEnergyNaNHandling) {
    TargetKnowledgeBase kb(0, 2);

    float nan_val = std::numeric_limits<float>::quiet_NaN();
    EXPECT_TRUE(kb.accumulate_grid_energies({-3.0f, nan_val}));

    auto means = kb.grid_mean_energies();
    EXPECT_FLOAT_EQ(means[0], -3.0f);
    EXPECT_FLOAT_EQ(means[1], 0.0f);  // No valid observations at point 1
}

TEST(TargetKnowledgeBase, GridEnergyInfHandling) {
    TargetKnowledgeBase kb(0, 2);

    float inf_val = std::numeric_limits<float>::infinity();
    EXPECT_TRUE(kb.accumulate_grid_energies({-5.0f, inf_val}));

    auto means = kb.grid_mean_energies();
    EXPECT_FLOAT_EQ(means[0], -5.0f);
    EXPECT_FLOAT_EQ(means[1], 0.0f);  // Inf skipped
}

TEST(TargetKnowledgeBase, GridDimensionMismatch) {
    TargetKnowledgeBase kb(0, 3);
    EXPECT_FALSE(kb.accumulate_grid_energies({1.0f, 2.0f}));  // wrong size
    EXPECT_FALSE(kb.accumulate_grid_energies({}));              // empty
    EXPECT_EQ(kb.grid_observation_count(), 0);
}

TEST(TargetKnowledgeBase, GridHotspotFraction) {
    TargetKnowledgeBase kb(0, 3);

    // Point 0: always favorable (< -1.0)
    // Point 1: never favorable
    // Point 2: favorable in 1 of 2 ligands
    kb.accumulate_grid_energies({-5.0f, 0.0f, -2.0f});
    kb.accumulate_grid_energies({-3.0f, 0.5f, 0.0f});

    auto fractions = kb.grid_hotspot_fraction(-1.0f);
    ASSERT_EQ(fractions.size(), 3u);
    EXPECT_FLOAT_EQ(fractions[0], 1.0f);   // 2/2
    EXPECT_FLOAT_EQ(fractions[1], 0.0f);   // 0/2
    EXPECT_FLOAT_EQ(fractions[2], 0.5f);   // 1/2
}

TEST(TargetKnowledgeBase, GridHotspotFractionEmpty) {
    TargetKnowledgeBase kb(0, 3);
    auto fractions = kb.grid_hotspot_fraction();
    ASSERT_EQ(fractions.size(), 3u);
    for (auto f : fractions) EXPECT_FLOAT_EQ(f, 0.0f);
}

// ============================================================================
// TargetKnowledgeBase — Binding subsite hits
// ============================================================================

TEST(TargetKnowledgeBase, SubsiteHitAccumulation) {
    TargetKnowledgeBase kb;

    kb.accumulate_binding_center(1.0f, 2.0f, 3.0f, -8.0, "ligandA");
    kb.accumulate_binding_center(4.0f, 5.0f, 6.0f, -6.0, "ligandB");

    EXPECT_EQ(kb.hit_count(), 2);

    auto hits = kb.all_hits();
    ASSERT_EQ(hits.size(), 2u);
    EXPECT_FLOAT_EQ(hits[0].center[0], 1.0f);
    EXPECT_NEAR(hits[0].energy, -8.0, 1e-10);
    EXPECT_EQ(hits[0].ligand_name, "ligandA");
    EXPECT_FLOAT_EQ(hits[1].center[0], 4.0f);
}

TEST(TargetKnowledgeBase, AllHitsReturnsCopy) {
    TargetKnowledgeBase kb;
    kb.accumulate_binding_center(0, 0, 0, -5.0, "lig1");

    auto snapshot = kb.all_hits();
    EXPECT_EQ(snapshot.size(), 1u);

    // Add another after snapshot
    kb.accumulate_binding_center(1, 1, 1, -3.0, "lig2");

    // Snapshot should NOT have the new hit
    EXPECT_EQ(snapshot.size(), 1u);

    // Fresh query should
    EXPECT_EQ(kb.all_hits().size(), 2u);
}

// ============================================================================
// TargetKnowledgeBase — Clear and invariants
// ============================================================================

TEST(TargetKnowledgeBase, ClearResetsAll) {
    TargetKnowledgeBase kb(3, 5);

    kb.accumulate_conformer_weights({0.5, 0.3, 0.2});
    kb.accumulate_grid_energies({-1.0f, -2.0f, -3.0f, -4.0f, -5.0f});
    kb.accumulate_binding_center(0, 0, 0, -5.0, "lig");

    kb.clear();

    EXPECT_EQ(kb.conformer_observation_count(), 0);
    EXPECT_EQ(kb.grid_observation_count(), 0);
    EXPECT_EQ(kb.hit_count(), 0);
    // Dimensions preserved
    EXPECT_EQ(kb.n_models(), 3);
    EXPECT_EQ(kb.n_grid_points(), 5);
    EXPECT_TRUE(kb.check_invariants());
}

// ============================================================================
// TargetKnowledgeBase — Serialization round-trip
// ============================================================================

TEST(TargetKnowledgeBase, SerializeDeserializeRoundTrip) {
    TargetKnowledgeBase original(3, 4);

    original.accumulate_conformer_weights({0.6, 0.3, 0.1});
    original.accumulate_conformer_weights({0.4, 0.4, 0.2});
    original.accumulate_grid_energies({-2.0f, -1.0f, 0.5f, -3.0f});
    original.accumulate_grid_energies({-1.0f, -0.5f, 0.0f, -2.0f});
    original.accumulate_binding_center(10.0f, 20.0f, 30.0f, -7.5, "aspirin");
    original.accumulate_binding_center(15.0f, 25.0f, 35.0f, -5.0, "ibuprofen");

    auto buf = original.serialize();
    EXPECT_GT(buf.size(), 0u);

    // Deserialize into a fresh KB
    TargetKnowledgeBase restored;
    restored.deserialize_merge(buf.data(), buf.size());

    // Verify conformer data
    EXPECT_EQ(restored.conformer_observation_count(), 2);
    auto orig_post = original.conformer_posterior();
    auto rest_post = restored.conformer_posterior();
    ASSERT_EQ(orig_post.size(), rest_post.size());
    for (size_t i = 0; i < orig_post.size(); ++i) {
        EXPECT_NEAR(orig_post[i], rest_post[i], 1e-10);
    }

    // Verify grid data
    EXPECT_EQ(restored.grid_observation_count(), 2);
    auto orig_means = original.grid_mean_energies();
    auto rest_means = restored.grid_mean_energies();
    ASSERT_EQ(orig_means.size(), rest_means.size());
    for (size_t i = 0; i < orig_means.size(); ++i) {
        EXPECT_FLOAT_EQ(orig_means[i], rest_means[i]);
    }

    // Verify hits
    auto orig_hits = original.all_hits();
    auto rest_hits = restored.all_hits();
    ASSERT_EQ(orig_hits.size(), rest_hits.size());
    for (size_t i = 0; i < orig_hits.size(); ++i) {
        EXPECT_FLOAT_EQ(orig_hits[i].center[0], rest_hits[i].center[0]);
        EXPECT_FLOAT_EQ(orig_hits[i].center[1], rest_hits[i].center[1]);
        EXPECT_FLOAT_EQ(orig_hits[i].center[2], rest_hits[i].center[2]);
        EXPECT_NEAR(orig_hits[i].energy, rest_hits[i].energy, 1e-10);
        EXPECT_EQ(orig_hits[i].ligand_name, rest_hits[i].ligand_name);
    }

    EXPECT_TRUE(restored.check_invariants());
}

TEST(TargetKnowledgeBase, SerializeMergeAddsData) {
    TargetKnowledgeBase kb_a(2, 3);
    TargetKnowledgeBase kb_b(2, 3);

    kb_a.accumulate_conformer_weights({0.7, 0.3});
    kb_a.accumulate_grid_energies({-1.0f, -2.0f, -3.0f});
    kb_a.accumulate_binding_center(1, 2, 3, -5.0, "lig_a");

    kb_b.accumulate_conformer_weights({0.4, 0.6});
    kb_b.accumulate_grid_energies({-2.0f, -1.0f, -1.0f});
    kb_b.accumulate_binding_center(4, 5, 6, -3.0, "lig_b");

    // Merge B into A
    auto buf_b = kb_b.serialize();
    kb_a.deserialize_merge(buf_b.data(), buf_b.size());

    // A should now have 2 conformer observations, 2 grid observations, 2 hits
    EXPECT_EQ(kb_a.conformer_observation_count(), 2);
    EXPECT_EQ(kb_a.grid_observation_count(), 2);
    EXPECT_EQ(kb_a.hit_count(), 2);
    EXPECT_TRUE(kb_a.check_invariants());
}

TEST(TargetKnowledgeBase, DeserializeEmptyBuffer) {
    TargetKnowledgeBase kb(3, 3);
    kb.deserialize_merge(nullptr, 0);
    kb.deserialize_merge("", 0);
    kb.deserialize_merge("short", 3);  // too short for header
    EXPECT_TRUE(kb.check_invariants());
}

TEST(TargetKnowledgeBase, SerializeEmptyKB) {
    TargetKnowledgeBase kb;
    auto buf = kb.serialize();
    EXPECT_GT(buf.size(), 0u);  // at least the header

    TargetKnowledgeBase restored;
    restored.deserialize_merge(buf.data(), buf.size());
    EXPECT_EQ(restored.n_models(), 0);
    EXPECT_EQ(restored.n_grid_points(), 0);
    EXPECT_EQ(restored.hit_count(), 0);
}

// ============================================================================
// TargetKnowledgeBase — Concurrent access
// ============================================================================

TEST(TargetKnowledgeBase, ConcurrentConformerAccumulation) {
    TargetKnowledgeBase kb(4, 0);

    const int n_threads = 8;
    const int ops_per_thread = 100;

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (int t = 0; t < n_threads; ++t) {
        threads.emplace_back([&kb, t]() {
            for (int i = 0; i < ops_per_thread; ++i) {
                double w = 1.0 / 4.0;
                kb.accumulate_conformer_weights({w, w, w, w});
                (void)kb.conformer_posterior();
            }
        });
    }
    for (auto& t : threads) t.join();

    EXPECT_EQ(kb.conformer_observation_count(), n_threads * ops_per_thread);
    EXPECT_TRUE(kb.check_invariants());
}

TEST(TargetKnowledgeBase, ConcurrentMixedAccumulation) {
    TargetKnowledgeBase kb(2, 3);

    const int n_threads = 8;
    const int ops_per_thread = 50;

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (int t = 0; t < n_threads; ++t) {
        threads.emplace_back([&kb, t]() {
            for (int i = 0; i < ops_per_thread; ++i) {
                kb.accumulate_conformer_weights({0.5, 0.5});
                kb.accumulate_grid_energies({-1.0f, 0.0f, 1.0f});
                kb.accumulate_binding_center(
                    static_cast<float>(t), static_cast<float>(i), 0.0f,
                    -static_cast<double>(i), "lig_" + std::to_string(t));
                (void)kb.all_hits();
                (void)kb.grid_mean_energies();
            }
        });
    }
    for (auto& t : threads) t.join();

    EXPECT_EQ(kb.conformer_observation_count(), n_threads * ops_per_thread);
    EXPECT_EQ(kb.grid_observation_count(), n_threads * ops_per_thread);
    EXPECT_EQ(kb.hit_count(), n_threads * ops_per_thread);
    EXPECT_TRUE(kb.check_invariants());
}

TEST(TargetKnowledgeBase, AllHitsRefReturnsReference) {
    TargetKnowledgeBase kb(1, 0);
    kb.accumulate_binding_center(1.0f, 2.0f, 3.0f, -5.0, "ref_lig");

    const auto& ref = kb.all_hits_ref();
    ASSERT_EQ(ref.size(), 1u);
    EXPECT_FLOAT_EQ(ref[0].center[0], 1.0f);
    EXPECT_FLOAT_EQ(ref[0].center[1], 2.0f);
    EXPECT_FLOAT_EQ(ref[0].center[2], 3.0f);
    EXPECT_EQ(ref[0].ligand_name, "ref_lig");

    // all_hits() returns a copy; all_hits_ref() returns the same data by ref
    auto copy = kb.all_hits();
    EXPECT_EQ(copy.size(), ref.size());
    EXPECT_EQ(copy[0].ligand_name, ref[0].ligand_name);
}

// ============================================================================
// SharedPosePool — Construction
// ============================================================================

TEST(SharedPosePool, Construction) {
    SharedPosePool pool(50);
    EXPECT_EQ(pool.capacity(), 50);
    EXPECT_EQ(pool.count(), 0);
    EXPECT_FALSE(pool.is_full());
    // Sentinel is DBL_MAX (not infinity — avoids UB under -ffast-math)
    EXPECT_EQ(pool.best_energy(), std::numeric_limits<double>::max());
    EXPECT_EQ(pool.worst_energy(), std::numeric_limits<double>::max());
}

TEST(SharedPosePool, ZeroCapacityThrows) {
    EXPECT_THROW(SharedPosePool(0), std::invalid_argument);
}

TEST(SharedPosePool, NegativeCapacityThrows) {
    EXPECT_THROW(SharedPosePool(-5), std::invalid_argument);
}

TEST(SharedPosePool, SingleCapacity) {
    SharedPosePool pool(1);

    SharedPose p1(-5.0, 1.0f, 2.0f, 3.0f, 0, 0);
    EXPECT_TRUE(pool.publish(p1));
    EXPECT_TRUE(pool.is_full());
    EXPECT_EQ(pool.count(), 1);
    EXPECT_DOUBLE_EQ(pool.best_energy(), -5.0);

    // Better pose replaces the only slot
    SharedPose p2(-10.0, 4.0f, 5.0f, 6.0f, 1, 1);
    EXPECT_TRUE(pool.publish(p2));
    EXPECT_DOUBLE_EQ(pool.best_energy(), -10.0);

    // Worse pose rejected
    SharedPose p3(-3.0, 7.0f, 8.0f, 9.0f, 2, 2);
    EXPECT_FALSE(pool.publish(p3));
    EXPECT_DOUBLE_EQ(pool.best_energy(), -10.0);
}

// ============================================================================
// SharedPosePool — Publish and retrieve
// ============================================================================

TEST(SharedPosePool, PublishAndGetTop) {
    SharedPosePool pool(10);

    SharedPose p1; p1.energy = -5.0; p1.source_region = 0;
    SharedPose p2; p2.energy = -10.0; p2.source_region = 1;
    SharedPose p3; p3.energy = -3.0; p3.source_region = 2;

    EXPECT_TRUE(pool.publish(p1));
    EXPECT_TRUE(pool.publish(p2));
    EXPECT_TRUE(pool.publish(p3));

    auto top = pool.get_top(2);
    ASSERT_EQ(top.size(), 2u);
    EXPECT_DOUBLE_EQ(top[0].energy, -10.0);
    EXPECT_DOUBLE_EQ(top[1].energy, -5.0);
}

TEST(SharedPosePool, EvictsWorst) {
    SharedPosePool pool(3);

    for (int i = 0; i < 5; ++i) {
        SharedPose p;
        p.energy = -static_cast<double>(i);
        p.source_region = i;
        pool.publish(p);
    }

    auto top = pool.get_top(3);
    ASSERT_EQ(top.size(), 3u);
    EXPECT_DOUBLE_EQ(top[0].energy, -4.0);
    EXPECT_DOUBLE_EQ(top[1].energy, -3.0);
    EXPECT_DOUBLE_EQ(top[2].energy, -2.0);
}

TEST(SharedPosePool, GetTopClampsToCount) {
    SharedPosePool pool(100);

    SharedPose p; p.energy = -1.0;
    pool.publish(p);

    auto top = pool.get_top(50);
    ASSERT_EQ(top.size(), 1u);  // only 1 pose in pool
}

TEST(SharedPosePool, GetTopZero) {
    SharedPosePool pool(10);
    SharedPose p; p.energy = -1.0;
    pool.publish(p);

    auto top = pool.get_top(0);
    EXPECT_TRUE(top.empty());
}

TEST(SharedPosePool, GetTopNegative) {
    SharedPosePool pool(10);
    SharedPose p; p.energy = -1.0;
    pool.publish(p);

    auto top = pool.get_top(-5);
    EXPECT_TRUE(top.empty());
}

// ============================================================================
// SharedPosePool — NaN rejection
// ============================================================================

TEST(SharedPosePool, RejectsNaN) {
    SharedPosePool pool(10);

    SharedPose nan_pose;
    nan_pose.energy = std::numeric_limits<double>::quiet_NaN();
    EXPECT_FALSE(pool.publish(nan_pose));
    EXPECT_EQ(pool.count(), 0);
}

TEST(SharedPosePool, RejectsInfinity) {
    SharedPosePool pool(10);

    SharedPose inf_pose;
    inf_pose.energy = std::numeric_limits<double>::infinity();
    EXPECT_FALSE(pool.publish(inf_pose));

    SharedPose neg_inf_pose;
    neg_inf_pose.energy = -std::numeric_limits<double>::infinity();
    EXPECT_FALSE(pool.publish(neg_inf_pose));

    EXPECT_EQ(pool.count(), 0);
}

// ============================================================================
// SharedPosePool — Energy accessors
// ============================================================================

TEST(SharedPosePool, BestWorstEnergy) {
    SharedPosePool pool(10);

    SharedPose p1; p1.energy = -5.0;
    SharedPose p2; p2.energy = -10.0;
    SharedPose p3; p3.energy = -3.0;

    pool.publish(p1);
    pool.publish(p2);
    pool.publish(p3);

    EXPECT_DOUBLE_EQ(pool.best_energy(), -10.0);
    EXPECT_DOUBLE_EQ(pool.worst_energy(), -3.0);
    EXPECT_EQ(pool.count(), 3);
}

// ============================================================================
// SharedPosePool — Serialization
// ============================================================================

TEST(SharedPosePool, SerializeDeserialize) {
    SharedPosePool pool(10);

    SharedPose p1(-8.0, 1.0f, 2.0f, 3.0f, 0, 10);
    SharedPose p2(-6.0, 4.0f, 5.0f, 6.0f, 1, 20);
    pool.publish(p1);
    pool.publish(p2);

    auto buf = pool.serialize();

    SharedPosePool pool2(10);
    pool2.deserialize_merge(buf.data(), buf.size());

    auto top = pool2.get_top(5);
    ASSERT_EQ(top.size(), 2u);
    EXPECT_DOUBLE_EQ(top[0].energy, -8.0);
    EXPECT_DOUBLE_EQ(top[1].energy, -6.0);
    EXPECT_EQ(top[0].source_region, 0);
    EXPECT_EQ(top[1].source_region, 1);
    EXPECT_FLOAT_EQ(top[0].grid_coor[0], 1.0f);
}

TEST(SharedPosePool, DeserializeTruncatedBuffer) {
    SharedPosePool pool(10);

    // Buffer too small for header
    pool.deserialize_merge("abc", 3);
    EXPECT_EQ(pool.count(), 0);

    // Buffer with valid header but truncated payload
    char buf[sizeof(int)] = {};
    int count = 5;
    std::memcpy(buf, &count, sizeof(int));
    pool.deserialize_merge(buf, sizeof(int));  // no payload
    EXPECT_EQ(pool.count(), 0);
}

TEST(SharedPosePool, DeserializeNullBuffer) {
    SharedPosePool pool(10);
    pool.deserialize_merge(nullptr, 100);
    EXPECT_EQ(pool.count(), 0);
}

TEST(SharedPosePool, DeserializeNegativeCount) {
    SharedPosePool pool(10);

    int neg = -1;
    char buf[sizeof(int)];
    std::memcpy(buf, &neg, sizeof(int));
    pool.deserialize_merge(buf, sizeof(buf));
    EXPECT_EQ(pool.count(), 0);
}

TEST(SharedPosePool, SerializeMergeFromMultiplePools) {
    SharedPosePool pool_a(10);
    SharedPosePool pool_b(10);
    SharedPosePool combined(5);

    // Pool A has poses at -8, -6
    pool_a.publish(SharedPose(-8.0, 0, 0, 0, 0, 0));
    pool_a.publish(SharedPose(-6.0, 0, 0, 0, 0, 0));

    // Pool B has poses at -10, -4
    pool_b.publish(SharedPose(-10.0, 0, 0, 0, 1, 0));
    pool_b.publish(SharedPose(-4.0, 0, 0, 0, 1, 0));

    // Merge both into combined (capacity 5)
    auto buf_a = pool_a.serialize();
    auto buf_b = pool_b.serialize();
    combined.deserialize_merge(buf_a.data(), buf_a.size());
    combined.deserialize_merge(buf_b.data(), buf_b.size());

    auto top = combined.get_top(5);
    ASSERT_EQ(top.size(), 4u);
    EXPECT_DOUBLE_EQ(top[0].energy, -10.0);
    EXPECT_DOUBLE_EQ(top[1].energy, -8.0);
    EXPECT_DOUBLE_EQ(top[2].energy, -6.0);
    EXPECT_DOUBLE_EQ(top[3].energy, -4.0);
}

// ============================================================================
// SharedPosePool — Clear
// ============================================================================

TEST(SharedPosePool, ClearResetsPool) {
    SharedPosePool pool(10);

    pool.publish(SharedPose(-5.0, 0, 0, 0, 0, 0));
    pool.publish(SharedPose(-3.0, 0, 0, 0, 0, 0));
    EXPECT_EQ(pool.count(), 2);

    pool.clear();
    EXPECT_EQ(pool.count(), 0);
    EXPECT_FALSE(pool.is_full());
    // Sentinel is DBL_MAX (not infinity — avoids UB under -ffast-math)
    EXPECT_EQ(pool.best_energy(), std::numeric_limits<double>::max());
}

// ============================================================================
// SharedPosePool — Move semantics
// ============================================================================

TEST(SharedPosePool, MoveConstruction) {
    SharedPosePool pool(10);
    pool.publish(SharedPose(-8.0, 1.0f, 2.0f, 3.0f, 0, 0));
    pool.publish(SharedPose(-5.0, 4.0f, 5.0f, 6.0f, 1, 0));

    SharedPosePool moved(std::move(pool));
    EXPECT_EQ(moved.count(), 2);
    EXPECT_DOUBLE_EQ(moved.best_energy(), -8.0);
    EXPECT_EQ(moved.capacity(), 10);
}

TEST(SharedPosePool, MoveAssignment) {
    SharedPosePool pool(10);
    pool.publish(SharedPose(-8.0, 1.0f, 2.0f, 3.0f, 0, 0));

    SharedPosePool other(5);
    other = std::move(pool);
    EXPECT_EQ(other.count(), 1);
    EXPECT_EQ(other.capacity(), 10);
}

// ============================================================================
// SharedPosePool — Concurrent publish (std::thread)
// ============================================================================

TEST(SharedPosePool, ConcurrentPublishStdThread) {
    SharedPosePool pool(100);
    const int n_threads = 8;
    const int poses_per_thread = 200;

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (int t = 0; t < n_threads; ++t) {
        threads.emplace_back([&pool, t]() {
            std::mt19937 rng(static_cast<unsigned>(t * 12345));
            std::uniform_real_distribution<double> dist(-100.0, 0.0);
            for (int i = 0; i < poses_per_thread; ++i) {
                SharedPose p;
                p.energy = dist(rng);
                p.source_region = t;
                p.generation = i;
                pool.publish(p);
            }
        });
    }
    for (auto& t : threads) t.join();

    auto top = pool.get_top(100);
    EXPECT_EQ(static_cast<int>(top.size()), 100);

    // Verify sorted order
    for (size_t i = 1; i < top.size(); ++i) {
        EXPECT_LE(top[i - 1].energy, top[i].energy)
            << "Sort violation at index " << i;
    }
}

TEST(SharedPosePool, ConcurrentPublishAndRead) {
    SharedPosePool pool(50);
    const int n_writers = 4;
    const int n_readers = 4;
    const int ops = 100;

    std::vector<std::thread> threads;
    threads.reserve(n_writers + n_readers);

    // Writers
    for (int t = 0; t < n_writers; ++t) {
        threads.emplace_back([&pool, t]() {
            for (int i = 0; i < ops; ++i) {
                SharedPose p;
                p.energy = -static_cast<double>(t * ops + i);
                pool.publish(p);
            }
        });
    }

    // Readers
    for (int t = 0; t < n_readers; ++t) {
        threads.emplace_back([&pool]() {
            for (int i = 0; i < ops; ++i) {
                auto top = pool.get_top(10);
                (void)pool.count();
                (void)pool.best_energy();
                (void)pool.worst_energy();
                (void)pool.is_full();
            }
        });
    }

    for (auto& t : threads) t.join();

    // Just verify no crashes and invariants hold
    auto top = pool.get_top(50);
    for (size_t i = 1; i < top.size(); ++i) {
        EXPECT_LE(top[i - 1].energy, top[i].energy);
    }
}

// ============================================================================
// SharedPosePool — Sorted order stress test
// ============================================================================

TEST(SharedPosePool, SortedOrderStress) {
    SharedPosePool pool(20);
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(-50.0, 50.0);

    for (int i = 0; i < 1000; ++i) {
        SharedPose p;
        p.energy = dist(rng);
        p.source_region = i % 10;
        pool.publish(p);
    }

    auto top = pool.get_top(20);
    ASSERT_EQ(top.size(), 20u);
    for (size_t i = 1; i < top.size(); ++i) {
        EXPECT_LE(top[i - 1].energy, top[i].energy)
            << "Sort violation at index " << i
            << ": " << top[i-1].energy << " > " << top[i].energy;
    }
}

// ============================================================================
// SharedPosePool — Duplicate energy handling
// ============================================================================

TEST(SharedPosePool, DuplicateEnergies) {
    SharedPosePool pool(5);

    for (int i = 0; i < 10; ++i) {
        SharedPose p;
        p.energy = -5.0;  // all same energy
        p.source_region = i;
        pool.publish(p);
    }

    // Should fill to capacity, then reject equals (not strictly less)
    EXPECT_EQ(pool.count(), 5);

    auto top = pool.get_top(5);
    for (const auto& p : top) {
        EXPECT_DOUBLE_EQ(p.energy, -5.0);
    }
}
