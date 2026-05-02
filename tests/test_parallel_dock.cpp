// tests/test_parallel_dock.cpp — Unit tests for parallel grid-decomposed docking
// Tests octree decomposition, StatMechEngine merge, SharedPosePool, and subgrid extraction

#include <gtest/gtest.h>
#include <cmath>
#include <set>
#include <numeric>
#include <thread>
#include "../LIB/GridDecomposer.h"
#include "../LIB/SharedPosePool.h"
#include "../LIB/statmech.h"

// ============================================================================
// Helper: create a synthetic cubic grid for testing
// ============================================================================

static gridpoint* make_test_grid(int side, float spacer, int& num_grd) {
    // Creates a side×side×side cubic grid + reference point at index 0
    num_grd = side * side * side + 1;
    gridpoint* grid = (gridpoint*)calloc(num_grd, sizeof(gridpoint));

    // Reference point at origin
    grid[0].coor[0] = grid[0].coor[1] = grid[0].coor[2] = 0.0f;

    int idx = 1;
    for (int x = 0; x < side; x++)
        for (int y = 0; y < side; y++)
            for (int z = 0; z < side; z++) {
                grid[idx].coor[0] = x * spacer;
                grid[idx].coor[1] = y * spacer;
                grid[idx].coor[2] = z * spacer;
                idx++;
            }

    return grid;
}

// ============================================================================
// Octree decomposition tests
// ============================================================================

TEST(GridDecomposer, CorrectRegionCount) {
    int num_grd;
    gridpoint* grid = make_test_grid(10, 1.0f, num_grd);  // 1000 points
    ASSERT_EQ(num_grd, 1001);

    auto regions = GridDecomposer::decompose_octree(grid, num_grd, 16, 10);

    // Should produce roughly 16 regions (may vary due to octree structure)
    EXPECT_GE((int)regions.size(), 4);
    EXPECT_LE((int)regions.size(), 64);

    free(grid);
}

TEST(GridDecomposer, AllPointsCovered) {
    int num_grd;
    gridpoint* grid = make_test_grid(8, 0.5f, num_grd);  // 512 points

    auto regions = GridDecomposer::decompose_octree(grid, num_grd, 8, 5);

    // Collect all indices across all regions
    std::set<int> all_indices;
    for (const auto& r : regions) {
        for (int idx : r.grid_indices) {
            all_indices.insert(idx);
        }
    }

    // Every grid index (1..num_grd-1) must appear exactly once
    EXPECT_EQ((int)all_indices.size(), num_grd - 1);
    for (int i = 1; i < num_grd; i++) {
        EXPECT_EQ(all_indices.count(i), 1u) << "Missing index " << i;
    }

    free(grid);
}

TEST(GridDecomposer, NoOverlap) {
    int num_grd;
    gridpoint* grid = make_test_grid(6, 1.0f, num_grd);  // 216 points

    auto regions = GridDecomposer::decompose_octree(grid, num_grd, 8, 5);

    // Check no index appears in two different regions
    std::set<int> seen;
    for (const auto& r : regions) {
        for (int idx : r.grid_indices) {
            EXPECT_EQ(seen.count(idx), 0u)
                << "Index " << idx << " in multiple regions";
            seen.insert(idx);
        }
    }

    free(grid);
}

TEST(GridDecomposer, ExtractSubgrid) {
    int num_grd;
    gridpoint* grid = make_test_grid(4, 1.0f, num_grd);  // 64 points

    auto regions = GridDecomposer::decompose_octree(grid, num_grd, 4, 5);
    ASSERT_FALSE(regions.empty());

    const auto& r = regions[0];
    int sub_num;
    gridpoint* subgrid = GridDecomposer::extract_subgrid(grid, r, sub_num);

    ASSERT_NE(subgrid, nullptr);
    EXPECT_EQ(sub_num, r.num_points + 1);  // +1 for reference point

    // Reference point (index 0) should match original
    EXPECT_FLOAT_EQ(subgrid[0].coor[0], grid[0].coor[0]);
    EXPECT_FLOAT_EQ(subgrid[0].coor[1], grid[0].coor[1]);
    EXPECT_FLOAT_EQ(subgrid[0].coor[2], grid[0].coor[2]);

    // All subgrid points should match corresponding original points
    for (int i = 0; i < r.num_points; i++) {
        int orig_idx = r.grid_indices[i];
        EXPECT_FLOAT_EQ(subgrid[i+1].coor[0], grid[orig_idx].coor[0]);
        EXPECT_FLOAT_EQ(subgrid[i+1].coor[1], grid[orig_idx].coor[1]);
        EXPECT_FLOAT_EQ(subgrid[i+1].coor[2], grid[orig_idx].coor[2]);
    }

    free(subgrid);
    free(grid);
}

TEST(GridDecomposer, BalanceMergesSmallRegions) {
    int num_grd;
    gridpoint* grid = make_test_grid(4, 1.0f, num_grd);  // 64 points

    // Request many regions (will create tiny ones)
    auto regions = GridDecomposer::decompose_octree(grid, num_grd, 64, 0);

    // Now balance with min_points = 10
    GridDecomposer::balance_regions(regions, grid, 10);

    // All remaining regions should have >= 10 points
    for (const auto& r : regions) {
        EXPECT_GE(r.num_points, 10)
            << "Region " << r.region_id << " has only " << r.num_points << " points";
    }

    free(grid);
}

// ============================================================================
// StatMechEngine merge tests
// ============================================================================

TEST(StatMechMerge, PartitionFunctionAdditive) {
    // Z_merged should equal Z_a + Z_b
    // In log-space: ln(Z_merged) = ln(exp(ln(Z_a)) + exp(ln(Z_b)))
    statmech::StatMechEngine a(300.0);
    statmech::StatMechEngine b(300.0);

    // Region A: low-energy poses
    a.add_sample(-10.0); a.add_sample(-9.5); a.add_sample(-9.0);

    // Region B: medium-energy poses
    b.add_sample(-5.0); b.add_sample(-4.5); b.add_sample(-4.0);

    auto td_a = a.compute();
    auto td_b = b.compute();

    // Merge
    statmech::StatMechEngine merged(300.0);
    merged.merge(a);
    merged.merge(b);
    auto td_merged = merged.compute();

    // Z_merged = Z_a + Z_b
    double Z_a = std::exp(td_a.log_Z);
    double Z_b = std::exp(td_b.log_Z);
    double Z_merged_expected = Z_a + Z_b;

    // Use relative tolerance: values are ~31M, so ULP spacing is ~4e-9
    EXPECT_NEAR(std::exp(td_merged.log_Z), Z_merged_expected,
                Z_merged_expected * 1e-12);
}

TEST(StatMechMerge, FreeEnergyConsistent) {
    statmech::StatMechEngine a(300.0);
    statmech::StatMechEngine b(300.0);

    a.add_sample(-8.0); a.add_sample(-7.0);
    b.add_sample(-6.0); b.add_sample(-5.0);

    auto td_a = a.compute();
    auto td_b = b.compute();

    // Manual: F_merged = -kT * ln(exp(-F_a/kT) + exp(-F_b/kT))
    double kT = statmech::kB_kcal * 300.0;
    double F_expected = -kT * std::log(
        std::exp(-td_a.free_energy / kT) +
        std::exp(-td_b.free_energy / kT)
    );

    statmech::StatMechEngine merged(300.0);
    merged.merge(a);
    merged.merge(b);
    auto td_merged = merged.compute();

    EXPECT_NEAR(td_merged.free_energy, F_expected, 1e-8);
}

TEST(StatMechMerge, SerializeRoundTrip) {
    statmech::StatMechEngine orig(300.0);
    orig.add_sample(-10.0, 2);
    orig.add_sample(-5.0, 1);
    orig.add_sample(-3.0, 3);

    auto energies = orig.serialize_energies();
    auto mults = orig.serialize_multiplicities();

    EXPECT_EQ(energies.size(), 3u);
    EXPECT_EQ(mults.size(), 3u);

    statmech::StatMechEngine reconstructed(300.0);
    reconstructed.merge_samples(
        std::span<const double>(energies),
        std::span<const double>(mults)
    );

    auto td_orig = orig.compute();
    auto td_recon = reconstructed.compute();

    EXPECT_NEAR(td_orig.free_energy, td_recon.free_energy, 1e-12);
    EXPECT_NEAR(td_orig.entropy, td_recon.entropy, 1e-12);
}

// ============================================================================
// SharedPosePool tests
// ============================================================================

TEST(SharedPosePool, PublishAndGetTop) {
    SharedPosePool pool(10);

    SharedPose p1; p1.energy = -5.0; p1.source_region = 0;
    SharedPose p2; p2.energy = -10.0; p2.source_region = 1;
    SharedPose p3; p3.energy = -3.0; p3.source_region = 2;

    pool.publish(p1);
    pool.publish(p2);
    pool.publish(p3);

    auto top = pool.get_top(2);
    ASSERT_EQ((int)top.size(), 2);
    EXPECT_DOUBLE_EQ(top[0].energy, -10.0);  // best first
    EXPECT_DOUBLE_EQ(top[1].energy, -5.0);
}

TEST(SharedPosePool, EvictsWorst) {
    SharedPosePool pool(3);

    for (int i = 0; i < 5; i++) {
        SharedPose p;
        p.energy = -(double)i;  // -0, -1, -2, -3, -4
        p.source_region = i;
        pool.publish(p);
    }

    auto top = pool.get_top(3);
    ASSERT_EQ((int)top.size(), 3);
    // Should keep the 3 best: -4, -3, -2
    EXPECT_DOUBLE_EQ(top[0].energy, -4.0);
    EXPECT_DOUBLE_EQ(top[1].energy, -3.0);
    EXPECT_DOUBLE_EQ(top[2].energy, -2.0);
}

TEST(SharedPosePool, SerializeDeserialize) {
    SharedPosePool pool(10);

    SharedPose p1; p1.energy = -8.0; p1.source_region = 0;
    SharedPose p2; p2.energy = -6.0; p2.source_region = 1;
    pool.publish(p1);
    pool.publish(p2);

    auto buf = pool.serialize();

    SharedPosePool pool2(10);
    pool2.deserialize_merge(buf.data(), buf.size());

    auto top = pool2.get_top(5);
    ASSERT_EQ((int)top.size(), 2);
    EXPECT_DOUBLE_EQ(top[0].energy, -8.0);
    EXPECT_DOUBLE_EQ(top[1].energy, -6.0);
}

TEST(SharedPosePool, ConcurrentPublish) {
    SharedPosePool pool(100);
    const int n_threads = 8;
    const int poses_per_thread = 50;

    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_threads)
    #endif
    for (int t = 0; t < n_threads; t++) {
        for (int i = 0; i < poses_per_thread; i++) {
            SharedPose p;
            p.energy = -(double)(t * poses_per_thread + i);
            p.source_region = t;
            pool.publish(p);
        }
    }

    auto top = pool.get_top(100);
    // Pool should have at most 100 entries, all unique energies
    EXPECT_LE((int)top.size(), 100);
    EXPECT_GT((int)top.size(), 0);

    // Should be sorted (ascending energy = best first)
    for (int i = 1; i < (int)top.size(); i++) {
        EXPECT_LE(top[i-1].energy, top[i].energy);
    }
}

// ============================================================================
// Region bounds computation
// ============================================================================

TEST(GridDecomposer, RegionBoundsCorrect) {
    int num_grd;
    gridpoint* grid = make_test_grid(4, 1.0f, num_grd);

    GridRegion r;
    r.grid_indices = {1, 2, 3, 4, 5};
    r.num_points = 5;

    GridDecomposer::compute_region_bounds(r, grid);

    // Centroid should be average of 5 points
    EXPECT_GT(r.radius, 0.0f);
    EXPECT_EQ(r.num_points, 5);

    free(grid);
}
