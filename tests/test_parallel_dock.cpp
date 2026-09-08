// tests/test_parallel_dock.cpp — Unit tests for parallel grid-decomposed docking
// Tests octree decomposition, StatMechEngine merge, SharedPosePool, and subgrid extraction

#include <gtest/gtest.h>
#include <cmath>
#include <set>
#include <numeric>
#include <latch>
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <memory>
#include <string>
#include <thread>
#include <vector>
#include <barrier>
#include "../LIB/AtomCopyExtent.h"
#include "../LIB/SerialRegionExecution.h"
#include "../LIB/RegionScoringWorkspace.h"
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

// ============================================================================
// 1-based atom[] copy (ParallelDock create_workspace)
// ============================================================================
// FlexAID stores live atoms at indices [1, atm_cnt]. Gaboom ParEvalWS copies
// `atoms + natm + 1`. The old ParallelDock line was:
//   ws.atoms_copy.assign(atoms_, atoms_ + FA_->atm_cnt);
// which is a half-open range of length atm_cnt and drops atoms[atm_cnt].
// flexaid_one_based_copy_n is what create_workspace uses; if it returned
// atm_cnt this test fails.

TEST(ParallelDockWorkspace, AtomCopyExtentIsAtmCntPlusOne) {
    EXPECT_EQ(flexaid_one_based_copy_n(0), 1);
    EXPECT_EQ(flexaid_one_based_copy_n(1), 2);
    EXPECT_EQ(flexaid_one_based_copy_n(7), 8);
    EXPECT_NE(flexaid_one_based_copy_n(7), 7);
}

TEST(ParallelDockWorkspace, AtomCopyIncludesLastOneBasedAtom) {
    const int atm_cnt = 3;
    atom atoms[4]{};
    atoms[1].number = 101;
    atoms[2].number = 102;
    atoms[3].number = 103;
    atoms[3].coor[0] = 42.0f;

    std::vector<atom> dropped(atoms, atoms + atm_cnt);
    std::vector<atom> kept(atoms, atoms + flexaid_one_based_copy_n(atm_cnt));

    ASSERT_EQ(dropped.size(), static_cast<size_t>(atm_cnt));
    ASSERT_EQ(kept.size(), static_cast<size_t>(atm_cnt) + 1);
    EXPECT_EQ(dropped.back().number, 102) << "assign(..., atm_cnt) stops at atoms[atm_cnt-1]";
    EXPECT_EQ(kept[atm_cnt].number, 103);
    EXPECT_FLOAT_EQ(kept[atm_cnt].coor[0], 42.0f);
}

// ============================================================================
// Production region scheduling and scratch ownership
// ============================================================================

namespace {
void set_region_test_seed_env(const char* value) {
#ifdef _WIN32
    _putenv_s("FLEXAID_SEED", value ? value : "");
#else
    if (value) setenv("FLEXAID_SEED", value, 1);
    else unsetenv("FLEXAID_SEED");
#endif
}

struct RegionSeedStateGuard {
    const bool had_env = std::getenv("FLEXAID_SEED") != nullptr;
    const std::string env = had_env ? std::getenv("FLEXAID_SEED") : "";
    const bool had_master = flexaids_rng::has_master_seed();
    const std::uint64_t master = flexaids_rng::master_seed();

    ~RegionSeedStateGuard() {
        set_region_test_seed_env(had_env ? env.c_str() : nullptr);
        flexaids_rng::set_master_seed(master);
        flexaids_rng::g_has_master_seed.store(had_master);
    }
};
}  // namespace

TEST(SerialRegionExecution, ExplicitParentSeedSurvivesPriorRegionEpochs) {
    RegionSeedStateGuard restore;
    set_region_test_seed_env("12345");
    std::uint64_t first = 0, repeated = 0;
    flexaids_rng::set_master_seed(98765);
    EXPECT_TRUE(flexaids::resolve_region_parent_seed(0, first));
    EXPECT_EQ(first, 12345u);

    // Exercise the same mutation each actual region GA performs.
    flexaids_rng::set_master_seed(flexaids::region_seed(first, 3));
    EXPECT_TRUE(flexaids::resolve_region_parent_seed(0, repeated));
    EXPECT_EQ(repeated, first);
    EXPECT_EQ(flexaids::region_seed(repeated, 0), flexaids::region_seed(first, 0));

    // A configured GB seed still wins; a master seed is only a fallback.
    EXPECT_TRUE(flexaids::resolve_region_parent_seed(444, repeated));
    EXPECT_EQ(repeated, 444u);
    set_region_test_seed_env(nullptr);
    flexaids_rng::set_master_seed(777);
    EXPECT_TRUE(flexaids::resolve_region_parent_seed(0, repeated));
    EXPECT_EQ(repeated, 777u);
}

TEST(SerialRegionExecution, ConcurrentDispatchersCannotOverlapRegionCallbacks) {
    // Same dispatcher as ParallelDockManager::run; no stand-in thread_local.
    std::atomic<int> active{0}, peak{0}, completed{0};
    std::latch ready(2);
    auto dispatch = [&] {
        ready.arrive_and_wait();
        flexaids::for_each_serial_region(0, 3, 1, [&](int region) {
            const int now = ++active;
            int prior = peak.load();
            while (prior < now && !peak.compare_exchange_weak(prior, now)) {}
            flexaids_rng::set_master_seed(flexaids::region_seed(12345, region));
            const auto epoch = flexaids_rng::g_seed_epoch.load();
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
            EXPECT_EQ(flexaids_rng::g_seed_epoch.load(), epoch);
            --active;
            ++completed;
        });
    };
    std::thread first(dispatch), second(dispatch);
    first.join();
    second.join();
    EXPECT_EQ(completed.load(), 6);
    EXPECT_EQ(peak.load(), 1);
}

TEST(SerialRegionExecution, RegionSeedsFollowGlobalIndexAcrossRankAssignments) {
    std::vector<int> local(8), distributed(8);
    flexaids::for_each_serial_region(0, 8, 1, [&](int region) {
        local[region] = flexaids::region_seed(12345, region);
        EXPECT_GT(local[region], 0);
    });
    for (int rank = 0; rank < 3; ++rank)
        flexaids::for_each_serial_region(rank, 8, 3, [&](int region) {
            distributed[region] = flexaids::region_seed(12345, region);
        });
    EXPECT_EQ(distributed, local);
    EXPECT_EQ(std::set<int>(local.begin(), local.end()).size(), local.size());
    EXPECT_NE(flexaids::region_seed(98765, 0), local[0]);
}

TEST(RegionScoringWorkspace, OwnsScratchRebindsOptresAndSupportsContactReallocation) {
    auto fa = std::make_unique<FA_Global>();
    GB_Global gb{};
    VC_Global vc{};
    atom atoms[3]{};
    resid residues[2]{};
    OptRes optres[1]{};
    fa->atm_cnt = fa->atm_cnt_real = 2;
    fa->res_cnt = fa->num_optres = fa->ntypes = 1;
    fa->optres = optres;
    atoms[2].optres = &optres[0];
    atoms[2].coor[0] = 17.0f;
    optres[0].cf.com = 42.0;
    ca_struct parent_contact{};
    parent_contact.area = 8.5;
    vc.ca_rec = &parent_contact; // deliberately non-heap: must never free/realloc parent
    vc.ca_recsize = 100;
    auto first = std::make_unique<flexaids::RegionScoringWorkspace>(*fa, gb, vc, atoms, residues);
    auto second = std::make_unique<flexaids::RegionScoringWorkspace>(*fa, gb, vc, atoms, residues);
    EXPECT_EQ(first->atoms_copy[2].optres, &first->optres[0]);
    EXPECT_EQ(second->atoms_copy[2].optres, &second->optres[0]);
    EXPECT_EQ(first->atoms_copy[1].optres, nullptr);
    EXPECT_NE(first->fa.optres, fa->optres);
    EXPECT_NE(first->fa.contacts, second->fa.contacts);
    EXPECT_NE(first->fa.contributions, second->fa.contributions);
    EXPECT_NE(first->vc.Calc, second->vc.Calc);
    EXPECT_NE(first->vc.Calclist, second->vc.Calclist);
    EXPECT_NE(first->vc.ca_rec, vc.ca_rec);
    EXPECT_NE(first->vc.ca_index, second->vc.ca_index);
    EXPECT_NE(first->vc.seed, second->vc.seed);
    EXPECT_NE(first->vc.contlist, second->vc.contlist);
    EXPECT_NE(first->vc.ptorder, second->vc.ptorder);
    EXPECT_NE(first->vc.centerpt, second->vc.centerpt);
    EXPECT_NE(first->vc.poly, second->vc.poly);
    EXPECT_NE(first->vc.cont, second->vc.cont);
    EXPECT_NE(first->vc.vedge, second->vc.vedge);
    EXPECT_NE(first->vc.scorable_list, second->vc.scorable_list);
    first->atoms_copy[2].optres->cf.com = -5.0;
    first->vc.Calc[0].atom = &first->atoms_copy[2];
    EXPECT_DOUBLE_EQ(optres[0].cf.com, 42.0);
    EXPECT_DOUBLE_EQ(second->optres[0].cf.com, 42.0);
    // save_areas() is permitted to realloc this exact allocation. Using a
    // vector's data() here is an allocator mismatch, even without concurrency.
    auto* grown = static_cast<ca_struct*>(std::realloc(first->vc.ca_rec, 10100 * sizeof(ca_struct)));
    ASSERT_NE(grown, nullptr);
    first->vc.ca_rec = grown;
    first->vc.ca_recsize = 10100;
    first.reset();
    EXPECT_FLOAT_EQ(parent_contact.area, 8.5f);
    EXPECT_EQ(atoms[2].optres, optres);
    EXPECT_FLOAT_EQ(atoms[2].coor[0], 17.0f);
}
