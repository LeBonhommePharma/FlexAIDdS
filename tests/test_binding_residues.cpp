// tests/test_binding_residues.cpp
// Unit tests for BindingResidues — key residue identification from MIF scores

#include <gtest/gtest.h>
#include "BindingResidues.h"
#include "MIFGrid.h"
#include <cstring>
#include <cmath>

// ===========================================================================
// HELPERS
// ===========================================================================

// Build a simple test setup: protein atoms around a binding site with grid points
struct TestSetup {
    std::vector<gridpoint> grid;
    std::vector<atom> atoms;
    std::vector<resid> residues;
    std::vector<float> mif_energies;
    cavity_detect::SpatialGrid spatial_grid;

    void build() {
        spatial_grid.build(atoms);
    }
};

static TestSetup make_binding_site() {
    TestSetup ts;

    // 3 residues: ASP (favorable), GLY (neutral), PHE (favorable)
    ts.residues.resize(3);
    memset(ts.residues.data(), 0, 3 * sizeof(resid));

    strncpy(ts.residues[0].name, "ASP", 3); ts.residues[0].name[3] = '\0';
    ts.residues[0].number = 42; ts.residues[0].chn = 'A';

    strncpy(ts.residues[1].name, "GLY", 3); ts.residues[1].name[3] = '\0';
    ts.residues[1].number = 43; ts.residues[1].chn = 'A';

    strncpy(ts.residues[2].name, "PHE", 3); ts.residues[2].name[3] = '\0';
    ts.residues[2].number = 44; ts.residues[2].chn = 'A';

    // 6 atoms: 2 per residue
    ts.atoms.resize(6);
    memset(ts.atoms.data(), 0, 6 * sizeof(atom));

    // ASP atoms near favorable grid points
    ts.atoms[0].coor[0] = 1.0f; ts.atoms[0].coor[1] = 0.0f; ts.atoms[0].coor[2] = 0.0f;
    ts.atoms[0].radius = 1.7f; ts.atoms[0].ofres = 0;
    ts.atoms[1].coor[0] = 1.5f; ts.atoms[1].coor[1] = 0.5f; ts.atoms[1].coor[2] = 0.0f;
    ts.atoms[1].radius = 1.5f; ts.atoms[1].ofres = 0;

    // GLY atoms far from favorable grid points
    ts.atoms[2].coor[0] = 20.0f; ts.atoms[2].coor[1] = 20.0f; ts.atoms[2].coor[2] = 20.0f;
    ts.atoms[2].radius = 1.7f; ts.atoms[2].ofres = 1;
    ts.atoms[3].coor[0] = 21.0f; ts.atoms[3].coor[1] = 20.0f; ts.atoms[3].coor[2] = 20.0f;
    ts.atoms[3].radius = 1.5f; ts.atoms[3].ofres = 1;

    // PHE atoms near moderately favorable grid points
    ts.atoms[4].coor[0] = -1.0f; ts.atoms[4].coor[1] = 0.0f; ts.atoms[4].coor[2] = 0.0f;
    ts.atoms[4].radius = 1.8f; ts.atoms[4].ofres = 2;
    ts.atoms[5].coor[0] = -1.5f; ts.atoms[5].coor[1] = 0.5f; ts.atoms[5].coor[2] = 0.0f;
    ts.atoms[5].radius = 1.8f; ts.atoms[5].ofres = 2;

    // Grid: 5 points (index 0 unused)
    ts.grid.resize(6);
    memset(ts.grid.data(), 0, 6 * sizeof(gridpoint));
    // Point 1: very favorable, near ASP
    ts.grid[1].coor[0] = 2.0f; ts.grid[1].coor[1] = 0.0f; ts.grid[1].coor[2] = 0.0f;
    ts.grid[1].index = 1;
    // Point 2: favorable, near ASP
    ts.grid[2].coor[0] = 2.5f; ts.grid[2].coor[1] = 0.5f; ts.grid[2].coor[2] = 0.0f;
    ts.grid[2].index = 2;
    // Point 3: moderately favorable, near PHE
    ts.grid[3].coor[0] = -2.0f; ts.grid[3].coor[1] = 0.0f; ts.grid[3].coor[2] = 0.0f;
    ts.grid[3].index = 3;
    // Point 4: neutral, near nothing
    ts.grid[4].coor[0] = 10.0f; ts.grid[4].coor[1] = 10.0f; ts.grid[4].coor[2] = 10.0f;
    ts.grid[4].index = 4;
    // Point 5: unfavorable
    ts.grid[5].coor[0] = 0.0f; ts.grid[5].coor[1] = 0.0f; ts.grid[5].coor[2] = 0.0f;
    ts.grid[5].index = 5;

    // MIF energies: [0]=unused, [1]=-5.0, [2]=-3.0, [3]=-1.5, [4]=0.0, [5]=+2.0
    ts.mif_energies = {0.0f, -5.0f, -3.0f, -1.5f, 0.0f, 2.0f};

    ts.build();
    return ts;
}

// ===========================================================================
// TESTS
// ===========================================================================

TEST(BindingResidues, EmptyGrid) {
    gridpoint grid[1];
    memset(grid, 0, sizeof(grid));
    float mif[] = {0.0f};

    atom atoms[1];
    memset(atoms, 0, sizeof(atoms));
    resid residues[1];
    memset(residues, 0, sizeof(residues));

    cavity_detect::SpatialGrid sg;
    std::vector<atom> av(atoms, atoms + 1);
    sg.build(av);

    auto results = binding_residues::identify_key_residues(
        grid, 1, mif, atoms, 1, residues, sg);
    EXPECT_TRUE(results.empty());
}

TEST(BindingResidues, NullMIF) {
    gridpoint grid[2];
    memset(grid, 0, sizeof(grid));

    auto results = binding_residues::identify_key_residues(
        grid, 2, nullptr, nullptr, 0, nullptr,
        cavity_detect::SpatialGrid{});
    EXPECT_TRUE(results.empty());
}

TEST(BindingResidues, IdentifiesASPAsMostFavorable) {
    auto ts = make_binding_site();

    auto results = binding_residues::identify_key_residues(
        ts.grid.data(), static_cast<int>(ts.grid.size()),
        ts.mif_energies.data(),
        ts.atoms.data(), static_cast<int>(ts.atoms.size()),
        ts.residues.data(), ts.spatial_grid,
        100.0f,   // top 100% — consider all grid points
        4.5f);

    // ASP should be first (most favorable MIF score)
    ASSERT_GE(results.size(), 1u);
    EXPECT_STREQ(results[0].name, "ASP");
    EXPECT_EQ(results[0].number, 42);
    EXPECT_EQ(results[0].chain, 'A');
    EXPECT_LT(results[0].mif_score, 0.0f);
}

TEST(BindingResidues, PHEIsSecondMostFavorable) {
    auto ts = make_binding_site();

    auto results = binding_residues::identify_key_residues(
        ts.grid.data(), static_cast<int>(ts.grid.size()),
        ts.mif_energies.data(),
        ts.atoms.data(), static_cast<int>(ts.atoms.size()),
        ts.residues.data(), ts.spatial_grid,
        100.0f, 4.5f);

    ASSERT_GE(results.size(), 2u);
    EXPECT_STREQ(results[0].name, "ASP");
    EXPECT_STREQ(results[1].name, "PHE");
}

TEST(BindingResidues, GLYNotInResults) {
    // GLY is far from all favorable grid points — should not appear
    auto ts = make_binding_site();

    auto results = binding_residues::identify_key_residues(
        ts.grid.data(), static_cast<int>(ts.grid.size()),
        ts.mif_energies.data(),
        ts.atoms.data(), static_cast<int>(ts.atoms.size()),
        ts.residues.data(), ts.spatial_grid,
        100.0f, 4.5f);

    for (const auto& r : results) {
        EXPECT_STRNE(r.name, "GLY") << "GLY should not be near favorable grid points";
    }
}

TEST(BindingResidues, ContactCountPositive) {
    auto ts = make_binding_site();

    auto results = binding_residues::identify_key_residues(
        ts.grid.data(), static_cast<int>(ts.grid.size()),
        ts.mif_energies.data(),
        ts.atoms.data(), static_cast<int>(ts.atoms.size()),
        ts.residues.data(), ts.spatial_grid,
        100.0f, 4.5f);

    for (const auto& r : results) {
        EXPECT_GT(r.contact_count, 0);
    }
}

TEST(BindingResidues, MinDistanceIsReasonable) {
    auto ts = make_binding_site();

    auto results = binding_residues::identify_key_residues(
        ts.grid.data(), static_cast<int>(ts.grid.size()),
        ts.mif_energies.data(),
        ts.atoms.data(), static_cast<int>(ts.atoms.size()),
        ts.residues.data(), ts.spatial_grid,
        100.0f, 4.5f);

    for (const auto& r : results) {
        EXPECT_GT(r.min_distance, 0.0f);
        EXPECT_LE(r.min_distance, 4.5f);
    }
}

TEST(BindingResidues, TopKFilterReducesResults) {
    auto ts = make_binding_site();

    // With 100% — all favorable grid points considered
    auto all_results = binding_residues::identify_key_residues(
        ts.grid.data(), static_cast<int>(ts.grid.size()),
        ts.mif_energies.data(),
        ts.atoms.data(), static_cast<int>(ts.atoms.size()),
        ts.residues.data(), ts.spatial_grid,
        100.0f, 4.5f);

    // With 20% — only the most favorable grid point
    auto top_results = binding_residues::identify_key_residues(
        ts.grid.data(), static_cast<int>(ts.grid.size()),
        ts.mif_energies.data(),
        ts.atoms.data(), static_cast<int>(ts.atoms.size()),
        ts.residues.data(), ts.spatial_grid,
        20.0f, 4.5f);

    // Top results should have fewer or equal residues
    EXPECT_LE(top_results.size(), all_results.size());
}

TEST(BindingResidues, SmallCutoffExcludesDistantAtoms) {
    auto ts = make_binding_site();

    // Very tight cutoff — only atoms very close to grid points
    auto results = binding_residues::identify_key_residues(
        ts.grid.data(), static_cast<int>(ts.grid.size()),
        ts.mif_energies.data(),
        ts.atoms.data(), static_cast<int>(ts.atoms.size()),
        ts.residues.data(), ts.spatial_grid,
        100.0f, 0.5f);  // 0.5 Å cutoff — very tight

    // With such a tight cutoff, we should get very few or no results
    // since atoms are at least 1.0 Å from grid points
    EXPECT_LE(results.size(), 1u);
}

TEST(BindingResidues, ScoreIsSortedDescending) {
    auto ts = make_binding_site();

    auto results = binding_residues::identify_key_residues(
        ts.grid.data(), static_cast<int>(ts.grid.size()),
        ts.mif_energies.data(),
        ts.atoms.data(), static_cast<int>(ts.atoms.size()),
        ts.residues.data(), ts.spatial_grid,
        100.0f, 4.5f);

    for (size_t i = 1; i < results.size(); ++i) {
        EXPECT_LE(results[i - 1].mif_score, results[i].mif_score)
            << "Results should be sorted by MIF score (ascending = most favorable first)";
    }
}

TEST(BindingResidues, PrintDoesNotCrash) {
    auto ts = make_binding_site();

    auto results = binding_residues::identify_key_residues(
        ts.grid.data(), static_cast<int>(ts.grid.size()),
        ts.mif_energies.data(),
        ts.atoms.data(), static_cast<int>(ts.atoms.size()),
        ts.residues.data(), ts.spatial_grid,
        100.0f, 4.5f);

    // Just verify it doesn't crash
    binding_residues::print_key_residues(results, 5);
    binding_residues::print_key_residues({}, 5);
}

// ===========================================================================
// AUTO-FLEX TESTS
// ===========================================================================

// Helper to create a minimal FA_Global for auto-flex testing
static FA_Global make_test_fa(float* mif_energies, int mif_count, int atm_cnt) {
    FA_Global fa;
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wnontrivial-memcall"
    std::memset(&fa, 0, sizeof(FA_Global));
#pragma clang diagnostic pop
    fa.mif_energies = mif_energies;
    fa.mif_count = mif_count;
    fa.atm_cnt_real = atm_cnt;
    fa.autoflex_enabled = 1;
    fa.autoflex_max = 5;
    fa.MIN_FLEX_RESIDUE = 10;
    fa.flex_res = nullptr;
    fa.nflxsc = 0;
    // num_grd must match grid size
    return fa;
}

TEST(AutoFlex, AddsASPAsFlexible) {
    auto ts = make_binding_site();
    auto fa = make_test_fa(ts.mif_energies.data(), 5,
                           static_cast<int>(ts.atoms.size()));
    fa.num_grd = static_cast<int>(ts.grid.size());

    // Mark all residues as protein type
    for (auto& r : ts.residues) r.type = 0;

    int added = binding_residues::add_key_residues_as_flexible(
        &fa, ts.grid.data(), ts.atoms.data(), ts.residues.data(),
        5, 100.0f, 4.5f);

    EXPECT_GE(added, 1);
    EXPECT_EQ(fa.nflxsc, added);
    // ASP should be first (most favorable)
    EXPECT_STREQ(fa.flex_res[0].name, "ASP");
    EXPECT_EQ(fa.flex_res[0].num, 42);
    EXPECT_EQ(fa.flex_res[0].chn, 'A');
    EXPECT_GT(fa.flex_res[0].prob, 0.0f);  // set_intprob should have set it

    free(fa.flex_res);
}

TEST(AutoFlex, SkipsGLYALAPRO) {
    auto ts = make_binding_site();
    // Change ASP to GLY
    strncpy(ts.residues[0].name, "GLY", 3);
    // Change PHE to ALA
    strncpy(ts.residues[2].name, "ALA", 3);

    auto fa = make_test_fa(ts.mif_energies.data(), 5,
                           static_cast<int>(ts.atoms.size()));
    fa.num_grd = static_cast<int>(ts.grid.size());
    for (auto& r : ts.residues) r.type = 0;

    int added = binding_residues::add_key_residues_as_flexible(
        &fa, ts.grid.data(), ts.atoms.data(), ts.residues.data(),
        5, 100.0f, 4.5f);

    // GLY and ALA should be skipped
    EXPECT_EQ(added, 0);
    free(fa.flex_res);
}

TEST(AutoFlex, SkipsAlreadyFlexible) {
    auto ts = make_binding_site();
    auto fa = make_test_fa(ts.mif_energies.data(), 5,
                           static_cast<int>(ts.atoms.size()));
    fa.num_grd = static_cast<int>(ts.grid.size());
    for (auto& r : ts.residues) r.type = 0;

    // Pre-add ASP as flexible
    fa.flex_res = static_cast<flxsc*>(calloc(10, sizeof(flxsc)));
    fa.nflxsc = 1;
    strncpy(fa.flex_res[0].name, "ASP", 3);
    fa.flex_res[0].inum = 0;  // residue index 0 = ASP

    int added = binding_residues::add_key_residues_as_flexible(
        &fa, ts.grid.data(), ts.atoms.data(), ts.residues.data(),
        5, 100.0f, 4.5f);

    // ASP should be skipped (already flexible), PHE should be added
    EXPECT_GE(added, 1);
    // PHE should be the newly added one
    bool found_phe = false;
    for (int i = 1; i < fa.nflxsc; ++i) {
        if (strcmp(fa.flex_res[i].name, "PHE") == 0) found_phe = true;
    }
    EXPECT_TRUE(found_phe);

    free(fa.flex_res);
}

TEST(AutoFlex, RespectsMaxLimit) {
    auto ts = make_binding_site();
    auto fa = make_test_fa(ts.mif_energies.data(), 5,
                           static_cast<int>(ts.atoms.size()));
    fa.num_grd = static_cast<int>(ts.grid.size());
    for (auto& r : ts.residues) r.type = 0;

    int added = binding_residues::add_key_residues_as_flexible(
        &fa, ts.grid.data(), ts.atoms.data(), ts.residues.data(),
        1,  // max 1
        100.0f, 4.5f);

    EXPECT_EQ(added, 1);
    EXPECT_EQ(fa.nflxsc, 1);

    free(fa.flex_res);
}

TEST(AutoFlex, DisabledWhenNoMIF) {
    auto ts = make_binding_site();
    auto fa = make_test_fa(nullptr, 0,
                           static_cast<int>(ts.atoms.size()));
    fa.num_grd = static_cast<int>(ts.grid.size());

    int added = binding_residues::add_key_residues_as_flexible(
        &fa, ts.grid.data(), ts.atoms.data(), ts.residues.data(), 5);

    EXPECT_EQ(added, 0);
}
