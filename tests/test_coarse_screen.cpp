// test_coarse_screen.cpp — Unit tests for NRGRank CoarseScreen + TwoStageScreen
//
// Tests cover:
//   1. Energy matrix integrity (nrgrank_matrix.h)
//   2. SYBYL type lookup
//   3. IndexCubeGrid construction + neighbor queries
//   4. CF precomputation (get_cf_list kernel)
//   5. Rotation generation (729 orientations, dedup)
//   6. Anchor point generation + clash cleaning
//   7. Single-ligand scoring
//   8. Batch screening + sorting
//   9. TwoStageScreener pipeline
//  10. MOL2/SDF file parsing round-trip
//
// Method (clean-room): DesCôteaux, Mailhot, Najmanovich,
// bioRxiv 2025.02.17.638675v2; Zenodo 10.5281/zenodo.16861024.
//
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>
#include "CoarseScreen.h"
#include "TwoStageScreen.h"
#include "nrgrank_matrix.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numeric>
#include <vector>

using namespace nrgrank;

// ═══════════════════════════════════════════════════════════════════════
//  Test helpers
// ═══════════════════════════════════════════════════════════════════════

/// Create a small synthetic "protein" with a few atoms in a known geometry.
static std::vector<TargetAtom> make_test_target() {
    std::vector<TargetAtom> atoms;
    // Approximate a binding pocket: atoms arranged in a cube from (0,0,0) to (20,20,20)
    // with different SYBYL types
    int type_cycle[] = {3, 4, 14, 10, 13}; // C.3, C.AR, O.3, N.AR, O.2
    float radii[]    = {1.88f, 1.76f, 1.46f, 1.64f, 1.42f};
    int idx = 0;
    for (float x = 0.f; x <= 20.f; x += 4.f) {
        for (float y = 0.f; y <= 20.f; y += 4.f) {
            for (float z = 0.f; z <= 20.f; z += 4.f) {
                TargetAtom ta;
                ta.pos = {x, y, z};
                ta.type   = type_cycle[idx % 5];
                ta.radius = radii[idx % 5];
                atoms.push_back(ta);
                idx++;
            }
        }
    }
    return atoms;
}

/// Create binding site spheres covering (5,5,5)–(15,15,15)
static std::vector<BindingSiteSphere> make_test_spheres() {
    std::vector<BindingSiteSphere> spheres;
    // 8 spheres at corners of (5,5,5)–(15,15,15) with radius 3.0
    for (float x : {5.f, 15.f})
        for (float y : {5.f, 15.f})
            for (float z : {5.f, 15.f})
                spheres.push_back({{x, y, z}, 3.0f});
    return spheres;
}

/// Create a simple test ligand (3-atom "molecule")
static ScreenLigand make_test_ligand(const std::string& name = "test_lig") {
    ScreenLigand lig;
    lig.name = name;
    // Triangle: C.3, N.AR, O.2 — ~2 Å spacing
    lig.atoms.push_back({{0.f, 0.f, 0.f}, 3});   // C.3
    lig.atoms.push_back({{2.f, 0.f, 0.f}, 10});  // N.AR
    lig.atoms.push_back({{1.f, 1.7f, 0.f}, 13}); // O.2
    return lig;
}

// ═══════════════════════════════════════════════════════════════════════
//  1. Energy matrix
// ═══════════════════════════════════════════════════════════════════════

TEST(NRGRankMatrix, Dimensions) {
    EXPECT_EQ(MATRIX_DIM, 41);
    EXPECT_EQ(NUM_ATOM_TYPES, 39);
}

TEST(NRGRankMatrix, SymmetrySpot) {
    // Check a few known symmetric entries
    EXPECT_DOUBLE_EQ(kEnergyMatrix[1][2], kEnergyMatrix[2][1]); // C.1↔C.2
    EXPECT_DOUBLE_EQ(kEnergyMatrix[3][13], kEnergyMatrix[13][3]); // C.3↔O.2
    EXPECT_DOUBLE_EQ(kEnergyMatrix[10][11], kEnergyMatrix[11][10]); // N.AR↔N.AM
}

TEST(NRGRankMatrix, KnownValues) {
    // C.1↔C.2 = -170.52
    EXPECT_NEAR(kEnergyMatrix[1][2], -170.52, 0.01);
    // C.AR↔C.AR = -397.40
    EXPECT_NEAR(kEnergyMatrix[4][4], -397.40, 0.01);
    // O.3↔O.3 = 133.30
    EXPECT_NEAR(kEnergyMatrix[14][14], 133.30, 0.01);
    // Row/col 0 should be all zeros
    for (int i = 0; i < MATRIX_DIM; ++i) {
        EXPECT_DOUBLE_EQ(kEnergyMatrix[0][i], 0.0);
        EXPECT_DOUBLE_EQ(kEnergyMatrix[i][0], 0.0);
    }
}

TEST(NRGRankMatrix, FloatVersion) {
    const auto& mf = EnergyMatrixF::instance();
    EXPECT_NEAR(mf.data[1][2], -170.52f, 0.01f);
    EXPECT_NEAR(mf.data[4][4], -397.40f, 0.01f);
    // Verify all entries match
    for (int i = 0; i < MATRIX_DIM; ++i)
        for (int j = 0; j < MATRIX_DIM; ++j)
            EXPECT_NEAR(mf.data[i][j], static_cast<float>(kEnergyMatrix[i][j]), 0.01f);
}

// ═══════════════════════════════════════════════════════════════════════
//  2. SYBYL type lookup
// ═══════════════════════════════════════════════════════════════════════

TEST(SybylLookup, KnownTypes) {
    float rad;
    EXPECT_EQ(sybyl_type_lookup("C.3",  rad), 3);  EXPECT_NEAR(rad, 1.88f, 0.01f);
    EXPECT_EQ(sybyl_type_lookup("N.AR", rad), 10); EXPECT_NEAR(rad, 1.64f, 0.01f);
    EXPECT_EQ(sybyl_type_lookup("O.2",  rad), 13); EXPECT_NEAR(rad, 1.42f, 0.01f);
    EXPECT_EQ(sybyl_type_lookup("ZN",   rad), 35); EXPECT_NEAR(rad, 0.74f, 0.01f);
}

TEST(SybylLookup, CaseInsensitive) {
    float rad;
    EXPECT_EQ(sybyl_type_lookup("c.3",  rad), 3);
    EXPECT_EQ(sybyl_type_lookup("n.ar", rad), 10);
    EXPECT_EQ(sybyl_type_lookup("C.Ar", rad), 4); // C.AR
}

TEST(SybylLookup, UnknownReturnsDummy) {
    float rad;
    EXPECT_EQ(sybyl_type_lookup("XYZ", rad), 39);
    EXPECT_NEAR(rad, 2.0f, 0.01f);
}

// ═══════════════════════════════════════════════════════════════════════
//  3. IndexCubeGrid
// ═══════════════════════════════════════════════════════════════════════

TEST(IndexCubeGrid, Construction) {
    auto atoms = make_test_target();
    IndexCubeGrid grid;
    grid.build(atoms, 6.56f);

    EXPECT_GT(grid.nx(), 0);
    EXPECT_GT(grid.ny(), 0);
    EXPECT_GT(grid.nz(), 0);
    EXPECT_GT(grid.max_atoms_per_cell(), 0);
    EXPECT_NEAR(grid.cell_width(), 6.56f, 0.001f);
}

TEST(IndexCubeGrid, AllAtomsIndexed) {
    auto atoms = make_test_target();
    IndexCubeGrid grid;
    grid.build(atoms, 6.56f);

    // Count all atoms across all cells
    int total = 0;
    for (int x = 0; x < grid.nx(); ++x)
        for (int y = 0; y < grid.ny(); ++y)
            for (int z = 0; z < grid.nz(); ++z) {
                const int* cell = grid.cell(x, y, z);
                for (int s = 0; s < grid.max_atoms_per_cell(); ++s) {
                    if (cell[s] == IndexCubeGrid::kPlaceholder) break;
                    total++;
                }
            }
    EXPECT_EQ(total, static_cast<int>(atoms.size()));
}

TEST(IndexCubeGrid, WorldToGrid) {
    auto atoms = make_test_target();
    IndexCubeGrid grid;
    grid.build(atoms, 6.56f);

    // Atom at (0,0,0) should be in a valid cell
    int ix, iy, iz;
    grid.world_to_grid({0.f, 0.f, 0.f}, ix, iy, iz);
    EXPECT_TRUE(grid.in_bounds(ix, iy, iz));

    // Check the atom is actually in that cell
    const int* cell = grid.cell(ix, iy, iz);
    bool found = false;
    for (int s = 0; s < grid.max_atoms_per_cell(); ++s) {
        int idx = cell[s];
        if (idx == IndexCubeGrid::kPlaceholder) break;
        if (atoms[idx].pos.x == 0.f && atoms[idx].pos.y == 0.f && atoms[idx].pos.z == 0.f) {
            found = true;
            break;
        }
    }
    EXPECT_TRUE(found);
}

// ═══════════════════════════════════════════════════════════════════════
//  4. CF precomputation
// ═══════════════════════════════════════════════════════════════════════

TEST(CoarseScreener, PrecomputeCF) {
    auto atoms   = make_test_target();
    auto spheres = make_test_spheres();

    CoarseScreener cs;
    CoarseScreenConfig cfg;
    cfg.use_clash = false; // simpler for testing
    cs.set_config(cfg);
    cs.prepare_target(atoms, spheres);

    EXPECT_TRUE(cs.is_prepared());
    EXPECT_GT(cs.cf_nx(), 0);
    EXPECT_GT(cs.cf_ny(), 0);
    EXPECT_GT(cs.cf_nz(), 0);
}

// ═══════════════════════════════════════════════════════════════════════
//  5. Rotation generation
// ═══════════════════════════════════════════════════════════════════════

TEST(Rotations, Count9PerAxis) {
    // With 9 per axis = 729 total, dedup removes ~27% typically
    // Python NRGRank gets ~576 unique from 729
    // We just verify we get a reasonable number in [400, 729]
    auto rots = CoarseScreener::generate_rotations(9);
    EXPECT_GE(static_cast<int>(rots.size()), 400);
    EXPECT_LE(static_cast<int>(rots.size()), 729);
}

TEST(Rotations, IdentityIncluded) {
    // With any per_axis >= 1, rotation (0,0,0) = identity should be included
    auto rots = CoarseScreener::generate_rotations(4);
    bool has_identity = false;
    for (const auto& m : rots) {
        // Identity: m[0]=1, m[4]=1, m[8]=1, rest ~0
        if (std::fabs(m[0]-1.f) < 0.01f && std::fabs(m[4]-1.f) < 0.01f &&
            std::fabs(m[8]-1.f) < 0.01f &&
            std::fabs(m[1]) < 0.01f && std::fabs(m[2]) < 0.01f &&
            std::fabs(m[3]) < 0.01f && std::fabs(m[5]) < 0.01f &&
            std::fabs(m[6]) < 0.01f && std::fabs(m[7]) < 0.01f) {
            has_identity = true;
            break;
        }
    }
    EXPECT_TRUE(has_identity);
}

TEST(Rotations, AreOrthonormal) {
    auto rots = CoarseScreener::generate_rotations(9);
    for (const auto& m : rots) {
        // Check each row is unit length
        float r0 = m[0]*m[0] + m[1]*m[1] + m[2]*m[2];
        float r1 = m[3]*m[3] + m[4]*m[4] + m[5]*m[5];
        float r2 = m[6]*m[6] + m[7]*m[7] + m[8]*m[8];
        EXPECT_NEAR(r0, 1.f, 0.01f);
        EXPECT_NEAR(r1, 1.f, 0.01f);
        EXPECT_NEAR(r2, 1.f, 0.01f);
        // Check rows are orthogonal
        float d01 = m[0]*m[3] + m[1]*m[4] + m[2]*m[5];
        float d02 = m[0]*m[6] + m[1]*m[7] + m[2]*m[8];
        float d12 = m[3]*m[6] + m[4]*m[7] + m[5]*m[8];
        EXPECT_NEAR(d01, 0.f, 0.01f);
        EXPECT_NEAR(d02, 0.f, 0.01f);
        EXPECT_NEAR(d12, 0.f, 0.01f);
    }
}

// ═══════════════════════════════════════════════════════════════════════
//  6. Anchor points
// ═══════════════════════════════════════════════════════════════════════

TEST(CoarseScreener, AnchorGeneration) {
    auto atoms   = make_test_target();
    auto spheres = make_test_spheres();

    CoarseScreener cs;
    CoarseScreenConfig cfg;
    cfg.use_clash = false;
    cs.set_config(cfg);
    cs.prepare_target(atoms, spheres);

    EXPECT_TRUE(cs.is_prepared());
    // With 1.5 Å spacing in a 10×10×10 Å box → ~296 points before cleaning
    // After removing those within 2 Å of target atoms, should be fewer
    EXPECT_GT(cs.num_anchors(), 0u);
    // Upper bound: (10/1.5)^3 ≈ 296
    EXPECT_LT(cs.num_anchors(), 500u);
}

// ═══════════════════════════════════════════════════════════════════════
//  7. Single-ligand scoring
// ═══════════════════════════════════════════════════════════════════════

TEST(CoarseScreener, ScoreOne) {
    auto atoms   = make_test_target();
    auto spheres = make_test_spheres();

    CoarseScreener cs;
    CoarseScreenConfig cfg;
    cfg.use_clash = false;
    cfg.rotations_per_axis = 4; // fewer rotations for speed in test
    cs.set_config(cfg);
    cs.prepare_target(atoms, spheres);

    auto lig = make_test_ligand("aspirin");
    auto result = cs.screen_one(lig);

    EXPECT_EQ(result.name, "aspirin");
    // Score should be a real number (not the default sentinel)
    // For a reasonable pocket, we expect a negative score
    EXPECT_LT(result.score, cfg.default_cf);
}

TEST(CoarseScreener, EmptyLigandReturnsDefault) {
    auto atoms   = make_test_target();
    auto spheres = make_test_spheres();

    CoarseScreener cs;
    CoarseScreenConfig cfg;
    cfg.use_clash = false;
    cs.set_config(cfg);
    cs.prepare_target(atoms, spheres);

    ScreenLigand empty;
    empty.name = "empty";
    auto result = cs.screen_one(empty);

    EXPECT_EQ(result.score, cfg.default_cf);
}

// ═══════════════════════════════════════════════════════════════════════
//  8. Batch screening
// ═══════════════════════════════════════════════════════════════════════

TEST(CoarseScreener, BatchScreen) {
    auto atoms   = make_test_target();
    auto spheres = make_test_spheres();

    CoarseScreener cs;
    CoarseScreenConfig cfg;
    cfg.use_clash = false;
    cfg.rotations_per_axis = 4;
    cs.set_config(cfg);
    cs.prepare_target(atoms, spheres);

    // Create 5 different ligands
    std::vector<ScreenLigand> ligs;
    for (int i = 0; i < 5; ++i) {
        auto lig = make_test_ligand("lig_" + std::to_string(i));
        // Vary atom positions slightly
        for (auto& a : lig.atoms) {
            a.pos.x += i * 0.5f;
            a.pos.y += i * 0.3f;
        }
        ligs.push_back(std::move(lig));
    }

    auto results = cs.screen(ligs);

    EXPECT_EQ(results.size(), 5u);
    // Results should be sorted by score
    for (size_t i = 1; i < results.size(); ++i)
        EXPECT_LE(results[i-1].score, results[i].score);
}

// ═══════════════════════════════════════════════════════════════════════
//  9. TwoStageScreener
// ═══════════════════════════════════════════════════════════════════════

TEST(TwoStageScreener, PipelineNoDocking) {
    auto atoms   = make_test_target();
    auto spheres = make_test_spheres();

    TwoStageScreener ts;
    TwoStageConfig cfg;
    cfg.coarse.use_clash = false;
    cfg.coarse.rotations_per_axis = 4;
    cfg.top_n = 3;
    cfg.write_coarse_csv = false;
    ts.set_config(cfg);
    ts.prepare_target(atoms, spheres);

    std::vector<ScreenLigand> ligs;
    for (int i = 0; i < 5; ++i)
        ligs.push_back(make_test_ligand("lig_" + std::to_string(i)));

    auto results = ts.run(ligs);
    EXPECT_EQ(results.size(), 5u);
    // Coarse rank should be assigned
    for (const auto& r : results) {
        EXPECT_GT(r.coarse_rank, 0);
        EXPECT_EQ(r.full_rank, 0); // no full dock callback set
    }
}

TEST(TwoStageScreener, PipelineWithMockDock) {
    auto atoms   = make_test_target();
    auto spheres = make_test_spheres();

    TwoStageScreener ts;
    TwoStageConfig cfg;
    cfg.coarse.use_clash = false;
    cfg.coarse.rotations_per_axis = 4;
    cfg.top_n = 3;
    cfg.write_coarse_csv = false;
    ts.set_config(cfg);
    ts.prepare_target(atoms, spheres);

    // Mock dock callback: return negative of coarse score (reverses ranking)
    ts.set_full_dock_callback([](const ScreenLigand&, const ScreenResult& cr) {
        return -cr.score;
    });

    std::vector<ScreenLigand> ligs;
    for (int i = 0; i < 5; ++i)
        ligs.push_back(make_test_ligand("lig_" + std::to_string(i)));

    auto results = ts.run(ligs);
    EXPECT_EQ(results.size(), 5u);

    // Top 3 should have full_rank assigned
    int docked_count = 0;
    for (const auto& r : results)
        if (r.full_rank > 0) docked_count++;
    EXPECT_EQ(docked_count, 3);
}

// ═══════════════════════════════════════════════════════════════════════
// 10. File I/O round-trip (write then parse MOL2)
// ═══════════════════════════════════════════════════════════════════════

TEST(FileIO, WriteParseMol2) {
    // Write a minimal MOL2 file
    const char* mol2_content = R"(
@<TRIPOS>MOLECULE
caffeine
   14 15  0  0  0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
      1 C1       2.000  0.000  0.000 C.ar      1 LIG    0.0000
      2 N1       3.400  0.000  0.000 N.ar      1 LIG    0.0000
      3 O1       1.000  1.700  0.000 O.2       1 LIG    0.0000
      4 H1       0.500  0.500  0.000 H         1 LIG    0.0000
@<TRIPOS>BOND
      1     1     2 ar
      2     1     3 2
)";

    std::string tmpfile = "/tmp/test_coarse_screen_caffeine.mol2";
    {
        std::ofstream fout(tmpfile);
        fout << mol2_content;
    }

    auto ligs = CoarseScreener::load_ligands_mol2(tmpfile);
    ASSERT_EQ(ligs.size(), 1u);
    EXPECT_EQ(ligs[0].name, "caffeine");
    // 3 heavy atoms (H1 is skipped)
    EXPECT_EQ(ligs[0].atoms.size(), 3u);

    // Check types
    float dummy_r;
    EXPECT_EQ(ligs[0].atoms[0].type, sybyl_type_lookup("C.AR", dummy_r));
    EXPECT_EQ(ligs[0].atoms[1].type, sybyl_type_lookup("N.AR", dummy_r));
    EXPECT_EQ(ligs[0].atoms[2].type, sybyl_type_lookup("O.2",  dummy_r));

    std::remove(tmpfile.c_str());
}

TEST(FileIO, ParseTargetMol2) {
    const char* mol2_content = R"(
@<TRIPOS>MOLECULE
protein_chain
   4 0 0 0 0
PROTEIN
NO_CHARGES

@<TRIPOS>ATOM
      1 CA       0.000  0.000  0.000 C.3       1 ALA    0.0000
      2 N        1.500  0.000  0.000 N.am      1 ALA    0.0000
      3 O        0.000  1.200  0.000 O.2       1 ALA    0.0000
      4 H        0.500  0.600  0.000 H         1 ALA    0.0000
@<TRIPOS>BOND
)";

    std::string tmpfile = "/tmp/test_target.mol2";
    {
        std::ofstream fout(tmpfile);
        fout << mol2_content;
    }

    auto atoms = parse_target_mol2(tmpfile);
    ASSERT_EQ(atoms.size(), 3u); // H skipped
    EXPECT_EQ(atoms[0].type, 3);  // C.3
    EXPECT_EQ(atoms[1].type, 11); // N.AM
    EXPECT_EQ(atoms[2].type, 13); // O.2

    std::remove(tmpfile.c_str());
}

TEST(FileIO, ParseBindingSitePDB) {
    const char* pdb_content =
        "ATOM      1  O   SPH X   1      10.000  20.000  30.000  1.00  5.00\n"
        "ATOM      2  O   SPH X   2      15.000  25.000  35.000  1.00  3.50\n";

    std::string tmpfile = "/tmp/test_cleft.pdb";
    {
        std::ofstream fout(tmpfile);
        fout << pdb_content;
    }

    auto spheres = parse_binding_site_pdb(tmpfile);
    ASSERT_EQ(spheres.size(), 2u);
    EXPECT_NEAR(spheres[0].center.x, 10.0f, 0.01f);
    EXPECT_NEAR(spheres[0].center.y, 20.0f, 0.01f);
    EXPECT_NEAR(spheres[0].center.z, 30.0f, 0.01f);
    EXPECT_NEAR(spheres[0].radius,    5.0f, 0.01f);

    std::remove(tmpfile.c_str());
}

// ═══════════════════════════════════════════════════════════════════════
// Additional edge cases
// ═══════════════════════════════════════════════════════════════════════

TEST(CoarseScreener, UnpreparedReturnsSentinel) {
    CoarseScreener cs;
    auto lig = make_test_ligand();
    auto result = cs.screen_one(lig);
    EXPECT_EQ(result.score, cs.config().default_cf);
}

TEST(CoarseScreener, WithClash) {
    auto atoms   = make_test_target();
    auto spheres = make_test_spheres();

    CoarseScreener cs;
    CoarseScreenConfig cfg;
    cfg.use_clash = true;
    cfg.rotations_per_axis = 4;
    cs.set_config(cfg);
    cs.prepare_target(atoms, spheres);

    EXPECT_TRUE(cs.is_prepared());

    auto lig = make_test_ligand();
    auto result = cs.screen_one(lig);
    // Should get a score (might be sentinel if all poses clash, which is ok)
    EXPECT_LE(result.score, cfg.default_cf);
}

TEST(CoarseScreener, DeterministicScoring) {
    auto atoms   = make_test_target();
    auto spheres = make_test_spheres();

    CoarseScreener cs;
    CoarseScreenConfig cfg;
    cfg.use_clash = false;
    cfg.rotations_per_axis = 4;
    cs.set_config(cfg);
    cs.prepare_target(atoms, spheres);

    auto lig = make_test_ligand();
    auto r1 = cs.screen_one(lig);
    auto r2 = cs.screen_one(lig);

    // Same input should give same output
    EXPECT_FLOAT_EQ(r1.score, r2.score);
    EXPECT_EQ(r1.best_rotation, r2.best_rotation);
    EXPECT_EQ(r1.best_anchor, r2.best_anchor);
}

TEST(TwoStageScreener, CSVOutput) {
    auto atoms   = make_test_target();
    auto spheres = make_test_spheres();

    TwoStageScreener ts;
    TwoStageConfig cfg;
    cfg.coarse.use_clash = false;
    cfg.coarse.rotations_per_axis = 4;
    cfg.top_n = 2;
    cfg.output_dir = "/tmp/test_screen_output";
    cfg.write_coarse_csv = true;
    ts.set_config(cfg);
    ts.prepare_target(atoms, spheres);

    std::vector<ScreenLigand> ligs;
    ligs.push_back(make_test_ligand("mol_A"));
    ligs.push_back(make_test_ligand("mol_B"));

    auto results = ts.run(ligs);
    EXPECT_EQ(results.size(), 2u);

    // Check CSV was written
    std::string csv_path = cfg.output_dir + "/coarse_screen.csv";
    std::ifstream fin(csv_path);
    EXPECT_TRUE(fin.is_open());
    std::string header;
    std::getline(fin, header);
    EXPECT_FALSE(header.empty());
    EXPECT_NE(header.find("Rank"), std::string::npos);

    // Clean up
    std::filesystem::remove_all(cfg.output_dir);
}

TEST(TwoStageScreener, UnifiedCsvAndReceipt) {
    auto atoms   = make_test_target();
    auto spheres = make_test_spheres();

    TwoStageScreener ts;
    TwoStageConfig cfg;
    cfg.coarse.use_clash = false;
    cfg.coarse.rotations_per_axis = 4;
    cfg.top_n = 2;
    cfg.output_dir = "/tmp/test_screen_unified";
    cfg.write_coarse_csv = false;
    cfg.stage2_kind_label = "surrogate";
    ts.set_config(cfg);
    ts.prepare_target(atoms, spheres);
    ts.set_full_dock_callback([](const ScreenLigand&, const ScreenResult& cr) {
        return cr.score - 1.0f;
    });

    std::vector<ScreenLigand> ligs;
    ligs.push_back(make_test_ligand("act_A"));
    ligs.push_back(make_test_ligand("act_B"));
    ligs.push_back(make_test_ligand("dec_C"));

    auto results = ts.run(ligs);
    ASSERT_EQ(results.size(), 3u);
    TwoStageScreener::write_unified_csv(cfg.output_dir + "/unified.csv", results);
    ASSERT_TRUE(TwoStageScreener::write_screen_receipt(
        cfg.output_dir, static_cast<int>(ligs.size()), cfg.top_n, true));

    std::ifstream csv(cfg.output_dir + "/unified.csv");
    ASSERT_TRUE(csv.is_open());
    std::string header;
    std::getline(csv, header);
    EXPECT_NE(header.find("coarse_rank"), std::string::npos);
    EXPECT_NE(header.find("stage2_score"), std::string::npos);
    EXPECT_NE(header.find("stage2_kind"), std::string::npos);

    std::ifstream rec(cfg.output_dir + "/RUN_RECEIPT.json");
    ASSERT_TRUE(rec.is_open());
    std::string body((std::istreambuf_iterator<char>(rec)),
                     std::istreambuf_iterator<char>());
    EXPECT_NE(body.find("10.1101/2025.02.17.638675"), std::string::npos);
    EXPECT_NE(body.find("10.5281/zenodo.16861024"), std::string::npos);
    EXPECT_NE(body.find("stage2_callback"), std::string::npos);

    auto names = TwoStageScreener::top_names(results, 2);
    EXPECT_EQ(names.size(), 2u);
    std::filesystem::remove_all(cfg.output_dir);
}

TEST(TwoStageScreener, Stage2ScoreLessPutsFiniteBeforeNan) {
    TwoStageResult a, b, c;
    a.full_dock_score = -2.f; a.coarse_rank = 3;
    b.full_dock_score = -1.f; b.coarse_rank = 1;
    c.full_dock_score = std::numeric_limits<float>::quiet_NaN(); c.coarse_rank = 2;
    EXPECT_TRUE(stage2_score_less(a, b));
    EXPECT_FALSE(stage2_score_less(b, a));
    EXPECT_TRUE(stage2_score_less(a, c));
    EXPECT_TRUE(stage2_score_less(b, c));
    EXPECT_FALSE(stage2_score_less(c, a));
    EXPECT_FALSE(stage2_score_less(c, c));
}

TEST(TwoStageScreener, PpropMixedKeepDropSortsFiniteFirst) {
    auto atoms   = make_test_target();
    auto spheres = make_test_spheres();

    TwoStageScreener ts;
    TwoStageConfig cfg;
    cfg.coarse.use_clash = false;
    cfg.coarse.rotations_per_axis = 4;
    cfg.top_n = 3;
    cfg.write_coarse_csv = false;
    ts.set_config(cfg);
    ts.prepare_target(atoms, spheres);
    ts.set_full_dock_callback([](const ScreenLigand&, const ScreenResult& cr) {
        return cr.score;
    });

#ifdef _WIN32
    _putenv_s("FLEXAIDDS_PPROP_MAX", "0.4");
#else
    setenv("FLEXAIDDS_PPROP_MAX", "0.4", 1);
#endif
    std::vector<ScreenLigand> ligs;
    for (int i = 0; i < 5; ++i)
        ligs.push_back(make_test_ligand("lig_" + std::to_string(i)));
    auto results = ts.run(ligs);
#ifdef _WIN32
    _putenv_s("FLEXAIDDS_PPROP_MAX", "");
#else
    unsetenv("FLEXAIDDS_PPROP_MAX");
#endif

    ASSERT_EQ(results.size(), 5u);
    int finite = 0, nan = 0;
    for (int i = 0; i < 3; ++i) {
        if (std::isfinite(results[static_cast<size_t>(i)].full_dock_score)) {
            ++finite;
            EXPECT_GT(results[static_cast<size_t>(i)].full_rank, 0);
        } else {
            ++nan;
            EXPECT_EQ(results[static_cast<size_t>(i)].full_rank, 0);
        }
    }
    EXPECT_GE(finite, 1);
    EXPECT_GE(nan, 1);
    // Total order: no finite after a NaN in the Stage-2 window.
    bool saw_nan = false;
    for (int i = 0; i < 3; ++i) {
        const bool is_nan = !std::isfinite(results[static_cast<size_t>(i)].full_dock_score);
        if (saw_nan) EXPECT_TRUE(is_nan);
        if (is_nan) saw_nan = true;
    }
}

TEST(TwoStageScreener, MiniEnrichmentRecoversActives) {
    auto atoms   = make_test_target();
    auto spheres = make_test_spheres();

    CoarseScreener cs;
    CoarseScreenConfig cfg;
    cfg.use_clash = false;
    cfg.rotations_per_axis = 4;
    cfg.top_n = 2;
    cs.set_config(cfg);
    cs.prepare_target(atoms, spheres);

    std::vector<ScreenLigand> ligs;
    ligs.push_back(make_test_ligand("ACT_1"));
    ligs.push_back(make_test_ligand("ACT_2"));
    for (int i = 0; i < 8; ++i) {
        ScreenLigand dec;
        dec.name = "DEC_" + std::to_string(i);
        // Far dummy atom: should be out of the pocket / sentinel-ish
        dec.atoms.push_back({{400.f + float(i), 400.f, 400.f}, 39});
        ligs.push_back(dec);
    }

    auto results = cs.screen(ligs);
    ASSERT_GE(results.size(), 2u);
    std::string n0 = results[0].name;
    std::string n1 = results[1].name;
    EXPECT_TRUE(n0.rfind("ACT_", 0) == 0) << "top1=" << n0;
    EXPECT_TRUE(n1.rfind("ACT_", 0) == 0) << "top2=" << n1;
}

TEST(TwoStageScreen, PrefilterKeepFailsClosedOnNameMismatch) {
    std::vector<std::string> keep{"SCREEN_A", "SCREEN_B"};
    std::vector<std::string> lib{"file_stem_1", "file_stem_2"};
    std::vector<std::string> resolved;
    EXPECT_FALSE(resolve_prefilter_keep(keep, lib, resolved));
    EXPECT_TRUE(resolved.empty());

    std::vector<std::string> mixed{"SCREEN_A", "file_stem_1"};
    EXPECT_FALSE(resolve_prefilter_keep(keep, mixed, resolved));
    EXPECT_TRUE(resolved.empty());

    std::vector<std::string> ok_lib{"SCREEN_B", "SCREEN_A", "OTHER"};
    ASSERT_TRUE(resolve_prefilter_keep(keep, ok_lib, resolved));
    ASSERT_EQ(resolved.size(), 2u);
    EXPECT_EQ(resolved[0], "SCREEN_A");
    EXPECT_EQ(resolved[1], "SCREEN_B");

    EXPECT_FALSE(resolve_prefilter_keep({}, lib, resolved));
}

TEST(NRGRankMatrix, ProvenancePinsPublishedName) {
    // Self-consistency pin. Do not fetch GPL upstream bytes to "verify".
    double sum = 0.0;
    for (int i = 0; i < MATRIX_DIM; ++i)
        for (int j = 0; j < MATRIX_DIM; ++j)
            sum += kEnergyMatrix[i][j];
    EXPECT_NEAR(kEnergyMatrix[1][2], -170.52, 0.01);
    EXPECT_LT(sum, 0.0);  // net attractive table
}
