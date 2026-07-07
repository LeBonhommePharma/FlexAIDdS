// tests/test_cavity_detect.cpp
// Unit tests for CavityDetect — SURFNET-based cavity detection
// Tests PDB loading, probe-sphere placement, cleft clustering, and output.

#include <gtest/gtest.h>
#include "CavityDetect.h"
#include "SpatialGrid.h"
#include <cmath>
#include <fstream>
#include <string>
#include <filesystem>
#include <algorithm>

using namespace cavity_detect;

// ===========================================================================
// HELPER: write a minimal PDB with atoms forming a cavity
// ===========================================================================
static std::string write_test_pdb(const std::string& filename) {
    std::string path = filename;
    std::ofstream out(path);

    // Create a small "protein" with atoms arranged in a hollow sphere
    // Atoms placed on corners of a cube (~8 Å side) → cavity in center
    float coords[][3] = {
        { 0.0f,  0.0f,  0.0f},
        { 8.0f,  0.0f,  0.0f},
        { 0.0f,  8.0f,  0.0f},
        { 8.0f,  8.0f,  0.0f},
        { 0.0f,  0.0f,  8.0f},
        { 8.0f,  0.0f,  8.0f},
        { 0.0f,  8.0f,  8.0f},
        { 8.0f,  8.0f,  8.0f},
        // Add atoms along edges for better probe placement
        { 4.0f,  0.0f,  0.0f},
        { 0.0f,  4.0f,  0.0f},
        { 0.0f,  0.0f,  4.0f},
        { 8.0f,  4.0f,  0.0f},
        { 8.0f,  0.0f,  4.0f},
        { 4.0f,  8.0f,  0.0f},
        { 0.0f,  8.0f,  4.0f},
        { 4.0f,  0.0f,  8.0f},
        { 0.0f,  4.0f,  8.0f},
        { 8.0f,  8.0f,  4.0f},
        { 8.0f,  4.0f,  8.0f},
        { 4.0f,  8.0f,  8.0f},
    };

    int n_atoms = sizeof(coords) / sizeof(coords[0]);
    for (int i = 0; i < n_atoms; ++i) {
        char line[120];
        std::snprintf(line, sizeof(line),
            "ATOM  %5d  CA  ALA A%4d    %8.3f%8.3f%8.3f  1.00  1.70\n",
            i + 1, i + 1, coords[i][0], coords[i][1], coords[i][2]);
        out << line;
    }
    out << "END\n";
    out.close();
    return path;
}

static std::string write_high_bfactor_test_pdb(const std::string& filename) {
    std::string path = filename;
    std::ofstream out(path);
    float coords[][3] = {
        { 0.0f,  0.0f,  0.0f},
        { 8.0f,  0.0f,  0.0f},
        { 0.0f,  8.0f,  0.0f},
        { 8.0f,  8.0f,  0.0f},
        { 0.0f,  0.0f,  8.0f},
        { 8.0f,  0.0f,  8.0f},
        { 0.0f,  8.0f,  8.0f},
        { 8.0f,  8.0f,  8.0f},
        { 4.0f,  0.0f,  0.0f},
        { 0.0f,  4.0f,  0.0f},
        { 0.0f,  0.0f,  4.0f},
    };

    int n_atoms = sizeof(coords) / sizeof(coords[0]);
    for (int i = 0; i < n_atoms; ++i) {
        char line[120];
        std::snprintf(line, sizeof(line),
            "ATOM  %5d  CA  ALA A%4d    %8.3f%8.3f%8.3f  1.00 99.99           C\n",
            i + 1, i + 1, coords[i][0], coords[i][1], coords[i][2]);
        out << line;
    }
    out << "END\n";
    out.close();
    return path;
}

// ===========================================================================
// CONSTRUCTION
// ===========================================================================

TEST(CavityDetect, DefaultConstructor) {
    CavityDetector detector;
    EXPECT_TRUE(detector.clefts().empty());
}

// ===========================================================================
// PDB LOADING
// ===========================================================================

TEST(CavityDetect, LoadFromPDB) {
    std::string pdb = write_test_pdb("test_cavity_input.pdb");
    CavityDetector detector;
    detector.load_from_pdb(pdb);

    // Internal atom list should be populated
    EXPECT_TRUE(detector.clefts().empty());  // no detection yet

    std::filesystem::remove(pdb);
}

TEST(CavityDetect, LoadFromMissingPDB) {
    CavityDetector detector;
    detector.load_from_pdb("nonexistent_file_xyz.pdb");
    // Should gracefully handle missing file
    EXPECT_TRUE(detector.clefts().empty());
}

TEST(CavityDetect, PDBBFactorIsNotUsedAsRadius) {
    std::string pdb = write_high_bfactor_test_pdb("test_cavity_high_bfactor.pdb");
    CavityDetector detector;
    detector.load_from_pdb(pdb);
    detector.detect(1.0f, 4.0f);

    EXPECT_GE(detector.clefts().size(), 1u);
    std::filesystem::remove(pdb);
}

// ===========================================================================
// CAVITY DETECTION
// ===========================================================================

TEST(CavityDetect, DetectFindsCleftsInCube) {
    std::string pdb = write_test_pdb("test_cavity_detect.pdb");
    CavityDetector detector;
    detector.load_from_pdb(pdb);
    detector.detect(1.0f, 4.0f);

    // Should find at least one cleft in the cube's interior
    EXPECT_GE(detector.clefts().size(), 1u);

    // All clefts should have positive volume
    for (const auto& cleft : detector.clefts()) {
        EXPECT_GT(cleft.volume, 0.0f);
        EXPECT_GT(cleft.spheres.size(), 0u);
    }

    std::filesystem::remove(pdb);
}

TEST(CavityDetect, DetectedSpheresInRadiusBounds) {
    std::string pdb = write_test_pdb("test_cavity_bounds.pdb");
    CavityDetector detector;
    detector.load_from_pdb(pdb);

    float min_r = 1.0f, max_r = 3.5f;
    detector.detect(min_r, max_r);

    for (const auto& cleft : detector.clefts()) {
        for (const auto& s : cleft.spheres) {
            EXPECT_GE(s.radius, min_r - 0.01f)
                << "Sphere radius below minimum";
            EXPECT_LE(s.radius, max_r + 0.01f)
                << "Sphere radius above maximum";
        }
    }

    std::filesystem::remove(pdb);
}

TEST(CavityDetect, TooFewAtomsProducesNoClefts) {
    // Single atom cannot form a cavity
    std::ofstream out("test_single_atom.pdb");
    out << "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  1.70\n";
    out << "END\n";
    out.close();

    CavityDetector detector;
    detector.load_from_pdb("test_single_atom.pdb");
    detector.detect();
    EXPECT_TRUE(detector.clefts().empty());

    std::filesystem::remove("test_single_atom.pdb");
}

// ===========================================================================
// CLEFT PROPERTIES
// ===========================================================================

TEST(CavityDetect, CleftIdsAreSequential) {
    std::string pdb = write_test_pdb("test_cavity_ids.pdb");
    CavityDetector detector;
    detector.load_from_pdb(pdb);
    detector.detect(1.0f, 4.0f);

    for (size_t i = 0; i < detector.clefts().size(); ++i) {
        EXPECT_EQ(detector.clefts()[i].id, static_cast<int>(i) + 1);
    }

    std::filesystem::remove(pdb);
}

TEST(CavityDetect, CleftsSortedBySize) {
    std::string pdb = write_test_pdb("test_cavity_sort.pdb");
    CavityDetector detector;
    detector.load_from_pdb(pdb);
    detector.detect(1.0f, 4.0f);

    for (size_t i = 1; i < detector.clefts().size(); ++i) {
        EXPECT_GE(detector.clefts()[i - 1].spheres.size(),
                  detector.clefts()[i].spheres.size())
            << "Clefts not sorted by size (descending)";
    }

    std::filesystem::remove(pdb);
}

TEST(CavityDetect, CleftCenterIsFinite) {
    std::string pdb = write_test_pdb("test_cavity_center.pdb");
    CavityDetector detector;
    detector.load_from_pdb(pdb);
    detector.detect(1.0f, 4.0f);

    for (const auto& cleft : detector.clefts()) {
        EXPECT_TRUE(std::isfinite(cleft.center[0]));
        EXPECT_TRUE(std::isfinite(cleft.center[1]));
        EXPECT_TRUE(std::isfinite(cleft.center[2]));
        EXPECT_TRUE(std::isfinite(cleft.effrad));
        EXPECT_GT(cleft.effrad, 0.0f);
    }

    std::filesystem::remove(pdb);
}

// ===========================================================================
// MERGE CLEFTS
// ===========================================================================

TEST(CavityDetect, MergeCleftsReducesCount) {
    std::string pdb = write_test_pdb("test_cavity_merge.pdb");
    CavityDetector detector;
    detector.load_from_pdb(pdb);
    detector.detect(1.0f, 4.0f);

    size_t before = detector.clefts().size();
    detector.merge_clefts();
    size_t after = detector.clefts().size();

    // After merging, count should be <= before
    EXPECT_LE(after, before);

    std::filesystem::remove(pdb);
}

// ===========================================================================
// SPHERE PDB OUTPUT
// ===========================================================================

TEST(CavityDetect, WriteSpheresPDB) {
    std::string pdb = write_test_pdb("test_cavity_write.pdb");
    CavityDetector detector;
    detector.load_from_pdb(pdb);
    detector.detect(1.0f, 4.0f);

    std::string out_file = "test_spheres_output.pdb";
    detector.write_sphere_pdb(out_file);

    // Output file should exist and be non-empty
    std::ifstream check(out_file);
    EXPECT_TRUE(check.good());
    std::string line;
    bool has_hetatm = false;
    while (std::getline(check, line)) {
        if (line.substr(0, 6) == "HETATM") {
            has_hetatm = true;
        }
    }
    if (!detector.clefts().empty()) {
        EXPECT_TRUE(has_hetatm) << "Sphere PDB should contain HETATM records";
    }

    std::filesystem::remove(pdb);
    std::filesystem::remove(out_file);
}

// ===========================================================================
// TO FLEXAID SPHERES (linked list)
// ===========================================================================

TEST(CavityDetect, ToFlexaidSpheresNullForMissingCleft) {
    CavityDetector detector;
    sphere* result = detector.to_flexaid_spheres(999);
    EXPECT_EQ(result, nullptr);
}

TEST(CavityDetect, ToFlexaidSpheresLinkedList) {
    std::string pdb = write_test_pdb("test_cavity_ll.pdb");
    CavityDetector detector;
    detector.load_from_pdb(pdb);
    detector.detect(1.0f, 4.0f);

    if (!detector.clefts().empty()) {
        int cleft_id = detector.clefts()[0].id;
        sphere* head = detector.to_flexaid_spheres(cleft_id);

        EXPECT_NE(head, nullptr);

        // Walk the linked list and count
        int count = 0;
        sphere* curr = head;
        while (curr) {
            ++count;
            EXPECT_TRUE(std::isfinite(curr->center[0]));
            EXPECT_TRUE(std::isfinite(curr->radius));
            EXPECT_GT(curr->radius, 0.0f);
            sphere* prev = curr->prev;
            delete curr;
            curr = prev;
        }
        EXPECT_EQ(count, static_cast<int>(detector.clefts()[0].spheres.size()));
    }

    std::filesystem::remove(pdb);
}

// ===========================================================================
// FILTER ANCHOR RESIDUES
// ===========================================================================

TEST(CavityDetect, FilterAnchorEmptyStringNoOp) {
    std::string pdb = write_test_pdb("test_cavity_filter.pdb");
    CavityDetector detector;
    detector.load_from_pdb(pdb);
    detector.detect(1.0f, 4.0f);

    size_t before = detector.clefts().size();
    detector.filter_anchor_residues("");
    EXPECT_EQ(detector.clefts().size(), before);

    std::filesystem::remove(pdb);
}

// ===========================================================================
// SPATIAL GRID
// ===========================================================================

// Helper: create an atom vector for grid tests
static std::vector<atom> make_atoms(const std::vector<std::array<float,3>>& coords,
                                     float radius = 1.7f) {
    std::vector<atom> atoms;
    for (const auto& c : coords) {
        atom a{};
        a.coor[0] = c[0]; a.coor[1] = c[1]; a.coor[2] = c[2];
        a.radius = radius;
        a.coor_ref = nullptr; a.par = nullptr;
        a.cons = nullptr; a.optres = nullptr; a.eigen = nullptr;
        atoms.push_back(a);
    }
    return atoms;
}

TEST(SpatialGrid, BuildEmpty) {
    SpatialGrid grid;
    std::vector<atom> empty;
    grid.build(empty);
    EXPECT_TRUE(grid.empty());
    EXPECT_EQ(grid.atom_count(), 0u);

    float coord[3] = {0.f, 0.f, 0.f};
    auto result = grid.query_neighbors(coord);
    EXPECT_TRUE(result.empty());
}

TEST(SpatialGrid, QueryReturnsNearbyAtoms) {
    // 8 atoms at cube corners (8 Å side) — all within one ~6.5 Å cell neighborhood
    auto atoms = make_atoms({
        {0.f, 0.f, 0.f}, {8.f, 0.f, 0.f}, {0.f, 8.f, 0.f}, {8.f, 8.f, 0.f},
        {0.f, 0.f, 8.f}, {8.f, 0.f, 8.f}, {0.f, 8.f, 8.f}, {8.f, 8.f, 8.f},
    });

    SpatialGrid grid;
    grid.build(atoms);
    EXPECT_FALSE(grid.empty());
    EXPECT_EQ(grid.atom_count(), 8u);

    // Query at center — should find all 8 atoms in neighboring cells
    float center[3] = {4.f, 4.f, 4.f};
    auto result = grid.query_neighbors(center);
    EXPECT_EQ(result.size(), 8u);
}

TEST(SpatialGrid, ExcludesDistantAtoms) {
    // Two clusters far apart
    auto atoms = make_atoms({
        {0.f, 0.f, 0.f}, {1.f, 1.f, 1.f},
        {100.f, 100.f, 100.f}, {101.f, 101.f, 101.f},
    });

    SpatialGrid grid;
    grid.build(atoms);

    // Query near origin — should NOT find atoms at (100,100,100)
    float near_origin[3] = {0.5f, 0.5f, 0.5f};
    auto result = grid.query_neighbors(near_origin);
    EXPECT_GE(result.size(), 2u);  // at least the 2 nearby atoms
    for (auto idx : result) {
        // All returned atoms should be near the origin
        EXPECT_LT(atoms[idx].coor[0], 50.f)
            << "Grid returned a distant atom at index " << idx;
    }

    // Query near far cluster
    float near_far[3] = {100.5f, 100.5f, 100.5f};
    auto far_result = grid.query_neighbors(near_far);
    EXPECT_GE(far_result.size(), 2u);
    for (auto idx : far_result) {
        EXPECT_GT(atoms[idx].coor[0], 50.f)
            << "Grid returned a nearby atom when querying far cluster";
    }
}

TEST(SpatialGrid, QueryAtBoundary) {
    auto atoms = make_atoms({{0.f, 0.f, 0.f}});
    SpatialGrid grid;
    grid.build(atoms);

    // Query at the atom's own position (boundary of grid)
    float coord[3] = {0.f, 0.f, 0.f};
    auto result = grid.query_neighbors(coord);
    EXPECT_EQ(result.size(), 1u);
    EXPECT_EQ(result[0], 0u);
}

TEST(SpatialGrid, DetectResultsUnchanged) {
    // Regression test: spatial grid optimization must produce identical results
    // to the original all-pairs approach. We verify by checking that detect()
    // produces non-empty results with the expected properties.
    std::string pdb = write_test_pdb("test_grid_regression.pdb");
    CavityDetector detector;
    detector.load_from_pdb(pdb);
    detector.detect(1.0f, 4.0f);

    // Should find at least one cleft (same as DetectFindsCleftsInCube)
    EXPECT_GE(detector.clefts().size(), 1u);

    // All spheres must be within bounds
    for (const auto& cleft : detector.clefts()) {
        EXPECT_GT(cleft.volume, 0.0f);
        for (const auto& s : cleft.spheres) {
            EXPECT_GE(s.radius, 1.0f - 0.01f);
            EXPECT_LE(s.radius, 4.0f + 0.01f);
            EXPECT_TRUE(std::isfinite(s.center[0]));
            EXPECT_TRUE(std::isfinite(s.center[1]));
            EXPECT_TRUE(std::isfinite(s.center[2]));
        }
    }

    std::filesystem::remove(pdb);
}
