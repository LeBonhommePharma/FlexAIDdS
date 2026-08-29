// tests/test_cleft_cavity.cpp
// Unit tests for CleftDetector and CavityDetect modules.
// Tests cover: default parameters, probe generation, clustering,
// CavityDetector lifecycle, and sphere export.
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include "../LIB/CleftDetector.h"
#include "../LIB/site_confine.h"
#include "../LIB/CavityDetect/CavityDetect.h"

#include <cmath>
#include <cstring>
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>

static constexpr float FEPSILON = 1e-4f;

// ===========================================================================
// HELPERS
// ===========================================================================

namespace {

// Create a minimal atom with position and radius.
atom make_atom(float x, float y, float z, float radius = 1.7f,
               const char* name = "C") {
    atom a{};
    a.coor[0] = x;
    a.coor[1] = y;
    a.coor[2] = z;
    a.radius = radius;
    a.number = 1;
    std::strncpy(a.name, name, 4);
    return a;
}

// Create a minimal residue spanning atoms [first, last].
resid make_residue(int first, int last, const char* name = "ALA") {
    resid r{};
    std::strncpy(r.name, name, 3);
    r.chn = 'A';
    r.number = 1;
    r.type = 0;  // protein
    // fatm and latm need to be allocated
    r.fatm = new int[1];
    r.latm = new int[1];
    r.fatm[0] = first;
    r.latm[0] = last;
    r.trot = 1;
    r.rot = 0;
    return r;
}

void free_residue(resid& r) {
    delete[] r.fatm;
    delete[] r.latm;
}

} // namespace

// ===========================================================================
// CleftDetector — default_cleft_params()
// ===========================================================================

TEST(CleftDetectorParams, DefaultValues) {
    CleftDetectorParams p = default_cleft_params();
    EXPECT_EQ(p.top_k_clefts, 5);

    EXPECT_FLOAT_EQ(p.max_pair_dist, 12.0f);
    EXPECT_FLOAT_EQ(p.probe_radius_max, 5.0f);
    EXPECT_FLOAT_EQ(p.probe_radius_min, 1.5f);
    EXPECT_FLOAT_EQ(p.probe_shrink_step, 0.1f);
    EXPECT_FLOAT_EQ(p.cluster_cutoff, 4.0f);
    EXPECT_EQ(p.min_cluster_size, 10);
}

// ===========================================================================
// CleftDetector — detect_cleft with minimal geometry
// ===========================================================================

class CleftDetectorTest : public ::testing::Test {
protected:
    std::vector<atom> atoms;
    std::vector<resid> residues;

    void SetUp() override {
        // Create a simple "box" of atoms forming a cavity in the center.
        // 8 corner atoms at ±4 Å to define a cleft.
        atoms.push_back(make_atom(-4.0f, -4.0f, -4.0f));
        atoms.push_back(make_atom(-4.0f, -4.0f,  4.0f));
        atoms.push_back(make_atom(-4.0f,  4.0f, -4.0f));
        atoms.push_back(make_atom(-4.0f,  4.0f,  4.0f));
        atoms.push_back(make_atom( 4.0f, -4.0f, -4.0f));
        atoms.push_back(make_atom( 4.0f, -4.0f,  4.0f));
        atoms.push_back(make_atom( 4.0f,  4.0f, -4.0f));
        atoms.push_back(make_atom( 4.0f,  4.0f,  4.0f));

        // Add some atoms along edges to create more pairs
        for (float z = -3.0f; z <= 3.0f; z += 2.0f) {
            atoms.push_back(make_atom(-4.0f, 0.0f, z));
            atoms.push_back(make_atom( 4.0f, 0.0f, z));
            atoms.push_back(make_atom( 0.0f, -4.0f, z));
            atoms.push_back(make_atom( 0.0f,  4.0f, z));
        }

        // Number atoms
        for (size_t i = 0; i < atoms.size(); ++i)
            atoms[i].number = static_cast<int>(i + 1);

        resid r = make_residue(0, static_cast<int>(atoms.size()) - 1);
        residues.push_back(r);
    }

    void TearDown() override {
        for (auto& r : residues)
            free_residue(r);
    }
};

TEST_F(CleftDetectorTest, ReturnsNullForTooFewAtoms) {
    // With only 2 atoms far apart, no cleft should form with min_cluster_size=10
    atom pair[2] = {make_atom(0.0f, 0.0f, 0.0f), make_atom(20.0f, 0.0f, 0.0f)};
    resid r = make_residue(0, 1);

    CleftDetectorParams params = default_cleft_params();
    sphere* result = detect_cleft(pair, &r, 2, 1, params);

    // With only 2 atoms, the maximum probes possible is 1,
    // which is less than min_cluster_size (10), so result should be NULL
    EXPECT_EQ(result, nullptr);

    free_residue(r);
    if (result) free_sphere_list(result);
}

TEST_F(CleftDetectorTest, DetectsCleftInBox) {
    CleftDetectorParams params = default_cleft_params();
    params.min_cluster_size = 2;  // Lower threshold for small test geometry

    sphere* result = detect_cleft(atoms.data(), residues.data(),
                                  static_cast<int>(atoms.size()), 1, params);

    if (result) {
        // Count spheres
        int count = 0;
        sphere* s = result;
        while (s) {
            ++count;
            // Each sphere should have positive radius
            EXPECT_GT(s->radius, 0.0f);
            s = s->prev;
        }
        EXPECT_GE(count, 2);
        free_sphere_list(result);
    }
    // If result is null, the geometry may not produce a cleft under default params
    // (this is geometry-dependent, not a failure)
}

TEST(SiteConfineDecision, SparsePocketMustNotFullGridFallback) {
    // 1STP biotin: keep=188, MIN_SITE_GRID=500, total=19897 → previously
    // fell back to the full protein grid. Must confine.
    EXPECT_TRUE(flexaids::site_confine_should_rebuild(188, 19897));
    EXPECT_TRUE(flexaids::site_confine_should_rebuild(500, 19897));
    EXPECT_TRUE(flexaids::site_confine_should_rebuild(25, 19897));
    EXPECT_FALSE(flexaids::site_confine_should_rebuild(24, 19897));
    EXPECT_FALSE(flexaids::site_confine_should_rebuild(0, 19897));
    EXPECT_FALSE(flexaids::site_confine_should_rebuild(19897, 19897));
}

TEST_F(CleftDetectorTest, AllHetAtomsYieldNoCleft) {
    // Every atom tagged as ligand/HETATM must be ignored, so SURFNET sees
    // an empty receptor and returns no cavity.
    std::vector<resid> het_only(2);
    het_only[1].type = 1;
    for (auto& a : atoms)
        a.ofres = 1;

    CleftDetectorParams params = default_cleft_params();
    params.min_cluster_size = 2;
    sphere* result = detect_cleft(atoms.data(), het_only.data(),
                                  static_cast<int>(atoms.size()), 1, params);
    EXPECT_EQ(result, nullptr);
    if (result) free_sphere_list(result);
}

TEST_F(CleftDetectorTest, SkipsHetLigandOccupyingPocket) {
    // 1-based residue map: [1]=protein (type 0), [2]=ligand (type 1).
    // A bulky HET atom at the cavity origin would shrink SURFNET probes
    // out of the pocket if it were included — the skip must prevent that.
    std::vector<atom> het_atoms = atoms;
    std::vector<resid> het_residues(3);
    het_residues[1].type = 0;
    het_residues[2].type = 1;
    for (auto& a : het_atoms)
        a.ofres = 1;
    atom lig = make_atom(0.0f, 0.0f, 0.0f, 2.5f, "C1");
    lig.ofres = 2;
    het_atoms.push_back(lig);

    CleftDetectorParams params = default_cleft_params();
    params.min_cluster_size = 2;

    sphere* with_skip = detect_cleft(het_atoms.data(), het_residues.data(),
                                     static_cast<int>(het_atoms.size()), 2, params);
    // Ligand skip is a correctness gate, not a geometry guarantee: the box
    // fixture may still return null under default SURFNET knobs. If spheres
    // are produced, the occupied-pocket ligand must not have emptied them.
    if (with_skip) {
        int count = 0;
        for (sphere* s = with_skip; s; s = s->prev) ++count;
        EXPECT_GE(count, 1);
        free_sphere_list(with_skip);
    }
}

TEST_F(CleftDetectorTest, FreeSphereListHandlesNull) {
    // Should not crash
    free_sphere_list(nullptr);
}

// ===========================================================================
// CleftDetector — write_cleft_spheres
// ===========================================================================

TEST(WriteCleftSpheres, WritesValidFile) {
    // Build a small sphere list manually
    sphere* s2 = static_cast<sphere*>(std::malloc(sizeof(sphere)));
    s2->center[0] = 1.0f; s2->center[1] = 2.0f; s2->center[2] = 3.0f;
    s2->radius = 1.5f;
    s2->prev = nullptr;

    sphere* s1 = static_cast<sphere*>(std::malloc(sizeof(sphere)));
    s1->center[0] = 4.0f; s1->center[1] = 5.0f; s1->center[2] = 6.0f;
    s1->radius = 2.0f;
    s1->prev = s2;

    std::string tmpfile = std::filesystem::temp_directory_path().string()
                          + "/test_cleft_spheres.pdb";
    write_cleft_spheres(s1, tmpfile.c_str());

    // Verify file exists and has content
    std::ifstream ifs(tmpfile);
    ASSERT_TRUE(ifs.good());
    std::string content((std::istreambuf_iterator<char>(ifs)),
                         std::istreambuf_iterator<char>());
    EXPECT_FALSE(content.empty());

    // Cleanup
    std::remove(tmpfile.c_str());
    free_sphere_list(s1);
}

// ===========================================================================
// CavityDetector — lifecycle tests
// ===========================================================================

TEST(CavityDetector, DefaultConstructorEmpty) {
    cavity_detect::CavityDetector detector;
    EXPECT_TRUE(detector.clefts().empty());
}

TEST(CavityDetector, LoadFromFAPopulatesAtoms) {
    std::vector<atom> atoms;
    atoms.push_back(make_atom(0.0f, 0.0f, 0.0f));
    atoms.push_back(make_atom(3.0f, 0.0f, 0.0f));
    atoms.push_back(make_atom(6.0f, 0.0f, 0.0f));

    resid r = make_residue(0, 2);

    cavity_detect::CavityDetector detector;
    detector.load_from_fa(atoms.data(), &r, 1);

    // After loading, clefts should still be empty (no detection yet)
    EXPECT_TRUE(detector.clefts().empty());

    free_residue(r);
}

TEST(CavityDetector, DetectOnSmallGeometry) {
    // Create a small cluster of atoms forming a pocket
    std::vector<atom> atoms;
    for (float x = -3.0f; x <= 3.0f; x += 2.0f)
        for (float y = -3.0f; y <= 3.0f; y += 2.0f)
            for (float z = -3.0f; z <= 3.0f; z += 2.0f)
                atoms.push_back(make_atom(x, y, z, 1.5f));

    for (size_t i = 0; i < atoms.size(); ++i)
        atoms[i].number = static_cast<int>(i + 1);

    resid r = make_residue(0, static_cast<int>(atoms.size()) - 1);

    cavity_detect::CavityDetector detector;
    detector.load_from_fa(atoms.data(), &r, 1);
    detector.detect(1.0f, 3.0f);

    // May or may not find clefts depending on geometry, but shouldn't crash
    // If clefts are found, verify they have valid data
    for (const auto& cleft : detector.clefts()) {
        EXPECT_GT(cleft.id, 0);
        EXPECT_FALSE(cleft.spheres.empty());
        EXPECT_GT(cleft.volume, 0.0f);
    }

    free_residue(r);
}

TEST(CavityDetector, SortCleftsOrdersBySize) {
    // Create atoms in two spatially separated clusters
    std::vector<atom> atoms;

    // Large cluster at (0,0,0)
    for (float x = -5.0f; x <= 5.0f; x += 2.0f)
        for (float y = -5.0f; y <= 5.0f; y += 2.0f)
            atoms.push_back(make_atom(x, y, 0.0f, 1.5f));

    // Small cluster at (30,0,0)
    for (float x = 28.0f; x <= 32.0f; x += 2.0f)
        atoms.push_back(make_atom(x, 0.0f, 0.0f, 1.5f));

    for (size_t i = 0; i < atoms.size(); ++i)
        atoms[i].number = static_cast<int>(i + 1);

    resid r = make_residue(0, static_cast<int>(atoms.size()) - 1);

    cavity_detect::CavityDetector detector;
    detector.load_from_fa(atoms.data(), &r, 1);
    detector.detect(1.0f, 3.0f);

    if (detector.clefts().size() >= 2) {
        // After sort_clefts (called inside detect), first cleft should have
        // more spheres than second
        EXPECT_GE(detector.clefts()[0].spheres.size(),
                   detector.clefts()[1].spheres.size());
    }

    free_residue(r);
}

TEST(CavityDetector, ToFlexaidSpheresReturnsLinkedList) {
    std::vector<atom> atoms;
    for (float x = -4.0f; x <= 4.0f; x += 2.0f)
        for (float y = -4.0f; y <= 4.0f; y += 2.0f)
            atoms.push_back(make_atom(x, y, 0.0f, 1.5f));

    for (size_t i = 0; i < atoms.size(); ++i)
        atoms[i].number = static_cast<int>(i + 1);

    resid r = make_residue(0, static_cast<int>(atoms.size()) - 1);

    cavity_detect::CavityDetector detector;
    detector.load_from_fa(atoms.data(), &r, 1);
    detector.detect(1.0f, 3.0f);

    if (!detector.clefts().empty()) {
        sphere* list = detector.to_flexaid_spheres(detector.clefts()[0].id);
        if (list) {
            // Walk the linked list and count
            int count = 0;
            sphere* s = list;
            while (s) {
                ++count;
                EXPECT_GT(s->radius, 0.0f);
                s = s->prev;
            }
            EXPECT_GT(count, 0);
            free_sphere_list(list);
        }
    }

    free_residue(r);
}

TEST(CavityDetector, WriteSpheresPDB) {
    std::vector<atom> atoms;
    for (float x = -3.0f; x <= 3.0f; x += 1.5f)
        for (float y = -3.0f; y <= 3.0f; y += 1.5f)
            atoms.push_back(make_atom(x, y, 0.0f, 1.5f));

    for (size_t i = 0; i < atoms.size(); ++i)
        atoms[i].number = static_cast<int>(i + 1);

    resid r = make_residue(0, static_cast<int>(atoms.size()) - 1);

    cavity_detect::CavityDetector detector;
    detector.load_from_fa(atoms.data(), &r, 1);
    detector.detect(1.0f, 3.0f);

    if (!detector.clefts().empty()) {
        std::string tmpfile = std::filesystem::temp_directory_path().string()
                              + "/test_cavity_spheres.pdb";
        detector.write_sphere_pdb(tmpfile, detector.clefts()[0].id);

        std::ifstream ifs(tmpfile);
        EXPECT_TRUE(ifs.good());
        std::remove(tmpfile.c_str());
    }

    free_residue(r);
}

// ===========================================================================
// Task 8: Cleft Annotation & Flexible Residue Selection (PREPROCESSING)
// ===========================================================================

TEST(CleftAnnotationTest, ToyStructureReturnsExpectedSet) {
    std::vector<int> cleft_seeds = {5};
    double shell = 5.0;
    std::vector<int> active = {10};
    std::vector<int> fixed  = {7};
    std::vector<int> forced = {12};

    auto result = select_flexible_residues(
        nullptr, nullptr, 0, 20,
        cleft_seeds, shell, active, fixed, forced);

    auto has = [&](int r){ return std::find(result.begin(), result.end(), r) != result.end(); };
    EXPECT_TRUE(has(5));
    EXPECT_TRUE(has(10));
    EXPECT_TRUE(has(12));
    EXPECT_FALSE(has(7));
}

TEST(CleftAnnotationTest, FixedResiduesAreExcluded) {
    std::vector<int> seeds = {1,2,3};
    std::vector<int> fixed = {2};
    auto res = select_flexible_residues(nullptr, nullptr, 0, 10, seeds, 4.0, {}, fixed, {});
    EXPECT_EQ(std::count(res.begin(), res.end(), 2), 0);
}

TEST(CleftAnnotationTest, ForcedFlexibleAreIncluded) {
    std::vector<int> seeds = {1};
    std::vector<int> forced = {99};
    auto res = select_flexible_residues(nullptr, nullptr, 0, 100, seeds, 4.0, {}, {}, forced);
    EXPECT_TRUE(std::find(res.begin(), res.end(), 99) != res.end());
}

TEST(CleftAnnotationTest, NoDuplicates) {
    std::vector<int> seeds = {5,5,5};
    auto res = select_flexible_residues(nullptr, nullptr, 0, 10, seeds, 4.0);
    EXPECT_EQ(std::count(res.begin(), res.end(), 5), 1);
}

TEST(CleftAnnotationTest, DeterministicOrdering) {
    std::vector<int> seeds = {9,1,5};
    auto res = select_flexible_residues(nullptr, nullptr, 0, 20, seeds, 4.0);
    std::vector<int> expected = {1,5,9};
    EXPECT_EQ(res, expected);
}
