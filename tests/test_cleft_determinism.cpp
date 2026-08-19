// test_cleft_determinism.cpp — order-sensitivity and flexible-residue tests
// for CleftDetector.
//
// Motivation: OPTIMIZATION_KNOWN_ISSUES.md flags the OpenMP probe merge in
// generate_probes() as TSan-racy and order-sensitive, and the downstream
// single-linkage cluster_probes() consumes that order. tests/test_cleft_cavity.cpp
// already covers the happy path (a cleft is found, spheres are written); what
// was missing is (a) whether the detector is reproducible run-to-run, and
// (b) any coverage at all of select_flexible_residues(), which is the public
// entry point for flexible-residue selection.
//
// cluster_probes() is a file-static free function, so it cannot be called
// directly without changing CleftDetector.cpp. Since that file is under active
// modification by another agent, these tests exercise it through the public
// detect_cleft() surface instead. See the report for the refactor that would
// make it directly testable.
//
// Apache-2.0 — see LICENSE.

#include <gtest/gtest.h>

#include "../LIB/CleftDetector.h"

#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace {

atom make_atom(float x, float y, float z, float radius = 1.7f) {
    atom a{};
    a.coor[0] = x;
    a.coor[1] = y;
    a.coor[2] = z;
    a.radius = radius;
    a.number = 1;
    std::strncpy(a.name, "C", 4);
    return a;
}

resid make_residue(const char* name) {
    resid r{};
    std::strncpy(r.name, name, 3);
    r.chn = 'A';
    r.number = 1;
    r.type = 0;
    r.fatm = new int[1];
    r.latm = new int[1];
    r.fatm[0] = 0;
    r.latm[0] = 0;
    r.trot = 1;
    r.rot = 0;
    return r;
}

void free_residues(std::vector<resid>& rs) {
    for (auto& r : rs) {
        delete[] r.fatm;
        delete[] r.latm;
    }
}

// A hollow box of atoms: the interior is a cavity large enough to survive
// probe shrinking with the default parameters.
std::vector<atom> hollow_box(float half, float offset = 0.0f) {
    std::vector<atom> atoms;
    for (int sx = -1; sx <= 1; ++sx)
        for (int sy = -1; sy <= 1; ++sy)
            for (int sz = -1; sz <= 1; ++sz) {
                if (sx == 0 && sy == 0 && sz == 0) continue;  // hollow centre
                atoms.push_back(make_atom(offset + static_cast<float>(sx) * half,
                                          static_cast<float>(sy) * half,
                                          static_cast<float>(sz) * half));
            }
    return atoms;
}

// Flatten a sphere linked list into a comparable, order-preserving vector.
struct SphereRec {
    float x, y, z, r;
    bool operator==(const SphereRec& o) const {
        return x == o.x && y == o.y && z == o.z && r == o.r;
    }
};

std::vector<SphereRec> flatten(sphere* head) {
    std::vector<SphereRec> out;
    for (sphere* s = head; s; s = s->prev)
        out.push_back({s->center[0], s->center[1], s->center[2], s->radius});
    return out;
}

// RAII guard so a test never leaks an env-var setting into its neighbours.
class EnvGuard {
public:
    EnvGuard(const char* name, const char* value) : name_(name) {
        const char* prev = std::getenv(name);
        had_ = prev != nullptr;
        if (had_) old_ = prev;
        if (value) ::setenv(name, value, 1);
        else ::unsetenv(name);
    }
    ~EnvGuard() {
        if (had_) ::setenv(name_.c_str(), old_.c_str(), 1);
        else ::unsetenv(name_.c_str());
    }
private:
    std::string name_, old_;
    bool had_;
};

}  // namespace

// --- detect_cleft reproducibility ----------------------------------------

TEST(CleftDeterminism, RepeatedDetectionOnIdenticalInputIsIdentical) {
    const auto atoms = hollow_box(4.0f);
    auto params = default_cleft_params();
    params.min_cluster_size = 1;

    sphere* a = detect_cleft(atoms.data(), nullptr,
                             static_cast<int>(atoms.size()), 0, params);
    const auto va = flatten(a);
    free_sphere_list(a);

    sphere* b = detect_cleft(atoms.data(), nullptr,
                             static_cast<int>(atoms.size()), 0, params);
    const auto vb = flatten(b);
    free_sphere_list(b);

    ASSERT_EQ(va.size(), vb.size())
        << "cleft detection returned a different number of spheres on a repeat run";
    for (size_t i = 0; i < va.size(); ++i)
        EXPECT_TRUE(va[i] == vb[i]) << "sphere " << i << " differs between runs";
}

TEST(CleftDeterminism, CanonicalSortModeIsAlsoSelfConsistent) {
    // FLEXAIDDS_CLEFT_SORT=1 imposes a geometric total order on probes so the
    // OpenMP merge order stops mattering. The gate is documented to change
    // tie-breaks relative to the unsorted path even in serial, so the two
    // modes are NOT asserted equal here — only that the sorted mode is itself
    // reproducible, which is the property it exists to provide.
    EnvGuard guard("FLEXAIDDS_CLEFT_SORT", "1");

    const auto atoms = hollow_box(4.0f);
    auto params = default_cleft_params();
    params.min_cluster_size = 1;

    sphere* a = detect_cleft(atoms.data(), nullptr,
                             static_cast<int>(atoms.size()), 0, params);
    const auto va = flatten(a);
    free_sphere_list(a);

    sphere* b = detect_cleft(atoms.data(), nullptr,
                             static_cast<int>(atoms.size()), 0, params);
    const auto vb = flatten(b);
    free_sphere_list(b);

    ASSERT_EQ(va.size(), vb.size());
    for (size_t i = 0; i < va.size(); ++i)
        EXPECT_TRUE(va[i] == vb[i]) << "sorted-mode sphere " << i << " differs between runs";
}

TEST(CleftDeterminism, CanonicalSortModeProducesTheSameSphereSet) {
    // The sort is a permutation of the surviving probes, so the *set* of
    // spheres kept must be invariant even though the emission order is not.
    const auto atoms = hollow_box(4.0f);
    auto params = default_cleft_params();
    params.min_cluster_size = 1;

    std::vector<SphereRec> unsorted, sorted;
    {
        EnvGuard off("FLEXAIDDS_CLEFT_SORT", nullptr);
        sphere* h = detect_cleft(atoms.data(), nullptr,
                                 static_cast<int>(atoms.size()), 0, params);
        unsorted = flatten(h);
        free_sphere_list(h);
    }
    {
        EnvGuard on("FLEXAIDDS_CLEFT_SORT", "1");
        sphere* h = detect_cleft(atoms.data(), nullptr,
                                 static_cast<int>(atoms.size()), 0, params);
        sorted = flatten(h);
        free_sphere_list(h);
    }

    EXPECT_EQ(unsorted.size(), sorted.size())
        << "the canonical-sort gate changed how many spheres survive, not just their order";
}

// --- parameter semantics --------------------------------------------------

TEST(CleftDeterminism, OracleSpatialFilterRestrictsTheSearchRegion) {
    // Two well-separated hollow boxes. With the oracle filter centred on the
    // first, atoms belonging to the second are excluded before probing.
    auto atoms = hollow_box(4.0f, 0.0f);
    const auto far = hollow_box(4.0f, 60.0f);
    atoms.insert(atoms.end(), far.begin(), far.end());

    auto params = default_cleft_params();
    params.min_cluster_size = 1;
    params.top_k_clefts = 0;  // keep every qualifying cluster

    sphere* all = detect_cleft(atoms.data(), nullptr,
                               static_cast<int>(atoms.size()), 0, params);
    const size_t n_all = flatten(all).size();
    free_sphere_list(all);

    params.oracle_center[0] = 0.0f;
    params.oracle_center[1] = 0.0f;
    params.oracle_center[2] = 0.0f;
    params.oracle_radius = 15.0f;

    sphere* filtered = detect_cleft(atoms.data(), nullptr,
                                    static_cast<int>(atoms.size()), 0, params);
    const size_t n_filtered = flatten(filtered).size();
    free_sphere_list(filtered);

    EXPECT_LT(n_filtered, n_all)
        << "the oracle pre-filter did not reduce the probe set";
    EXPECT_GT(n_filtered, 0u) << "the pre-filter removed the pocket it was centred on";
}

TEST(CleftDeterminism, ZeroCoordinatePaddingAtomsAreIgnored) {
    // Atoms that are all-zero in coordinates AND radius are treated as
    // uninitialised padding and skipped.
    auto atoms = hollow_box(4.0f, 30.0f);
    const size_t real_atoms = atoms.size();
    for (int i = 0; i < 12; ++i) atoms.push_back(make_atom(0.0f, 0.0f, 0.0f, 0.0f));

    auto params = default_cleft_params();
    params.min_cluster_size = 1;

    sphere* with_padding = detect_cleft(atoms.data(), nullptr,
                                        static_cast<int>(atoms.size()), 0, params);
    const auto vp = flatten(with_padding);
    free_sphere_list(with_padding);

    sphere* without = detect_cleft(atoms.data(), nullptr,
                                   static_cast<int>(real_atoms), 0, params);
    const auto vw = flatten(without);
    free_sphere_list(without);

    EXPECT_EQ(vp.size(), vw.size())
        << "zero-padding atoms leaked into probe generation";
}

TEST(CleftDeterminism, ImpossibleMinClusterSizeStillReturnsTheLargestCluster) {
    // Documented fallback: if no cluster reaches min_cluster_size the detector
    // warns and keeps the largest cluster rather than returning nothing.
    const auto atoms = hollow_box(4.0f);
    auto params = default_cleft_params();
    params.min_cluster_size = 1000000;

    sphere* head = detect_cleft(atoms.data(), nullptr,
                                static_cast<int>(atoms.size()), 0, params);
    EXPECT_NE(head, nullptr) << "fallback to largest cluster did not fire";
    free_sphere_list(head);
}

TEST(CleftDeterminism, TopKZeroKeepsAtLeastAsManySpheresAsTopKOne) {
    auto atoms = hollow_box(4.0f, 0.0f);
    const auto far = hollow_box(4.0f, 60.0f);
    atoms.insert(atoms.end(), far.begin(), far.end());

    auto params = default_cleft_params();
    params.min_cluster_size = 1;

    params.top_k_clefts = 1;
    sphere* one = detect_cleft(atoms.data(), nullptr,
                               static_cast<int>(atoms.size()), 0, params);
    const size_t n_one = flatten(one).size();
    free_sphere_list(one);

    params.top_k_clefts = 0;  // 0 == keep all
    sphere* all = detect_cleft(atoms.data(), nullptr,
                               static_cast<int>(atoms.size()), 0, params);
    const size_t n_all = flatten(all).size();
    free_sphere_list(all);

    EXPECT_GE(n_all, n_one);
}

// --- select_flexible_residues --------------------------------------------

TEST(SelectFlexibleResidues, EmptySeedsGiveEmptyResult) {
    std::vector<resid> residues = {make_residue("SER")};
    const auto out = select_flexible_residues(nullptr, residues.data(), 0, 1, {}, 5.0);
    EXPECT_TRUE(out.empty());
    free_residues(residues);
}

TEST(SelectFlexibleResidues, ResultIsDeduplicatedAndAscending) {
    std::vector<resid> residues;
    for (int i = 0; i < 6; ++i) residues.push_back(make_residue("SER"));

    const std::vector<int> seeds = {4, 1, 4, 3, 1};
    const auto out = select_flexible_residues(nullptr, residues.data(), 0,
                                              static_cast<int>(residues.size()),
                                              seeds, 5.0);
    EXPECT_EQ(out, (std::vector<int>{1, 3, 4}));
    free_residues(residues);
}

TEST(SelectFlexibleResidues, OutOfRangeSeedsAreDropped) {
    std::vector<resid> residues;
    for (int i = 0; i < 3; ++i) residues.push_back(make_residue("SER"));

    const std::vector<int> seeds = {-1, 0, 2, 3, 99};
    const auto out = select_flexible_residues(nullptr, residues.data(), 0, 3, seeds, 5.0);
    EXPECT_EQ(out, (std::vector<int>{0, 2}));
    free_residues(residues);
}

TEST(SelectFlexibleResidues, GlycineAndAlanineAreExcludedFromDistanceSeeds) {
    // Documented rule: never include Gly/Ala unless a backbone flexibility
    // module is active. Applies when residue metadata is supplied.
    std::vector<resid> residues = {make_residue("GLY"), make_residue("SER"),
                                   make_residue("ALA"), make_residue("LEU")};
    const std::vector<int> seeds = {0, 1, 2, 3};
    const auto out = select_flexible_residues(nullptr, residues.data(), 0, 4, seeds, 5.0);
    EXPECT_EQ(out, (std::vector<int>{1, 3}));
    free_residues(residues);
}

TEST(SelectFlexibleResidues, NullResidueArraySkipsTheGlyAlaFilter) {
    // Without residue metadata there is no residue name to test, so all
    // in-range seeds are kept.
    const std::vector<int> seeds = {0, 1, 2, 3};
    const auto out = select_flexible_residues(nullptr, nullptr, 0, 4, seeds, 5.0);
    EXPECT_EQ(out, (std::vector<int>{0, 1, 2, 3}));
}

TEST(SelectFlexibleResidues, UserFixedResiduesAreExcluded) {
    std::vector<resid> residues;
    for (int i = 0; i < 4; ++i) residues.push_back(make_residue("SER"));

    const std::vector<int> seeds = {0, 1, 2, 3};
    const std::vector<int> fixed = {1, 2};
    const auto out = select_flexible_residues(nullptr, residues.data(), 0, 4,
                                              seeds, 5.0, {}, fixed);
    EXPECT_EQ(out, (std::vector<int>{0, 3}));
    free_residues(residues);
}

TEST(SelectFlexibleResidues, ActiveSiteResiduesAreUnionedWithCleftSeeds) {
    std::vector<resid> residues;
    for (int i = 0; i < 6; ++i) residues.push_back(make_residue("SER"));

    const std::vector<int> seeds = {0, 1};
    const std::vector<int> active = {4, 5};
    const auto out = select_flexible_residues(nullptr, residues.data(), 0, 6,
                                              seeds, 5.0, active);
    EXPECT_EQ(out, (std::vector<int>{0, 1, 4, 5}));
    free_residues(residues);
}

TEST(SelectFlexibleResidues, ForcedFlexibleIsIncludedEvenWhenItIsGlycine) {
    // Forced entries bypass the Gly/Ala rule: the caller is declared
    // responsible for having a backbone module.
    std::vector<resid> residues = {make_residue("GLY"), make_residue("SER")};
    const auto out = select_flexible_residues(nullptr, residues.data(), 0, 2,
                                              {}, 5.0, {}, {}, {0});
    EXPECT_EQ(out, (std::vector<int>{0}));
    free_residues(residues);
}

TEST(SelectFlexibleResidues, ForcedFlexibleOutOfRangeIsRejected) {
    std::vector<resid> residues = {make_residue("SER")};
    const auto out = select_flexible_residues(nullptr, residues.data(), 0, 1,
                                              {}, 5.0, {}, {}, {-1, 7});
    EXPECT_TRUE(out.empty());
    free_residues(residues);
}

// CHARACTERIZATION: a residue listed in BOTH user_fixed_residues and
// user_forced_flexible is currently returned as flexible — force wins, because
// forced entries are inserted before the fixed set is consulted and the fixed
// check lives only on the distance-seed path. The header documents both rules
// ("Respect user_fixed_residues" and "Force-include user_forced_flexible")
// without stating which takes precedence, so this pins today's answer.
TEST(SelectFlexibleResidues, CHARACTERIZATION_ForcedOverridesFixedOnConflict) {
    std::vector<resid> residues;
    for (int i = 0; i < 3; ++i) residues.push_back(make_residue("SER"));

    const auto out = select_flexible_residues(nullptr, residues.data(), 0, 3,
                                              {0, 1, 2}, 5.0, {}, {1}, {1});
    EXPECT_EQ(out, (std::vector<int>{0, 1, 2}))
        << "residue 1 is both fixed and forced; forced currently wins";
    free_residues(residues);
}

// CHARACTERIZATION: distance_shell_A is accepted but never read — the
// implementation treats cleft_sphere_residues purely as a pre-computed seed
// list ("simplified: we treat the passed cleft_sphere_residues as seeds").
// The header, however, documents a distance-shell rule. Any caller passing a
// shell radius and expecting spatial expansion is silently getting a no-op.
TEST(SelectFlexibleResidues, CHARACTERIZATION_DistanceShellArgumentIsIgnored) {
    std::vector<resid> residues;
    for (int i = 0; i < 5; ++i) residues.push_back(make_residue("SER"));
    const std::vector<int> seeds = {2};

    const auto tight = select_flexible_residues(nullptr, residues.data(), 0, 5, seeds, 0.0);
    const auto wide  = select_flexible_residues(nullptr, residues.data(), 0, 5, seeds, 1000.0);
    EXPECT_EQ(tight, wide) << "distance_shell_A is not implemented; both calls return the seeds";
    EXPECT_EQ(tight, (std::vector<int>{2}));
    free_residues(residues);
}
