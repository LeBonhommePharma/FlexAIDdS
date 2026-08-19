// test_instream_clustering.cpp — unit tests for flexaids::InStreamCluster
//
// InStreamCluster maintains a bounded medoid set online during GA evolution.
// It sits on the determinism-critical path: the medoid set it produces feeds
// cluster-representative selection, and it contains two order-sensitive
// constructs of the same family as the sorts flagged in
// OPTIMIZATION_KNOWN_ISSUES.md:
//
//   * finalize() sorts medoids with std::sort on a bare double best_score
//     (unstable, no secondary key)
//   * merge_closest_pair() removes the merged medoid by swap-with-last,
//     which permutes the surviving medoid order
//
// All tests are deterministic and fixture-free. Tests that pin behaviour the
// code does not explicitly promise are labelled CHARACTERIZATION.
//
// Apache-2.0 — see LICENSE.

#include "InStreamClustering.h"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <vector>

using flexaids::ClusterMedoid;
using flexaids::InStreamCluster;

namespace {

std::vector<float> flat(std::initializer_list<std::initializer_list<float>> rows) {
    std::vector<float> out;
    for (const auto& r : rows)
        for (float v : r) out.push_back(v);
    return out;
}

}  // namespace

// --- construction / reset -------------------------------------------------

TEST(InStreamCluster, StartsEmpty) {
    InStreamCluster isc(2.0f, 100, 6);
    EXPECT_EQ(isc.cluster_count(), 0);
    EXPECT_EQ(isc.total_merged(), 0);
    EXPECT_FLOAT_EQ(isc.rmsd_threshold(), 2.0f);
    EXPECT_TRUE(isc.snapshot().empty());
}

TEST(InStreamCluster, NonPositiveMaxRepresentativesFallsBackToDefault) {
    InStreamCluster isc(2.0f, 0, 3);
    const std::vector<float> genes = flat({{0, 0, 0}});
    const double score = -1.0;
    isc.merge_elites(genes.data(), &score, 1, 0, 3);
    EXPECT_EQ(isc.cluster_count(), 1);
}

TEST(InStreamCluster, ResetClearsMedoidsAndCounter) {
    InStreamCluster isc(2.0f, 100, 3);
    const std::vector<float> genes = flat({{0, 0, 0}, {50, 50, 50}});
    const double scores[2] = {-1.0, -2.0};
    isc.merge_elites(genes.data(), scores, 2, 1, 3);
    ASSERT_EQ(isc.cluster_count(), 2);

    isc.reset();
    EXPECT_EQ(isc.cluster_count(), 0);
    EXPECT_EQ(isc.total_merged(), 0);
}

// --- input guards ---------------------------------------------------------

TEST(InStreamCluster, RejectsEmptyAndNullInput) {
    InStreamCluster isc(2.0f, 100, 3);
    const std::vector<float> genes = flat({{0, 0, 0}});
    const double score = -1.0;

    isc.merge_elites(genes.data(), &score, 0, 0, 3);
    isc.merge_elites(genes.data(), &score, 1, 0, 0);
    isc.merge_elites(nullptr, &score, 1, 0, 3);
    isc.merge_elites(genes.data(), nullptr, 1, 0, 3);

    EXPECT_EQ(isc.cluster_count(), 0);
    EXPECT_EQ(isc.total_merged(), 0);
}

TEST(InStreamCluster, SkipsNonFiniteScoresWithoutCountingThem) {
    InStreamCluster isc(2.0f, 100, 3);
    const std::vector<float> genes = flat({{0, 0, 0}, {50, 50, 50}, {99, 99, 99}});
    const double scores[3] = {std::numeric_limits<double>::quiet_NaN(),
                              std::numeric_limits<double>::infinity(),
                              -3.0};
    isc.merge_elites(genes.data(), scores, 3, 7, 3);

    EXPECT_EQ(isc.cluster_count(), 1);
    EXPECT_EQ(isc.total_merged(), 1);
    ASSERT_EQ(isc.snapshot().size(), 1u);
    EXPECT_DOUBLE_EQ(isc.snapshot()[0].best_score, -3.0);
}

TEST(InStreamCluster, NumGenesInferredOnFirstCallWhenNotSetAtConstruction) {
    InStreamCluster isc(2.0f, 100, 0);
    const std::vector<float> genes = flat({{1, 2, 3, 4}});
    const double score = -1.0;
    isc.merge_elites(genes.data(), &score, 1, 0, 4);

    ASSERT_EQ(isc.cluster_count(), 1);
    EXPECT_EQ(isc.snapshot()[0].genes_ic.size(), 4u);
}

// --- threshold semantics --------------------------------------------------

TEST(InStreamCluster, NearbyElitesCollapseIntoOneCluster) {
    InStreamCluster isc(2.0f, 100, 3);
    const std::vector<float> genes = flat({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}});
    const double scores[3] = {-1.0, -2.0, -0.5};
    isc.merge_elites(genes.data(), scores, 3, 5, 3);

    ASSERT_EQ(isc.cluster_count(), 1);
    EXPECT_EQ(isc.total_merged(), 3);
    const ClusterMedoid& m = isc.snapshot()[0];
    EXPECT_EQ(m.member_count, 3);
    EXPECT_DOUBLE_EQ(m.best_score, -2.0);
}

TEST(InStreamCluster, DistantElitesFormSeparateClusters) {
    InStreamCluster isc(2.0f, 100, 3);
    const std::vector<float> genes = flat({{0, 0, 0}, {100, 0, 0}, {0, 100, 0}});
    const double scores[3] = {-1.0, -2.0, -3.0};
    isc.merge_elites(genes.data(), scores, 3, 0, 3);
    EXPECT_EQ(isc.cluster_count(), 3);
    EXPECT_EQ(isc.total_merged(), 3);
}

TEST(InStreamCluster, MedoidCentreMovesToBetterScoringMember) {
    InStreamCluster isc(5.0f, 100, 3);
    const std::vector<float> genes = flat({{0, 0, 0}, {1, 1, 1}});
    const double scores[2] = {-1.0, -9.0};
    isc.merge_elites(genes.data(), scores, 2, 3, 3);

    ASSERT_EQ(isc.cluster_count(), 1);
    const ClusterMedoid& m = isc.snapshot()[0];
    EXPECT_DOUBLE_EQ(m.best_score, -9.0);
    EXPECT_FLOAT_EQ(m.genes_ic[0], 1.0f);
    EXPECT_FLOAT_EQ(m.genes_ic[1], 1.0f);
    EXPECT_FLOAT_EQ(m.genes_ic[2], 1.0f);
}

TEST(InStreamCluster, MedoidCentreStaysWhenNewMemberIsWorse) {
    InStreamCluster isc(5.0f, 100, 3);
    const std::vector<float> genes = flat({{0, 0, 0}, {1, 1, 1}});
    const double scores[2] = {-9.0, -1.0};
    isc.merge_elites(genes.data(), scores, 2, 3, 3);

    ASSERT_EQ(isc.cluster_count(), 1);
    const ClusterMedoid& m = isc.snapshot()[0];
    EXPECT_DOUBLE_EQ(m.best_score, -9.0);
    EXPECT_FLOAT_EQ(m.genes_ic[0], 0.0f);
    EXPECT_EQ(m.member_count, 2);
}

TEST(InStreamCluster, GenerationBookkeeping) {
    InStreamCluster isc(5.0f, 100, 3);
    const std::vector<float> a = flat({{0, 0, 0}});
    const double sa = -1.0;
    isc.merge_elites(a.data(), &sa, 1, 10, 3);

    const std::vector<float> b = flat({{1, 0, 0}});
    const double sb = -2.0;
    isc.merge_elites(b.data(), &sb, 1, 25, 3);

    ASSERT_EQ(isc.cluster_count(), 1);
    const ClusterMedoid& m = isc.snapshot()[0];
    EXPECT_EQ(m.first_seen_gen, 10);
    EXPECT_EQ(m.last_updated_gen, 25);
}

// --- angular gene handling ------------------------------------------------

TEST(InStreamCluster, AngularGenesWrapAcrossThePeriodicBoundary) {
    InStreamCluster isc(0.5f, 100, 4);
    const std::vector<float> wrapped = flat({{0, 0, 0, 179.0f}, {0, 0, 0, -179.0f}});
    const double scores[2] = {-1.0, -2.0};
    isc.merge_elites(wrapped.data(), scores, 2, 0, 4);
    EXPECT_EQ(isc.cluster_count(), 1)
        << "179 and -179 degrees should wrap to a 2-degree gap";

    InStreamCluster far(0.5f, 100, 4);
    const std::vector<float> apart = flat({{0, 0, 0, 0.0f}, {0, 0, 0, 90.0f}});
    far.merge_elites(apart.data(), scores, 2, 0, 4);
    EXPECT_EQ(far.cluster_count(), 2);
}

TEST(InStreamCluster, TranslationDominatesAngularContribution) {
    InStreamCluster torsion(1.0f, 100, 4);
    const std::vector<float> t = flat({{0, 0, 0, 0.0f}, {0, 0, 0, 10.0f}});
    const double scores[2] = {-1.0, -2.0};
    torsion.merge_elites(t.data(), scores, 2, 0, 4);
    EXPECT_EQ(torsion.cluster_count(), 1);

    InStreamCluster shift(1.0f, 100, 4);
    const std::vector<float> s = flat({{0, 0, 0, 0.0f}, {10, 0, 0, 0.0f}});
    shift.merge_elites(s.data(), scores, 2, 0, 4);
    EXPECT_EQ(shift.cluster_count(), 2);
}

TEST(InStreamCluster, FewerThanThreeGenesUsesTranslationPathOnly) {
    InStreamCluster isc(1.0f, 100, 2);
    const std::vector<float> genes = flat({{0, 0}, {5, 0}});
    const double scores[2] = {-1.0, -2.0};
    isc.merge_elites(genes.data(), scores, 2, 0, 2);
    EXPECT_EQ(isc.cluster_count(), 2);
}

// --- capacity bound -------------------------------------------------------

TEST(InStreamCluster, CapacityBoundIsNeverExceeded) {
    const int kCap = 4;
    InStreamCluster isc(0.01f, kCap, 3);

    std::vector<float> genes;
    std::vector<double> scores;
    for (int i = 0; i < 40; ++i) {
        genes.push_back(static_cast<float>(i) * 10.0f);
        genes.push_back(0.0f);
        genes.push_back(0.0f);
        scores.push_back(-static_cast<double>(i));
    }
    isc.merge_elites(genes.data(), scores.data(), 40, 2, 3);

    EXPECT_LE(isc.cluster_count(), kCap);
    EXPECT_EQ(isc.total_merged(), 40);
}

TEST(InStreamCluster, CapacityMergeKeepsTheBestScoreSeen) {
    const int kCap = 2;
    InStreamCluster isc(0.01f, kCap, 3);

    std::vector<float> genes;
    std::vector<double> scores;
    for (int i = 0; i < 12; ++i) {
        genes.push_back(static_cast<float>(i) * 10.0f);
        genes.push_back(0.0f);
        genes.push_back(0.0f);
        scores.push_back(-static_cast<double>(i));
    }
    isc.merge_elites(genes.data(), scores.data(), 12, 0, 3);

    const auto medoids = isc.finalize();
    ASSERT_FALSE(medoids.empty());
    EXPECT_DOUBLE_EQ(medoids.front().best_score, -11.0);

    int total_members = 0;
    for (const auto& m : medoids) total_members += m.member_count;
    EXPECT_EQ(total_members, 12);
}

// --- finalize() ordering --------------------------------------------------

TEST(InStreamCluster, FinalizeSortsByBestScoreAscending) {
    InStreamCluster isc(0.01f, 100, 3);
    const std::vector<float> genes =
        flat({{0, 0, 0}, {100, 0, 0}, {200, 0, 0}, {300, 0, 0}});
    const double scores[4] = {-1.0, -9.0, -5.0, -3.0};
    isc.merge_elites(genes.data(), scores, 4, 0, 3);

    const auto medoids = isc.finalize();
    ASSERT_EQ(medoids.size(), 4u);
    for (size_t i = 1; i < medoids.size(); ++i)
        EXPECT_LE(medoids[i - 1].best_score, medoids[i].best_score);
    EXPECT_DOUBLE_EQ(medoids.front().best_score, -9.0);
}

TEST(InStreamCluster, FinalizeIsIdempotent) {
    InStreamCluster isc(0.01f, 100, 3);
    const std::vector<float> genes = flat({{0, 0, 0}, {100, 0, 0}, {200, 0, 0}});
    const double scores[3] = {-1.0, -3.0, -2.0};
    isc.merge_elites(genes.data(), scores, 3, 0, 3);

    const auto first = isc.finalize();
    const auto second = isc.finalize();
    ASSERT_EQ(first.size(), second.size());
    for (size_t i = 0; i < first.size(); ++i)
        EXPECT_DOUBLE_EQ(first[i].best_score, second[i].best_score);
}

// CHARACTERIZATION: finalize() uses std::sort (unstable) on a bare double with
// no secondary key. When medoids share a best_score their relative order is
// not defined by the algorithm. This pins only the property the pipeline
// actually relies on -- repeated calls on the same object agree -- and does
// not assert a tie order the code never promised.
TEST(InStreamCluster, CHARACTERIZATION_TiedScoresOnlyGuaranteeRunStability) {
    InStreamCluster isc(0.01f, 100, 3);
    const std::vector<float> genes =
        flat({{0, 0, 0}, {100, 0, 0}, {200, 0, 0}, {300, 0, 0}});
    const double scores[4] = {-4.0, -4.0, -4.0, -4.0};
    isc.merge_elites(genes.data(), scores, 4, 0, 3);

    const auto a = isc.finalize();
    const auto b = isc.finalize();
    ASSERT_EQ(a.size(), 4u);
    for (size_t i = 0; i < a.size(); ++i)
        EXPECT_FLOAT_EQ(a[i].genes_ic[0], b[i].genes_ic[0]);
}

// --- whole-pipeline determinism -------------------------------------------

TEST(InStreamCluster, IdenticalInputProducesIdenticalMedoidSet) {
    auto run = [] {
        InStreamCluster isc(1.5f, 8, 5);
        std::vector<float> genes;
        std::vector<double> scores;
        for (int i = 0; i < 60; ++i) {
            genes.push_back(static_cast<float>((i * 7) % 23));
            genes.push_back(static_cast<float>((i * 13) % 19));
            genes.push_back(static_cast<float>((i * 3) % 11));
            genes.push_back(static_cast<float>((i * 31) % 360));
            genes.push_back(static_cast<float>((i * 17) % 360));
            scores.push_back(-static_cast<double>((i * 11) % 37));
        }
        isc.merge_elites(genes.data(), scores.data(), 60, 4, 5);
        return isc.finalize();
    };

    const auto a = run();
    const auto b = run();
    ASSERT_EQ(a.size(), b.size());
    for (size_t i = 0; i < a.size(); ++i) {
        EXPECT_DOUBLE_EQ(a[i].best_score, b[i].best_score);
        EXPECT_EQ(a[i].member_count, b[i].member_count);
        ASSERT_EQ(a[i].genes_ic.size(), b[i].genes_ic.size());
        for (size_t g = 0; g < a[i].genes_ic.size(); ++g)
            EXPECT_FLOAT_EQ(a[i].genes_ic[g], b[i].genes_ic[g]);
    }
}

TEST(InStreamCluster, BatchingDoesNotChangeTheClusterCount) {
    const std::vector<float> genes = flat({{0, 0, 0},
                                           {1, 0, 0},
                                           {100, 0, 0},
                                           {101, 0, 0},
                                           {200, 0, 0},
                                           {201, 0, 0}});
    const double scores[6] = {-1.0, -2.0, -3.0, -4.0, -5.0, -6.0};

    InStreamCluster one(2.0f, 100, 3);
    one.merge_elites(genes.data(), scores, 6, 0, 3);

    InStreamCluster many(2.0f, 100, 3);
    for (int b = 0; b < 3; ++b)
        many.merge_elites(genes.data() + b * 6, scores + b * 2, 2, b, 3);

    EXPECT_EQ(one.cluster_count(), many.cluster_count());
    EXPECT_EQ(one.total_merged(), many.total_merged());
}

// CHARACTERIZATION: a medoid whose stored gene vector length differs from the
// num_genes of the incoming batch is skipped by the nearest-medoid scan, so a
// new medoid is always created rather than the call being rejected. Mixing
// gene counts on one instance silently partitions the medoid set.
TEST(InStreamCluster, CHARACTERIZATION_MismatchedGeneCountCreatesDisjointMedoids) {
    InStreamCluster isc(50.0f, 100, 0);
    const std::vector<float> three = flat({{0, 0, 0}});
    const double s1 = -1.0;
    isc.merge_elites(three.data(), &s1, 1, 0, 3);

    const std::vector<float> four = flat({{0, 0, 0, 0}});
    const double s2 = -2.0;
    isc.merge_elites(four.data(), &s2, 1, 1, 4);

    EXPECT_EQ(isc.cluster_count(), 2)
        << "gene-count mismatch is silently tolerated, not merged or rejected";
}
