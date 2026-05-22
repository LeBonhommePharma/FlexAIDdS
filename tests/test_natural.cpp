// tests/test_natural.cpp
// Unit tests for the NATURaL module:
//   RibosomeElongation (Zhao 2011 master equation)
//   TransloconInsertion (Hessa 2007 scale)
//   NucleationSeedDetector (RNA hairpins, G-quads, protein helix/hydrophobic)
//   NATURaLConfig auto_configure helpers
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>

#include "../LIB/NATURaL/RibosomeElongation.h"
#include "../LIB/NATURaL/TransloconInsertion.h"
#include "../LIB/NATURaL/NucleationDetector.h"
#include "../LIB/NATURaL/NATURaLDualAssembly.h"

#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <limits>

using namespace ribosome;
using namespace translocon;
using namespace natural;

static constexpr double TOL = 1e-6;

// ===========================================================================
// CodonRateTable
// ===========================================================================

TEST(CodonRateTable, BuildEcoliSucceeds) {
    EXPECT_NO_THROW(CodonRateTable::build_ecoli());
}

TEST(CodonRateTable, BuildHumanSucceeds) {
    EXPECT_NO_THROW(CodonRateTable::build_human());
}

TEST(CodonRateTable, EcoliMeanRateInPhysiologicalRange) {
    auto tbl = CodonRateTable::build_ecoli();
    // Wohlgemuth 2008: E. coli mean ≈ 10–20 aa/s; our constant is 16.5
    EXPECT_GT(tbl.mean_rate_aa_per_s, 5.0);
    EXPECT_LT(tbl.mean_rate_aa_per_s, 30.0);
}

TEST(CodonRateTable, HumanMeanRateLowerThanEcoli) {
    auto ecoli = CodonRateTable::build_ecoli();
    auto human = CodonRateTable::build_human();
    EXPECT_LT(human.mean_rate_aa_per_s, ecoli.mean_rate_aa_per_s);
}

TEST(CodonRateTable, UnknownCodonReturnsMeanRate) {
    auto tbl = CodonRateTable::build_ecoli();
    double r = tbl.rate("XYZ");
    // Should fall back to mean rather than crash
    EXPECT_GT(r, 0.0);
    EXPECT_TRUE(std::isfinite(r));
}

TEST(CodonRateTable, MeanRateOverSequence) {
    auto tbl = CodonRateTable::build_ecoli();
    std::vector<std::string> codons(10, "AAA"); // all lysine (Lys)
    double mean = tbl.mean_rate(codons);
    EXPECT_GT(mean, 0.0);
    EXPECT_TRUE(std::isfinite(mean));
}

TEST(CodonRateTable, PauseSitesReturnedAreSubset) {
    auto tbl = CodonRateTable::build_ecoli();
    std::vector<std::string> codons = {"AAA","CGA","AAA","CGA","AAA"};
    auto pauses = tbl.pause_sites(codons);
    for (int idx : pauses) {
        EXPECT_GE(idx, 0);
        EXPECT_LT(idx, static_cast<int>(codons.size()));
    }
}

TEST(CodonRateTable, MFCMapNonEmpty) {
    auto mfc = CodonRateTable::aa_to_mfc(Organism::EcoliK12);
    EXPECT_FALSE(mfc.empty());
    // Common amino acids should be present
    EXPECT_NE(mfc.find('A'), mfc.end()); // Ala
    EXPECT_NE(mfc.find('G'), mfc.end()); // Gly
}

// ===========================================================================
// RibosomeElongation
// ===========================================================================

class RibosomeElongationTest : public ::testing::Test {
protected:
    static constexpr int N = 60; // short protein for speed
    CodonRateTable tbl = CodonRateTable::build_ecoli();
    std::string seq = std::string(N, 'A'); // all-Ala
    std::vector<std::string> codons; // empty → use mean rates

    RibosomeElongation make_engine() {
        return RibosomeElongation(seq, codons, tbl);
    }
};

TEST_F(RibosomeElongationTest, ConstructionSucceeds) {
    EXPECT_NO_THROW(make_engine());
}

TEST_F(RibosomeElongationTest, NResiduesMatchesInput) {
    auto eng = make_engine();
    EXPECT_EQ(eng.n_residues(), N);
}

TEST_F(RibosomeElongationTest, ElongationRatesPositive) {
    auto eng = make_engine();
    for (double r : eng.elongation_rates()) {
        EXPECT_GT(r, 0.0);
        EXPECT_TRUE(std::isfinite(r));
    }
}

TEST_F(RibosomeElongationTest, MeanArrivalTimeMonotonicIncreasing) {
    auto eng = make_engine();
    for (int i = 1; i < N; ++i)
        EXPECT_GE(eng.mean_arrival_time(i), eng.mean_arrival_time(i - 1));
}

TEST_F(RibosomeElongationTest, MeanTotalTimeEqualsAnalyticSum) {
    auto eng = make_engine();
    // Analytic: T = 1/k_ini + Σ 1/k_n + 1/k_ter
    double analytic = 1.0 / K_INI_DEFAULT;
    for (double k : eng.elongation_rates())
        analytic += 1.0 / k;
    analytic += 1.0 / K_TERM_DEFAULT;

    EXPECT_NEAR(eng.mean_total_time(), analytic, 1e-6);
}

TEST_F(RibosomeElongationTest, ValidateMasterEquationPasses) {
    auto vr = validate_master_equation(30, Organism::EcoliK12);
    EXPECT_TRUE(vr.passed) << vr.message;
    EXPECT_LT(vr.relative_error, 0.15); // within 15%
}

TEST_F(RibosomeElongationTest, FoldingWindowsNotEmpty) {
    auto eng = make_engine();
    auto windows = eng.folding_windows();
    // For a 60-aa protein with a 34-aa tunnel, there should be folding windows
    EXPECT_FALSE(windows.empty());
}

TEST_F(RibosomeElongationTest, FoldingWindowProbabilitiesInRange) {
    auto eng = make_engine();
    for (const auto& w : eng.folding_windows()) {
        EXPECT_GE(w.p_folded_cotrans, 0.0);
        EXPECT_LE(w.p_folded_cotrans, 1.0);
        EXPECT_GT(w.t_available, 0.0);
    }
}

TEST_F(RibosomeElongationTest, TimeWeightedScoreFinite) {
    auto eng = make_engine();
    double score = eng.time_weighted_score([](int){ return -10.0; });
    EXPECT_TRUE(std::isfinite(score));
}

TEST_F(RibosomeElongationTest, TimeWeightedScoreConstantEqualsConstant) {
    auto eng = make_engine();
    double C = -5.5;
    // ∫ C dw = C (time-weighted mean of constant is the constant)
    double score = eng.time_weighted_score([C](int){ return C; });
    EXPECT_NEAR(score, C, 0.01);
}

TEST_F(RibosomeElongationTest, PauseSitesRatesBelowThreshold) {
    auto eng = make_engine();
    double mean_r = 0.0;
    for (double r : eng.elongation_rates()) mean_r += r;
    mean_r /= eng.n_residues();

    for (int idx : eng.pause_sites()) {
        EXPECT_LT(eng.elongation_rates()[idx],
                  mean_r * RIBOSOME_PAUSE_THRESHOLD * 1.01); // 1% tolerance
    }
}

TEST(NATURaLGrowthEntropy, ReturnsNatsNotBits) {
    const double h2 = growth_entropy_nats({0.0, 1.0});
    const double h4 = growth_entropy_nats({0.0, 1.0, 2.0, 3.0});

    EXPECT_NEAR(h2, std::log(2.0), 1e-9);
    EXPECT_NEAR(h4, std::log(4.0), 1e-9);
}

TEST(NATURaLGrowthEntropy, RepeatedEqualOccupancyBinsReturnAnalyticNats) {
    const std::vector<double> two_bins = {
        -1.0, -1.0, -1.0, -1.0,
         1.0,  1.0,  1.0,  1.0
    };
    const std::vector<double> four_bins = {
        0.0, 0.0,
        1.0, 1.0,
        2.0, 2.0,
        3.0, 3.0
    };

    EXPECT_NEAR(growth_entropy_nats(two_bins), std::log(2.0), 1e-9);
    EXPECT_NEAR(growth_entropy_nats(four_bins), std::log(4.0), 1e-9);
}

TEST(NATURaLGrowthEntropy, BinBoundaryAtMaximumDoesNotCreateExtraEntropy) {
    const std::vector<double> samples = {
        0.0, 0.0, 0.0,
        31.0, 31.0, 31.0
    };

    EXPECT_NEAR(growth_entropy_nats(samples), std::log(2.0), 1e-9);
}

TEST(NATURaLGrowthEntropy, DegenerateTrajectoryHasZeroEntropy) {
    EXPECT_DOUBLE_EQ(growth_entropy_nats({}), 0.0);
    EXPECT_DOUBLE_EQ(growth_entropy_nats({-3.0, -3.0, -3.0}), 0.0);
}

// ===========================================================================
// TransloconInsertion
// ===========================================================================

class TransloconTest : public ::testing::Test {
protected:
    TransloconInsertion ti{310.0, 0.5, 34};

    // Hydrophobic TM window (should insert spontaneously)
    static std::string hydrophobic_tm() { return std::string(19, 'L'); }
    // Hydrophilic window (should not insert)
    static std::string hydrophilic() { return std::string(19, 'D'); }
};

TEST_F(TransloconTest, HydrophobicWindowInserts) {
    auto window = ti.check_window(hydrophobic_tm(), 0);
    EXPECT_LT(window.deltaG_insert, 0.0);  // negative ΔG → spontaneous
    EXPECT_GT(window.p_insert, 0.5);
    EXPECT_TRUE(window.is_inserted);
}

TEST_F(TransloconTest, HydrophilicWindowDoesNotInsert) {
    auto window = ti.check_window(hydrophilic(), 0);
    EXPECT_GT(window.deltaG_insert, 0.0);  // positive ΔG → disfavored
    EXPECT_LT(window.p_insert, 0.5);
    EXPECT_FALSE(window.is_inserted);
}

TEST_F(TransloconTest, InsertionProbabilityInRange) {
    for (const auto& seq : {std::string(19, 'L'), std::string(19, 'A'), std::string(19, 'D')}) {
        auto w = ti.check_window(seq, 0);
        EXPECT_GE(w.p_insert, 0.0);
        EXPECT_LE(w.p_insert, 1.0);
    }
}

TEST_F(TransloconTest, ScoreWindowFinite) {
    std::string seq = "LLLLLAAAAALLLLLAAAAALLLL";
    double score = ti.score_window(seq, 0, TM_WINDOW_LEN);
    EXPECT_TRUE(std::isfinite(score));
}

TEST_F(TransloconTest, ScanLongSequence) {
    // Should return a window for each valid position
    std::string seq = "MAAALLLLLLLLLLLLLLLLLAAA"; // 1 TM helix in middle
    EXPECT_NO_THROW(ti.scan(seq));
    auto windows = ti.scan(seq);
    EXPECT_FALSE(windows.empty());
}

TEST_F(TransloconTest, PositionWeightSymmetric) {
    // position_weight(0, 19) should equal position_weight(18, 19)
    double w0 = position_weight(0, TM_WINDOW_LEN);
    double w18 = position_weight(TM_WINDOW_LEN - 1, TM_WINDOW_LEN);
    EXPECT_NEAR(w0, w18, 1e-9);
}

TEST_F(TransloconTest, PositionWeightPeaksAtCenter) {
    double wc = position_weight(TM_WINDOW_LEN / 2, TM_WINDOW_LEN);
    double we = position_weight(0, TM_WINDOW_LEN);
    EXPECT_GT(wc, we);
}

TEST_F(TransloconTest, AccessorsMatchConstructorArgs) {
    EXPECT_NEAR(ti.temperature_K(), 310.0, TOL);
    EXPECT_NEAR(ti.insertion_threshold(), 0.5, TOL);
    EXPECT_EQ(ti.tunnel_length(), 34);
}

TEST_F(TransloconTest, HessaScaleLeuNegative) {
    // L (Leu) has one of the most negative ΔG values → should be < 0
    EXPECT_LT(HESSA_SCALE['L'], 0.0);
}

TEST_F(TransloconTest, HessaScaleAspNegativePositive) {
    // D (Asp) is highly hydrophilic → positive ΔG
    EXPECT_GT(HESSA_SCALE['D'], 0.0);
}

// ===========================================================================
// NucleationSeedDetector — protein detectors
// ===========================================================================

TEST(NucleationDetector, ProteinHydrophobicClusterDetected) {
    // ILVFMW run of 5
    auto seeds = NucleationSeedDetector::detect_protein_hydrophobic("AAILLLLLAAA");
    EXPECT_FALSE(seeds.empty());
    for (const auto& s : seeds) {
        EXPECT_EQ(s.type, NucleationSeed::Type::PROTEIN_HYDROPHOBIC);
        EXPECT_GE(s.folding_rate_boost, 1.0);
    }
}

TEST(NucleationDetector, ProteinHydrophobicShortRunNotDetected) {
    // Only 2 hydrophobic residues → below min_run=4
    auto seeds = NucleationSeedDetector::detect_protein_hydrophobic("AAILLAA");
    bool all_short = std::all_of(seeds.begin(), seeds.end(),
        [](const NucleationSeed& s){
            return s.type != NucleationSeed::Type::PROTEIN_HYDROPHOBIC
                || (s.end_pos - s.start_pos + 1) < 4;
        });
    // No 4+ hydrophobic run should be detected
    EXPECT_TRUE(seeds.empty() || all_short);
}

TEST(NucleationDetector, ProteinHelixDetected) {
    // AELKAAED — Ala and Glu/Lys are decent helix formers
    // Use a sequence with known good helix propensity: all Ala
    std::string helical(12, 'A');
    auto seeds = NucleationSeedDetector::detect_protein_helix(helical);
    // Ala P_α = ~1.45 (high), so all windows should be seeds
    EXPECT_FALSE(seeds.empty());
    for (const auto& s : seeds)
        EXPECT_EQ(s.type, NucleationSeed::Type::PROTEIN_HELIX);
}

TEST(NucleationDetector, ProteinHelixLowPropensityNotDetected) {
    // Proline is a helix breaker (P_α ≈ 0.57)
    std::string breaker(12, 'P');
    auto seeds = NucleationSeedDetector::detect_protein_helix(breaker);
    EXPECT_TRUE(seeds.empty());
}

TEST(NucleationDetector, HelixBoostAboveOne) {
    std::string helical(12, 'A');
    auto seeds = NucleationSeedDetector::detect_protein_helix(helical);
    for (const auto& s : seeds)
        EXPECT_GE(s.folding_rate_boost, 1.0);
}

// ===========================================================================
// NucleationSeedDetector — RNA detectors
// ===========================================================================

TEST(NucleationDetector, GQuadrupletDetected) {
    // Classic G3+ pattern: GGGNGGGNGGGNGGG (where N is spacer)
    std::string gquad = "GGGUGGGUGGGUGGG";
    auto seeds = NucleationSeedDetector::detect_rna_gquads(gquad);
    EXPECT_FALSE(seeds.empty());
    for (const auto& s : seeds)
        EXPECT_EQ(s.type, NucleationSeed::Type::RNA_GQUADRUPLEX);
}

TEST(NucleationDetector, GQuadBoostIsLarge) {
    std::string gquad = "GGGUGGGUGGGUGGG";
    auto seeds = NucleationSeedDetector::detect_rna_gquads(gquad);
    if (!seeds.empty())
        EXPECT_GT(seeds[0].folding_rate_boost, 2.0);
}

TEST(NucleationDetector, RNAHairpinDetected) {
    // Simple palindrome: GCGCAAAGCGC (stem GCGC + AAAA loop + complement GCGC)
    std::string hairpin = "GCGCAAAAGCGC";
    auto seeds = NucleationSeedDetector::detect_rna_hairpins(hairpin);
    EXPECT_FALSE(seeds.empty());
    for (const auto& s : seeds)
        EXPECT_EQ(s.type, NucleationSeed::Type::RNA_HAIRPIN);
}

TEST(NucleationDetector, PositionBoostMapLength) {
    std::string seq = "GGGUGGGUGGGUGGG";
    auto seeds = NucleationSeedDetector::detect_rna_gquads(seq);
    auto boost_map = NucleationSeedDetector::position_boost_map(seeds, static_cast<int>(seq.size()));
    EXPECT_EQ(boost_map.size(), seq.size());
}

TEST(NucleationDetector, PositionBoostMapBaselineOne) {
    std::string seq(20, 'A'); // no seeds
    auto seeds = NucleationSeedDetector::detect(seq, false);
    auto boost_map = NucleationSeedDetector::position_boost_map(seeds, static_cast<int>(seq.size()));
    for (double b : boost_map)
        EXPECT_GE(b, 1.0);
}

// ===========================================================================
// NATURaLDualAssembly — config helpers (no receptor/ligand atoms needed)
// ===========================================================================

TEST(NATURaLConfig, DefaultsAreReasonable) {
    NATURaLConfig cfg;
    EXPECT_FALSE(cfg.enabled);
    EXPECT_GT(cfg.temperature_K, 200.0);
    EXPECT_LT(cfg.temperature_K, 400.0);
    EXPECT_GT(cfg.mg_concentration_mM, 0.0);
}

TEST(NATURaLConfig, IsNucleotideLigandReturnsFalseForEmpty) {
    EXPECT_FALSE(is_nucleotide_ligand(nullptr, 0));
}

TEST(NATURaLConfig, IsNucleicAcidReceptorReturnsFalseForEmpty) {
    EXPECT_FALSE(is_nucleic_acid_receptor(nullptr, 0));
}

// ===========================================================================
// Human protofibril DualAssembly planning
// ===========================================================================

TEST(ProtofibrilDualAssemblyPlan, BuildsCanonicalAndReciprocalTracks) {
    auto tracks = make_human_protofibril_tracks(90, 45, true);
    ASSERT_EQ(tracks.size(), 4u);

    int canonical = 0;
    int reciprocal = 0;
    int transcription = 0;
    int translation = 0;

    for (const auto& tr : tracks) {
        if (tr.role_policy == DockingRolePolicy::ProtofibrilAsTarget) {
            ++canonical;
            EXPECT_TRUE(tr.protofibril_is_target);
        } else {
            ++reciprocal;
            EXPECT_FALSE(tr.protofibril_is_target);
        }

        if (tr.process == GrowthProcess::Transcription) ++transcription;
        if (tr.process == GrowthProcess::Translation)   ++translation;
    }

    EXPECT_EQ(canonical, 2);
    EXPECT_EQ(reciprocal, 2);
    EXPECT_EQ(transcription, 2);
    EXPECT_EQ(translation, 2);
}

TEST(ProtofibrilDualAssemblyPlan, UsesHumanCellGrowthClocks) {
    auto tracks = make_human_protofibril_tracks(100, 50, false);
    ASSERT_EQ(tracks.size(), 2u);

    const auto& tx = tracks[0];
    const auto& tl = tracks[1];

    ASSERT_EQ(tx.process, GrowthProcess::Transcription);
    ASSERT_EQ(tl.process, GrowthProcess::Translation);

    EXPECT_DOUBLE_EQ(tx.mean_elongation_rate, MEAN_NT_RATE_HUMAN);
    EXPECT_DOUBLE_EQ(tl.mean_elongation_rate, MEAN_EL_RATE_HUMAN);
    EXPECT_DOUBLE_EQ(tx.tunnel_length, RNAP_TUNNEL_NT);
    EXPECT_DOUBLE_EQ(tl.tunnel_length, TUNNEL_LENGTH_AA);

    // Human transcription is nuclear; keep it in the synchronized clock but do
    // not silently claim direct protofibril encounter in the default protocol.
    EXPECT_FALSE(tx.direct_encounter_allowed);
    EXPECT_TRUE(tl.direct_encounter_allowed);
}

TEST(ProtofibrilDualAssemblyPlan, CompletionTimesIncludeInitiationAndDwellTime) {
    auto tracks = make_human_protofibril_tracks(100, 50, false);
    ASSERT_EQ(tracks.size(), 2u);

    const auto& tx = tracks[0];
    const auto& tl = tracks[1];

    EXPECT_NEAR(tx.dwell_time_s(), 1.0 / MEAN_NT_RATE_HUMAN, 1e-12);
    EXPECT_NEAR(tl.dwell_time_s(), 1.0 / MEAN_EL_RATE_HUMAN, 1e-12);

    EXPECT_NEAR(tx.completion_time_s(),
                1.0 / K_RNAP_INI_DEFAULT + 100.0 / MEAN_NT_RATE_HUMAN,
                1e-12);
    EXPECT_NEAR(tl.completion_time_s(),
                1.0 / K_INI_DEFAULT + 50.0 / MEAN_EL_RATE_HUMAN,
                1e-12);

    EXPECT_NEAR(tx.first_exposed_time_s(),
                1.0 / K_RNAP_INI_DEFAULT
                    + (std::floor(RNAP_TUNNEL_NT) + 1.0) / MEAN_NT_RATE_HUMAN,
                1e-12);
    EXPECT_NEAR(tl.first_exposed_time_s(),
                1.0 / K_INI_DEFAULT
                    + (std::floor(TUNNEL_LENGTH_AA) + 1.0) / MEAN_EL_RATE_HUMAN,
                1e-12);
}

TEST(ProtofibrilDualAssemblyPlan, ParallelScheduleIsTimeSortedAndTracksExposure) {
    auto tracks = make_human_protofibril_tracks(12, 36, false);
    auto events = build_parallel_growth_schedule(tracks);
    ASSERT_EQ(events.size(), 48u);

    for (size_t i = 1; i < events.size(); ++i)
        EXPECT_LE(events[i - 1].t_arrival, events[i].t_arrival);

    auto tx_first_exposed = std::find_if(events.begin(), events.end(),
        [](const ParallelGrowthEvent& e) {
            return e.process == GrowthProcess::Transcription && e.exposed_units > 0;
        });
    ASSERT_NE(tx_first_exposed, events.end());
    EXPECT_EQ(tx_first_exposed->unit_index, static_cast<int>(RNAP_TUNNEL_NT) + 1);

    auto tl_first_exposed = std::find_if(events.begin(), events.end(),
        [](const ParallelGrowthEvent& e) {
            return e.process == GrowthProcess::Translation && e.exposed_units > 0;
        });
    ASSERT_NE(tl_first_exposed, events.end());
    EXPECT_EQ(tl_first_exposed->unit_index, static_cast<int>(TUNNEL_LENGTH_AA) + 1);
}

TEST(ProtofibrilDualAssemblyPlan, RejectsNegativeChainLengths) {
    EXPECT_THROW(make_human_protofibril_tracks(-1, 10), std::invalid_argument);
    EXPECT_THROW(make_human_protofibril_tracks(10, -1), std::invalid_argument);
}

TEST(ProtofibrilDualAssemblyPlan, ZeroLengthInputsProduceNoTracksOrEvents) {
    auto tracks = make_human_protofibril_tracks(0, 0, true);
    EXPECT_TRUE(tracks.empty());

    auto events = build_parallel_growth_schedule(tracks);
    EXPECT_TRUE(events.empty());
}

TEST(ProtofibrilDualAssemblyPlan, RejectsNonPositiveOrNonFiniteElongationRates) {
    InVivoAssemblyTrack bad_rate;
    bad_rate.chain_length = 5;
    bad_rate.mean_elongation_rate = 0.0;

    EXPECT_THROW(build_parallel_growth_schedule({bad_rate}), std::invalid_argument);

    bad_rate.mean_elongation_rate = -1.0;
    EXPECT_THROW(build_parallel_growth_schedule({bad_rate}), std::invalid_argument);

    bad_rate.mean_elongation_rate = std::numeric_limits<double>::infinity();
    EXPECT_THROW(build_parallel_growth_schedule({bad_rate}), std::invalid_argument);
}

// ===========================================================================
// MAIN
// ===========================================================================

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
