// test_nascent_chain_scheduler.cpp — GoogleTest for NascentChainScheduler.
//
// Covers: schedule ordering, checkpoint counts, tunnel/chaperone gating, T/L
// discriminator at the three entropy regimes, and FibrilGrowthOracle Ξ math.
//
// Copyright 2026 Le Bonhomme Pharma. SPDX-License-Identifier: Apache-2.0
#include <gtest/gtest.h>

#include "NATURaL/NascentChainScheduler.h"
#include "NATURaL/FibrilGrowthOracle.h"
#include "NATURaL/DualAssemblyRunner.h"
#include "NATURaL/NATURaLDualAssembly.h"
#include "NATURaL/RibosomeElongation.h"
#include "ShannonThermoStack/ShannonThermoStack.h"
#include "statmech.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ─── Scheduling ──────────────────────────────────────────────────────────────
TEST(NascentChainScheduler, EmitsMonotoneCheckpoints) {
    auto tracks = natural::make_human_protofibril_tracks(
        /*transcript_nt=*/300, /*peptide_aa=*/100,
        /*include_reciprocal_controls=*/true);
    natural::NascentChainScheduler sched(tracks, /*interval=*/10);

    EXPECT_GE(sched.schedule().size(), 4u);  // ≥ 2 tracks × ≥ 2 checkpoints

    // Schedule must be monotone in (t_arrival_s, track_index).
    for (size_t i = 1; i < sched.schedule().size(); ++i) {
        const auto& prev = sched.schedule()[i - 1];
        const auto& cur  = sched.schedule()[i];
        EXPECT_LE(prev.t_arrival_s, cur.t_arrival_s + 1e-12);
    }
}

TEST(NascentChainScheduler, FirstCheckpointPastTunnelPlusChaperoneSkip) {
    auto tracks = natural::make_human_protofibril_tracks(
        /*transcript_nt=*/0, /*peptide_aa=*/100,
        /*include_reciprocal_controls=*/false);
    natural::NascentChainScheduler sched(tracks, /*interval=*/10);

    // Tunnel = 34 aa, chaperone skip = 6 → first translation checkpoint at L = 40.
    bool found_translation = false;
    for (const auto& ck : sched.schedule()) {
        if (ck.process == natural::GrowthProcess::Translation) {
            found_translation = true;
            EXPECT_GE(ck.L_k, 40);
            EXPECT_FALSE(ck.in_tunnel);
            EXPECT_FALSE(ck.chaperone_shielded);
            break;
        }
    }
    EXPECT_TRUE(found_translation);
}

TEST(NascentChainScheduler, IteratorExhaustsAfterNCalls) {
    auto tracks = natural::make_human_protofibril_tracks(0, 100, false);
    natural::NascentChainScheduler sched(tracks, 10);
    const size_t n = sched.schedule().size();

    size_t emitted = 0;
    while (sched.has_next()) { sched.next(); ++emitted; }
    EXPECT_EQ(emitted, n);
    EXPECT_FALSE(sched.has_next());
    EXPECT_THROW((void)sched.next(), std::out_of_range);
}

TEST(NascentChainScheduler, RejectsBadInterval) {
    auto tracks = natural::make_human_protofibril_tracks(0, 100, false);
    EXPECT_THROW(
        (natural::NascentChainScheduler(tracks, /*interval=*/0)),
        std::invalid_argument);
    EXPECT_THROW(
        (natural::NascentChainScheduler(tracks, /*interval=*/-1)),
        std::invalid_argument);
}

// ─── T/L discriminator ──────────────────────────────────────────────────────
TEST(TLDiscriminator, HardRegimeAssignsLowerEntropySystem) {
    natural::CheckpointOutcome out;
    // H_A = 0.3 nats (< ln 2 ≈ 0.693) → A wins. H_B = 1.5 (> soft).
    natural::DualAssemblyRunner::assign_tl(/*H_A=*/0.3, /*H_B=*/1.5, '?', out);
    EXPECT_EQ(out.tl_primary, 'A');
    EXPECT_FALSE(out.tl_deferred);
    EXPECT_DOUBLE_EQ(out.tl_weight, 1.0);

    natural::CheckpointOutcome out2;
    natural::DualAssemblyRunner::assign_tl(/*H_A=*/1.5, /*H_B=*/0.3, '?', out2);
    EXPECT_EQ(out2.tl_primary, 'B');
    EXPECT_FALSE(out2.tl_deferred);
    EXPECT_DOUBLE_EQ(out2.tl_weight, 1.0);
}

TEST(TLDiscriminator, DeferredRegimeInheritsPrevAssignment) {
    natural::CheckpointOutcome out;
    // Both above soft threshold (2 ln 2 ≈ 1.386).
    natural::DualAssemblyRunner::assign_tl(/*H_A=*/1.5, /*H_B=*/1.7, /*prev=*/'B', out);
    EXPECT_EQ(out.tl_primary, 'B');
    EXPECT_TRUE(out.tl_deferred);
    EXPECT_DOUBLE_EQ(out.tl_weight, 0.0);
}

TEST(TLDiscriminator, SoftRegimeProducesIntermediateWeight) {
    natural::CheckpointOutcome out;
    // H_lower in (ln 2, 2 ln 2). Pick H_lower = 1.0 → w = (1.386 − 1.0)/(1.386 − 0.693).
    constexpr double H_lower = 1.0;
    natural::DualAssemblyRunner::assign_tl(H_lower, /*H_B=*/2.0, /*prev=*/'?', out);
    EXPECT_EQ(out.tl_primary, 'A');
    EXPECT_FALSE(out.tl_deferred);
    constexpr double soft = shannon_thermo::kHSC_soft_nats;
    constexpr double hard = shannon_thermo::kHSC_hard_nats;
    const double expected_w = (soft - H_lower) / (soft - hard);
    EXPECT_NEAR(out.tl_weight, expected_w, 1e-9);
    EXPECT_GT(out.tl_weight, 0.0);
    EXPECT_LT(out.tl_weight, 1.0);
}

TEST(TLDiscriminator, BothInfiniteInheritsPrev) {
    natural::CheckpointOutcome out;
    const double inf = std::numeric_limits<double>::infinity();
    natural::DualAssemblyRunner::assign_tl(inf, inf, /*prev=*/'A', out);
    EXPECT_EQ(out.tl_primary, 'A');
    EXPECT_TRUE(out.tl_deferred);
    EXPECT_DOUBLE_EQ(out.tl_weight, 0.0);
}

// ─── FibrilGrowthOracle ─────────────────────────────────────────────────────
TEST(FibrilGrowthOracle, AnalyticXiForKnownZ) {
    natural::FibrilGrowthOracle oracle(/*T_K=*/310.15);

    // Inject a synthetic monomer ensemble with a known F (and hence ln Z = -F/kT).
    statmech::StatMechEngine engine(310.15);
    // ΔG = -2 kcal/mol (favourable). Add a single sample at this energy.
    engine.add_sample(-2.0, 1.0);

    // At c = 1 µM, z = 1e-6. With F ≈ -2 kcal/mol and kT ≈ 0.616 kcal/mol at 310 K,
    // ln Z ≈ 2 / 0.616 ≈ 3.247. z·Z ≈ 1e-6 · e^3.247 ≈ 2.57e-5.
    // p ≈ z·Z / (1 + z·Z) ≈ 2.57e-5.
    auto decision = oracle.gate(engine, /*c_monomer_M=*/1.0e-6);
    EXPECT_GT(decision.p_elong, 0.0);
    EXPECT_LT(decision.p_elong, 1.0);
    EXPECT_TRUE(std::isfinite(decision.dG_elong));

    // At c = 1 M (z = 1), p should be much higher.
    auto decision2 = oracle.gate(engine, /*c_monomer_M=*/1.0);
    EXPECT_GT(decision2.p_elong, decision.p_elong);

    const double kT = statmech::kB_kcal * 310.15;
    EXPECT_NEAR(decision.dG_elong - decision2.dG_elong,
                kT * std::log(1.0e6),
                1e-9);
}

TEST(FibrilGrowthOracle, RejectsNonPositiveConcentration) {
    natural::FibrilGrowthOracle oracle(310.15);
    statmech::StatMechEngine engine(310.15);
    engine.add_sample(-1.0, 1.0);
    EXPECT_THROW(oracle.gate(engine, 0.0), std::invalid_argument);
    EXPECT_THROW(oracle.gate(engine, -1e-9), std::invalid_argument);
}

TEST(FibrilGrowthOracle, GateDecisionTracksConcentrationCorrectedFreeEnergy) {
    natural::FibrilGrowthOracle oracle(/*T_K=*/310.15, /*thresh=*/0.5);
    statmech::StatMechEngine engine(310.15);
    engine.add_sample(-20.0, 1.0);

    const auto low_conc = oracle.gate(engine, 1.0e-30);
    const auto high_conc = oracle.gate(engine, 1.0);

    EXPECT_FALSE(low_conc.gated_in);
    EXPECT_TRUE(high_conc.gated_in);
    EXPECT_LT(low_conc.p_elong, high_conc.p_elong);
    EXPECT_GT(low_conc.dG_elong, high_conc.dG_elong);
}

TEST(FibrilGrowthOracle, RejectsBadTemperature) {
    EXPECT_THROW((natural::FibrilGrowthOracle(0.0)), std::invalid_argument);
    EXPECT_THROW((natural::FibrilGrowthOracle(310.15, /*thresh=*/-0.1)), std::invalid_argument);
    EXPECT_THROW((natural::FibrilGrowthOracle(310.15, /*thresh=*/1.5)), std::invalid_argument);
}

// ─── DualAssemblyRunner end-to-end with synthetic backends ──────────────────
namespace {
natural::GAResult synthetic_engine(double T, int n_poses, double width, double dG, unsigned seed) {
    statmech::StatMechEngine eng(T);
    std::mt19937 rng(seed);
    std::normal_distribution<double> rd(0.0, std::max(1e-3, width));
    std::normal_distribution<double> ed(dG, 1.0);
    std::vector<double> rmsds;
    rmsds.reserve(n_poses);
    for (int i = 0; i < n_poses; ++i) {
        eng.add_sample(ed(rng), 1.0);
        rmsds.push_back(std::abs(rd(rng)));
    }
    return {std::move(eng), std::move(rmsds)};
}

std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string field;
    while (std::getline(ss, field, ','))
        fields.push_back(field);
    return fields;
}
} // namespace

TEST(DualAssemblyRunner, RunsEndToEndWithSyntheticBackend) {
    natural::DualAssemblyConfig cfg;
    cfg.protofibril_pdb = "fake.pdb";
    cfg.monomer_pdb     = "fake_monomer.pdb";
    cfg.sequence_fasta  = std::string(80, 'A');  // 80 aa
    cfg.checkpoint_interval = 10;
    cfg.sim_c_interval = 3;
    cfg.include_reciprocal_controls = false;     // keep small
    cfg.sim_c_enabled = true;
    cfg.output_csv = "/tmp/dual_assembly_smoke.csv";
    cfg.nascent_pdb_dir = "/tmp/dual_assembly_smoke_pdbs";

    auto sim_a = [](const std::string&, const std::string&, int L_k, double T) {
        // Narrow as L_k grows so H crosses thresholds.
        return synthetic_engine(T, 64, std::max(0.5, 8.0 - 0.1 * L_k), -1.0 - 0.05 * L_k, 1u + L_k);
    };
    auto sim_b = [](const std::string&, const std::string&, int L_k, double T) {
        return synthetic_engine(T, 32, std::max(0.5, 4.0 - 0.02 * L_k), -0.5 - 0.02 * L_k, 2u + L_k);
    };
    auto sim_c = [](const std::string&, const std::string&, double T) {
        return synthetic_engine(T, 32, 1.5, -3.0, 3u);
    };
    auto trunc = [](const std::string& seq, int L_k, const std::string&) {
        return std::string("/tmp/nascent_L") + std::to_string(L_k) + ".pdb";
    };

    natural::DualAssemblyRunner runner(std::move(cfg), sim_a, sim_b, sim_c, trunc);
    auto history = runner.run();
    EXPECT_GT(history.size(), 0u);
    // First eligible checkpoint must be at L_k >= 40 (tunnel + chaperone skip).
    bool found_eligible = false;
    for (const auto& [ck, out] : history) {
        if (!ck.in_tunnel && !ck.chaperone_shielded && ck.direct_encounter_allowed) {
            EXPECT_GE(ck.L_k, 40);
            EXPECT_TRUE(std::isfinite(out.H_A_nats));
            EXPECT_TRUE(std::isfinite(out.dG_A_kcal));
            found_eligible = true;
            break;
        }
    }
    EXPECT_TRUE(found_eligible);

    std::ifstream csv("/tmp/dual_assembly_smoke.csv");
    ASSERT_TRUE(csv.good());
    std::string header;
    ASSERT_TRUE(std::getline(csv, header));
    EXPECT_NE(header.find("tl_basis"), std::string::npos);
    const std::vector<std::string> expected_header = {
        "checkpoint_idx",
        "L_k",
        "process",
        "track_name",
        "role_policy",
        "t_arrival_s",
        "in_tunnel",
        "H_A_nats",
        "H_B_nats",
        "dG_A_kcal",
        "dG_B_kcal",
        "tl_primary",
        "tl_weight",
        "tl_deferred",
        "tl_basis",
        "p_elong",
        "dG_elong_kcal",
        "sim_c_evaluated",
        "sim_c_gated_in",
        "protofibril_state_idx",
        "protofibril_structure_updated"
    };
    EXPECT_EQ(split_csv_line(header), expected_header);

    std::string first_row;
    ASSERT_TRUE(std::getline(csv, first_row));
    EXPECT_NE(first_row.find("pose_entropy_heuristic"), std::string::npos);
    EXPECT_EQ(split_csv_line(first_row).size(), expected_header.size());
}

TEST(DualAssemblyRunner, ReciprocalControlSwapsPrimaryTargetLigandOrder) {
    natural::DualAssemblyConfig cfg;
    cfg.protofibril_pdb = "proto.pdb";
    cfg.sequence_fasta  = std::string(45, 'A');
    cfg.checkpoint_interval = 100;
    cfg.include_reciprocal_controls = true;
    cfg.sim_c_enabled = false;
    cfg.output_csv = "/tmp/dual_assembly_reciprocal.csv";
    cfg.nascent_pdb_dir = "/tmp/dual_assembly_reciprocal_pdbs";

    std::vector<std::pair<std::string, std::string>> sim_a_calls;
    auto sim_a = [&](const std::string& target,
                     const std::string& ligand,
                     int L_k,
                     double T) {
        sim_a_calls.emplace_back(target, ligand);
        return synthetic_engine(T, 16, 1.0, -1.0, 10u + L_k);
    };
    auto trunc = [](const std::string&, int L_k, const std::string&) {
        return std::string("/tmp/dual_assembly_reciprocal_L") + std::to_string(L_k) + ".pdb";
    };

    natural::DualAssemblyRunner runner(std::move(cfg), sim_a, nullptr, nullptr, trunc);
    auto history = runner.run();

    ASSERT_GT(history.size(), 0u);
    ASSERT_EQ(sim_a_calls.size(), 2u);
    EXPECT_EQ(sim_a_calls[0].first, "proto.pdb");
    EXPECT_EQ(sim_a_calls[0].second, "/tmp/dual_assembly_reciprocal_L40.pdb");
    EXPECT_EQ(sim_a_calls[1].first, "/tmp/dual_assembly_reciprocal_L40.pdb");
    EXPECT_EQ(sim_a_calls[1].second, "proto.pdb");

    std::vector<natural::DockingRolePolicy> translation_roles;
    for (const auto& [ck, out] : history) {
        (void)out;
        if (ck.process == natural::GrowthProcess::Translation)
            translation_roles.push_back(ck.role_policy);
    }
    ASSERT_EQ(translation_roles.size(), 2u);
    EXPECT_EQ(translation_roles[0], natural::DockingRolePolicy::ProtofibrilAsTarget);
    EXPECT_EQ(translation_roles[1], natural::DockingRolePolicy::ReciprocalControl);
}

TEST(DualAssemblyRunner, SimBGateUsesPriorEntropyPerTrack) {
    natural::DualAssemblyConfig cfg;
    cfg.protofibril_pdb = "proto.pdb";
    cfg.monomer_pdb = "monomer.pdb";
    cfg.sequence_fasta = std::string(60, 'A');
    cfg.checkpoint_interval = 20;
    cfg.include_reciprocal_controls = true;
    cfg.sim_c_enabled = false;
    cfg.output_csv = "/tmp/dual_assembly_track_state.csv";
    cfg.nascent_pdb_dir = "/tmp/dual_assembly_track_state_pdbs";

    auto sim_a = [](const std::string& target,
                    const std::string&,
                    int L_k,
                    double T) {
        statmech::StatMechEngine eng(T);
        eng.add_sample(-1.0, 1.0);
        std::vector<double> coords;
        if (target == "proto.pdb" && L_k == 40) {
            coords = {0.0, 0.0, 0.0, 0.0}; // canonical track collapses
        } else {
            coords = {0.0, 1.0, 2.0, 3.0, 4.0}; // high entropy, above soft threshold
        }
        return natural::GAResult{std::move(eng), std::move(coords)};
    };

    std::vector<std::string> sim_b_targets;
    auto sim_b = [&](const std::string& target,
                     const std::string&,
                     int,
                     double T) {
        sim_b_targets.push_back(target);
        statmech::StatMechEngine eng(T);
        eng.add_sample(-0.5, 1.0);
        return natural::GAResult{std::move(eng), std::vector<double>{0.0, 0.0}};
    };

    auto trunc = [](const std::string&, int L_k, const std::string&) {
        return std::string("/tmp/dual_assembly_track_state_L") + std::to_string(L_k) + ".pdb";
    };

    natural::DualAssemblyRunner runner(std::move(cfg), sim_a, sim_b, nullptr, trunc);
    (void)runner.run();

    ASSERT_EQ(sim_b_targets.size(), 1u);
    EXPECT_EQ(sim_b_targets[0], "/tmp/dual_assembly_track_state_L60.pdb");
}

TEST(DualAssemblyRunner, AcceptedSimCAdvanceChangesSubsequentProtofibrilTarget) {
    natural::DualAssemblyConfig cfg;
    cfg.protofibril_pdb = "proto_state_0.pdb";
    cfg.monomer_pdb = "monomer.pdb";
    cfg.sequence_fasta = std::string(60, 'A');
    cfg.checkpoint_interval = 20;
    cfg.include_reciprocal_controls = false;
    cfg.sim_c_enabled = true;
    cfg.sim_c_interval = 1;
    cfg.monomer_conc_M = 1.0;
    cfg.output_csv = "/tmp/dual_assembly_advance.csv";
    cfg.nascent_pdb_dir = "/tmp/dual_assembly_advance_pdbs";

    std::vector<std::string> sim_a_targets;
    auto sim_a = [&](const std::string& target,
                     const std::string&,
                     int L_k,
                     double T) {
        sim_a_targets.push_back(target);
        return synthetic_engine(T, 16, 1.0, -1.0, 100u + L_k);
    };

    auto sim_c = [](const std::string&, const std::string&, double T) {
        statmech::StatMechEngine eng(T);
        eng.add_sample(-20.0, 1.0); // overwhelming acceptance at 1 M
        return natural::GAResult{std::move(eng), std::vector<double>{0.0, 0.0}};
    };

    auto trunc = [](const std::string&, int L_k, const std::string&) {
        return std::string("/tmp/dual_assembly_advance_L") + std::to_string(L_k) + ".pdb";
    };

    auto advance = [](const std::string&,
                      const std::string&,
                      int next_state_index) {
        return std::string("proto_state_") + std::to_string(next_state_index) + ".pdb";
    };

    natural::DualAssemblyRunner runner(std::move(cfg), sim_a, nullptr, sim_c, trunc, advance);
    auto history = runner.run();

    ASSERT_EQ(sim_a_targets.size(), 2u);
    EXPECT_EQ(sim_a_targets[0], "proto_state_0.pdb");
    EXPECT_EQ(sim_a_targets[1], "proto_state_1.pdb");

    auto first_translation = std::find_if(history.begin(), history.end(),
        [](const auto& item) {
            return item.first.process == natural::GrowthProcess::Translation;
        });
    ASSERT_NE(first_translation, history.end());
    EXPECT_TRUE(first_translation->second.sim_c_evaluated);
    EXPECT_TRUE(first_translation->second.sim_c_gated_in);
    EXPECT_EQ(first_translation->second.protofibril_state_index, 1);
    EXPECT_TRUE(first_translation->second.protofibril_structure_updated);
}

TEST(DualAssemblyRunner, RejectedSimCDoesNotAdvanceProtofibrilState) {
    natural::DualAssemblyConfig cfg;
    cfg.protofibril_pdb = "proto_state_0.pdb";
    cfg.monomer_pdb = "monomer.pdb";
    cfg.sequence_fasta = std::string(60, 'A');
    cfg.checkpoint_interval = 20;
    cfg.include_reciprocal_controls = false;
    cfg.sim_c_enabled = true;
    cfg.sim_c_interval = 1;
    cfg.monomer_conc_M = 1.0e-12;
    cfg.output_csv = "/tmp/dual_assembly_reject.csv";
    cfg.nascent_pdb_dir = "/tmp/dual_assembly_reject_pdbs";

    std::vector<std::string> sim_a_targets;
    auto sim_a = [&](const std::string& target,
                     const std::string&,
                     int L_k,
                     double T) {
        sim_a_targets.push_back(target);
        return synthetic_engine(T, 16, 1.0, -1.0, 200u + L_k);
    };

    auto sim_c = [](const std::string&, const std::string&, double T) {
        statmech::StatMechEngine eng(T);
        eng.add_sample(20.0, 1.0);
        return natural::GAResult{std::move(eng), std::vector<double>{0.0, 0.0}};
    };

    auto trunc = [](const std::string&, int L_k, const std::string&) {
        return std::string("/tmp/dual_assembly_reject_L") + std::to_string(L_k) + ".pdb";
    };

    natural::DualAssemblyRunner runner(std::move(cfg), sim_a, nullptr, sim_c, trunc);
    auto history = runner.run();

    ASSERT_EQ(sim_a_targets.size(), 2u);
    EXPECT_EQ(sim_a_targets[0], "proto_state_0.pdb");
    EXPECT_EQ(sim_a_targets[1], "proto_state_0.pdb");

    auto first_translation = std::find_if(history.begin(), history.end(),
        [](const auto& item) {
            return item.first.process == natural::GrowthProcess::Translation;
        });
    ASSERT_NE(first_translation, history.end());
    EXPECT_TRUE(first_translation->second.sim_c_evaluated);
    EXPECT_FALSE(first_translation->second.sim_c_gated_in);
    EXPECT_EQ(first_translation->second.protofibril_state_index, 0);
    EXPECT_FALSE(first_translation->second.protofibril_structure_updated);
}

TEST(DualAssemblyRunner, AcceptedSimCWithoutAdvanceCallbackDoesNotClaimStructureUpdate) {
    natural::DualAssemblyConfig cfg;
    cfg.protofibril_pdb = "proto_state_0.pdb";
    cfg.monomer_pdb = "monomer.pdb";
    cfg.sequence_fasta = std::string(60, 'A');
    cfg.checkpoint_interval = 20;
    cfg.include_reciprocal_controls = false;
    cfg.sim_c_enabled = true;
    cfg.sim_c_interval = 1;
    cfg.monomer_conc_M = 1.0;
    cfg.output_csv = "/tmp/dual_assembly_accept_no_update.csv";
    cfg.nascent_pdb_dir = "/tmp/dual_assembly_accept_no_update_pdbs";

    std::vector<std::string> sim_a_targets;
    auto sim_a = [&](const std::string& target,
                     const std::string&,
                     int L_k,
                     double T) {
        sim_a_targets.push_back(target);
        return synthetic_engine(T, 16, 1.0, -1.0, 300u + L_k);
    };

    auto sim_c = [](const std::string&, const std::string&, double T) {
        statmech::StatMechEngine eng(T);
        eng.add_sample(-20.0, 1.0);
        return natural::GAResult{std::move(eng), std::vector<double>{0.0, 0.0}};
    };

    auto trunc = [](const std::string&, int L_k, const std::string&) {
        return std::string("/tmp/dual_assembly_accept_no_update_L") + std::to_string(L_k) + ".pdb";
    };

    natural::DualAssemblyRunner runner(std::move(cfg), sim_a, nullptr, sim_c, trunc);
    auto history = runner.run();

    ASSERT_EQ(sim_a_targets.size(), 2u);
    EXPECT_EQ(sim_a_targets[0], "proto_state_0.pdb");
    EXPECT_EQ(sim_a_targets[1], "proto_state_0.pdb");

    auto first_translation = std::find_if(history.begin(), history.end(),
        [](const auto& item) {
            return item.first.process == natural::GrowthProcess::Translation;
        });
    ASSERT_NE(first_translation, history.end());
    EXPECT_TRUE(first_translation->second.sim_c_evaluated);
    EXPECT_TRUE(first_translation->second.sim_c_gated_in);
    EXPECT_EQ(first_translation->second.protofibril_state_index, 1);
    EXPECT_FALSE(first_translation->second.protofibril_structure_updated);
}

TEST(DualAssemblyRunner, SimCCadenceIsPerTrackNotGlobal) {
    natural::DualAssemblyConfig cfg;
    cfg.protofibril_pdb = "proto.pdb";
    cfg.monomer_pdb = "monomer.pdb";
    cfg.sequence_fasta = std::string(60, 'A');
    cfg.checkpoint_interval = 20;
    cfg.include_reciprocal_controls = true;
    cfg.sim_c_enabled = true;
    cfg.sim_c_interval = 2;
    cfg.monomer_conc_M = 1.0;
    cfg.output_csv = "/tmp/dual_assembly_sim_c_cadence.csv";
    cfg.nascent_pdb_dir = "/tmp/dual_assembly_sim_c_cadence_pdbs";

    auto sim_a = [](const std::string&,
                    const std::string&,
                    int L_k,
                    double T) {
        return synthetic_engine(T, 16, 1.0, -1.0, 400u + L_k);
    };

    int sim_c_calls = 0;
    auto sim_c = [&](const std::string&, const std::string&, double T) {
        ++sim_c_calls;
        statmech::StatMechEngine eng(T);
        eng.add_sample(-20.0, 1.0);
        return natural::GAResult{std::move(eng), std::vector<double>{0.0, 0.0}};
    };

    auto trunc = [](const std::string&, int L_k, const std::string&) {
        return std::string("/tmp/dual_assembly_sim_c_cadence_L") + std::to_string(L_k) + ".pdb";
    };

    natural::DualAssemblyRunner runner(std::move(cfg), sim_a, nullptr, sim_c, trunc);
    auto history = runner.run();

    std::vector<bool> translation_sim_c_flags;
    for (const auto& [ck, out] : history) {
        if (ck.process == natural::GrowthProcess::Translation)
            translation_sim_c_flags.push_back(out.sim_c_evaluated);
    }

    ASSERT_EQ(translation_sim_c_flags.size(), 4u);
    EXPECT_FALSE(translation_sim_c_flags[0]);
    EXPECT_FALSE(translation_sim_c_flags[1]);
    EXPECT_TRUE(translation_sim_c_flags[2]);
    EXPECT_TRUE(translation_sim_c_flags[3]);
    EXPECT_EQ(sim_c_calls, 2);
}
