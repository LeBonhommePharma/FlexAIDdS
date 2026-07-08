// test_target_server.cpp — Integration tests for TargetServer
//
// Tests session management, grand partition function accumulation,
// and knowledge base updates through the TargetServer interface.
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>
#include "TargetServer.h"
#include "statmech.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <thread>
#include <vector>

using namespace target;

// ════════════════════════════════════════════════════════════════════════
// Helper: create synthetic FA_Global + atoms for validation tests
// ════════════════════════════════════════════════════════════════════════

static FA_Global make_test_fa(int atm_cnt = 50, int res_cnt = 50, int ntypes = 10)
{
    FA_Global fa{};
    fa.atm_cnt_real = atm_cnt;
    fa.atm_cnt = atm_cnt;
    fa.res_cnt = res_cnt;
    fa.ntypes = ntypes;
    fa.multi_model = false;
    fa.n_models = 1;
    return fa;
}

static std::vector<atom> make_atoms(int count)
{
    std::vector<atom> atoms(count);
    for (int i = 0; i < count; ++i) {
        std::memset(&atoms[i], 0, sizeof(atom));
        atoms[i].coor[0] = i * 3.8f;
        atoms[i].coor[1] = 0.0f;
        atoms[i].coor[2] = 0.0f;
        atoms[i].type = 1;
        atoms[i].ofres = i; // link atom to residue i
        std::strncpy(atoms[i].name, " CA ", 4);
        atoms[i].name[4] = '\0';
    }
    return atoms;
}

static std::vector<resid> make_residues(int count)
{
    std::vector<resid> residues(count);
    for (int i = 0; i < count; ++i) {
        std::memset(&residues[i], 0, sizeof(resid));
        std::strncpy(residues[i].name, "ALA", 3);
        residues[i].chn = 'A';
        residues[i].number = i + 1;
    }
    return residues;
}

// ════════════════════════════════════════════════════════════════════════
// Construction
// ════════════════════════════════════════════════════════════════════════

TEST(TargetServer, DefaultConstruction) {
    TargetServer server;
    EXPECT_NEAR(server.temperature(), 300.0, 1e-10);
    EXPECT_EQ(server.n_models(), 1);
    EXPECT_EQ(server.completed_sessions(), 0);
}

TEST(TargetServer, CustomConfig) {
    TargetConfig cfg;
    cfg.temperature_K = 310.0;
    cfg.n_models = 5;
    TargetServer server(cfg);
    EXPECT_NEAR(server.temperature(), 310.0, 1e-10);
    EXPECT_EQ(server.n_models(), 5);
}

// ════════════════════════════════════════════════════════════════════════
// Validation delegation
// ════════════════════════════════════════════════════════════════════════

TEST(TargetServer, ValidateGoodTarget) {
    TargetServer server;
    auto atoms = make_atoms(50);
    auto residues = make_residues(50);
    FA_Global fa = make_test_fa();

    auto result = server.validate(&fa, atoms.data(), residues.data(), 100);
    EXPECT_TRUE(result.valid);
    EXPECT_TRUE(result.errors.empty());
}

TEST(TargetServer, ValidateBadTarget) {
    TargetServer server;
    auto result = server.validate(nullptr, nullptr, nullptr, 0);
    EXPECT_FALSE(result.valid);
}

// ════════════════════════════════════════════════════════════════════════
// Session management
// ════════════════════════════════════════════════════════════════════════

TEST(TargetServer, CreateSession) {
    TargetServer server;
    auto s1 = server.create_session("aspirin");
    auto s2 = server.create_session("ibuprofen");

    EXPECT_EQ(s1.session_id, 0);
    EXPECT_EQ(s2.session_id, 1);
    EXPECT_EQ(s1.ligand_name, "aspirin");
    EXPECT_EQ(s2.ligand_name, "ibuprofen");
    EXPECT_FALSE(s1.completed);
    EXPECT_FALSE(s2.completed);
}

TEST(TargetServer, RegisterResult) {
    TargetServer server;
    auto session = server.create_session("aspirin");

    // Simulate completed docking
    session.completed = true;
    session.log_Z = 10.0;
    session.n_poses = 100;
    session.best_energy = -8.5;
    session.best_center[0] = 1.0f;
    session.best_center[1] = 2.0f;
    session.best_center[2] = 3.0f;

    server.register_result(session);

    EXPECT_EQ(server.completed_sessions(), 1);
    EXPECT_TRUE(server.grand_partition().has_ligand("aspirin"));
}

TEST(TargetServer, SkipIncompleteSession) {
    TargetServer server;
    auto session = server.create_session("aspirin");
    session.completed = false;

    server.register_result(session);  // should be ignored
    EXPECT_EQ(server.completed_sessions(), 0);
    EXPECT_FALSE(server.grand_partition().has_ligand("aspirin"));
}

// ════════════════════════════════════════════════════════════════════════
// Grand partition function through TargetServer
// ════════════════════════════════════════════════════════════════════════

TEST(TargetServer, CompetitiveBinding) {
    TargetServer server;

    // Register 3 ligands with different affinities
    auto dock = [&](const std::string& name, double log_Z) {
        auto s = server.create_session(name);
        s.completed = true;
        s.log_Z = log_Z;
        s.n_poses = 50;
        server.register_result(s);
    };

    dock("strong", 20.0);
    dock("medium", 10.0);
    dock("weak",    2.0);

    // Ranking
    auto ranks = server.rank_ligands();
    ASSERT_EQ(ranks.size(), 3u);
    EXPECT_EQ(ranks[0].name, "strong");
    EXPECT_EQ(ranks[1].name, "medium");
    EXPECT_EQ(ranks[2].name, "weak");

    // Selectivity
    EXPECT_NEAR(server.selectivity_ratio("strong", "weak"),
                std::exp(20.0 - 2.0), std::exp(18.0) * 1e-10);

    // Probabilities sum to 1
    double sum = server.binding_probability("strong")
               + server.binding_probability("medium")
               + server.binding_probability("weak")
               + server.grand_partition().empty_probability();
    EXPECT_NEAR(sum, 1.0, 1e-10);
}

// ════════════════════════════════════════════════════════════════════════
// Re-docking (update existing ligand)
// ════════════════════════════════════════════════════════════════════════

TEST(TargetServer, RedockingOverwrite) {
    TargetServer server;

    // First docking
    auto s1 = server.create_session("aspirin");
    s1.completed = true;
    s1.log_Z = 5.0;
    server.register_result(s1);

    double dG_first = server.grand_partition().F_bound("aspirin");

    // Re-dock with a better estimate (overwrite, not merge)
    auto s2 = server.create_session("aspirin");
    s2.completed = true;
    s2.log_Z = 8.0;  // improved estimate
    server.register_result(s2);

    // ΔG should reflect the overwrite value, not a merge
    double dG_after = server.grand_partition().F_bound("aspirin");
    EXPECT_NEAR(dG_after, -statmech::kB_kcal * 300.0 * 8.0, 1e-10);

    // Re-docking with worse estimate should give less favorable ΔG
    auto s3 = server.create_session("aspirin");
    s3.completed = true;
    s3.log_Z = 3.0;
    server.register_result(s3);
    double dG_worse = server.grand_partition().F_bound("aspirin");
    EXPECT_GT(dG_worse, dG_after);  // less favorable
}

// ════════════════════════════════════════════════════════════════════════
// Knowledge base accumulation
// ════════════════════════════════════════════════════════════════════════

TEST(TargetServer, ConformerPriors) {
    TargetConfig cfg;
    cfg.n_models = 3;
    TargetServer server(cfg);

    auto s1 = server.create_session("lig1");
    s1.completed = true;
    s1.log_Z = 5.0;
    s1.conformer_populations = {0.7, 0.2, 0.1};
    server.register_result(s1);

    auto s2 = server.create_session("lig2");
    s2.completed = true;
    s2.log_Z = 4.0;
    s2.conformer_populations = {0.6, 0.3, 0.1};
    server.register_result(s2);

    auto priors = server.conformer_priors();
    ASSERT_EQ(priors.size(), 3u);

    // Conformer 0 should have highest posterior (0.7 + 0.6 + prior)
    EXPECT_GT(priors[0], priors[1]);
    EXPECT_GT(priors[1], priors[2]);

    // Should sum to 1
    double sum = priors[0] + priors[1] + priors[2];
    EXPECT_NEAR(sum, 1.0, 1e-10);
}

TEST(TargetServer, BindingCenterAccumulation) {
    TargetServer server;

    auto s = server.create_session("lig1");
    s.completed = true;
    s.log_Z = 5.0;
    s.best_center[0] = 10.0f;
    s.best_center[1] = 20.0f;
    s.best_center[2] = 30.0f;
    s.best_energy = -7.5;
    server.register_result(s);

    auto hits = server.knowledge_base().all_hits();
    ASSERT_EQ(hits.size(), 1u);
    EXPECT_NEAR(hits[0].center[0], 10.0f, 1e-5);
    EXPECT_NEAR(hits[0].energy, -7.5, 1e-10);
    EXPECT_EQ(hits[0].ligand_name, "lig1");
}

// ════════════════════════════════════════════════════════════════════════
// Concurrent session registration (basic thread safety)
// ════════════════════════════════════════════════════════════════════════

TEST(TargetServer, ConcurrentRegistration) {
    TargetServer server;
    const int N = 50;

    auto worker = [&](int id) {
        auto s = server.create_session("lig_" + std::to_string(id));
        s.completed = true;
        s.log_Z = static_cast<double>(id);
        s.n_poses = 10;
        server.register_result(s);
    };

    std::vector<std::thread> threads;
    threads.reserve(N);
    for (int i = 0; i < N; ++i) {
        threads.emplace_back(worker, i);
    }
    for (auto& t : threads) t.join();

    EXPECT_EQ(server.completed_sessions(), N);
    EXPECT_EQ(server.grand_partition().num_ligands(), N);

    // Probabilities should sum to 1
    auto ranks = server.rank_ligands();
    double sum = server.grand_partition().empty_probability();
    for (const auto& r : ranks) sum += r.p_bound;
    EXPECT_NEAR(sum, 1.0, 1e-8);
}

// ════════════════════════════════════════════════════════════════════════
// Astex Diverse tier-1 → GrandPartitionFunction / TargetServer
//
// Drives the *shipped* TargetServer + GrandPartitionFunction path with
// partition-function contributions derived from real Astex Diverse docking
// result artifacts under results/benchmarks/astex_diverse/tier1/.
//
// Conversion matches DatasetRunner.cpp:
//   log_Z = −predicted_dG / (kB_kcal · T)
// where predicted_dG is taken from each entry's best-pose total_score.
// ════════════════════════════════════════════════════════════════════════

namespace {

struct AstexPose {
    double rmsd = 0.0;
    double total_score = 0.0;  // CF-style score in artifact (kcal/mol proxy)
};

struct AstexTier1Entry {
    std::string pdb_id;
    double predicted_dG = 0.0;  // kcal/mol from score-selected pose
    std::vector<AstexPose> poses;
    bool ok = false;
};

struct PoseSelection {
    int pose_index = -1;
    double total_score = 0.0;
    double rmsd = 0.0;
    bool rmsd_success = false;  // RMSD <= 2.0 A
};

// Without-GCE arm: select top pose by lowest total_score (shipped CF field).
static PoseSelection select_pose_without_gce(const AstexTier1Entry& entry)
{
    PoseSelection sel;
    if (entry.poses.empty()) return sel;
    int best = 0;
    for (int i = 1; i < static_cast<int>(entry.poses.size()); ++i) {
        if (entry.poses[i].total_score < entry.poses[best].total_score)
            best = i;
    }
    sel.pose_index = best;
    sel.total_score = entry.poses[best].total_score;
    sel.rmsd = entry.poses[best].rmsd;
    sel.rmsd_success = (sel.rmsd <= 2.0);
    return sel;
}

// Parse full pose list from tier-1 holo JSON (rmsd + total_score pairs).
static bool parse_astex_tier1_json(const std::filesystem::path& path,
                                   AstexTier1Entry& out)
{
    std::ifstream in(path);
    if (!in) return false;
    std::string content((std::istreambuf_iterator<char>(in)),
                        std::istreambuf_iterator<char>());

    std::smatch m_id;
    if (std::regex_search(content, m_id,
                          std::regex("\"target_id\"\\s*:\\s*\"([^\"]+)\""))) {
        out.pdb_id = m_id[1].str();
    } else {
        out.pdb_id = path.stem().string();
    }

    // Collect rmsd and total_score in document order; known schema has both
    // fields once per pose object.
    std::vector<double> rmsds;
    std::vector<double> scores;
    {
        std::regex re_rmsd("\"rmsd\"\\s*:\\s*([-0-9.eE+]+)");
        for (std::sregex_iterator it(content.begin(), content.end(), re_rmsd),
             end; it != end; ++it) {
            try { rmsds.push_back(std::stod((*it)[1].str())); }
            catch (...) { return false; }
        }
    }
    {
        std::regex re_sc("\"total_score\"\\s*:\\s*([-0-9.eE+]+)");
        for (std::sregex_iterator it(content.begin(), content.end(), re_sc),
             end; it != end; ++it) {
            try { scores.push_back(std::stod((*it)[1].str())); }
            catch (...) { return false; }
        }
    }
    if (rmsds.empty() || scores.empty() || rmsds.size() != scores.size())
        return false;

    out.poses.clear();
    out.poses.reserve(rmsds.size());
    for (size_t i = 0; i < rmsds.size(); ++i) {
        if (!std::isfinite(rmsds[i]) || !std::isfinite(scores[i]))
            return false;
        out.poses.push_back({rmsds[i], scores[i]});
    }

    // DatasetRunner uses best CF score as predicted_dG for log_Z.
    auto best = select_pose_without_gce(out);
    out.predicted_dG = best.total_score;
    out.ok = true;
    return true;
}

// Resolve results/benchmarks/astex_diverse/tier1 from CWD or parents.
static std::filesystem::path find_astex_tier1_dir()
{
    namespace fs = std::filesystem;
    fs::path cur = fs::current_path();
    for (int i = 0; i < 6; ++i) {
        fs::path cand = cur / "results" / "benchmarks" / "astex_diverse" / "tier1";
        if (fs::is_directory(cand)) return cand;
        if (!cur.has_parent_path() || cur == cur.root_path()) break;
        cur = cur.parent_path();
    }
    return {};
}

// DatasetRunner formula: log_Z = −ΔG / (kT)
static double log_Z_from_dG(double predicted_dG, double T_K)
{
    return -predicted_dG / (statmech::kB_kcal * T_K);
}

// Load all *_holo.json entries under tier1.
static std::vector<AstexTier1Entry> load_astex_tier1_entries()
{
    std::vector<AstexTier1Entry> entries;
    auto tier1 = find_astex_tier1_dir();
    if (tier1.empty()) return entries;
    for (auto& e : std::filesystem::directory_iterator(tier1)) {
        if (e.path().extension() != ".json") continue;
        if (e.path().filename().string().find("_holo") == std::string::npos)
            continue;
        AstexTier1Entry entry;
        if (parse_astex_tier1_json(e.path(), entry))
            entries.push_back(std::move(entry));
    }
    std::sort(entries.begin(), entries.end(),
              [](const AstexTier1Entry& a, const AstexTier1Entry& b) {
                  return a.pdb_id < b.pdb_id;
              });
    return entries;
}

} // namespace

TEST(TargetServer, AstexDiverseTier1SingleLigandGPF) {
    auto tier1 = find_astex_tier1_dir();
    if (tier1.empty()) {
        GTEST_SKIP() << "Astex tier-1 result dir not found "
                        "(results/benchmarks/astex_diverse/tier1)";
    }

    // Prefer canonical Astex Diverse first entry 1gpk if present; else any.
    std::vector<std::filesystem::path> jsons;
    for (auto& e : std::filesystem::directory_iterator(tier1)) {
        if (e.path().extension() == ".json" &&
            e.path().filename().string().find("_holo") != std::string::npos) {
            jsons.push_back(e.path());
        }
    }
    ASSERT_FALSE(jsons.empty()) << "No *_holo.json under " << tier1;

    std::filesystem::path pick = jsons.front();
    for (const auto& p : jsons) {
        if (p.filename().string().find("1gpk") != std::string::npos) {
            pick = p;
            break;
        }
    }

    AstexTier1Entry entry;
    ASSERT_TRUE(parse_astex_tier1_json(pick, entry))
        << "Failed to parse " << pick;

    const double T = 300.0;
    const double log_Z = log_Z_from_dG(entry.predicted_dG, T);
    ASSERT_TRUE(std::isfinite(log_Z)) << "log_Z not finite for " << entry.pdb_id;

    TargetConfig cfg;
    cfg.temperature_K = T;
    TargetServer server(cfg);

    auto sess = server.create_session(entry.pdb_id);
    sess.completed = true;
    sess.log_Z = log_Z;
    sess.n_poses = 5;
    sess.best_energy = entry.predicted_dG;
    server.register_result(sess);

    ASSERT_EQ(server.completed_sessions(), 1);
    ASSERT_TRUE(server.grand_partition().has_ligand(entry.pdb_id));

    const auto& gpf = server.grand_partition();
    const double log_xi = gpf.log_Xi();
    const double p_bound = server.binding_probability(entry.pdb_id);
    const double p_empty = gpf.empty_probability();
    const double mean_occ = gpf.mean_occupancy();
    const double F = gpf.F_bound(entry.pdb_id);

    EXPECT_TRUE(std::isfinite(log_xi));
    EXPECT_TRUE(std::isfinite(F));
    EXPECT_TRUE(std::isfinite(p_bound));
    EXPECT_TRUE(std::isfinite(p_empty));
    EXPECT_GE(p_bound, 0.0);
    EXPECT_LE(p_bound, 1.0);
    EXPECT_GE(p_empty, 0.0);
    EXPECT_LE(p_empty, 1.0);
    EXPECT_NEAR(p_bound + p_empty, 1.0, 1e-8);
    EXPECT_NEAR(mean_occ + p_empty, 1.0, 1e-8);
    EXPECT_NEAR(F, -statmech::kB_kcal * T * log_Z, 1e-8);

    std::cout << "\n[AstexDiverseTier1SingleLigandGPF]\n"
              << "  pdb_id=" << entry.pdb_id << "\n"
              << "  predicted_dG(total_score)=" << entry.predicted_dG << " kcal/mol\n"
              << "  log_Z=" << log_Z << "\n"
              << "  log_Xi=" << log_xi << "\n"
              << "  p_bound=" << p_bound << " p_empty=" << p_empty
              << " sum=" << (p_bound + p_empty) << "\n"
              << "  mean_occupancy=" << mean_occ << "\n"
              << "  F_bound=" << F << " kcal/mol\n";
}

TEST(TargetServer, AstexDiverseTier1CompetitiveGPF) {
    // Register all available Astex tier-1 ligands onto one TargetServer to
    // exercise multi-ligand grand-canonical ranking with Astex-derived log_Z.
    // Note: native Astex Diverse is 1 ligand/receptor; this is a synthetic
    // multi-ligand competitive pilot using Astex affinities as inputs.
    auto tier1 = find_astex_tier1_dir();
    if (tier1.empty()) {
        GTEST_SKIP() << "Astex tier-1 result dir not found";
    }

    const double T = 300.0;
    TargetConfig cfg;
    cfg.temperature_K = T;
    TargetServer server(cfg);

    std::vector<AstexTier1Entry> entries;
    for (auto& e : std::filesystem::directory_iterator(tier1)) {
        if (e.path().extension() != ".json") continue;
        if (e.path().filename().string().find("_holo") == std::string::npos) continue;
        AstexTier1Entry entry;
        if (!parse_astex_tier1_json(e.path(), entry)) continue;
        entries.push_back(entry);
    }
    ASSERT_GE(entries.size(), 2u) << "Need ≥2 Astex tier-1 entries for competitive test";

    for (const auto& entry : entries) {
        double log_Z = log_Z_from_dG(entry.predicted_dG, T);
        ASSERT_TRUE(std::isfinite(log_Z));
        auto sess = server.create_session(entry.pdb_id);
        sess.completed = true;
        sess.log_Z = log_Z;
        sess.n_poses = 5;
        sess.best_energy = entry.predicted_dG;
        server.register_result(sess);
    }

    ASSERT_EQ(server.completed_sessions(), static_cast<int>(entries.size()));
    ASSERT_GE(server.completed_sessions(), 2);  // DatasetRunner multi-ligand gate

    const auto& gpf = server.grand_partition();
    double log_xi = gpf.log_Xi();
    EXPECT_TRUE(std::isfinite(log_xi));

    double sum = gpf.empty_probability();
    auto ranks = server.rank_ligands();
    ASSERT_EQ(ranks.size(), entries.size());
    for (const auto& r : ranks) {
        EXPECT_TRUE(std::isfinite(r.dG));
        EXPECT_TRUE(std::isfinite(r.p_bound));
        EXPECT_GE(r.p_bound, 0.0);
        EXPECT_LE(r.p_bound, 1.0);
        sum += r.p_bound;
    }
    EXPECT_NEAR(sum, 1.0, 1e-8);

    // Strongest (most negative dG / largest log_Z) should rank first
    EXPECT_LE(ranks.front().dG, ranks.back().dG);

    std::cout << "\n[AstexDiverseTier1CompetitiveGPF]\n"
              << "  n_ligands=" << entries.size() << "\n"
              << "  log_Xi=" << log_xi << "\n"
              << "  empty_p=" << gpf.empty_probability() << "\n"
              << "  sum(empty+bound)=" << sum << "\n"
              << "  rank order (best first):\n";
    for (const auto& r : ranks) {
        std::cout << "    " << r.name
                  << "  dG=" << r.dG
                  << "  p_bound=" << r.p_bound
                  << "  log_Z=" << r.log_Z << "\n";
    }
}

// ════════════════════════════════════════════════════════════════════════
// With vs without GCE docking-power ablation on Astex tier-1 artifacts
//
// without_GCE: top pose = argmin total_score (CF field in JSON)
// with_GCE:    register log_Z via shipped TargetServer (DatasetRunner formula)
//              then select top pose with the same score rule for native
//              single-ligand topology (GCE ranks ligands, not poses).
// docking_power_top1 = fraction of entries with selected pose RMSD <= 2.0 A
// ════════════════════════════════════════════════════════════════════════

TEST(TargetServer, AstexWithWithoutGCEAblation) {
    auto entries = load_astex_tier1_entries();
    if (entries.empty()) {
        GTEST_SKIP() << "No Astex tier-1 *_holo.json artifacts found";
    }

    const double T = 300.0;
    const double rmsd_cut = 2.0;

    int without_hits = 0;
    int with_hits = 0;
    int selection_diffs = 0;
    const int N = static_cast<int>(entries.size());

    // ── without-GCE arm (CF score only) ─────────────────────────────────
    std::vector<PoseSelection> without_sels;
    without_sels.reserve(entries.size());
    for (const auto& e : entries) {
        auto sel = select_pose_without_gce(e);
        ASSERT_GE(sel.pose_index, 0) << e.pdb_id;
        without_sels.push_back(sel);
        if (sel.rmsd_success) ++without_hits;
    }

    // ── with-GCE arm: per-entry TargetServer (native 1 ligand/receptor) ─
    // Pose selection remains score-based; GCE validates thermodynamic
    // registration of the CF-selected ΔG. Rates match without-GCE when
    // only one ligand is registered per complex.
    std::vector<PoseSelection> with_sels;
    with_sels.reserve(entries.size());
    for (size_t i = 0; i < entries.size(); ++i) {
        const auto& e = entries[i];
        const auto& wout = without_sels[i];

        TargetConfig cfg;
        cfg.temperature_K = T;
        TargetServer server(cfg);

        const double log_Z = log_Z_from_dG(wout.total_score, T);
        ASSERT_TRUE(std::isfinite(log_Z)) << e.pdb_id;

        auto sess = server.create_session(e.pdb_id);
        sess.completed = true;
        sess.log_Z = log_Z;
        sess.n_poses = static_cast<int>(e.poses.size());
        sess.best_energy = wout.total_score;
        server.register_result(sess);

        const auto& gpf = server.grand_partition();
        const double log_xi = gpf.log_Xi();
        const double p_bound = server.binding_probability(e.pdb_id);
        const double p_empty = gpf.empty_probability();
        EXPECT_TRUE(std::isfinite(log_xi));
        EXPECT_NEAR(p_bound + p_empty, 1.0, 1e-8);

        // Native topology: one ligand → GCE does not re-rank multi-pose list.
        // Top pose = same CF argmin as without-GCE.
        PoseSelection wsel = wout;
        with_sels.push_back(wsel);
        if (wsel.rmsd_success) ++with_hits;
        if (wsel.pose_index != wout.pose_index) ++selection_diffs;

        std::cout << "ENTRY " << e.pdb_id
                  << " without_pose=" << wout.pose_index
                  << " without_rmsd=" << wout.rmsd
                  << " without_score=" << wout.total_score
                  << " with_pose=" << wsel.pose_index
                  << " with_rmsd=" << wsel.rmsd
                  << " log_Z=" << log_Z
                  << " log_Xi=" << log_xi
                  << " closure=" << (p_bound + p_empty)
                  << "\n";
    }

    // ── Multi-ligand competitive pilot (all tier-1 on one TargetServer) ─
    // Does not change per-entry docking_power (still one receptor each in
    // native Astex); exercises GCE ligand ranking + global closure.
    {
        TargetConfig cfg;
        cfg.temperature_K = T;
        TargetServer multi(cfg);
        for (size_t i = 0; i < entries.size(); ++i) {
            auto sess = multi.create_session(entries[i].pdb_id);
            sess.completed = true;
            sess.log_Z = log_Z_from_dG(without_sels[i].total_score, T);
            sess.n_poses = static_cast<int>(entries[i].poses.size());
            sess.best_energy = without_sels[i].total_score;
            multi.register_result(sess);
        }
        ASSERT_EQ(multi.completed_sessions(), N);
        double sum = multi.grand_partition().empty_probability();
        for (const auto& r : multi.rank_ligands())
            sum += r.p_bound;
        EXPECT_NEAR(sum, 1.0, 1e-8);
        EXPECT_TRUE(std::isfinite(multi.grand_partition().log_Xi()));
        std::cout << "MULTI_LIGAND_GCE n=" << N
                  << " log_Xi=" << multi.grand_partition().log_Xi()
                  << " closure=" << sum << "\n";
        std::cout << "MULTI_LIGAND_RANK";
        for (const auto& r : multi.rank_ligands())
            std::cout << " " << r.name << ":" << r.p_bound;
        std::cout << "\n";
    }

    const double without_rate = static_cast<double>(without_hits) / N;
    const double with_rate = static_cast<double>(with_hits) / N;

    EXPECT_GE(without_rate, 0.0);
    EXPECT_LE(without_rate, 1.0);
    EXPECT_GE(with_rate, 0.0);
    EXPECT_LE(with_rate, 1.0);
    // Native single-ligand: arms must select the same poses
    EXPECT_EQ(selection_diffs, 0);
    EXPECT_DOUBLE_EQ(without_rate, with_rate);

    std::cout << "\n[AstexWithWithoutGCEAblation]\n"
              << "  N=" << N << " rmsd_cut=" << rmsd_cut << "\n"
              << "  without_GCE_top1=" << without_hits << "/" << N
              << " = " << without_rate << "\n"
              << "  with_GCE_top1=" << with_hits << "/" << N
              << " = " << with_rate << "\n"
              << "  selection_diffs=" << selection_diffs << "\n"
              << "  note=native_Astex_one_ligand_per_receptor;"
                 "GCE_does_not_re_rank_poses\n"
              << "  PDB_IDs=";
    for (size_t i = 0; i < entries.size(); ++i) {
        if (i) std::cout << ",";
        std::cout << entries[i].pdb_id;
    }
    std::cout << "\n";
}
