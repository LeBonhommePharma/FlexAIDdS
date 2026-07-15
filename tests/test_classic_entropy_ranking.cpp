// tests/test_classic_entropy_ranking.cpp
// Classic FlexAID entropy ranking (soft-β SoftBeta G̃ ≡ ACF elects rank-0).
// Rollback: FA->force_cf_rank_emission = true (P3b CF emission / physical F).
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include "../LIB/BindingMode.h"
#include "../LIB/statmech.h"
#include "../LIB/gaboom.h"
#include "../LIB/SoftBetaFreeEnergy.h"
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include <limits>

namespace {

// Minimal ACF-vs-CF election policy (mirrors cluster.cpp gate, pure logic).
// Returns index of rank-0 under classic entropy vs force_cf.
int elect_rank0(const std::vector<double>& acf,
                const std::vector<double>& cf,
                bool force_cf_rank_emission,
                unsigned temperature)
{
    const int n = static_cast<int>(acf.size());
    if (n == 0) return -1;
    std::vector<int> order(n);
    for (int i = 0; i < n; ++i) order[i] = i;

    if (temperature > 0) {
        std::stable_sort(order.begin(), order.end(),
            [&](int a, int b) { return acf[a] < acf[b]; });
    }
    const bool classic = (temperature > 0) && !force_cf_rank_emission;
    if (!classic) {
        std::stable_sort(order.begin(), order.end(),
            [&](int a, int b) { return cf[a] < cf[b]; });
    }
    return order[0];
}

}  // namespace

class ClassicEntropyRankingTest : public ::testing::Test {
protected:
    FA_Global* mock_fa = nullptr;
    GB_Global* mock_gb = nullptr;
    VC_Global* mock_vc = nullptr;
    chromosome* mock_chroms = nullptr;
    genlim* mock_gene_lim = nullptr;
    atom* mock_atoms = nullptr;
    resid* mock_residue = nullptr;
    gridpoint* mock_cleftgrid = nullptr;
    BindingPopulation* pop = nullptr;
    static constexpr double T = 300.0;
    static constexpr double EPS = 1e-6;
    static constexpr int N_CHROM = 40;

    void SetUp() override {
        mock_fa = new FA_Global();
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wnontrivial-memcall"
        std::memset(mock_fa, 0, sizeof(FA_Global));
#pragma clang diagnostic pop
        mock_fa->temperature = static_cast<uint>(T);
        mock_fa->force_cf_rank_emission = false;  // classic product default
        mock_fa->beta = 1.0 / T;

        mock_gb = new GB_Global();
        std::memset(mock_gb, 0, sizeof(GB_Global));
        mock_gb->num_genes = 6;

        mock_vc = new VC_Global();
        std::memset(mock_vc, 0, sizeof(VC_Global));

        mock_chroms = new chromosome[N_CHROM];
        for (int i = 0; i < N_CHROM; ++i) {
            mock_chroms[i].genes = new gene[mock_gb->num_genes];
            std::memset(mock_chroms[i].genes, 0, sizeof(gene) * mock_gb->num_genes);
            mock_chroms[i].evalue = 0.0;
            mock_chroms[i].app_evalue = 0.0;
            mock_chroms[i].fitnes = 0.0;
            mock_chroms[i].status = 'n';
        }
        mock_gene_lim = new genlim[mock_gb->num_genes];
        std::memset(mock_gene_lim, 0, sizeof(genlim) * mock_gb->num_genes);
        mock_atoms = new atom[10];
        std::memset(mock_atoms, 0, sizeof(atom) * 10);
        mock_residue = new resid[2];
        std::memset(mock_residue, 0, sizeof(resid) * 2);
        mock_cleftgrid = new gridpoint[100];
        std::memset(mock_cleftgrid, 0, sizeof(gridpoint) * 100);

        pop = new BindingPopulation(
            mock_fa, mock_gb, mock_vc,
            mock_chroms, mock_gene_lim,
            mock_atoms, mock_residue, mock_cleftgrid,
            N_CHROM);
    }

    void TearDown() override {
        delete pop;
        delete[] mock_cleftgrid;
        delete[] mock_residue;
        delete[] mock_atoms;
        delete[] mock_gene_lim;
        for (int i = 0; i < N_CHROM; ++i) delete[] mock_chroms[i].genes;
        delete[] mock_chroms;
        delete mock_vc;
        delete mock_gb;
        delete mock_fa;
    }

    Pose make_pose(double cf, int index) {
        mock_chroms[index].app_evalue = cf;
        mock_chroms[index].evalue = cf;
        std::vector<float> empty;
        Pose p(&mock_chroms[index], index, 0, 0.0f, static_cast<uint>(T), empty);
        p.CF = static_cast<float>(cf);
        return p;
    }
};

// ─── 1HNN-class emission policy (pure logic, no full GA) ───────────────────
// Live .cad: cluster ACF=-263.4 (freq 29, CF weak) vs ACF=-49.3 (CF=-189.9).
TEST(ClassicEntropyEmissionPolicy, DefaultElectsLowestACF) {
    // indices: 0 = CF champion, 3 = ACF champion (1HNN-like)
    std::vector<double> acf = {-49.3, -83.4, -48.9, -263.4, -221.2};
    std::vector<double> cf  = {-189.9, -120.0, -100.0, -72.1, -80.0};
    EXPECT_EQ(elect_rank0(acf, cf, /*force_cf=*/false, 300), 3);
}

TEST(ClassicEntropyEmissionPolicy, ForceCFRestoresP3bLowestCF) {
    std::vector<double> acf = {-49.3, -83.4, -48.9, -263.4, -221.2};
    std::vector<double> cf  = {-189.9, -120.0, -100.0, -72.1, -80.0};
    EXPECT_EQ(elect_rank0(acf, cf, /*force_cf=*/true, 300), 0);
}

TEST(ClassicEntropyEmissionPolicy, TemperatureZeroForcesCF) {
    std::vector<double> acf = {-49.3, -263.4};
    std::vector<double> cf  = {-189.9, -72.1};
    EXPECT_EQ(elect_rank0(acf, cf, /*force_cf=*/false, 0), 0);
}

// ─── Dense basin beats singleton best CF (classic BindingMode F) ───────────
TEST_F(ClassicEntropyRankingTest, DenseModeBeatsSingletonBestCF) {
    ASSERT_FALSE(mock_fa->force_cf_rank_emission);

    BindingMode dense(pop);
    // 25 poses at middling CF ≈ -100
    for (int i = 0; i < 25; ++i) {
        Pose p = make_pose(-100.0 - 0.1 * (i % 5), i);
        dense.add_Pose(p);
    }

    BindingMode singleton(pop);
    // One deep false minimum at CF = -190
    Pose best = make_pose(-190.0, 30);
    singleton.add_Pose(best);

    // Register both into global Z (classic contract)
    pop->add_BindingMode(dense);
    pop->add_BindingMode(singleton);

    // Re-fetch after Entropize path via set_energy on add
    const BindingMode& m0 = pop->get_binding_mode(0);
    const BindingMode& m1 = pop->get_binding_mode(1);

    // Rank-0 after Entropize must be the dense mode (classic entropy).
    // get_binding_mode sorts by classic F.
    EXPECT_GT(m0.get_BindingMode_size(), m1.get_BindingMode_size());
    EXPECT_LT(m0.get_cached_energy(), m1.get_cached_energy());
    EXPECT_GT(m0.compute_entropy(), 0.0);
}

TEST_F(ClassicEntropyRankingTest, ForceCFPrefersLowestCFSingleton) {
    mock_fa->force_cf_rank_emission = true;

    BindingMode dense(pop);
    for (int i = 0; i < 25; ++i) {
        Pose p = make_pose(-100.0, i);
        dense.add_Pose(p);
    }
    BindingMode singleton(pop);
    Pose best = make_pose(-190.0, 30);
    singleton.add_Pose(best);

    pop->add_BindingMode(dense);
    pop->add_BindingMode(singleton);

    const BindingMode& m0 = pop->get_binding_mode(0);
    // Physical per-mode F collapses toward min energy → singleton wins.
    EXPECT_EQ(m0.get_BindingMode_size(), 1);
}

TEST_F(ClassicEntropyRankingTest, SoftBetaPoseWeightNotPhysicalKB) {
    Pose p = make_pose(-90.0, 0);
    // Classic: exp(90/300) = exp(0.3) ≈ 1.34986
    const double classic_w = std::exp(90.0 / 300.0);
    EXPECT_NEAR(p.boltzmann_weight, classic_w, 1e-9);
    // Physical would be exp(90/(kB*300)) with kB≈0.001987 → huge
    const double physical_w = std::exp(90.0 / (statmech::kB_kcal * 300.0));
    EXPECT_GT(physical_w / classic_w, 10.0);
}

TEST_F(ClassicEntropyRankingTest, PhysicalLedgerStillCallableUnderClassicRanking) {
    BindingMode mode(pop);
    for (int i = 0; i < 5; ++i) {
        Pose p = make_pose(-50.0 - i, i);
        mode.add_Pose(p);
    }
    pop->add_BindingMode(mode);

    // Ranking energy is classic; diagnostic ledger remains physical StatMech.
    auto thermo = pop->get_binding_mode(0).get_thermodynamics();
    EXPECT_TRUE(std::isfinite(thermo.free_energy));
    EXPECT_TRUE(std::isfinite(thermo.entropy));
}

// Vibrational correction stays on classic ranking F (FlexAIDdS, not stripped).
// With normal_modes == 0 the correction is 0 — API still composes F_conf + vib.
TEST_F(ClassicEntropyRankingTest, ClassicRankingIncludesVibCorrectionTerm) {
    ASSERT_FALSE(mock_fa->force_cf_rank_emission);
    mock_fa->normal_modes = 0;  // vib term evaluates to 0, but path is exercised

    BindingMode mode(pop);
    for (int i = 0; i < 4; ++i) {
        Pose p = make_pose(-80.0 - i, i);
        mode.add_Pose(p);
    }
    pop->add_BindingMode(mode);

    const BindingMode& ranked = pop->get_binding_mode(0);
    const double F = ranked.compute_energy();
    const double H = ranked.compute_enthalpy();
    const double S = ranked.compute_entropy();
    const double F_conf = H - T * S;
    // No modes → vib = 0; ranking F must match classic conf free energy.
    EXPECT_NEAR(F, F_conf, EPS);
    EXPECT_GT(S, 0.0);
}

// ─── SoftBeta free-energy identity (shared ranking objective) ───────────────
// G̃ = H̃ − T·S̃ ≡ E_min − T ln Z  (cluster ACF form). Used by cluster.cpp,
// BindingMode classic ranking, and DatasetRunner S1 election.

TEST(SoftBetaIdentity, EmptyIsInf) {
    auto fe = flexaids::soft_beta::free_energy({}, 300.0);
    EXPECT_TRUE(std::isinf(fe.G));
}

TEST(SoftBetaIdentity, SingletonEqualsEnergy) {
    auto fe = flexaids::soft_beta::free_energy({-42.5}, 300.0);
    EXPECT_NEAR(fe.G, -42.5, 1e-12);
    EXPECT_NEAR(fe.H, -42.5, 1e-12);
    EXPECT_NEAR(fe.S, 0.0, 1e-12);
    EXPECT_EQ(fe.n, 1);
}

TEST(SoftBetaIdentity, GEqualsHminusTS) {
    const std::vector<double> energies = {-100.0, -98.0, -95.0, -90.0, -88.0};
    const double T = 300.0;
    auto fe = flexaids::soft_beta::free_energy(energies, T);
    EXPECT_NEAR(fe.G, fe.H - T * fe.S, 1e-9);
}

TEST(SoftBetaIdentity, GEqualsEminMinusTlnZ) {
    const std::vector<double> energies = {-100.0, -98.0, -95.0, -90.0, -88.0};
    const double T = 300.0;
    auto fe = flexaids::soft_beta::free_energy(energies, T);
    const double acf_form = fe.Emin - T * std::log(fe.Z);
    EXPECT_NEAR(fe.G, acf_form, 1e-9);
    EXPECT_NEAR(flexaids::soft_beta::acf(energies, T), fe.G, 1e-12);
}

TEST(SoftBetaIdentity, DenseBasinBeatsDeepSingleton) {
    // 1HNN-class: dense middling basin vs sparse deep false minimum
    std::vector<double> dense(29, -72.0);
    for (int i = 0; i < 29; ++i)
        dense[static_cast<size_t>(i)] = -72.0 - 0.05 * (i % 7);
    const std::vector<double> deep = {-189.9};
    const double T = 300.0;
    const double G_dense = flexaids::soft_beta::free_energy_G(dense, T);
    const double G_deep  = flexaids::soft_beta::free_energy_G(deep, T);
    EXPECT_LT(G_dense, G_deep);  // elect lowest G → dense wins
}

TEST(SoftBetaIdentity, OffsetInvariance) {
    // Adding constant Δ to all CF must shift G by Δ (soft-β property)
    const std::vector<double> energies = {-50.0, -48.0, -45.0};
    const double T = 250.0;
    const double delta = 12.3;
    std::vector<double> energies2 = energies;
    for (double& e : energies2) e += delta;
    auto a = flexaids::soft_beta::free_energy(energies, T);
    auto b = flexaids::soft_beta::free_energy(energies2, T);
    EXPECT_NEAR(b.G, a.G + delta, 1e-9);
    EXPECT_NEAR(b.S, a.S, 1e-9);
}

TEST_F(ClassicEntropyRankingTest, BindingModeMatchesSoftBetaLocal) {
    ASSERT_FALSE(mock_fa->force_cf_rank_emission);
    BindingMode mode(pop);
    std::vector<double> cfs;
    for (int i = 0; i < 8; ++i) {
        const double cf = -80.0 - 0.5 * i;
        cfs.push_back(cf);
        Pose p = make_pose(cf, i);
        mode.add_Pose(p);
    }
    pop->add_BindingMode(mode);
    const BindingMode& ranked = pop->get_binding_mode(0);
    const double T = 300.0;
    auto fe = flexaids::soft_beta::free_energy(cfs, T);
    EXPECT_NEAR(ranked.compute_enthalpy(), fe.H, 1e-9);
    EXPECT_NEAR(ranked.compute_entropy(), fe.S, 1e-9);
    // vib=0, nat=0 → ranking energy == SoftBeta G
    EXPECT_NEAR(ranked.compute_energy(), fe.G, 1e-9);
}
