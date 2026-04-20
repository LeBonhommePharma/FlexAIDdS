// tests/test_gist_hbond.cpp
// Unit tests for GIST blurry trilinear displacement and directional H-bond scoring
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include "../LIB/flexaid.h"
#include "../LIB/GISTEvaluator.h"
#include "../LIB/HBondEvaluator.h"
#include <cmath>
#include <fstream>
#include <filesystem>

static constexpr double EPS = 1e-6;

// ===========================================================================
// GIST EVALUATOR TESTS
// ===========================================================================

class GISTEvaluatorTest : public ::testing::Test {
protected:
    GISTEvaluator gist;
    std::string tmp_dir;

    void SetUp() override {
        tmp_dir = std::filesystem::temp_directory_path().string() + "/gist_test_" +
                  std::to_string(::testing::UnitTest::GetInstance()->random_seed());
        std::filesystem::create_directories(tmp_dir);
    }

    void TearDown() override {
        std::filesystem::remove_all(tmp_dir);
    }

    // Write a minimal 3x3x3 DX file with given data
    void write_dx(const std::string& path, const std::vector<double>& data,
                  double ox = 0.0, double oy = 0.0, double oz = 0.0,
                  double sp = 1.0, int nx = 3, int ny = 3, int nz = 3) {
        std::ofstream out(path);
        out << "# Test DX file\n";
        out << "object 1 class gridpositions counts " << nx << " " << ny << " " << nz << "\n";
        out << "origin " << ox << " " << oy << " " << oz << "\n";
        out << "delta " << sp << " 0 0\n";
        out << "delta 0 " << sp << " 0\n";
        out << "delta 0 0 " << sp << "\n";
        out << "object 2 class gridconnections counts " << nx << " " << ny << " " << nz << "\n";
        out << "object 3 class array type double rank 0 items " << data.size() << " data follows\n";
        for (size_t i = 0; i < data.size(); ++i) {
            out << data[i];
            if ((i + 1) % 3 == 0) out << "\n";
            else out << " ";
        }
        out << "\nattribute \"dep\" string \"positions\"\n";
    }
};

TEST_F(GISTEvaluatorTest, LoadValidDXFiles) {
    int n = 3 * 3 * 3;
    std::vector<double> dg_data(n, 2.0);   // all unfavorable
    std::vector<double> rho_data(n, 5.0);  // above density cutoff

    std::string dg_path = tmp_dir + "/dg.dx";
    std::string rho_path = tmp_dir + "/rho.dx";
    write_dx(dg_path, dg_data);
    write_dx(rho_path, rho_data);

    ASSERT_TRUE(gist.load_dx(dg_path, rho_path));
    EXPECT_TRUE(gist.loaded);
    EXPECT_EQ(gist.nx, 3);
    EXPECT_EQ(gist.ny, 3);
    EXPECT_EQ(gist.nz, 3);
}

TEST_F(GISTEvaluatorTest, LoadMissingFileReturnsFalse) {
    EXPECT_FALSE(gist.load_dx("/nonexistent/dg.dx", "/nonexistent/rho.dx"));
    EXPECT_FALSE(gist.loaded);
}

TEST_F(GISTEvaluatorTest, ScoreAtomReturnsZeroWhenNotLoaded) {
    EXPECT_DOUBLE_EQ(gist.score_atom(0.0, 0.0, 0.0, 1.5), 0.0);
}

TEST_F(GISTEvaluatorTest, ScoreAtomAtGridCenter) {
    // 3x3x3 grid: origin=(0,0,0), spacing=1.0
    // Set center voxel (1,1,1) to high DG and density
    int n = 27;
    std::vector<double> dg_data(n, 0.0);
    std::vector<double> rho_data(n, 0.0);

    // Index (1,1,1) = 1*9 + 1*3 + 1 = 13
    dg_data[13] = 3.0;
    rho_data[13] = 6.0;

    std::string dg_path = tmp_dir + "/dg2.dx";
    std::string rho_path = tmp_dir + "/rho2.dx";
    write_dx(dg_path, dg_data);
    write_dx(rho_path, rho_data);

    ASSERT_TRUE(gist.load_dx(dg_path, rho_path));

    // Score atom exactly at center voxel (1.0, 1.0, 1.0)
    gist.divisor = 2.0;
    double score = gist.score_atom(1.0, 1.0, 1.0, 2.0);
    // sigma = 2.0/2.0 = 1.0, at center dist_sq=0 → weight=exp(0)=1.0
    // score = 1.0 * 3.0 = 3.0
    EXPECT_NEAR(score, 3.0, EPS);
}

TEST_F(GISTEvaluatorTest, ScoreDecaysWithDistance) {
    int n = 27;
    std::vector<double> dg_data(n, 0.0);
    std::vector<double> rho_data(n, 0.0);
    dg_data[13] = 3.0;
    rho_data[13] = 6.0;

    std::string dg_path = tmp_dir + "/dg3.dx";
    std::string rho_path = tmp_dir + "/rho3.dx";
    write_dx(dg_path, dg_data);
    write_dx(rho_path, rho_data);

    ASSERT_TRUE(gist.load_dx(dg_path, rho_path));
    gist.divisor = 2.0;

    // Score at center
    double score_center = gist.score_atom(1.0, 1.0, 1.0, 2.0);
    // Score offset by 0.5 Å
    double score_offset = gist.score_atom(1.5, 1.0, 1.0, 2.0);
    // Score should decay with distance
    EXPECT_GT(score_center, score_offset);
    EXPECT_GT(score_offset, 0.0);
}

TEST_F(GISTEvaluatorTest, CutoffsFilterFavorableWater) {
    int n = 27;
    std::vector<double> dg_data(n, 0.5);   // below DG cutoff (1.0)
    std::vector<double> rho_data(n, 6.0);  // above density cutoff

    std::string dg_path = tmp_dir + "/dg4.dx";
    std::string rho_path = tmp_dir + "/rho4.dx";
    write_dx(dg_path, dg_data);
    write_dx(rho_path, rho_data);

    ASSERT_TRUE(gist.load_dx(dg_path, rho_path));

    // No voxels pass the free-energy cutoff → score = 0
    double score = gist.score_atom(1.0, 1.0, 1.0, 2.0);
    EXPECT_DOUBLE_EQ(score, 0.0);
}

// ===========================================================================
// H-BOND EVALUATOR TESTS
// ===========================================================================

class HBondEvaluatorTest : public ::testing::Test {
};

TEST_F(HBondEvaluatorTest, DonorTypeClassification) {
    // N types (6–12) are donors
    EXPECT_TRUE(hbond::is_donor_type(6));   // N.1
    EXPECT_TRUE(hbond::is_donor_type(8));   // N.3
    EXPECT_TRUE(hbond::is_donor_type(12));  // N.PL3
    EXPECT_TRUE(hbond::is_donor_type(14));  // O.3

    // C, S, P types are not donors
    EXPECT_FALSE(hbond::is_donor_type(1));  // C.1
    EXPECT_FALSE(hbond::is_donor_type(3));  // C.3
    EXPECT_FALSE(hbond::is_donor_type(17)); // S.2
    EXPECT_FALSE(hbond::is_donor_type(22)); // P.3
}

TEST_F(HBondEvaluatorTest, AcceptorTypeClassification) {
    // All N and O types (6–16) are acceptors
    EXPECT_TRUE(hbond::is_acceptor_type(6));   // N.1
    EXPECT_TRUE(hbond::is_acceptor_type(13));  // O.2
    EXPECT_TRUE(hbond::is_acceptor_type(15));  // O.CO2
    EXPECT_TRUE(hbond::is_acceptor_type(16));  // O.AR

    // Non-N/O types are not acceptors
    EXPECT_FALSE(hbond::is_acceptor_type(3));  // C.3
    EXPECT_FALSE(hbond::is_acceptor_type(17)); // S.2
}

TEST_F(HBondEvaluatorTest, SaltBridgePairDetection) {
    // N.4 (9) + O.CO2 (15) = salt bridge
    EXPECT_TRUE(hbond::is_salt_bridge_pair(9, 15));
    EXPECT_TRUE(hbond::is_salt_bridge_pair(15, 9));  // symmetric
    // N.3 (8) + O.CO2 (15)
    EXPECT_TRUE(hbond::is_salt_bridge_pair(8, 15));

    // N.AM (11) + O.CO2 (15) = NOT salt bridge (N.AM is not cationic)
    EXPECT_FALSE(hbond::is_salt_bridge_pair(11, 15));
    // Two donors = NOT salt bridge
    EXPECT_FALSE(hbond::is_salt_bridge_pair(8, 9));
}

TEST_F(HBondEvaluatorTest, LinearHBondGivesMaxMultiplier) {
    // D at (0,0,0), H at (1,0,0), A at (2,0,0) → theta = 180° (linear)
    double mult = hbond::angular_multiplier(
        0.0, 0.0, 0.0,   // donor
        1.0, 0.0, 0.0,   // hydrogen
        2.0, 0.0, 0.0,   // acceptor
        false);

    // At theta=180°: cos(180°) = -1, cos²=1, (pi - pi)^6 = 0 → exp(0)=1
    // g(180°) = 1.0 * 1.0 = 1.0
    EXPECT_NEAR(mult, 1.0, EPS);
}

TEST_F(HBondEvaluatorTest, RightAngleHBondGivesZero) {
    // D at (0,0,0), H at (1,0,0), A at (1,1,0) → theta = 90°
    double mult = hbond::angular_multiplier(
        0.0, 0.0, 0.0,   // donor
        1.0, 0.0, 0.0,   // hydrogen
        1.0, 1.0, 0.0,   // acceptor
        false);

    // At theta=90°: cos(90°) = 0 → cos² = 0 → result ≈ 0
    EXPECT_NEAR(mult, 0.0, EPS);
}

TEST_F(HBondEvaluatorTest, ObtuseAngleGivesSmallMultiplier) {
    // D at (0,0,0), H at (1,0,0), A at (2,1,0) → theta between 90° and 180°
    double mult = hbond::angular_multiplier(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        2.0, 1.0, 0.0,
        false);

    // Should be less than linear (1.0) but greater than right-angle (0.0)
    EXPECT_GT(mult, 0.0);
    EXPECT_LT(mult, 1.0);
}

TEST_F(HBondEvaluatorTest, SaltBridgeBroaderTolerance) {
    // Same geometry for standard vs salt bridge
    double mult_std = hbond::angular_multiplier(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        2.0, 2.0, 0.0,  // fairly bent
        false);

    double mult_salt = hbond::angular_multiplier(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        2.0, 2.0, 0.0,  // same geometry
        true);

    // Salt bridge should be more tolerant (higher multiplier for bent angles)
    EXPECT_GT(mult_salt, mult_std);
}

TEST_F(HBondEvaluatorTest, SaltBridgeLinearIsMax) {
    // Linear salt bridge → multiplier should be 1.0
    double mult = hbond::angular_multiplier(
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        2.0, 0.0, 0.0,
        true);

    EXPECT_NEAR(mult, 1.0, EPS);
}

TEST_F(HBondEvaluatorTest, HeavyAtomProxyOptimalDistance) {
    // Optimal D···A distance (2.9 Å) should give maximum multiplier
    double mult = hbond::angular_multiplier_heavy_atom(
        0.0, 0.0, 0.0,
        2.9, 0.0, 0.0,
        2.9, false);

    EXPECT_NEAR(mult, 1.0, 0.01);
}

TEST_F(HBondEvaluatorTest, HeavyAtomProxyLongDistance) {
    // Long D···A distance should give low multiplier
    double mult = hbond::angular_multiplier_heavy_atom(
        0.0, 0.0, 0.0,
        5.0, 0.0, 0.0,
        5.0, false);

    EXPECT_LT(mult, 0.01);
}

TEST_F(HBondEvaluatorTest, HeavyAtomSaltBridgeBypass) {
    // Salt bridge always returns 1.0 in heavy-atom proxy
    double mult = hbond::angular_multiplier_heavy_atom(
        0.0, 0.0, 0.0,
        5.0, 0.0, 0.0,
        5.0, true);

    EXPECT_NEAR(mult, 1.0, EPS);
}

TEST_F(HBondEvaluatorTest, ZeroLengthVectorsReturnZero) {
    // H coincident with D → degenerate
    double mult = hbond::angular_multiplier(
        1.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        2.0, 0.0, 0.0,
        false);

    EXPECT_NEAR(mult, 0.0, EPS);
}

// ===========================================================================
// cfstr FIELD VERIFICATION
// ===========================================================================

// Verify the cfstr struct has the new gist and hbond fields
TEST(CfstrFieldTest, NewFieldsExist) {
    // This is a compile-time test — if it compiles, the fields exist
    cfstr cf;
    cf.gist = 1.234;
    cf.hbond = -0.567;
    EXPECT_DOUBLE_EQ(cf.gist, 1.234);
    EXPECT_DOUBLE_EQ(cf.hbond, -0.567);
}

// ===========================================================================
// GIST DX PARSING EDGE CASES
// ===========================================================================

TEST_F(GISTEvaluatorTest, DXWithComments) {
    int n = 8;  // 2x2x2
    std::vector<double> data(n, 1.5);

    std::string path = tmp_dir + "/commented.dx";
    {
        std::ofstream out(path);
        out << "# This is a comment\n";
        out << "# Another comment\n";
        out << "object 1 class gridpositions counts 2 2 2\n";
        out << "origin 0 0 0\n";
        out << "delta 0.5 0 0\n";
        out << "delta 0 0.5 0\n";
        out << "delta 0 0 0.5\n";
        out << "object 2 class gridconnections counts 2 2 2\n";
        out << "object 3 class array type double rank 0 items 8 data follows\n";
        for (int i = 0; i < n; ++i) out << data[i] << " ";
        out << "\n";
    }

    GISTEvaluator g;
    std::vector<double> rho_data(n, 5.0);
    std::string rho_path = tmp_dir + "/rho_commented.dx";
    write_dx(rho_path, rho_data, 0, 0, 0, 0.5, 2, 2, 2);

    ASSERT_TRUE(g.load_dx(path, rho_path));
    EXPECT_EQ(g.nx, 2);
    EXPECT_NEAR(g.spacing, 0.5, 0.01);
}
