// test_matrix_v2_science.cpp — verify Priority-1 VCT matrix v2 science corrections
//
// Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

#include <gtest/gtest.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

TEST(MatrixV2Science, PriorityOneEntries) {
    const char* path = "MC_st0r5.2_6_v2_science.dat";
    std::ifstream probe(path);
    if (!probe.good()) {
        GTEST_SKIP() << "MC_st0r5.2_6_v2_science.dat not found — run make_fa_matrix_v2_science.py";
    }

    // Parse upper triangle in file order (820 lines).
    std::ifstream in(path);
    std::vector<std::pair<std::pair<int, int>, double>> entries;
    for (int ii = 1; ii <= 40; ++ii) {
        for (int jj = ii; jj <= 40; ++jj) {
            std::string line;
            ASSERT_TRUE(std::getline(in, line));
            auto eq = line.find('=');
            ASSERT_NE(eq, std::string::npos);
            entries.push_back({{ii, jj}, std::stod(line.substr(eq + 1))});
        }
    }
    ASSERT_EQ(entries.size(), 820u);

    auto get = [&](int a, int b) {
        if (a > b) std::swap(a, b);
        for (const auto& e : entries) {
            if (e.first.first == a && e.first.second == b) return e.second;
        }
        return 0.0;
    };

    EXPECT_NEAR(get(2, 4), -65.0, 0.01);
    EXPECT_NEAR(get(40, 13), 90.0, 0.01);
    EXPECT_NEAR(get(40, 14), 90.0, 0.01);
    EXPECT_NEAR(get(40, 15), 90.0, 0.01);

    // Unchanged diagonal sanity: C.ar self should match canonical.
    std::ifstream canon("MC_st0r5.2_6.dat");
    ASSERT_TRUE(canon.good());
    std::vector<double> canon_vals;
    for (int k = 0; k < 820; ++k) {
        std::string line;
        ASSERT_TRUE(std::getline(canon, line));
        auto eq = line.find('=');
        canon_vals.push_back(std::stod(line.substr(eq + 1)));
    }
    int idx_4_4 = 0;
    for (int ii = 1; ii <= 40; ++ii) {
        for (int jj = ii; jj <= 40; ++jj) {
            if (ii == 4 && jj == 4) {
                EXPECT_NEAR(get(4, 4), canon_vals[idx_4_4], 0.01);
                return;
            }
            ++idx_4_4;
        }
    }
}