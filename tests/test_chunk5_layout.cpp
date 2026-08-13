// tests/test_chunk5_layout.cpp — niche hash, ca_rec flatten, pProp (all default OFF)
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>
#include "niche_hash.h"
#include "ca_rec_flat.h"
#include "pprop.h"

#include <cstdlib>
#include <vector>

TEST(Chunk5Gates, DefaultOff) {
#ifdef _WIN32
    _putenv_s("FLEXAIDDS_NICHE_HASH", "");
    _putenv_s("FLEXAIDDS_CA_REC_FLAT", "");
    _putenv_s("FLEXAIDDS_PPROP_MAX", "");
#else
    unsetenv("FLEXAIDDS_NICHE_HASH");
    unsetenv("FLEXAIDDS_CA_REC_FLAT");
    unsetenv("FLEXAIDDS_PPROP_MAX");
#endif
    EXPECT_FALSE(flexaids::niche_hash_enabled());
    EXPECT_FALSE(flexaids::ca_rec_flat_enabled());
    EXPECT_LT(flexaids::pprop_max_cap(), 0.0);
#ifndef FLEXAIDS_USE_SOA_DISTANCES
    SUCCEED() << "SoA distances compile flag default OFF";
#endif
}

TEST(NicheHash, NeighborsIncludeSelfCellAndAdjacent) {
    float cents[] = {
        0.1f, 0.1f, 0.1f,   // 0 cell (0,0,0)
        0.2f, 0.2f, 0.2f,   // 1 same cell
        1.5f, 0.1f, 0.1f,   // 2 cell (1,0,0)
        10.f, 10.f, 10.f    // 3 far
    };
    auto map = flexaids::niche_hash_build(cents, 4, 1.0f);
    std::vector<int> neigh;
    flexaids::niche_hash_neighbors(map, flexaids::niche_cell_of(0.1f, 0.1f, 0.1f, 1.0f), neigh);
    bool saw0 = false, saw1 = false, saw2 = false, saw3 = false;
    for (int i : neigh) {
        if (i == 0) saw0 = true;
        if (i == 1) saw1 = true;
        if (i == 2) saw2 = true;
        if (i == 3) saw3 = true;
    }
    EXPECT_TRUE(saw0);
    EXPECT_TRUE(saw1);
    EXPECT_TRUE(saw2);
    EXPECT_FALSE(saw3);
}

TEST(CaRecFlat, WalkOrderMatchesPrevChain) {
    ca_struct rec[3]{};
    rec[0].prev = -1;
    rec[1].prev = 0;
    rec[2].prev = 1;
    int ca_index[2] = {2, -1};
    int out[8];
    const int n = flexaids::flatten_ca_rec(ca_index, rec, 0, out, 8);
    ASSERT_EQ(n, 3);
    EXPECT_EQ(out[0], 2);
    EXPECT_EQ(out[1], 1);
    EXPECT_EQ(out[2], 0);
}

TEST(PProp, RankFractionAndCap) {
    EXPECT_DOUBLE_EQ(flexaids::pprop(1, 10), 0.1);
    EXPECT_DOUBLE_EQ(flexaids::delta_pprop(0.2, 0.5), -0.3);
    EXPECT_TRUE(flexaids::pprop_keep(1, 10));  // no cap
#ifdef _WIN32
    _putenv_s("FLEXAIDDS_PPROP_MAX", "0.2");
#else
    setenv("FLEXAIDDS_PPROP_MAX", "0.2", 1);
#endif
    EXPECT_TRUE(flexaids::pprop_keep(1, 10));
    EXPECT_TRUE(flexaids::pprop_keep(2, 10));
    EXPECT_FALSE(flexaids::pprop_keep(3, 10));
#ifdef _WIN32
    _putenv_s("FLEXAIDDS_PPROP_MAX", "");
#else
    unsetenv("FLEXAIDDS_PPROP_MAX");
#endif
}
