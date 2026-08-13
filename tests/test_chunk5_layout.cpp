// tests/test_chunk5_layout.cpp — niche hash, ca_rec flatten, pProp (all default OFF)
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>
#include "niche_hash.h"
#include "ca_rec_flat.h"
#include "pprop.h"
#include "get_yval.h"

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

TEST(CaRecFlat, SkipFirstThenEachIndexOnce) {
    // Same walk vcfunction uses: flatten, take first, then ca_rec_next on skip.
    ca_struct rec[3]{};
    rec[0].prev = -1;
    rec[1].prev = 0;
    rec[2].prev = 1;
    int ca_index[1] = {2};
    int flat[8];
    const int nflat = flexaids::flatten_ca_rec(ca_index, rec, 0, flat, 8);
    ASSERT_EQ(nflat, 3);
    int flat_k = 1;
    int curr = flat[0];
    EXPECT_EQ(curr, 2);
    // Skip the first node (bonded / already-seen continue path).
    curr = flexaids::ca_rec_next(true, curr, rec, flat, nflat, flat_k);
    std::vector<int> seen;
    while (curr != -1) {
        seen.push_back(curr);
        curr = flexaids::ca_rec_next(true, curr, rec, flat, nflat, flat_k);
    }
    ASSERT_EQ(seen.size(), 2u);
    EXPECT_EQ(seen[0], 1);
    EXPECT_EQ(seen[1], 0);
    EXPECT_NE(seen[0], seen[1]);
}

TEST(CaRecFlat, MoreThan100ContactsMatchesUncappedWalkAndCf) {
    // S5: MAX_CONT is 100. A 120-contact chain must not drop a prefix.
    constexpr int kN = 120;
    std::vector<ca_struct> rec(kN);
    for (int i = 0; i < kN; ++i) {
        rec[i].prev = (i == 0) ? -1 : i - 1;
        rec[i].atom = 1000 + i;
        rec[i].area = 0.1 * (i + 1);
        rec[i].dist = 3.0;
    }
    int ca_index[1] = {kN - 1};

    std::vector<int> legacy;
    int curr = ca_index[0];
    while (curr != -1) {
        legacy.push_back(curr);
        curr = rec[curr].prev;
    }
    ASSERT_EQ(legacy.size(), static_cast<size_t>(kN));

    int too_small[100];
    EXPECT_EQ(flexaids::flatten_ca_rec(ca_index, rec.data(), 0, too_small, 100),
              -kN);

    std::vector<int> flat;
    const int n = flexaids::flatten_ca_rec_all(ca_index, rec.data(), 0, flat);
    ASSERT_EQ(n, kN);
    ASSERT_EQ(flat.size(), legacy.size());
    for (int i = 0; i < kN; ++i) EXPECT_EQ(flat[i], legacy[i]);

#ifdef _WIN32
    _putenv_s("FLEXAIDDS_CA_REC_FLAT", "1");
#else
    setenv("FLEXAIDDS_CA_REC_FLAT", "1", 1);
#endif
    EXPECT_TRUE(flexaids::ca_rec_flat_enabled());
#ifdef _WIN32
    _putenv_s("FLEXAIDDS_GET_YVAL_LUT", "");
#else
    unsetenv("FLEXAIDDS_GET_YVAL_LUT");
#endif

    energy_values ev{};
    ev.next_value = nullptr;
    std::vector<float> xs{0.f, 1.f};
    std::vector<float> ys{0.f, 10.f};
    std::vector<float> sl{(10.f - 0.f) / 1.f};
    energy_matrix em{};
    em.weight = 0;
    em.energy_values = &ev;
    em.flat_n = 2;
    em.flat_x = xs.data();
    em.flat_y = ys.data();
    em.flat_slope = sl.data();
    constexpr double kSurfA = 50.0;

    auto cf_of = [&](const std::vector<int>& idx) {
        double cf = 0.0;
        for (int id : idx)
            cf += get_yval(&em, rec[id].area / kSurfA);
        return cf;
    };
    EXPECT_DOUBLE_EQ(cf_of(flat), cf_of(legacy));
    EXPECT_NE(cf_of(std::vector<int>(legacy.begin(), legacy.begin() + 100)),
              cf_of(legacy));

#ifdef _WIN32
    _putenv_s("FLEXAIDDS_CA_REC_FLAT", "");
#else
    unsetenv("FLEXAIDDS_CA_REC_FLAT");
#endif
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
