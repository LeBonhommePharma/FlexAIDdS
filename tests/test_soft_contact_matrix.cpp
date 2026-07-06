// test_soft_contact_matrix.cpp — GoogleTest for the 256×256 soft contact matrix
//
// Tests: encoding, SYBYL projection, binary I/O, FastOPTICS clustering,
// scoring, and AVX batch correctness.

#include <gtest/gtest.h>
#include "atom_typing_256.h"
#include "soft_contact_matrix.h"
#include "shannon_matrix_scorer.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>
#include <algorithm>
#include <chrono>
#include <filesystem>

// ─── atom_typing_256 tests ──────────────────────────────────────────────────

TEST(AtomTyping256, EncodeDecodeRoundtrip) {
    for (int base = 0; base < 64; ++base) {
        for (int donor = 0; donor < 2; ++donor) {
            for (int acceptor = 0; acceptor < 2; ++acceptor) {
                uint8_t code = atom256::encode(
                    static_cast<uint8_t>(base), donor != 0, acceptor != 0);
                EXPECT_EQ(atom256::get_base(code), base);
                EXPECT_EQ(atom256::get_hbond_donor(code), donor != 0);
                EXPECT_EQ(atom256::get_hbond_acceptor(code), acceptor != 0);
                EXPECT_EQ(atom256::get_hbond(code), donor || acceptor);
            }
        }
    }
}

TEST(AtomTyping256, All256CodesUnique) {
    std::vector<uint8_t> codes;
    for (int base = 0; base < 64; ++base)
        for (int donor = 0; donor < 2; ++donor)
            for (int acceptor = 0; acceptor < 2; ++acceptor)
                codes.push_back(atom256::encode(
                    static_cast<uint8_t>(base), donor != 0, acceptor != 0));
    EXPECT_EQ(codes.size(), 256u);
    std::sort(codes.begin(), codes.end());
    auto it = std::unique(codes.begin(), codes.end());
    EXPECT_EQ(it, codes.end()) << "Duplicate codes found";
}

TEST(AtomTyping256, SybylToBaseMapping) {
    // All 40 SYBYL types should map to valid base types
    for (int s = 1; s <= 40; ++s) {
        uint8_t base = atom256::sybyl_to_base(s);
        EXPECT_GE(base, 0);
        EXPECT_LT(base, 64) << "SYBYL type " << s << " → base " << (int)base;
    }
}

TEST(AtomTyping256, NoSolventFallback) {
    // Previously collapsed types should now have distinct base types
    // SE(27), MG(28), SR(29), CU(30), MN(31), HG(32), CD(33), NI(34),
    // CO.OH(38), DUMMY(39) should NOT map to Solvent (31)
    int collapsed[] = {27, 28, 29, 30, 31, 32, 33, 34, 38, 39};
    for (int s : collapsed) {
        uint8_t base = atom256::sybyl_to_base(s);
        EXPECT_NE(base, atom256::Solvent)
            << "SYBYL type " << s << " should not collapse to Solvent";
    }
    // Only SYBYL 40 (actual solvent) should map to Solvent
    EXPECT_EQ(atom256::sybyl_to_base(40), atom256::Solvent);
}

TEST(AtomTyping256, BaseToSybylParent) {
    // First 22 base types should round-trip cleanly
    for (int b = 0; b < 22; ++b) {
        int sybyl = atom256::base_to_sybyl_parent(b);
        EXPECT_GE(sybyl, 1);
        EXPECT_LE(sybyl, 40);
        uint8_t back = atom256::sybyl_to_base(sybyl);
        EXPECT_EQ(back, b) << "base=" << b << " → sybyl=" << sybyl
                           << " → base=" << (int)back;
    }
    // Extended types (32-41) should round-trip through SYBYL
    for (int b = 32; b <= 41; ++b) {
        int sybyl = atom256::base_to_sybyl_parent(b);
        EXPECT_GE(sybyl, 1);
        EXPECT_LE(sybyl, 40);
        uint8_t back = atom256::sybyl_to_base(sybyl);
        EXPECT_EQ(back, b) << "base=" << b << " → sybyl=" << sybyl
                           << " → base=" << (int)back;
    }
}

TEST(AtomTyping256, ContextRefinement) {
    // C_ar with heteroatom neighbor → C_ar_hetadj
    EXPECT_EQ(atom256::refine_base_type(atom256::C_ar, true, true, false),
              atom256::C_ar_hetadj);
    // C_ar at bridgehead → C_pi_bridge
    EXPECT_EQ(atom256::refine_base_type(atom256::C_ar, true, false, true),
              atom256::C_pi_bridge);
    // Non-aromatic carbon unchanged
    EXPECT_EQ(atom256::refine_base_type(atom256::C_sp3, false, true, false),
              atom256::C_sp3);
}

TEST(AtomTyping256, HBondRoles) {
    EXPECT_FALSE(atom256::classify_hbond_donor(atom256::N_sp3, 0.0f, 1));
    EXPECT_TRUE(atom256::classify_hbond_acceptor(atom256::N_sp3, 0.0f, 0));
    EXPECT_FALSE(atom256::classify_hbond_donor(atom256::O_sp2, -0.5f, 0));
    EXPECT_TRUE(atom256::classify_hbond_acceptor(atom256::O_sp2, -0.5f, 0));
    EXPECT_FALSE(atom256::classify_hbond_donor(atom256::N_sp, 0.0f, 0));
    EXPECT_TRUE(atom256::classify_hbond_acceptor(atom256::N_sp, 0.0f, 0));
    EXPECT_TRUE(atom256::classify_hbond_donor(atom256::N_am, 0.0f, 1));
    EXPECT_FALSE(atom256::classify_hbond_acceptor(atom256::N_am, 0.0f, 1));
    EXPECT_FALSE(atom256::classify_hbond_donor(atom256::N_quat, 0.0f, 0));
    EXPECT_FALSE(atom256::classify_hbond_acceptor(atom256::N_quat, 0.0f, 0));
}

TEST(AtomTyping256, EncodeFromSybyl) {
    uint8_t code = atom256::encode_from_sybyl(4, 0.1f, 0);
    EXPECT_EQ(atom256::get_base(code), atom256::C_ar);
    EXPECT_FALSE(atom256::get_hbond_donor(code));
    EXPECT_FALSE(atom256::get_hbond_acceptor(code));

    code = atom256::encode_from_sybyl(11, 0.0f, 1);
    EXPECT_EQ(atom256::get_base(code), atom256::N_am);
    EXPECT_TRUE(atom256::get_hbond_donor(code));
    EXPECT_FALSE(atom256::get_hbond_acceptor(code));

    code = atom256::encode_from_sybyl(8, 0.0f, 0);
    EXPECT_EQ(atom256::get_base(code), atom256::N_sp3);
    EXPECT_FALSE(atom256::get_hbond_donor(code));
    EXPECT_TRUE(atom256::get_hbond_acceptor(code));

    // MG (28) → Metal_Mg, not Solvent
    code = atom256::encode_from_sybyl(28, 0.5f, 0);
    EXPECT_EQ(atom256::get_base(code), atom256::Metal_Mg);
    EXPECT_NE(atom256::get_base(code), atom256::Solvent);
}

// ─── SoftContactMatrix tests ────────────────────────────────────────────────

TEST(SoftContactMatrix, ZeroInitialisation) {
    scm::SoftContactMatrix mat;
    mat.zero();
    for (int i = 0; i < scm::MATRIX_SIZE; ++i)
        EXPECT_FLOAT_EQ(mat.data[i], 0.0f);
}

TEST(SoftContactMatrix, LookupSetRoundtrip) {
    scm::SoftContactMatrix mat;
    mat.zero();
    mat.set(10, 20, 3.14f);
    EXPECT_FLOAT_EQ(mat.lookup(10, 20), 3.14f);
    EXPECT_FLOAT_EQ(mat.lookup(20, 10), 0.0f);  // not symmetric until symmetrise
}

TEST(SoftContactMatrix, Symmetrise) {
    scm::SoftContactMatrix mat;
    mat.zero();
    mat.set(10, 20, 4.0f);
    mat.set(20, 10, 2.0f);
    mat.symmetrise();
    EXPECT_FLOAT_EQ(mat.lookup(10, 20), 3.0f);
    EXPECT_FLOAT_EQ(mat.lookup(20, 10), 3.0f);
}

TEST(SoftContactMatrix, ScoreContacts) {
    scm::SoftContactMatrix mat;
    mat.zero();
    mat.set(1, 2, 0.5f);
    mat.set(3, 4, 1.5f);

    uint8_t ta[] = {1, 3};
    uint8_t tb[] = {2, 4};
    float areas[] = {2.0f, 3.0f};

    float score = mat.score_contacts(ta, tb, areas, 2);
    EXPECT_FLOAT_EQ(score, 0.5f * 2.0f + 1.5f * 3.0f);
}

TEST(SoftContactMatrix, BinaryIORoundtrip) {
    scm::SoftContactMatrix mat;
    std::mt19937 rng(42);
    std::uniform_real_distribution<float> dist(-5.0f, 5.0f);
    for (int i = 0; i < scm::MATRIX_SIZE; ++i)
        mat.data[i] = dist(rng);

    const char* path = "/tmp/test_scm_roundtrip.shnn";
    ASSERT_TRUE(mat.save(path));

    scm::SoftContactMatrix loaded;
    ASSERT_TRUE(loaded.load(path));

    for (int i = 0; i < scm::MATRIX_SIZE; ++i)
        EXPECT_FLOAT_EQ(loaded.data[i], mat.data[i]);

    std::remove(path);
}

TEST(SoftContactMatrix, InvalidMagicFails) {
    const char* path = "/tmp/test_scm_badmagic.bin";
    FILE* fp = fopen(path, "wb");
    uint32_t bad = 0xDEADBEEF;
    fwrite(&bad, 4, 1, fp);
    fclose(fp);

    scm::SoftContactMatrix mat;
    EXPECT_FALSE(mat.load(path));
    std::remove(path);
}

TEST(SoftContactMatrix, ProjectTo40x40) {
    scm::SoftContactMatrix mat;
    mat.zero();

    // Set all C_sp × C_sp interactions (base=0) to 5.0
    for (int ci = 0; ci < 256; ++ci) {
        if (atom256::get_base(ci) == atom256::C_sp) {
            for (int cj = 0; cj < 256; ++cj) {
                if (atom256::get_base(cj) == atom256::C_sp)
                    mat.set(ci, cj, 5.0f);
            }
        }
    }

    auto proj = mat.project_to_40x40();
    // SYBYL C.1 = type 1, 0-indexed = 0
    EXPECT_FLOAT_EQ(proj[0 * 40 + 0], 5.0f);
    // Other types should be 0 (or whatever average of zeros)
    EXPECT_FLOAT_EQ(proj[1 * 40 + 1], 0.0f);
}

// ─── FastOPTICS tests ───────────────────────────────────────────────────────

TEST(FastOPTICS, DetectsIdenticalRows) {
    scm::SoftContactMatrix mat;
    mat.zero();

    // Make rows 0-49 identical, rows 50-99 identical (different from 0-49)
    std::mt19937 rng(123);
    std::uniform_real_distribution<float> dist(-1.0f, 1.0f);

    float pattern_a[scm::MATRIX_DIM], pattern_b[scm::MATRIX_DIM];
    for (int d = 0; d < scm::MATRIX_DIM; ++d) {
        pattern_a[d] = dist(rng);
        pattern_b[d] = dist(rng) + 10.0f;  // well separated
    }

    for (int i = 0; i < 50; ++i)
        std::memcpy(mat.row(i), pattern_a, sizeof(pattern_a));
    for (int i = 50; i < 100; ++i)
        std::memcpy(mat.row(i), pattern_b, sizeof(pattern_b));
    // Rows 100-255 are all zero (third cluster or noise)

    auto result = scm::find_super_clusters(mat, 3, 15, 42);
    EXPECT_GE(result.n_clusters, 2) << "Should detect at least 2 super-clusters";
    EXPECT_EQ(result.order.size(), 256u);
    EXPECT_EQ(result.cluster_labels.size(), 256u);

    // Rows 0-49 should be in the same cluster
    int label_a = result.cluster_labels[0];
    for (int i = 1; i < 50; ++i)
        EXPECT_EQ(result.cluster_labels[i], label_a);
}

TEST(FastOPTICS, ReachabilityNonNegative) {
    scm::SoftContactMatrix mat;
    std::mt19937 rng(99);
    std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
    for (int i = 0; i < scm::MATRIX_SIZE; ++i) mat.data[i] = dist(rng);

    auto result = scm::find_super_clusters(mat, 5, 10, 99);
    for (float r : result.reachability) {
        EXPECT_GE(r, 0.0f);
    }
}

// ─── ShannonMatrixScorer tests ──────────────────────────────────────────────

TEST(ShannonMatrixScorer, SinglePoseScoring) {
    scm::SoftContactMatrix mat;
    mat.zero();
    mat.set(1, 2, -2.0f);  // favourable contact

    scorer::ShannonMatrixScorer sc(mat, 300.0, 100.0f);

    std::vector<scorer::Contact> contacts;
    contacts.push_back({1, 2, 3.0f, 3.5f, 1.7f, 1.7f, 0.0f, 0.0f});

    auto result = sc.score_pose(contacts);
    EXPECT_TRUE(result.survived_filter);
    EXPECT_LT(result.matrix_score, 0.0f);  // favourable
}

TEST(ShannonMatrixScorer, FilterRejectsHighScores) {
    scm::SoftContactMatrix mat;
    mat.zero();
    mat.set(1, 2, 50.0f);  // unfavourable

    scorer::ShannonMatrixScorer sc(mat, 300.0, 10.0f);

    std::vector<scorer::Contact> contacts;
    contacts.push_back({1, 2, 1.0f, 3.5f, 1.7f, 1.7f, 0.0f, 0.0f});

    auto result = sc.score_pose(contacts);
    EXPECT_FALSE(result.survived_filter);
}

TEST(ShannonMatrixScorer, EnsembleEntropy) {
    scm::SoftContactMatrix mat;
    std::mt19937 rng(42);
    std::uniform_real_distribution<float> dist(-2.0f, 2.0f);
    for (int i = 0; i < scm::MATRIX_SIZE; ++i) mat.data[i] = dist(rng);

    scorer::ShannonMatrixScorer sc(mat, 300.0, 1000.0f);  // permissive filter

    // Create 10 poses with random contacts
    std::vector<std::vector<scorer::Contact>> poses(10);
    for (auto& p : poses) {
        for (int c = 0; c < 5; ++c) {
            uint8_t ta = rng() % 256;
            uint8_t tb = rng() % 256;
            p.push_back({ta, tb, 1.0f, 3.0f, 1.7f, 1.7f, 0.0f, 0.0f});
        }
    }

    auto result = sc.score_ensemble(poses);
    EXPECT_GT(result.n_survivors, 0);
    EXPECT_GE(result.shannonEntropy, 0.0);
    EXPECT_FALSE(std::isnan(result.deltaG));
}

// ─── performance benchmark ──────────────────────────────────────────────────

TEST(SoftContactMatrix, LookupPerformance) {
    scm::SoftContactMatrix mat;
    std::mt19937 rng(42);
    std::uniform_real_distribution<float> dist(-5.0f, 5.0f);
    for (int i = 0; i < scm::MATRIX_SIZE; ++i) mat.data[i] = dist(rng);

    const int N = 1'000'000;
    std::vector<uint8_t> ta(N), tb(N);
    std::vector<float> areas(N);
    for (int i = 0; i < N; ++i) {
        ta[i] = rng() % 256;
        tb[i] = rng() % 256;
        areas[i] = 1.0f;
    }

    auto start = std::chrono::high_resolution_clock::now();
    volatile float result = mat.score_contacts(ta.data(), tb.data(),
                                                areas.data(), N);
    auto end = std::chrono::high_resolution_clock::now();
    double ns = std::chrono::duration<double, std::nano>(end - start).count() / N;

    (void)result;
    printf("Lookup performance: %.2f ns/lookup (%d lookups)\n", ns, N);
    EXPECT_LT(ns, 50.0) << "Lookup should be < 50 ns even without cache warmth";
}
