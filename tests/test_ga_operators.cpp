// tests/test_ga_operators.cpp — untested GA helper functions + calc_rmsp
// Tests: calc_poss, set_bins, validate_dups, bin_print, deelig_search, calc_rmsp
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <map>

#include "flexaid.h"
#include "gaboom.h"
#include "ga_constants.h"

// calc_rmsp.cpp defines a float overload not declared in any header
float calc_rmsp(int npar, float* icv, float* icu);

// ═══════════════════════════════════════════════════════════════════════
// calc_poss — computes total search-space size from gene limits
// ═══════════════════════════════════════════════════════════════════════

class CalcPossTest : public ::testing::Test {
protected:
    void make_genlim(genlim& g, double mn, double mx, double dl, double bn, double nbn, int mp) {
        g.min = mn; g.max = mx; g.del = dl;
        g.bin = bn; g.nbin = nbn; g.map = mp;
    }
};

TEST_F(CalcPossTest, SingleGene) {
    genlim gl[1];
    make_genlim(gl[0], 0, 10, 1, 0, 10, 0);
    EXPECT_DOUBLE_EQ(calc_poss(gl, 1), 10.0);
}

TEST_F(CalcPossTest, TwoGenes) {
    genlim gl[2];
    make_genlim(gl[0], 0, 10, 1, 0, 10, 0);
    make_genlim(gl[1], 0, 5, 1, 0, 5, 0);
    EXPECT_DOUBLE_EQ(calc_poss(gl, 2), 50.0);
}

TEST_F(CalcPossTest, ThreeGenesProduct) {
    genlim gl[3];
    make_genlim(gl[0], 0, 0, 0, 0, 4, 0);
    make_genlim(gl[1], 0, 0, 0, 0, 3, 0);
    make_genlim(gl[2], 0, 0, 0, 0, 5, 0);
    EXPECT_DOUBLE_EQ(calc_poss(gl, 3), 60.0);
}

TEST_F(CalcPossTest, ZeroGenesReturnsZero) {
    genlim gl[1];
    make_genlim(gl[0], 0, 0, 0, 0, 10, 0);
    EXPECT_DOUBLE_EQ(calc_poss(gl, 0), 0.0);
}

TEST_F(CalcPossTest, SingleNbinOne) {
    genlim gl[1];
    make_genlim(gl[0], 0, 0, 0, 0, 1, 0);
    EXPECT_DOUBLE_EQ(calc_poss(gl, 1), 1.0);
}

TEST_F(CalcPossTest, LargeSearchSpace) {
    genlim gl[4];
    make_genlim(gl[0], 0, 0, 0, 0, 360, 0);
    make_genlim(gl[1], 0, 0, 0, 0, 360, 0);
    make_genlim(gl[2], 0, 0, 0, 0, 360, 0);
    make_genlim(gl[3], 0, 0, 0, 0, 100, 0);
    double expected = 360.0 * 360.0 * 360.0 * 100.0;
    EXPECT_DOUBLE_EQ(calc_poss(gl, 4), expected);
}

// ═══════════════════════════════════════════════════════════════════════
// set_bins — sets bin size and nbin from min/max/del
// Multi-gene overload sets both bin and nbin. Single-gene overload sets only bin.
// ═══════════════════════════════════════════════════════════════════════

TEST(SetBinsMultiTest, ExactDivision) {
    genlim arr[2];
    std::memset(arr, 0, sizeof(arr));
    arr[0].min = 0.0; arr[0].max = 10.0; arr[0].del = 2.0; arr[0].map = 0;
    arr[1].min = -5.0; arr[1].max = 5.0; arr[1].del = 1.0; arr[1].map = 0;

    set_bins(arr, 2);

    EXPECT_DOUBLE_EQ(arr[0].nbin, 5.0);
    EXPECT_DOUBLE_EQ(arr[0].bin, 1.0 / 5.0);
    EXPECT_DOUBLE_EQ(arr[1].nbin, 10.0);
    EXPECT_DOUBLE_EQ(arr[1].bin, 1.0 / 10.0);
}

TEST(SetBinsMultiTest, RemainderRoundsUp) {
    genlim arr[1];
    std::memset(arr, 0, sizeof(arr));
    arr[0].min = 0.0; arr[0].max = 10.0; arr[0].del = 3.0; arr[0].map = 0;

    set_bins(arr, 1);

    // 10/3 = 3.333, remainder > 0 → nbin = 3.333 + 1 = 4.333
    EXPECT_DOUBLE_EQ(arr[0].nbin, 10.0/3.0 + 1.0);
    EXPECT_DOUBLE_EQ(arr[0].bin, 1.0 / (10.0/3.0 + 1.0));
}

TEST(SetBinsMultiTest, ExactDivisionNoRoundUp) {
    genlim arr[1];
    std::memset(arr, 0, sizeof(arr));
    arr[0].min = 0.0; arr[0].max = 10.0; arr[0].del = 5.0; arr[0].map = 0;

    set_bins(arr, 1);

    EXPECT_DOUBLE_EQ(arr[0].nbin, 2.0);
    EXPECT_DOUBLE_EQ(arr[0].bin, 1.0 / 2.0);
}

TEST(SetBinsMultiTest, MapAddsOne) {
    genlim arr[1];
    std::memset(arr, 0, sizeof(arr));
    arr[0].min = 0.0; arr[0].max = 10.0; arr[0].del = 5.0; arr[0].map = 1;

    set_bins(arr, 1);

    EXPECT_DOUBLE_EQ(arr[0].nbin, 3.0);
    EXPECT_DOUBLE_EQ(arr[0].bin, 1.0 / 3.0);
}

TEST(SetBinsMultiTest, RemainderAndMapBothAdd) {
    genlim arr[1];
    std::memset(arr, 0, sizeof(arr));
    arr[0].min = 0.0; arr[0].max = 10.0; arr[0].del = 3.0; arr[0].map = 1;

    set_bins(arr, 1);

    // 10/3 = 3.333 + 1 (remainder) + 1 (map) = 5.333
    EXPECT_DOUBLE_EQ(arr[0].nbin, 10.0/3.0 + 1.0 + 1.0);
    EXPECT_DOUBLE_EQ(arr[0].bin, 1.0 / (10.0/3.0 + 2.0));
}

// Single-gene overload: sets bin but NOT nbin
TEST(SetBinsSingleTest, NoMap) {
    genlim gl;
    std::memset(&gl, 0, sizeof(genlim));
    gl.min = -1.0; gl.max = 1.0; gl.del = 0.5; gl.map = 0;
    set_bins(&gl);

    // 2.0/0.5 = 4.0 exact → no round-up, no map
    EXPECT_DOUBLE_EQ(gl.bin, 1.0 / 4.0);
}

TEST(SetBinsSingleTest, WithMap) {
    genlim gl;
    std::memset(&gl, 0, sizeof(genlim));
    gl.min = 0.0; gl.max = 5.0; gl.del = 1.0; gl.map = 1;
    set_bins(&gl);

    // 5.0/1.0 = 5.0 exact + 1 for map = 6
    EXPECT_DOUBLE_EQ(gl.bin, 1.0 / 6.0);
}

// ═══════════════════════════════════════════════════════════════════════
// validate_dups — auto-enables duplicates when pop > search space
// ═══════════════════════════════════════════════════════════════════════

class ValidateDupsTest : public ::testing::Test {
protected:
    GB_Global GB;
    genlim gl[2];
    void SetUp() override {
        std::memset(&GB, 0, sizeof(GB_Global));
        std::memset(gl, 0, sizeof(gl));
    }
};

TEST_F(ValidateDupsTest, PopSmallerThanPoss_NoAutoEnable) {
    GB.num_chrom = 10;
    GB.duplicates = 0;
    gl[0].nbin = 100; gl[1].nbin = 100;  // poss = 10000

    validate_dups(&GB, gl, 2);
    EXPECT_EQ(GB.duplicates, 0);
}

TEST_F(ValidateDupsTest, PopLargerThanPoss_AutoEnables) {
    GB.num_chrom = 200;
    GB.duplicates = 0;
    gl[0].nbin = 10; gl[1].nbin = 10;  // poss = 100

    validate_dups(&GB, gl, 2);
    EXPECT_EQ(GB.duplicates, 1);
}

TEST_F(ValidateDupsTest, AlreadyEnabledStaysEnabled) {
    GB.num_chrom = 200;
    GB.duplicates = 1;
    gl[0].nbin = 10; gl[1].nbin = 10;

    validate_dups(&GB, gl, 2);
    EXPECT_EQ(GB.duplicates, 1);
}

TEST_F(ValidateDupsTest, PopExactlyEqualsPoss_NoAutoEnable) {
    GB.num_chrom = 100;
    GB.duplicates = 0;
    gl[0].nbin = 10; gl[1].nbin = 10;  // poss = 100

    validate_dups(&GB, gl, 2);
    EXPECT_EQ(GB.duplicates, 0);
}

// ═══════════════════════════════════════════════════════════════════════
// bin_print — prints binary representation to stdout
// ═══════════════════════════════════════════════════════════════════════

class BinPrintTest : public ::testing::Test {
protected:
    void SetUp() override { fflush(stdout); }
    // Capture stdout via pipe
    std::string capture_bin_print(int dec, int len) {
        fflush(stdout);
        // Redirect stdout to a temp file
        char tmpfile[] = "/tmp/bintest_XXXXXX";
        int fd = mkstemp(tmpfile);
        FILE* saved = stdout;
        stdout = fdopen(fd, "w");

        bin_print(dec, len);

        fflush(stdout);
        fclose(stdout);
        stdout = saved;

        // Read back
        FILE* f = fopen(tmpfile, "r");
        char buf[256] = {};
        fgets(buf, sizeof(buf), f);
        fclose(f);
        unlink(tmpfile);
        return std::string(buf);
    }
};

TEST_F(BinPrintTest, ZeroLen4) {
    EXPECT_EQ(capture_bin_print(0, 4), "0000");
}

TEST_F(BinPrintTest, FiveLen3) {
    // 5 = 101 binary
    EXPECT_EQ(capture_bin_print(5, 3), "101");
}

TEST_F(BinPrintTest, SevenLen4) {
    // 7 = 0111 with len=4
    EXPECT_EQ(capture_bin_print(7, 4), "0111");
}

TEST_F(BinPrintTest, FifteenLen4) {
    // 15 = 1111
    EXPECT_EQ(capture_bin_print(15, 4), "1111");
}

TEST_F(BinPrintTest, OneLen8) {
    // 1 = 00000001 with len=8
    EXPECT_EQ(capture_bin_print(1, 8), "00000001");
}

// ═══════════════════════════════════════════════════════════════════════
// deelig_search — tree traversal for dead-end elimination
// ═══════════════════════════════════════════════════════════════════════

class DeeligSearchTest : public ::testing::Test {
protected:
    deelig_node_struct root;

    void SetUp() override {
        root.parent = nullptr;
        root.childs.clear();
    }

    // Build a simple tree: root -> level1[value] -> level2[value]
    void build_tree(const std::vector<std::pair<int,int>>& path_values) {
        // path_values: {(level, value)} — builds linear chain from root
        // For multi-branch at same level, call manually
        deelig_node_struct* node = &root;
        for (auto& [level, value] : path_values) {
            if (node->childs.find(value) == node->childs.end()) {
                auto* child = new deelig_node_struct();
                child->parent = node;
                node->childs[value] = child;
            }
            node = node->childs[value];
        }
    }

    void TearDown() override {
        // Recursively free tree (root is stack-allocated)
        free_node(&root);
        root.childs.clear();
    }

    void free_node(deelig_node_struct* n) {
        for (auto& [k, child] : n->childs) {
            free_node(child);
            delete child;
        }
        n->childs.clear();
    }
};

TEST_F(DeeligSearchTest, EmptyTreeNotFound) {
    int list[] = {0, 1, 2};  // fdih=2, list[1]=1, list[2]=2
    EXPECT_EQ(deelig_search(&root, list, 2), 0);
}

TEST_F(DeeligSearchTest, ExactPathFound) {
    // Build: root -> 1 -> 3 -> 5
    auto* n1 = new deelig_node_struct(); n1->parent = &root;
    auto* n2 = new deelig_node_struct(); n2->parent = n1;
    root.childs[1] = n1;
    n1->childs[3] = n2;

    int list[] = {0, 1, 3};  // fdih=2
    EXPECT_EQ(deelig_search(&root, list, 2), 1);

    // Cleanup
    delete n2; delete n1;
    root.childs.clear();
}

TEST_F(DeeligSearchTest, PartialPathNotFound) {
    auto* n1 = new deelig_node_struct(); n1->parent = &root;
    root.childs[1] = n1;
    n1->childs[3] = new deelig_node_struct();

    int list[] = {0, 1, 9};  // list[2]=9 not in n1->childs
    EXPECT_EQ(deelig_search(&root, list, 2), 0);

    free_node(n1); delete n1;
    root.childs.clear();
}

TEST_F(DeeligSearchTest, SentinelWildcard) {
    // root -> 1 -> SENTINEL -> 5
    auto* n1 = new deelig_node_struct(); n1->parent = &root;
    auto* n2 = new deelig_node_struct(); n2->parent = n1;
    root.childs[1] = n1;
    n1->childs[GA_DEELIG_SENTINEL] = n2;

    // Searching for 1->7 (not in tree) but SENTINEL wildcard matches
    int list[] = {0, 1, 7};
    EXPECT_EQ(deelig_search(&root, list, 2), 1);

    free_node(&root);
}

TEST_F(DeeligSearchTest, SentinelNotUsedWhenValueMatches) {
    // root -> 1 -> {3: nodeA, SENTINEL: nodeB}
    auto* n1 = new deelig_node_struct(); n1->parent = &root;
    root.childs[1] = n1;
    n1->childs[3] = new deelig_node_struct();
    n1->childs[GA_DEELIG_SENTINEL] = new deelig_node_struct();

    int list[] = {0, 1, 3};  // exact match found first
    EXPECT_EQ(deelig_search(&root, list, 2), 1);

    free_node(&root);
}

TEST_F(DeeligSearchTest, FdihZeroReturnsOne) {
    // fdih=0 means the for loop doesn't execute → falls through to return(1)
    int list[] = {0};
    EXPECT_EQ(deelig_search(&root, list, 0), 1);
}

TEST_F(DeeligSearchTest, DeepPathFound) {
    // Build a 4-level path: root -> 2 -> 4 -> 6 -> 8
    auto* n1 = new deelig_node_struct(); n1->parent = &root;
    auto* n2 = new deelig_node_struct(); n2->parent = n1;
    auto* n3 = new deelig_node_struct(); n3->parent = n2;
    root.childs[2] = n1;
    n1->childs[4] = n2;
    n2->childs[6] = n3;

    int list[] = {0, 2, 4, 6};  // fdih=3
    EXPECT_EQ(deelig_search(&root, list, 3), 1);

    free_node(&root);
}

// ═══════════════════════════════════════════════════════════════════════
// calc_rmsp — float RMSD between two arrays
// ═══════════════════════════════════════════════════════════════════════

class CalcRmspTest : public ::testing::Test {
protected:
    static constexpr float EPS = 1e-5f;
};

TEST_F(CalcRmspTest, IdenticalArrays) {
    float a[] = {1.0f, 2.0f, 3.0f};
    float b[] = {1.0f, 2.0f, 3.0f};
    EXPECT_NEAR(calc_rmsp(3, a, b), 0.0f, EPS);
}

TEST_F(CalcRmspTest, UnitOffset) {
    float a[] = {0.0f, 0.0f, 0.0f};
    float b[] = {1.0f, 0.0f, 0.0f};
    // RMSD = sqrt((1^2 + 0 + 0)/3) = sqrt(1/3)
    EXPECT_NEAR(calc_rmsp(3, a, b), std::sqrt(1.0f / 3.0f), EPS);
}

TEST_F(CalcRmspTest, AllOffset) {
    float a[] = {0.0f, 0.0f};
    float b[] = {3.0f, 4.0f};
    // RMSD = sqrt((9+16)/2) = sqrt(12.5)
    EXPECT_NEAR(calc_rmsp(2, a, b), std::sqrt(12.5f), EPS);
}

TEST_F(CalcRmspTest, SingleElement) {
    float a[] = {5.0f};
    float b[] = {8.0f};
    // RMSD = sqrt((3^2)/1) = 3.0
    EXPECT_NEAR(calc_rmsp(1, a, b), 3.0f, EPS);
}

TEST_F(CalcRmspTest, NegativeValues) {
    float a[] = {-1.0f, -2.0f};
    float b[] = {1.0f, 2.0f};
    // diffs: 2, 4 → RMSD = sqrt((4+16)/2) = sqrt(10)
    EXPECT_NEAR(calc_rmsp(2, a, b), std::sqrt(10.0f), EPS);
}

TEST_F(CalcRmspTest, MixedPositiveNegative) {
    float a[] = {1.0f, -1.0f, 0.0f};
    float b[] = {-1.0f, 1.0f, 0.0f};
    // diffs: 2, -2, 0 → RMSD = sqrt((4+4+0)/3) = sqrt(8/3)
    EXPECT_NEAR(calc_rmsp(3, a, b), std::sqrt(8.0f / 3.0f), EPS);
}

// ═══════════════════════════════════════════════════════════════════════
// swap_chrom — trivial chromosome swap
// ═══════════════════════════════════════════════════════════════════════

TEST(SwapChromTest, SwapsValues) {
    chromosome a, b;
    std::memset(&a, 0, sizeof(chromosome));
    std::memset(&b, 0, sizeof(chromosome));
    a.cf.com = 1.0;
    b.cf.com = 2.0;

    swap_chrom(&a, &b);

    EXPECT_DOUBLE_EQ(a.cf.com, 2.0);
    EXPECT_DOUBLE_EQ(b.cf.com, 1.0);
}

// ═══════════════════════════════════════════════════════════════════════
// hash_genes — gene hash function
// ═══════════════════════════════════════════════════════════════════════

TEST(HashGenesTest, SameGenesSameHash) {
    gene g1[3], g2[3];
    std::memset(g1, 0, sizeof(gene) * 3);
    std::memset(g2, 0, sizeof(gene) * 3);
    g1[0].to_ic = 1.0; g1[1].to_ic = 2.0; g1[2].to_ic = 3.0;
    g2[0].to_ic = 1.0; g2[1].to_ic = 2.0; g2[2].to_ic = 3.0;
    EXPECT_EQ(hash_genes(g1, 3), hash_genes(g2, 3));
}

TEST(HashGenesTest, DifferentGenesDifferentHash) {
    gene g1[2], g2[2];
    std::memset(g1, 0, sizeof(gene) * 2);
    std::memset(g2, 0, sizeof(gene) * 2);
    g1[0].to_ic = 1.0; g1[1].to_ic = 2.0;
    g2[0].to_ic = 9.0; g2[1].to_ic = 8.0;
    EXPECT_NE(hash_genes(g1, 2), hash_genes(g2, 2));
}
