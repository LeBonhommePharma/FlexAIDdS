#include <gtest/gtest.h>

#include "PoseProvenance.h"
#include "flexaid.h"

#include <array>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <memory>

namespace {

class PoseProvenanceTest : public ::testing::Test {
protected:
    void SetUp() override {
        old_seed_ = flexaids_rng::g_master_seed.load();
        old_has_seed_ = flexaids_rng::g_has_master_seed.load();
        old_epoch_ = flexaids_rng::g_seed_epoch.load();
    }
    void TearDown() override {
        flexaids_rng::g_master_seed.store(old_seed_);
        flexaids_rng::g_has_master_seed.store(old_has_seed_);
        flexaids_rng::g_seed_epoch.store(old_epoch_);
    }
private:
    std::uint64_t old_seed_ = 0;
    bool old_has_seed_ = false;
    std::uint64_t old_epoch_ = 0;
};

TEST_F(PoseProvenanceTest, FormatterMatchesLegacyBytesIncludingFullSeedRange) {
    for (const auto seed : {UINT64_C(0), UINT64_C(12345), UINT64_MAX}) {
        for (const int dirty : {0, 1, 2}) {
            char legacy[256];
            std::snprintf(legacy, sizeof(legacy),
                "REMARK FLEXAID.commit=%s FLEXAID.dirty=%d FLEXAID.seed=%llu\n",
                "012345ab", dirty, static_cast<unsigned long long>(seed));
            EXPECT_EQ(flexaids::pose_provenance::format_remark("012345ab", dirty, seed),
                      legacy);
        }
    }
    EXPECT_EQ(flexaids::pose_provenance::format_remark("unknown", 2, 0),
              "REMARK FLEXAID.commit=unknown FLEXAID.dirty=2 FLEXAID.seed=0\n");
}

TEST_F(PoseProvenanceTest, ReadsEffectiveSeedAndRetainsUninitializedFallback) {
    flexaids_rng::g_has_master_seed.store(false);
    flexaids_rng::g_master_seed.store(9876);
    EXPECT_EQ(flexaids::pose_provenance::remark(),
              flexaids::pose_provenance::format_remark(
                  FLEXAIDS_GIT_COMMIT, FLEXAIDS_GIT_DIRTY, 0));
    flexaids_rng::set_master_seed(UINT64_MAX);
    EXPECT_EQ(flexaids::pose_provenance::remark(),
              flexaids::pose_provenance::format_remark(
                  FLEXAIDS_GIT_COMMIT, FLEXAIDS_GIT_DIRTY, UINT64_MAX));
}

TEST_F(PoseProvenanceTest, ReadingProvenanceDoesNotInitializeOrAdvanceRng) {
    flexaids_rng::g_has_master_seed.store(false);
    const auto unseeded_epoch = flexaids_rng::g_seed_epoch.load();
    (void)flexaids::pose_provenance::remark();
    EXPECT_FALSE(flexaids_rng::has_master_seed());
    EXPECT_EQ(flexaids_rng::g_seed_epoch.load(), unseeded_epoch);

    flexaids_rng::set_master_seed(12345);
    auto& rng = flexaids_rng::lazy_thread_rng(0x9A800D);
    auto expected = rng;
    const auto epoch = flexaids_rng::g_seed_epoch.load();
    for (int i = 0; i < 8; ++i) {
        (void)flexaids::pose_provenance::add_to_remarks("REMARK CF=-1.25\n");
    }
    EXPECT_EQ(flexaids_rng::master_seed(), 12345u);
    EXPECT_EQ(flexaids_rng::g_seed_epoch.load(), epoch);
    EXPECT_TRUE(rng == expected);
    EXPECT_EQ(rng(), expected());
}

TEST_F(PoseProvenanceTest, FullBoundedScientificBufferIsPreservedByteForByte) {
    std::string body = "REMARK optimized structure\nREMARK CF=-123.45678\n";
    body.resize(MAX_REMARK - 1, 'x'); // includes a deliberately incomplete tail
    const auto before = body;
    const auto metadata = flexaids::pose_provenance::remark();
    auto expanded = flexaids::pose_provenance::add_to_remarks(body);
    const auto position = body.find('\n') + 1;
    EXPECT_EQ(expanded.substr(position, metadata.size()), metadata);
    expanded.erase(position, metadata.size());
    EXPECT_EQ(expanded, before);
    EXPECT_EQ(body, before);
}

TEST_F(PoseProvenanceTest, EmptyAndUnterminatedBodiesRetainTheirBytes) {
    const auto metadata = flexaids::pose_provenance::remark();
    EXPECT_EQ(flexaids::pose_provenance::add_to_remarks(""), metadata);
    EXPECT_EQ(flexaids::pose_provenance::add_to_remarks("REMARK CF=-1.0"),
              metadata + "REMARK CF=-1.0");
}

class WriterFixture {
public:
    WriterFixture() : fa(std::make_unique<FA_Global>()) {
        const auto tick = std::chrono::steady_clock::now().time_since_epoch().count();
        for (int attempt = 0; attempt < 100; ++attempt) {
            directory = std::filesystem::temp_directory_path() /
                ("pose-provenance-test-" + std::to_string(tick) + "-" +
                 std::to_string(attempt));
            if (std::filesystem::create_directory(directory)) break;
            if (attempt == 99) throw std::runtime_error("Cannot create test directory");
        }
        fa->res_cnt = 2;
        fa->num_het = 1;
        fa->het_res[1] = 2;
        for (int r = 1; r <= 2; ++r) {
            residues[r].fatm = &first[r - 1];
            residues[r].latm = &last[r - 1];
            residues[r].type = r - 1;
            residues[r].chn = 'A';
            residues[r].number = r;
            std::strcpy(residues[r].name, r == 1 ? "ALA" : "LIG");
        }
        for (int i = 1; i <= 3; ++i) {
            atoms[i].ofres = i == 1 ? 1 : 2;
            atoms[i].number = i == 1 ? 1 : 89999 + i;
            std::strcpy(atoms[i].name, i == 1 ? "CA" : i == 2 ? "CL1" : "BR1");
            std::strcpy(atoms[i].element, i == 1 ? "C" : i == 2 ? "Cl" : "Br");
            atoms[i].coor[0] = i * 1.25f;
            atoms[i].coor[1] = i * -2.5f;
            atoms[i].coor[2] = i * 0.5f;
        }
        atoms[2].bond[0] = atoms[3].bond[0] = 1;
        atoms[2].bond[1] = 3;
        atoms[3].bond[1] = 2;
    }
    ~WriterFixture() {
        std::error_code ignored;
        std::filesystem::remove_all(directory, ignored);
    }
    std::string write(const char* name, std::string remarks, int models = 0) {
        auto path = (directory / name).string();
        if (models == 0) {
            EXPECT_EQ(write_pdb(fa.get(), atoms.data(), residues.data(),
                                path.data(), remarks.data()), 0);
        } else {
            for (int i = 1; i <= models; ++i)
                EXPECT_EQ(write_MODEL_pdb(i == 1, i == models, i, fa.get(),
                    atoms.data(), residues.data(), path.data(), remarks.data()), 0);
        }
        std::ifstream input(path, std::ios::binary);
        if (!input) throw std::runtime_error("Missing PDB test output");
        return std::string(std::istreambuf_iterator<char>(input), {});
    }
    std::unique_ptr<FA_Global> fa;
    std::array<atom, 4> atoms{};
    std::array<resid, 3> residues{};
private:
    std::array<int, 2> first{1, 2};
    std::array<int, 2> last{1, 3};
    std::filesystem::path directory;
};

void expect_only_metadata_added(const std::string& before, std::string after) {
    const auto metadata = flexaids::pose_provenance::remark();
    const auto position = after.find(metadata);
    ASSERT_NE(position, std::string::npos);
    EXPECT_EQ(after.find(metadata, position + metadata.size()), std::string::npos);
    after.erase(position, metadata.size());
    EXPECT_EQ(after, before); // all scientific REMARKs, ATOM/HETATM, CONECT, END
    EXPECT_NE(before.find("ATOM  "), std::string::npos);
    EXPECT_NE(before.find("HETATM"), std::string::npos);
    EXPECT_NE(before.find("CONECT"), std::string::npos);
    EXPECT_NE(before.find("REMARK CF=-12.50000"), std::string::npos);
}

const std::string scientific_remarks =
    "REMARK optimized structure\nREMARK CF=-12.50000\nREMARK CF.com=-14.00000\n"
    "REMARK CF.wal= 1.50000\nREMARK rmsd_raw = 3.00000\n";

TEST_F(PoseProvenanceTest, RealPdbWriterAddsOnlyOneMetadataLine) {
    WriterFixture fixture;
    auto full_body = scientific_remarks;
    full_body.resize(MAX_REMARK - 1, 'x');
    expect_only_metadata_added(fixture.write("before.pdb", full_body),
        fixture.write("after.pdb", flexaids::pose_provenance::add_to_remarks(full_body)));
}

TEST_F(PoseProvenanceTest, RealSingleModelWriterPreservesScientificBytes) {
    WriterFixture fixture;
    expect_only_metadata_added(fixture.write("before.pdb", scientific_remarks, 1),
        fixture.write("after.pdb",
            flexaids::pose_provenance::add_to_remarks(scientific_remarks), 1));
}

TEST_F(PoseProvenanceTest, RealMultiModelWriterAddsMetadataOnlyToFirstHeader) {
    WriterFixture fixture;
    expect_only_metadata_added(fixture.write("before.pdb", scientific_remarks, 3),
        fixture.write("after.pdb",
            flexaids::pose_provenance::add_to_remarks(scientific_remarks), 3));
}

} // namespace
