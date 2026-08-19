// test_binary_snapshot.cpp — unit tests for the FASN binary snapshot format
//
// SnapshotWriter/SnapshotReader are the persistence layer for GA pose
// trajectories. A silent format or endianness regression here corrupts every
// downstream analysis without any visible error, so these tests pin the
// on-disk layout explicitly (magic, version, record stride) rather than only
// checking that a roundtrip happens to agree with itself.
//
// All tests use fixed inputs and a per-test temporary directory.
//
// Apache-2.0 — see LICENSE.

#include "BinarySnapshot.h"

#include <gtest/gtest.h>

#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

using flexaids::BinaryFileHeader;
using flexaids::SnapshotReader;
using flexaids::SnapshotWriter;

namespace {

class BinarySnapshotTest : public ::testing::Test {
protected:
    void SetUp() override {
        const ::testing::TestInfo* info =
            ::testing::UnitTest::GetInstance()->current_test_info();
        dir_ = fs::temp_directory_path() /
               (std::string("flexaidds_snap_") + info->name());
        fs::remove_all(dir_);
        fs::create_directories(dir_);
    }
    void TearDown() override {
        std::error_code ec;
        fs::remove_all(dir_, ec);
    }
    std::string path(const std::string& name) const {
        return (dir_ / name).string();
    }
    fs::path dir_;
};

std::vector<float> ramp(uint32_t n_atoms, float base) {
    std::vector<float> c(static_cast<size_t>(n_atoms) * 3);
    for (size_t i = 0; i < c.size(); ++i)
        c[i] = base + static_cast<float>(i) * 0.25f;
    return c;
}

}  // namespace

// --- size helpers ---------------------------------------------------------

TEST(BinarySnapshotLayout, HeaderIsSixteenPackedBytes) {
    EXPECT_EQ(flexaids::file_header_bytes(), 16u);
    EXPECT_EQ(sizeof(BinaryFileHeader), 16u);
}

TEST(BinarySnapshotLayout, RecordStrideIsScorePlusGenerationPlusCoords) {
    EXPECT_EQ(flexaids::snapshot_record_bytes(0), 8u);
    EXPECT_EQ(flexaids::snapshot_record_bytes(1), 8u + 12u);
    EXPECT_EQ(flexaids::snapshot_record_bytes(10), 8u + 120u);
}

TEST(BinarySnapshotLayout, EndianHelpersRoundTrip) {
    EXPECT_EQ(flexaids::le_to_host_u16(flexaids::host_to_le_u16(0xBEEF)), 0xBEEF);
    EXPECT_EQ(flexaids::le_to_host_u32(flexaids::host_to_le_u32(0xDEADBEEFu)),
              0xDEADBEEFu);
    EXPECT_FLOAT_EQ(flexaids::le_to_host_f32(flexaids::host_to_le_f32(-12.5f)),
                    -12.5f);
}

TEST(BinarySnapshotLayout, LittleEndianHostRoundTripIsIdentity) {
    // On a little-endian host the conversions must be pure identity, which is
    // what makes the mmap zero-copy path legitimate.
    if (!flexaids::host_is_little_endian()) GTEST_SKIP() << "big-endian host";
    EXPECT_EQ(flexaids::host_to_le_u32(0x01020304u), 0x01020304u);
}

// --- write / read roundtrip ----------------------------------------------

TEST_F(BinarySnapshotTest, RoundTripPreservesScoreGenerationAndCoordinates) {
    const uint32_t kAtoms = 5;
    const std::string p = path("round.fasn");

    {
        SnapshotWriter w(p, kAtoms);
        EXPECT_EQ(w.n_atoms(), kAtoms);
        EXPECT_EQ(w.count(), 0u);
        w.write_snapshot(-11.5f, 3, ramp(kAtoms, 1.0f));
        w.write_snapshot(-7.25f, 9, ramp(kAtoms, 100.0f));
        EXPECT_EQ(w.count(), 2u);
    }

    SnapshotReader r(p);
    ASSERT_EQ(r.n_snapshots(), 2u);
    EXPECT_EQ(r.n_atoms(), kAtoms);

    const auto all = r.read_all();
    ASSERT_EQ(all.size(), 2u);

    EXPECT_FLOAT_EQ(all[0].score, -11.5f);
    EXPECT_EQ(all[0].generation, 3u);
    EXPECT_EQ(all[0].coords, ramp(kAtoms, 1.0f));

    EXPECT_FLOAT_EQ(all[1].score, -7.25f);
    EXPECT_EQ(all[1].generation, 9u);
    EXPECT_EQ(all[1].coords, ramp(kAtoms, 100.0f));
}

TEST_F(BinarySnapshotTest, RandomAccessMatchesSequentialRead) {
    const uint32_t kAtoms = 3;
    const std::string p = path("random.fasn");
    {
        SnapshotWriter w(p, kAtoms);
        for (uint32_t i = 0; i < 8; ++i)
            w.write_snapshot(-static_cast<float>(i), i, ramp(kAtoms, static_cast<float>(i)));
    }

    SnapshotReader r(p);
    const auto all = r.read_all();
    ASSERT_EQ(all.size(), 8u);
    for (uint32_t i = 0; i < 8; ++i) {
        const auto s = r.read_snapshot(i);
        EXPECT_FLOAT_EQ(s.score, all[i].score);
        EXPECT_EQ(s.generation, all[i].generation);
        EXPECT_EQ(s.coords, all[i].coords);
    }
}

TEST_F(BinarySnapshotTest, RawPointerOverloadMatchesVectorOverload) {
    const uint32_t kAtoms = 4;
    const auto coords = ramp(kAtoms, 2.0f);

    const std::string pv = path("vec.fasn");
    { SnapshotWriter w(pv, kAtoms); w.write_snapshot(-1.0f, 1, coords); }

    const std::string pp = path("ptr.fasn");
    { SnapshotWriter w(pp, kAtoms); w.write_snapshot(-1.0f, 1, coords.data(), coords.size()); }

    SnapshotReader rv(pv), rp(pp);
    EXPECT_EQ(rv.read_snapshot(0).coords, rp.read_snapshot(0).coords);
    EXPECT_EQ(fs::file_size(pv), fs::file_size(pp));
}

TEST_F(BinarySnapshotTest, FileSizeIsExactlyHeaderPlusNRecords) {
    const uint32_t kAtoms = 7;
    const uint32_t kN = 5;
    const std::string p = path("size.fasn");
    {
        SnapshotWriter w(p, kAtoms);
        for (uint32_t i = 0; i < kN; ++i)
            w.write_snapshot(0.0f, i, ramp(kAtoms, 0.0f));
    }
    const auto expected =
        flexaids::file_header_bytes() + kN * flexaids::snapshot_record_bytes(kAtoms);
    EXPECT_EQ(fs::file_size(p), expected);
}

TEST_F(BinarySnapshotTest, OnDiskHeaderCarriesMagicAndVersion) {
    const std::string p = path("magic.fasn");
    { SnapshotWriter w(p, 2); w.write_snapshot(-1.0f, 0, ramp(2, 0.0f)); }

    std::ifstream in(p, std::ios::binary);
    ASSERT_TRUE(in.good());
    BinaryFileHeader h{};
    in.read(reinterpret_cast<char*>(&h), sizeof(h));
    EXPECT_EQ(std::memcmp(h.magic, "FASN", 4), 0);
    EXPECT_EQ(flexaids::le_to_host_u16(h.version), 1u);
    EXPECT_EQ(flexaids::le_to_host_u16(h.flags), 0u);
    EXPECT_EQ(flexaids::le_to_host_u32(h.n_snapshots), 1u);
    EXPECT_EQ(flexaids::le_to_host_u32(h.n_atoms), 2u);
}

TEST_F(BinarySnapshotTest, FlushUpdatesTheHeaderCountBeforeClose) {
    const std::string p = path("flush.fasn");
    SnapshotWriter w(p, 2);
    w.write_snapshot(-1.0f, 0, ramp(2, 0.0f));
    w.write_snapshot(-2.0f, 1, ramp(2, 1.0f));
    w.flush();

    // Reader must see both records while the writer is still open.
    SnapshotReader r(p);
    EXPECT_EQ(r.n_snapshots(), 2u);
}

TEST_F(BinarySnapshotTest, ZeroSnapshotsProducesAValidEmptyFile) {
    const std::string p = path("empty.fasn");
    { SnapshotWriter w(p, 4); }

    SnapshotReader r(p);
    EXPECT_EQ(r.n_snapshots(), 0u);
    EXPECT_EQ(r.n_atoms(), 4u);
    EXPECT_TRUE(r.read_all().empty());
}

// --- error paths ----------------------------------------------------------

TEST_F(BinarySnapshotTest, WriterRejectsCoordinateLengthMismatch) {
    const std::string p = path("mismatch.fasn");
    SnapshotWriter w(p, 5);
    const std::vector<float> too_short(9, 0.0f);  // 3 atoms, not 5
    EXPECT_THROW(w.write_snapshot(-1.0f, 0, too_short), std::runtime_error);
    EXPECT_EQ(w.count(), 0u) << "a rejected write must not be counted";
}

TEST_F(BinarySnapshotTest, WriterThrowsOnUnopenablePath) {
    const std::string bad = (dir_ / "no_such_subdir" / "x.fasn").string();
    EXPECT_THROW(SnapshotWriter(bad, 1), std::runtime_error);
}

TEST_F(BinarySnapshotTest, ReaderThrowsOnMissingFile) {
    EXPECT_THROW(SnapshotReader(path("does_not_exist.fasn")), std::runtime_error);
}

TEST_F(BinarySnapshotTest, ReaderRejectsWrongMagic) {
    const std::string p = path("badmagic.fasn");
    {
        std::ofstream out(p, std::ios::binary);
        BinaryFileHeader h{};
        std::memcpy(h.magic, "XXXX", 4);
        h.version = flexaids::host_to_le_u16(1);
        h.flags = 0;
        h.n_snapshots = 0;
        h.n_atoms = flexaids::host_to_le_u32(1);
        out.write(reinterpret_cast<const char*>(&h), sizeof(h));
    }
    EXPECT_THROW(SnapshotReader(p), std::runtime_error);
}

TEST_F(BinarySnapshotTest, ReaderRejectsTruncatedHeader) {
    const std::string p = path("trunc.fasn");
    {
        std::ofstream out(p, std::ios::binary);
        out.write("FAS", 3);  // shorter than the 16-byte header
    }
    EXPECT_THROW(SnapshotReader(p), std::runtime_error);
}

TEST_F(BinarySnapshotTest, ReadSnapshotOutOfRangeThrows) {
    const std::string p = path("oor.fasn");
    { SnapshotWriter w(p, 2); w.write_snapshot(-1.0f, 0, ramp(2, 0.0f)); }

    SnapshotReader r(p);
    EXPECT_THROW(r.read_snapshot(1), std::out_of_range);
    EXPECT_THROW(r.read_snapshot(99), std::out_of_range);
    EXPECT_NO_THROW(r.read_snapshot(0));
}

// --- format sniffing ------------------------------------------------------

TEST_F(BinarySnapshotTest, IsBinarySnapshotDetectsFormat) {
    const std::string good = path("good.fasn");
    { SnapshotWriter w(good, 1); w.write_snapshot(0.0f, 0, ramp(1, 0.0f)); }
    EXPECT_TRUE(SnapshotReader::is_binary_snapshot(good));

    const std::string text = path("pose.pdb");
    { std::ofstream out(text); out << "ATOM      1  CA  ALA A   1\n"; }
    EXPECT_FALSE(SnapshotReader::is_binary_snapshot(text));

    EXPECT_FALSE(SnapshotReader::is_binary_snapshot(path("absent.fasn")));
}

// --- mmap view ------------------------------------------------------------

TEST_F(BinarySnapshotTest, MmapCoordinatesAgreeWithParsedCoordinates) {
    const uint32_t kAtoms = 6;
    const std::string p = path("mmap.fasn");
    {
        SnapshotWriter w(p, kAtoms);
        w.write_snapshot(-3.0f, 0, ramp(kAtoms, 5.0f));
        w.write_snapshot(-4.0f, 1, ramp(kAtoms, 50.0f));
    }

    SnapshotReader r(p);
    for (uint32_t i = 0; i < 2; ++i) {
        const float* mapped = r.mmap_coordinates(i);
        if (mapped == nullptr) GTEST_SKIP() << "mmap unavailable on this platform";
        const auto parsed = r.read_snapshot(i);
        for (size_t k = 0; k < parsed.coords.size(); ++k)
            EXPECT_FLOAT_EQ(mapped[k], parsed.coords[k]) << "snapshot " << i << " elem " << k;
    }
}

TEST_F(BinarySnapshotTest, MmapCoordinatesOutOfRangeReturnsNull) {
    const std::string p = path("mmap_oor.fasn");
    { SnapshotWriter w(p, 2); w.write_snapshot(0.0f, 0, ramp(2, 0.0f)); }
    SnapshotReader r(p);
    EXPECT_EQ(r.mmap_coordinates(5), nullptr);
}

// --- PDB projection -------------------------------------------------------

TEST_F(BinarySnapshotTest, SnapshotToPdbEmitsOneAtomLinePerAtom) {
    SnapshotReader::Snapshot s;
    s.score = -9.5f;
    s.generation = 2;
    s.coords = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f};

    const std::vector<std::string> names = {" CA ", " CB "};
    const std::vector<std::string> res = {"ALA", "ALA"};
    const std::vector<int> resnum = {1, 1};
    const std::vector<char> chains = {'A', 'A'};
    const std::vector<std::string> elems = {" C", " C"};
    const std::vector<int> atomnum = {1, 2};

    const std::string out_path = path("out.pdb");
    FILE* f = std::fopen(out_path.c_str(), "w");
    ASSERT_NE(f, nullptr);
    flexaids::snapshot_to_pdb(s, names, res, resnum, chains, elems, atomnum,
                              "TEST REMARK", f);
    std::fclose(f);

    std::ifstream in(out_path);
    std::string line;
    int atom_lines = 0;
    bool saw_remark = false;
    while (std::getline(in, line)) {
        if (line.rfind("ATOM", 0) == 0) ++atom_lines;
        if (line.find("TEST REMARK") != std::string::npos) saw_remark = true;
    }
    EXPECT_EQ(atom_lines, 2);
    EXPECT_TRUE(saw_remark);
}
