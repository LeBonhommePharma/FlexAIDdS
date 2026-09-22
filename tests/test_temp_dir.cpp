// tests/test_temp_dir.cpp — behaviour of the never-throwing temp-dir resolver
// Apache-2.0 © 2026 Le Bonhomme Pharma
//
// THE INCIDENT THESE TESTS ENCODE
// -------------------------------
// std::filesystem::temp_directory_path() aborts the process when $TMPDIR names
// a directory that has been deleted. On 2026-09-20 nine cells of an 85-target
// Astex campaign died that way: $TMPDIR pointed into an agent session
// workspace which its harness swept after a few hours idle, and every engine
// invocation started afterwards exited at startup with zero poses. The driver
// reported it only as "Incomplete docking", so it read as nine bad targets.
// It was one bad assumption: that a path which resolved at launch still exists
// when the work runs.
//
// The control test below is the important one. It asserts the OLD behaviour
// still fails, because a regression suite that only exercises the fix cannot
// tell you whether the scenario it guards against is still reachable.

#include <gtest/gtest.h>

#include "../LIB/temp_dir.h"

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>

namespace fs = std::filesystem;

namespace {

class TempDirEnv : public ::testing::Test {
protected:
    void SetUp() override {
        const char* t = std::getenv("TMPDIR");
        old_tmpdir_ = t ? t : "";
        const char* o = std::getenv("FLEXAIDDS_TMPDIR");
        old_override_ = o ? o : "";
    }
    void TearDown() override {
        if (old_tmpdir_.empty()) unsetenv("TMPDIR");
        else setenv("TMPDIR", old_tmpdir_.c_str(), 1);
        if (old_override_.empty()) unsetenv("FLEXAIDDS_TMPDIR");
        else setenv("FLEXAIDDS_TMPDIR", old_override_.c_str(), 1);
        flexaids::temp_dir_reset_for_testing();
    }
    std::string old_tmpdir_, old_override_;
};

fs::path unique_under_tmp(const char* stem) {
    return fs::path("/tmp") / (std::string(stem) + "_" + std::to_string(::getpid()));
}

}  // namespace

// ---------------------------------------------------------------- control --
// If this ever starts passing without a throw, the hazard has gone away at the
// library level and these guards can be revisited. Until then it documents
// that the danger is live.
TEST_F(TempDirEnv, ThrowingOverloadStillAbortsOnADeletedTmpdir) {
    const fs::path ghost = unique_under_tmp("flexaidds_ghost");
    fs::remove_all(ghost);
    setenv("TMPDIR", ghost.string().c_str(), 1);

    bool threw = false;
    try { (void)fs::temp_directory_path(); }
    catch (const fs::filesystem_error&) { threw = true; }
    EXPECT_TRUE(threw) << "the failure mode this file guards against is no longer reproducible";
}

// ------------------------------------------------------------------- fix ---
TEST_F(TempDirEnv, ResolverSurvivesADeletedTmpdir) {
    const fs::path ghost = unique_under_tmp("flexaidds_ghost");
    fs::remove_all(ghost);
    setenv("TMPDIR", ghost.string().c_str(), 1);
    unsetenv("FLEXAIDDS_TMPDIR");
    flexaids::temp_dir_reset_for_testing();

    EXPECT_TRUE(flexaids::temp_dir_ok());
    const std::string got = flexaids::temp_dir();
    EXPECT_NE(got, ghost.string());
    EXPECT_TRUE(fs::exists(got));
    EXPECT_TRUE(fs::is_directory(got));
}

TEST_F(TempDirEnv, ResolvedDirectoryActuallyAcceptsAWrite) {
    unsetenv("FLEXAIDDS_TMPDIR");
    setenv("TMPDIR", "/tmp", 1);
    flexaids::temp_dir_reset_for_testing();

    const fs::path f = fs::path(flexaids::temp_dir()) / "flexaidds_gtest_write_probe";
    std::FILE* fp = std::fopen(f.string().c_str(), "wb");
    ASSERT_NE(fp, nullptr) << "resolver returned a directory that rejects writes";
    std::fclose(fp);
    fs::remove(f);
}

TEST_F(TempDirEnv, ExplicitOverrideWinsAndIsRecordedInProvenance) {
    const fs::path durable = unique_under_tmp("flexaidds_durable");
    fs::create_directories(durable);
    setenv("FLEXAIDDS_TMPDIR", durable.string().c_str(), 1);
    flexaids::temp_dir_reset_for_testing();

    EXPECT_EQ(flexaids::temp_dir(), durable.string());
    EXPECT_EQ(flexaids::temp_dir_source(), "FLEXAIDDS_TMPDIR");
    fs::remove_all(durable);
}

TEST_F(TempDirEnv, NonexistentOverrideFallsThroughRatherThanFailing) {
    setenv("FLEXAIDDS_TMPDIR", "/no/such/path/anywhere", 1);
    setenv("TMPDIR", "/tmp", 1);
    flexaids::temp_dir_reset_for_testing();

    EXPECT_TRUE(flexaids::temp_dir_ok());
    EXPECT_NE(flexaids::temp_dir(), "/no/such/path/anywhere");
}

// The distinguishing property versus an access()/stat() check: a directory can
// exist, report as a directory, and still reject every write.
TEST_F(TempDirEnv, ExistingButUnwritableDirectoryIsRejectedByTheWriteProbe) {
    if (::geteuid() == 0) GTEST_SKIP() << "root ignores mode bits";
    const fs::path ro = unique_under_tmp("flexaidds_ro");
    fs::remove_all(ro);
    fs::create_directories(ro);
    fs::permissions(ro, fs::perms::owner_read | fs::perms::owner_exec,
                    fs::perm_options::replace);

    setenv("FLEXAIDDS_TMPDIR", ro.string().c_str(), 1);
    setenv("TMPDIR", "/tmp", 1);
    flexaids::temp_dir_reset_for_testing();

    EXPECT_NE(flexaids::temp_dir(), ro.string())
        << "an unwritable directory was accepted; the probe is not probing";
    EXPECT_TRUE(flexaids::temp_dir_ok());

    fs::permissions(ro, fs::perms::owner_all, fs::perm_options::replace);
    fs::remove_all(ro);
}

// A long docking run must not silently migrate its scratch between restarts.
TEST_F(TempDirEnv, ResolutionIsCachedAndStableMidRun) {
    unsetenv("FLEXAIDDS_TMPDIR");
    setenv("TMPDIR", "/tmp", 1);
    flexaids::temp_dir_reset_for_testing();
    const std::string first = flexaids::temp_dir();

    setenv("TMPDIR", "/var/tmp", 1);          // change AFTER resolution
    EXPECT_EQ(flexaids::temp_dir(), first);
}
