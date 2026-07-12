// test_buffer_safety.cpp — regression tests for bounded legacy config parsing.

#include <gtest/gtest.h>

#include "read_input_utils.h"

#include <cstring>
#include <string>

TEST(BufferSafety, LongConfigValueIsTruncatedAndTerminated) {
    char dest[16];
    std::memset(dest, 'X', sizeof(dest));

    std::string line = "GISTDG " + std::string(256, 'a');
    flexaids_copy_config_value(line.c_str(), dest, sizeof(dest));

    EXPECT_EQ(std::strlen(dest), sizeof(dest) - 1);
    EXPECT_EQ(dest[sizeof(dest) - 1], '\0');
}

TEST(BufferSafety, ConfigValueTrimsTrailingNewline) {
    char dest[32];
    flexaids_copy_config_value("GISTDN /tmp/rho.dx\n", dest, sizeof(dest));

    EXPECT_STREQ(dest, "/tmp/rho.dx");
}

TEST(BufferSafety, NullOrEmptyInputsStayTerminated) {
    char dest[8];
    std::memset(dest, 'X', sizeof(dest));

    flexaids_copy_config_value(nullptr, dest, sizeof(dest));
    EXPECT_EQ(dest[0], '\0');

    flexaids_copy_config_value("GISTDG", dest, sizeof(dest));
    EXPECT_EQ(dest[0], '\0');
}

TEST(BufferSafety, ZeroDestSizeIsNoOp) {
    char dest[4] = {'A', 'B', 'C', 'D'};
    flexaids_copy_config_value("GISTDG /tmp/x.dx", dest, 0);
    // With dest_size == 0 the helper must not write.
    EXPECT_EQ(dest[0], 'A');
    EXPECT_EQ(dest[3], 'D');
}

TEST(BufferSafety, TrimsCarriageReturnAndNewline) {
    char dest[64];
    flexaids_copy_config_value("GISTDG /data/gist.dx\r\n", dest, sizeof(dest));
    EXPECT_STREQ(dest, "/data/gist.dx");
}

TEST(BufferSafety, MultipleSpacesBetweenKeyAndValue) {
    char dest[64];
    flexaids_copy_config_value("GISTDN    /tmp/rho.dx", dest, sizeof(dest));
    EXPECT_STREQ(dest, "/tmp/rho.dx");
}

TEST(BufferSafety, ExactFitValueIsTerminated) {
    // dest can hold 7 chars + NUL; value is exactly 7 chars.
    char dest[8];
    flexaids_copy_config_value("KEY abcdefg", dest, sizeof(dest));
    EXPECT_STREQ(dest, "abcdefg");
    EXPECT_EQ(dest[7], '\0');
}
