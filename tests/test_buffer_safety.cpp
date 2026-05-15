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
