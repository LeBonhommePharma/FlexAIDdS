// Unit tests for Wave 3.4 memetic env gate (drives shipped LIB/memetic_gate.h).
#include "memetic_gate.h"

#include <cstdlib>
#include <gtest/gtest.h>
#include <string>

namespace {

void clear_env(const char* key) {
#if defined(_WIN32)
    _putenv_s(key, "");
#else
    unsetenv(key);
#endif
}

void set_env(const char* key, const char* val) {
#if defined(_WIN32)
    _putenv_s(key, val);
#else
    setenv(key, val, 1);
#endif
}

class MemeticGateEnv : public ::testing::Test {
protected:
    void SetUp() override {
        clear_env("FLEXAIDDS_MEMETIC");
        clear_env("FLEXAIDDS_WALL_PILOT_PASS");
    }
    void TearDown() override {
        clear_env("FLEXAIDDS_MEMETIC");
        clear_env("FLEXAIDDS_WALL_PILOT_PASS");
    }
};

}  // namespace

TEST_F(MemeticGateEnv, DefaultOff) {
    EXPECT_EQ(flexaids::resolve_use_memetic_from_env(), 0);
}

TEST_F(MemeticGateEnv, MemeticAloneDoesNotEnable) {
    set_env("FLEXAIDDS_MEMETIC", "1");
    EXPECT_EQ(flexaids::resolve_use_memetic_from_env(), 0);
}

TEST_F(MemeticGateEnv, WallPassAloneDoesNotEnable) {
    set_env("FLEXAIDDS_WALL_PILOT_PASS", "1");
    EXPECT_EQ(flexaids::resolve_use_memetic_from_env(), 0);
}

TEST_F(MemeticGateEnv, BothEnables) {
    set_env("FLEXAIDDS_MEMETIC", "1");
    set_env("FLEXAIDDS_WALL_PILOT_PASS", "1");
    EXPECT_EQ(flexaids::resolve_use_memetic_from_env(), 1);
}

TEST_F(MemeticGateEnv, ZeroMeansOff) {
    set_env("FLEXAIDDS_MEMETIC", "0");
    set_env("FLEXAIDDS_WALL_PILOT_PASS", "1");
    EXPECT_EQ(flexaids::resolve_use_memetic_from_env(), 0);
}
