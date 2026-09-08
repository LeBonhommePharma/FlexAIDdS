// CORE-2 / CORE-6: exercise the production voronoi_poly2 retry and exit paths.
// The test-only hook raises the edge-limit trigger on a small bounded hull;
// plane generation, jitter, RESTART, successful completion and error return are
// the actual LIB/Vcontacts.cpp implementation, not a mirror of its guard logic.
// Apache-2.0

#include <gtest/gtest.h>
#include "Vcontacts.h"
#include "RngSeed.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <tuple>
#include <vector>

#ifndef FLEXAIDS_VCONTACTS_TEST_HOOKS
#error "This regression must link production Vcontacts.cpp with its test hook enabled"
#endif

namespace {

class ScopedEnv {
public:
    ScopedEnv(const char* name, const char* value) : name_(name) {
        const char* previous = std::getenv(name);
        had_previous_ = previous != nullptr;
        if (previous) previous_ = previous;
        set(value);
    }
    ~ScopedEnv() { set(had_previous_ ? previous_.c_str() : nullptr); }
    ScopedEnv(const ScopedEnv&) = delete;
    ScopedEnv& operator=(const ScopedEnv&) = delete;
    void set(const char* value) {
#ifdef _WIN32
        _putenv_s(name_.c_str(), value ? value : "");
#else
        if (value) setenv(name_.c_str(), value, 1);
        else unsetenv(name_.c_str());
#endif
    }
private:
    std::string name_;
    std::string previous_;
    bool had_previous_ = false;
};

struct Hull {
    static constexpr int contact_count = 6;
    std::array<atom, 7> atoms{};
    std::array<atomsas, 7> calc{};
    std::array<contactlist, contact_count> contacts{};
    std::array<plane, contact_count + 4> planes{};
    std::array<vertex, contact_count> centers{};
    std::vector<vertex> poly = std::vector<vertex>(MAX_POLY);
    std::vector<edgevector> edges = std::vector<edgevector>(MAX_POLY);
    std::array<int, 21> seeds{};
    std::array<std::array<float, 3>, 7> pristine{};
    VC_Global vc{};

    Hull() {
        // Deliberately away from zero: the old restore-to-origin defect cannot
        // pass accidentally. Six neighbours bound a cube with half-width 1.
        const float origin[3] = {17.25f, -9.5f, 31.75f};
        const float offsets[7][3] = {
            {0, 0, 0}, {2, 0, 0}, {-2, 0, 0}, {0, 2, 0},
            {0, -2, 0}, {0, 0, 2}, {0, 0, -2}
        };
        for (int i = 0; i < 7; ++i) {
            atoms[i].number = 90000 + i;
            atoms[i].radius = 1.0f;
            calc[i].atom = &atoms[i];
            for (int axis = 0; axis < 3; ++axis) {
                atoms[i].coor[axis] = origin[axis] + offsets[i][axis];
                pristine[i][axis] = atoms[i].coor[axis];
            }
        }
        for (int i = 0; i < contact_count; ++i) {
            contacts[i].index = i + 1;
            contacts[i].dist = 2.0;
        }
        seeds.fill(-1);
        vc.Calc = calc.data();
        vc.poly = poly.data();
        vc.vedge = edges.data();
        vc.centerpt = centers.data();
        vc.seed = seeds.data();
        vc.planedef = 'B';
        vc.recalc = 1;
    }

    int run(int forced_failures, bool allow_retry = true) {
        auto& trace = flexaids_vct_test::trace();
        trace = {};
        trace.force_failsafe_passes = forced_failures;
        vc.recalc = allow_retry ? 1 : 0;
        return voronoi_poly2(&vc, 0, planes.data(), 2.0f,
                             contact_count, contacts.data());
    }

    void expect_input_preserved() const {
        for (int i = 0; i < 7; ++i) {
            EXPECT_EQ(std::memcmp(atoms[i].coor, pristine[i].data(),
                                  sizeof(atoms[i].coor)), 0) << "atom " << i;
        }
    }
};

using Gates = std::tuple<bool, bool, bool>;  // local, eager guard, keyed jitter
class VcontactsFailsafe : public ::testing::TestWithParam<Gates> {
protected:
    ScopedEnv local{"FLEXAIDDS_VCT_LOCAL_PERTURB", "0"};
    ScopedEnv guard{"FLEXAIDDS_VCT_COORD_GUARD", "0"};
    ScopedEnv diag{"FLEXAIDDS_VCT_FAILSAFE_DIAG", "0"};
    ScopedEnv keyed{"FLEXAIDDS_VORONOI_KEYED_JITTER", "0"};
    ScopedEnv stream_fix{"FLEXAIDDS_RNG_STREAM_FIX", "0"};

    void SetUp() override {
        local.set(std::get<0>(GetParam()) ? "1" : "0");
        guard.set(std::get<1>(GetParam()) ? "1" : "0");
        keyed.set(std::get<2>(GetParam()) ? "1" : "0");
        flexaids_rng::set_master_seed(12345);
    }

    void expect_real_perturbation(const Hull& hull) {
        const auto& trace = flexaids_vct_test::trace();
        EXPECT_EQ(trace.local_mode, std::get<0>(GetParam()));
        EXPECT_EQ(trace.coord_guard_requested, std::get<1>(GetParam()));
        EXPECT_NE(std::memcmp(trace.working_after_perturb,
                              hull.pristine[0].data(), 3 * sizeof(float)), 0);
        if (std::get<0>(GetParam())) {
            // Check inside the failsafe as well as after return: local mode
            // must not rely on an eventual guard to undo an illegal write.
            EXPECT_EQ(std::memcmp(trace.input_after_perturb,
                                  hull.pristine[0].data(), 3 * sizeof(float)), 0);
        } else {
            EXPECT_EQ(std::memcmp(trace.input_after_perturb,
                                  trace.working_after_perturb, 3 * sizeof(float)), 0);
        }
    }
};

TEST_P(VcontactsFailsafe, PristineCoordinatesAfterOneAndRepeatedRetries) {
    for (int failures : {1, 3}) {
        SCOPED_TRACE(failures);
        Hull hull;
        ASSERT_GT(hull.run(failures), 0);
        const auto& trace = flexaids_vct_test::trace();
        EXPECT_EQ(trace.force_failsafe_passes, 0);
        EXPECT_EQ(trace.failsafe_entries, failures);
        EXPECT_EQ(trace.hull_passes, failures + 1);
        EXPECT_EQ(trace.success_exits, 1);
        EXPECT_EQ(trace.error_exits, 0);
        expect_real_perturbation(hull);
        hull.expect_input_preserved();
    }
}

TEST_P(VcontactsFailsafe, PristineCoordinatesOnFailsafeErrorExit) {
    Hull hull;
    ASSERT_EQ(hull.run(1, false), -1);
    const auto& trace = flexaids_vct_test::trace();
    EXPECT_EQ(trace.force_failsafe_passes, 0);
    EXPECT_EQ(trace.failsafe_entries, 1);
    EXPECT_EQ(trace.hull_passes, 1);
    EXPECT_EQ(trace.success_exits, 0);
    EXPECT_EQ(trace.error_exits, 1);
    expect_real_perturbation(hull);
    hull.expect_input_preserved();
}

TEST_P(VcontactsFailsafe, ValidHullHasAnalyticalCubeVerticesWithoutRetries) {
    Hull hull;
    ASSERT_EQ(hull.run(0), 8);
    EXPECT_EQ(flexaids_vct_test::trace().failsafe_entries, 0);
    EXPECT_EQ(flexaids_vct_test::trace().hull_passes, 1);
    unsigned corners = 0;
    for (int i = 0; i < 8; ++i) {
        int corner = 0;
        for (int axis = 0; axis < 3; ++axis) {
            EXPECT_DOUBLE_EQ(std::abs(hull.poly[i].xi[axis]), 1.0);
            if (hull.poly[i].xi[axis] > 0.0) corner |= 1 << axis;
        }
        corners |= 1u << corner;
    }
    EXPECT_EQ(corners, 0xffu);
    hull.expect_input_preserved();
}

INSTANTIATE_TEST_SUITE_P(AllGateCombinations, VcontactsFailsafe,
                        ::testing::Combine(::testing::Bool(), ::testing::Bool(),
                                           ::testing::Bool()));

TEST(VcontactsFailsafeFlags, CanonicalTrueFalseAndFallbackSpellingsReachProduction) {
    ScopedEnv local("FLEXAIDDS_VCT_LOCAL_PERTURB", nullptr);
    ScopedEnv guard("FLEXAIDDS_VCT_COORD_GUARD", nullptr);
    ScopedEnv diag("FLEXAIDDS_VCT_FAILSAFE_DIAG", nullptr);
    ScopedEnv keyed("FLEXAIDDS_VORONOI_KEYED_JITTER", "1");
    const std::pair<const char*, bool> values[] = {
        {"1", true}, {"true", true}, {"yes", true}, {"on", true},
        {"TRUE", true}, {"YeS", true}, {"ON", true}, {"  true  ", true},
        {"0", false}, {"false", false}, {"no", false}, {"off", false},
        {"FALSE", false}, {"No", false}, {"OFF", false}, {"  off  ", false},
        {nullptr, false}, {"", false}, {"  ", false}, {"garbage", false},
        {"2", false}
    };
    for (const auto& value : values) {
        SCOPED_TRACE(value.first ? value.first : "unset");
        local.set(value.first);
        guard.set(value.first);
        diag.set(value.first);
        flexaids_rng::set_master_seed(12345);
        Hull hull;
        ::testing::internal::CaptureStderr();
        const int result = hull.run(1, false);
        const auto output = ::testing::internal::GetCapturedStderr();
        ASSERT_EQ(result, -1);
        const auto& trace = flexaids_vct_test::trace();
        EXPECT_EQ(trace.local_mode, value.second);
        EXPECT_EQ(trace.coord_guard_requested, value.second);
        EXPECT_EQ(trace.diagnostics_enabled, value.second);
        EXPECT_EQ(output.find("[VCT-FAILSAFE]") != std::string::npos, value.second);
        hull.expect_input_preserved();
    }
}

TEST(VcontactsFailsafeFlags, LocalDiagnosticsReportWorkingCopyDisplacement) {
    ScopedEnv local("FLEXAIDDS_VCT_LOCAL_PERTURB", "on");
    ScopedEnv guard("FLEXAIDDS_VCT_COORD_GUARD", "off");
    ScopedEnv diag("FLEXAIDDS_VCT_FAILSAFE_DIAG", "yes");
    ScopedEnv keyed("FLEXAIDDS_VORONOI_KEYED_JITTER", "1");
    flexaids_rng::set_master_seed(12345);
    Hull hull;
    ::testing::internal::CaptureStderr();
    const int result = hull.run(1, false);
    const auto output = ::testing::internal::GetCapturedStderr();
    ASSERT_EQ(result, -1);
    const auto& trace = flexaids_vct_test::trace();
    char expected[128];
    std::snprintf(expected, sizeof(expected), "d=(%.6f,%.6f,%.6f)",
                  trace.working_after_perturb[0] - hull.pristine[0][0],
                  trace.working_after_perturb[1] - hull.pristine[0][1],
                  trace.working_after_perturb[2] - hull.pristine[0][2]);
    EXPECT_NE(output.find(expected), std::string::npos);
    EXPECT_EQ(output.find("d=(0.000000,0.000000,0.000000)"), std::string::npos);
    hull.expect_input_preserved();
}

}  // namespace
