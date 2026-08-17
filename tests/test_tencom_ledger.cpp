// tests/test_tencom_ledger.cpp — ledger-only tENCoM λ (Wave 4)
//
// Flag default OFF. When ON, λ comes from TorsionalENM (real channel) and is
// tagged inert_on_election. Ranking/fitness APIs are not in this TU.

#include <gtest/gtest.h>

#include "tencom_ledger.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

class ScopedEnv {
public:
    ScopedEnv(const char* name, const char* value) : name_(name) {
        const char* prev = std::getenv(name);
        if (prev) {
            had_prev_ = true;
            prev_ = prev;
        }
#if defined(_WIN32)
        if (value) _putenv_s(name, value);
        else _putenv_s(name, "");
#else
        if (value) setenv(name, value, 1);
        else unsetenv(name);
#endif
    }
    ~ScopedEnv() {
#if defined(_WIN32)
        if (had_prev_) _putenv_s(name_.c_str(), prev_.c_str());
        else _putenv_s(name_.c_str(), "");
#else
        if (had_prev_) setenv(name_.c_str(), prev_.c_str(), 1);
        else unsetenv(name_.c_str());
#endif
    }
    ScopedEnv(const ScopedEnv&) = delete;
    ScopedEnv& operator=(const ScopedEnv&) = delete;

private:
    std::string name_;
    std::string prev_;
    bool had_prev_ = false;
};

std::vector<std::array<float, 3>> helix_ca(int n, float radius = 2.3f, float rise = 1.5f) {
    std::vector<std::array<float, 3>> ca;
    ca.reserve(static_cast<size_t>(n));
    const float turn = 100.0f * 3.14159265f / 180.0f;
    for (int r = 0; r < n; ++r) {
        ca.push_back({
            radius * std::cos(static_cast<float>(r) * turn),
            radius * std::sin(static_cast<float>(r) * turn),
            static_cast<float>(r) * rise,
        });
    }
    return ca;
}

}  // namespace

TEST(TencomLedger, FlagDefaultOff) {
    ScopedEnv off("FLEXAIDDS_LEDGER_TENCOM_LAMBDA", nullptr);
    EXPECT_FALSE(flexaids::ledger_tencom_lambda_enabled());
    const auto rec = flexaids::collect_tencom_lambda_from_ca(helix_ca(20));
    EXPECT_FALSE(rec.flag_enabled);
    EXPECT_EQ(rec.status, "flag_off");
    EXPECT_TRUE(rec.lambdas.empty());
    EXPECT_TRUE(rec.inert_on_election);
    EXPECT_TRUE(flexaids::format_tencom_lambda_remark(rec).empty());
}

TEST(TencomLedger, FlagOnWritesRealSpectrumAndInertTag) {
    ScopedEnv on("FLEXAIDDS_LEDGER_TENCOM_LAMBDA", "1");
    ASSERT_TRUE(flexaids::ledger_tencom_lambda_enabled());
    const auto rec = flexaids::collect_tencom_lambda_from_ca(helix_ca(20));
    EXPECT_TRUE(rec.flag_enabled);
    EXPECT_TRUE(rec.inert_on_election);
    EXPECT_TRUE(flexaids::kTencomLambdaInertOnElection);
    ASSERT_EQ(rec.status, "ok") << "tENCoM channel must supply real eigenvalues";
    ASSERT_GT(rec.n, 0);
    ASSERT_EQ(rec.n, static_cast<int>(rec.lambdas.size()));
    for (double lam : rec.lambdas) {
        EXPECT_TRUE(std::isfinite(lam));
    }
    EXPECT_LE(rec.lambda_min, rec.lambda_max);
    const std::string remark = flexaids::format_tencom_lambda_remark(rec);
    EXPECT_NE(remark.find("tencom_lambda_ledger"), std::string::npos);
    EXPECT_NE(remark.find("inert_on_election=1"), std::string::npos);
    EXPECT_NE(remark.find("status=ok"), std::string::npos);
}

TEST(TencomLedger, FlagOnTooFewCaIsUnavailableStillInert) {
    ScopedEnv on("FLEXAIDDS_LEDGER_TENCOM_LAMBDA", "1");
    const auto rec = flexaids::collect_tencom_lambda_from_ca(helix_ca(2));
    EXPECT_TRUE(rec.flag_enabled);
    EXPECT_EQ(rec.status, "unavailable");
    EXPECT_TRUE(rec.lambdas.empty());
    EXPECT_TRUE(rec.inert_on_election);
    const std::string remark = flexaids::format_tencom_lambda_remark(rec);
    EXPECT_NE(remark.find("inert_on_election=1"), std::string::npos);
    EXPECT_NE(remark.find("status=unavailable"), std::string::npos);
}

TEST(TencomLedger, FlagOnVersusOffDoesNotInventEigenvaluesFromAtomEigen) {
    // The forbidden path: treating atom::eigen x-components as λ. This TU never
    // reads atom::eigen. Flag off yields an empty spectrum even for a helix.
    ScopedEnv off("FLEXAIDDS_LEDGER_TENCOM_LAMBDA", "0");
    const auto rec = flexaids::collect_tencom_lambda_from_ca(helix_ca(20));
    EXPECT_TRUE(rec.lambdas.empty());
    EXPECT_EQ(rec.status, "flag_off");
}

TEST(TencomLedger, InertOnElectionIsACompileTimeContract) {
    EXPECT_TRUE(flexaids::kTencomLambdaInertOnElection);
}
