// tests/test_svib_zero_mode_stability.cpp — the zero-mode classification must not
// depend on the SIGN of a numerically-zero eigenvalue.
//
// THE DEFECT THIS GUARDS. tencm::TorsionalENM::vibrational_eigenvalues() used
// `eigenvalue > cutoff` with cutoff == 0.0 on the TORSIONAL basis, i.e. a bare
// sign test. A mathematically exact zero mode is emitted by an iterative
// self-adjoint solver as +/-epsilon with a sign set by the rounding path, so the
// verdict was noise. The exposure is not theoretical: the dS_vib path
// (ic2cf.cpp) passes eigenvalue_cutoff = 0.0 to compute_vibrational_entropy
// precisely so that vibrational_eigenvalues() is the ONLY filter, and
// cf->minus_T_dsvib enters get_cf_evalue with unit coefficient.
//
// WHY THIS TEST CAN RUN ON ONE MACHINE. Reproducing the failure honestly would
// need a second architecture / a different Eigen build. Instead these tests
// SIMULATE it deterministically: build a real spectrum, then perturb every
// eigenvalue that is zero to within the threshold (sign flip, ULP nudge, and a
// magnitude rescale that stays under the threshold) and assert that the dropped
// count, the kept spectrum and S_vib are unchanged. classify_spectrum() is a
// static function over a plain vector<double> so a perturbed spectrum can be fed
// in directly -- no second machine, no mocking of the solver.
//
// ABLE TO FAIL: with the pre-fix rule restored (which the LegacyGate test below
// exercises on purpose), DroppedCountIsInvariantUnderZeroModeSignFlip and
// SvibIsInvariantUnderZeroModeSignFlip both go red on the torsional basis.
//
// Apache-2.0 (c) 2026 Le Bonhomme Pharma.

#include <gtest/gtest.h>

#include "tENCoM/tencm.h"
#include "encom.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace {

using Basis = tencm::TorsionalENM::Basis;

std::vector<std::array<float,3>> helix_ca(int n) {
    std::vector<std::array<float,3>> ca;
    ca.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        const double t = i * 100.0 * M_PI / 180.0;
        ca.push_back({static_cast<float>(2.3 * std::cos(t)),
                      static_cast<float>(2.3 * std::sin(t)),
                      static_cast<float>(1.5 * i)});
    }
    return ca;
}

std::vector<atom> carbon_atoms(const std::vector<std::array<float,3>>& xyz) {
    std::vector<atom> v(xyz.size());
    for (size_t i = 0; i < xyz.size(); ++i) {
        std::memset(&v[i], 0, sizeof(atom));
        v[i].coor[0] = xyz[i][0];
        v[i].coor[1] = xyz[i][1];
        v[i].coor[2] = xyz[i][2];
        std::strncpy(v[i].element, "C", sizeof(v[i].element) - 1);
    }
    return v;
}

constexpr float kCut = 9.0f;
constexpr float kK0  = 1.0f;

struct Built {
    tencm::TorsionalENM enm;
    std::vector<atom>   backing;
};

void build_path(Built& b, Basis basis,
                const std::vector<std::array<float,3>>& xyz, float k0 = kK0) {
    if (basis == Basis::Torsional) {
        b.enm.build_from_ca(xyz, kCut, k0);
    } else {
        b.backing = carbon_atoms(xyz);
        b.enm.build_from_ligand(b.backing.data(), 0,
                                static_cast<int>(b.backing.size()), kCut, k0);
    }
}

std::vector<double> spectrum_of(const tencm::TorsionalENM& enm) {
    std::vector<double> ev;
    for (const auto& m : enm.modes()) ev.push_back(m.eigenvalue);
    std::sort(ev.begin(), ev.end());
    return ev;
}

const char* basis_name(Basis b) {
    return b == Basis::Torsional ? "torsional" : "cartesian";
}

// The simulated other architecture. Every eigenvalue that is zero to within
// `cutoff` is replaced by a DIFFERENT representation of the same zero: opposite
// sign, one ULP away, and (for `scale`) a different magnitude that is still
// under the threshold. A correct classifier cannot tell these apart.
std::vector<double> perturb_zero_modes(const std::vector<double>& ev,
                                       double cutoff, double scale = 1.0) {
    std::vector<double> out = ev;
    for (double& x : out) {
        if (std::fabs(x) > cutoff) continue;
        double flipped = -x * scale;
        // The rescale must stay INSIDE the zero band, otherwise the mode really
        // has become a different mode and reclassifying it would be correct --
        // the test would then be asserting something false.
        if (std::fabs(flipped) > cutoff)
            flipped = std::copysign(cutoff, flipped);
        x = std::nextafter(flipped, flipped >= 0.0 ? 1.0 : -1.0);
        if (std::fabs(x) > cutoff) x = flipped;   // nextafter must not escape either
    }
    std::sort(out.begin(), out.end());
    return out;
}

// The perturbation band must NOT be read back from classify_spectrum(), because
// under the legacy gate that cutoff is 0.0 and the perturbation would collapse to
// a no-op -- making these tests pass VACUOUSLY on the very code path they exist
// to condemn. Compute the scale-aware band directly so the perturbation is real
// under either rule, and the tests go genuinely red with the fix reverted.
double zero_band(const std::vector<double>& ev) {
    double lam_max = 0.0;
    for (double x : ev) if (std::isfinite(x)) lam_max = std::max(lam_max, x);
    return tencm::TorsionalENM::zero_mode_cutoff(lam_max);
}

double svib_production(const std::vector<double>& kept) {
    // Mirrors the dS_vib path: eigenvalue_cutoff = 0.0, so the classification
    // performed upstream is the only filter that applies.
    std::vector<encom::NormalMode> modes;
    modes.reserve(kept.size());
    int idx = 1;
    for (double x : kept) {
        encom::NormalMode m;
        m.index = idx++;
        m.eigenvalue = x;
        m.frequency = x > 0.0 ? std::sqrt(x) : 0.0;
        modes.push_back(m);
    }
    return encom::ENCoMEngine::compute_vibrational_entropy(modes, 300.0, 0.0)
               .S_vib_kcal_mol_K;
}

// RAII for the legacy gate so a failing assertion cannot leak it into the rest
// of the suite (getenv state is process-global).
struct ScopedEnv {
    const char* name;
    explicit ScopedEnv(const char* n, const char* v) : name(n) { ::setenv(n, v, 1); }
    ~ScopedEnv() { ::unsetenv(name); }
};

} // namespace

TEST(SvibZeroModeStability, ZeroModesExistOnBothBases) {
    // Guards the premise: if no spectrum carried a numerically-zero mode the
    // invariance tests below would be vacuously true.
    for (Basis basis : {Basis::Torsional, Basis::Cartesian}) {
        Built b;
        build_path(b, basis, helix_ca(24));
        ASSERT_TRUE(b.enm.is_built()) << basis_name(basis);
        const auto ev = spectrum_of(b.enm);
        const auto cls = tencm::TorsionalENM::classify_spectrum(ev, basis);
        EXPECT_GT(cls.n_zero, 0)
            << basis_name(basis)
            << ": no sub-threshold mode, so the perturbation tests are vacuous";
        EXPECT_GT(cls.cutoff, 0.0)
            << basis_name(basis) << ": a cutoff of 0.0 IS the sign test";
    }
}

TEST(SvibZeroModeStability, DroppedCountIsInvariantUnderZeroModeSignFlip) {
    for (Basis basis : {Basis::Torsional, Basis::Cartesian}) {
        for (int n : {12, 18, 20, 24, 32}) {
            Built b;
            build_path(b, basis, helix_ca(n));
            ASSERT_TRUE(b.enm.is_built()) << basis_name(basis) << " n=" << n;
            const auto ev  = spectrum_of(b.enm);
            const auto ref = tencm::TorsionalENM::classify_spectrum(ev, basis);

            for (double scale : {1.0, 10.0, 0.1}) {
                const auto pert = perturb_zero_modes(ev, zero_band(ev), scale);
                const auto got  = tencm::TorsionalENM::classify_spectrum(pert, basis);
                EXPECT_EQ(got.n_zero, ref.n_zero)
                    << basis_name(basis) << " n=" << n << " scale=" << scale
                    << ": zero-mode count changed when a zero changed SIGN";
                EXPECT_EQ(got.vibrational.size(), ref.vibrational.size())
                    << basis_name(basis) << " n=" << n << " scale=" << scale
                    << ": kept-mode count depends on the sign of a zero";
            }
        }
    }
}

TEST(SvibZeroModeStability, KeptSpectrumIsInvariantUnderZeroModeSignFlip) {
    for (Basis basis : {Basis::Torsional, Basis::Cartesian}) {
        Built b;
        build_path(b, basis, helix_ca(24));
        const auto ev   = spectrum_of(b.enm);
        const auto ref  = tencm::TorsionalENM::classify_spectrum(ev, basis);
        const auto pert = perturb_zero_modes(ev, zero_band(ev));
        const auto got  = tencm::TorsionalENM::classify_spectrum(pert, basis);

        ASSERT_EQ(got.vibrational.size(), ref.vibrational.size()) << basis_name(basis);
        for (size_t i = 0; i < ref.vibrational.size(); ++i)
            EXPECT_DOUBLE_EQ(got.vibrational[i], ref.vibrational[i])
                << basis_name(basis) << " kept mode " << i << " changed";
    }
}

TEST(SvibZeroModeStability, SvibIsInvariantUnderZeroModeSignFlip) {
    for (Basis basis : {Basis::Torsional, Basis::Cartesian}) {
        for (int n : {12, 18, 24}) {
            Built b;
            build_path(b, basis, helix_ca(n));
            const auto ev   = spectrum_of(b.enm);
            const auto ref  = tencm::TorsionalENM::classify_spectrum(ev, basis);
            const auto pert = perturb_zero_modes(ev, zero_band(ev));
            const auto got  = tencm::TorsionalENM::classify_spectrum(pert, basis);

            EXPECT_DOUBLE_EQ(svib_production(got.vibrational),
                             svib_production(ref.vibrational))
                << basis_name(basis) << " n=" << n
                << ": S_vib moved because a numerically-zero mode changed sign";
        }
    }
}

// The real-input version of the same statement, and the one that reproduces the
// observed cross-platform failure mode without a second platform: perturb the
// INPUT by one float ULP and re-diagonalise. Pre-fix this changed the dropped
// count for 11 of 24 single-coordinate nudges on helix24.
TEST(SvibZeroModeStability, SingleUlpInputPerturbationDoesNotChangeDroppedCount) {
    const int n = 24;
    Built ref;
    build_path(ref, Basis::Torsional, helix_ca(n));
    ASSERT_TRUE(ref.enm.is_built());
    int ref_dropped = -1, ref_expected = -1;
    (void)ref.enm.vibrational_eigenvalues(&ref_dropped, &ref_expected);

    for (int which = 0; which < n; ++which) {
        auto ca = helix_ca(n);
        ca[static_cast<size_t>(which)][0] =
            std::nextafterf(ca[static_cast<size_t>(which)][0], 1e30f);
        Built b;
        build_path(b, Basis::Torsional, ca);
        ASSERT_TRUE(b.enm.is_built()) << "atom " << which;
        int dropped = -1, expected = -1;
        (void)b.enm.vibrational_eigenvalues(&dropped, &expected);
        EXPECT_EQ(dropped, ref_dropped)
            << "a ONE-ULP nudge of atom " << which
            << " reclassified a mode: the zero-mode verdict is rounding noise";
    }
}

// The cutoff must not be so wide that real physics is discarded. A genuinely
// soft mode sits DECADES above the threshold and must survive.
TEST(SvibZeroModeStability, GenuinelySoftModeIsNotDiscarded) {
    for (Basis basis : {Basis::Torsional, Basis::Cartesian}) {
        Built b;
        build_path(b, basis, helix_ca(24));
        auto ev = spectrum_of(b.enm);
        const auto base = tencm::TorsionalENM::classify_spectrum(ev, basis);

        // 1000x the zero threshold is still ~5 decades below the softest real
        // mode, so it is unambiguously "soft but real".
        const double soft = zero_band(ev) * 1000.0;
        ev.push_back(soft);
        std::sort(ev.begin(), ev.end());
        const auto got = tencm::TorsionalENM::classify_spectrum(ev, basis);

        EXPECT_EQ(got.vibrational.size(), base.vibrational.size() + 1)
            << basis_name(basis) << ": a soft mode at 1000x cutoff was dropped";
        EXPECT_EQ(got.n_zero, base.n_zero) << basis_name(basis);
    }
}

// A mode below -cutoff is an imaginary frequency, not round-off. It must be
// counted separately and must never reach log().
TEST(SvibZeroModeStability, GenuinelyNegativeModeIsCountedNotKept) {
    Built b;
    build_path(b, Basis::Cartesian, helix_ca(18));
    auto ev = spectrum_of(b.enm);
    const auto base = tencm::TorsionalENM::classify_spectrum(ev, Basis::Cartesian);

    ev.push_back(-1000.0 * zero_band(ev));
    std::sort(ev.begin(), ev.end());
    const auto got = tencm::TorsionalENM::classify_spectrum(ev, Basis::Cartesian);

    EXPECT_EQ(got.n_negative, base.n_negative + 1);
    EXPECT_EQ(got.vibrational.size(), base.vibrational.size());
    for (double x : got.vibrational) EXPECT_GT(x, 0.0);
}

// Documents what the gate restores, and pins the asymmetry that WAS the defect:
// with the gate ON the torsional cutoff is 0.0 (sign test) while the Cartesian
// cutoff is unchanged. This is also the able-to-fail demonstration in-suite.
TEST(SvibZeroModeStability, LegacyGateRestoresTheSignTestOnTheTorsionalBasisOnly) {
    Built tors, cart;
    build_path(tors, Basis::Torsional, helix_ca(24));
    build_path(cart, Basis::Cartesian, helix_ca(18));
    const auto ev_t = spectrum_of(tors.enm);
    const auto ev_c = spectrum_of(cart.enm);

    const auto def_t = tencm::TorsionalENM::classify_spectrum(ev_t, Basis::Torsional);
    const auto def_c = tencm::TorsionalENM::classify_spectrum(ev_c, Basis::Cartesian);
    ASSERT_GT(def_t.cutoff, 0.0);

    ScopedEnv gate("FLEXAIDDS_TENCOM_LEGACY_ZERO_SIGN", "1");
    const auto leg_t = tencm::TorsionalENM::classify_spectrum(ev_t, Basis::Torsional);
    const auto leg_c = tencm::TorsionalENM::classify_spectrum(ev_c, Basis::Cartesian);

    EXPECT_DOUBLE_EQ(leg_t.cutoff, 0.0)
        << "legacy gate must reproduce the bare sign test on the torsional basis";

    // The Cartesian branch is bit-identical under either rule, which is what
    // keeps the recorded BU72 anchor valid.
    EXPECT_DOUBLE_EQ(leg_c.cutoff, def_c.cutoff);
    ASSERT_EQ(leg_c.vibrational.size(), def_c.vibrational.size());
    for (size_t i = 0; i < def_c.vibrational.size(); ++i)
        EXPECT_DOUBLE_EQ(leg_c.vibrational[i], def_c.vibrational[i])
            << "Cartesian kept mode " << i << " moved under the gate";

    // And the gate is what a sign flip can still move -- the defect, on demand.
    // `gate` is still in scope, so the legacy rule is still active here.
    //
    // The legacy rule must be demonstrably sign-sensitive, otherwise this file
    // is guarding nothing. Compare the kept SET, not its size: when the zero
    // modes carry OPPOSITE signs the sign test keeps exactly one of them, so a
    // flip swaps WHICH mode is kept and leaves the count unchanged. MEASURED on
    // an OpenMP build of helix24: zero=2, one positive and one negative, 22 of
    // 23 kept either way. That is the more dangerous form of the defect --
    // ic2cf.cpp's ev_free.size() != ev_field.size() conservation check compares
    // COUNTS and would pass while the geometric mean silently changed.
    const auto pert = perturb_zero_modes(ev_t, zero_band(ev_t));
    const auto leg_pert = tencm::TorsionalENM::classify_spectrum(pert, Basis::Torsional);
    const bool legacy_is_sign_sensitive =
        leg_pert.vibrational.size() != leg_t.vibrational.size() ||
        !std::equal(leg_pert.vibrational.begin(), leg_pert.vibrational.end(),
                    leg_t.vibrational.begin());
    EXPECT_TRUE(legacy_is_sign_sensitive)
        << "legacy rule was expected to be sign-sensitive, but the kept spectrum "
           "was identical after flipping every zero mode ("
        << leg_t.vibrational.size() << " kept before, "
        << leg_pert.vibrational.size() << " after). If this fails the spectrum "
           "carried no sub-threshold mode and the demonstration is vacuous.";

    // And the fixed rule must be invariant on the SAME perturbation, which is
    // the contrast that makes the fix meaningful rather than merely different.
    ::unsetenv("FLEXAIDDS_TENCOM_LEGACY_ZERO_SIGN");
    const auto fix_pert = tencm::TorsionalENM::classify_spectrum(pert, Basis::Torsional);
    ASSERT_EQ(fix_pert.vibrational.size(), def_t.vibrational.size());
    for (size_t i = 0; i < def_t.vibrational.size(); ++i)
        EXPECT_DOUBLE_EQ(fix_pert.vibrational[i], def_t.vibrational[i])
            << "fixed rule moved kept mode " << i << " under a zero-mode flip";
}

// An EMPTY value must mean OFF. A raw getenv first-char test makes it mean ON.
TEST(SvibZeroModeStability, EmptyGateValueMeansDisabled) {
    Built b;
    build_path(b, Basis::Torsional, helix_ca(24));
    const auto ev = spectrum_of(b.enm);
    const auto def = tencm::TorsionalENM::classify_spectrum(ev, Basis::Torsional);

    ScopedEnv gate("FLEXAIDDS_TENCOM_LEGACY_ZERO_SIGN", "");
    const auto got = tencm::TorsionalENM::classify_spectrum(ev, Basis::Torsional);
    EXPECT_DOUBLE_EQ(got.cutoff, def.cutoff)
        << "an empty gate value enabled the legacy rule";
}
