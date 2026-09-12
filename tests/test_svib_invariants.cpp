// tests/test_svib_invariants.cpp — physical invariants of the vibrational-entropy path.
//
// WHY INVARIANTS AND NOT GOLDEN NUMBERS
// A recorded expected value asserts only that today's output equals today's
// output. It pins behaviour without encoding understanding and goes stale the
// moment anything legitimately changes. Every assertion below instead states a
// property that MUST hold of a correct implementation, so each catches a CLASS
// of defect rather than one instance. Where an invariant is not applicable the
// reason is stated at the test site rather than the test being omitted silently.
//
// Each test in this file was verified to FAIL under a deliberate mutation of the
// code or data it guards; the mutation and the resulting failure are recorded in
// the PR. A test with no demonstrated failing case is a coverage line, not a test.
//
// Scope: this file exercises tencm::TorsionalENM and encom::ENCoMEngine directly.
// It contains no stubs and no substitutes -- every call reaches the real
// implementation. It does not touch LIB/ic2cf.cpp.
//
// Apache-2.0 (c) 2026 Le Bonhomme Pharma.

#include <gtest/gtest.h>

#include "tENCoM/tencm.h"
#include "encom.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <vector>

namespace {

// A synthetic alpha-helix-like CA trace: deterministic, connected, and dense
// enough that the contact graph is a single component at the default cutoff.
std::vector<std::array<float,3>> helix_ca(int n) {
    std::vector<std::array<float,3>> ca;
    ca.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        const double t = i * 100.0 * M_PI / 180.0;   // 100 deg per residue
        ca.push_back({static_cast<float>(2.3 * std::cos(t)),
                      static_cast<float>(2.3 * std::sin(t)),
                      static_cast<float>(1.5 * i)});
    }
    return ca;
}

// Minimal heavy-atom array for the Cartesian ANM entry point. Only coor[] and
// element are read by build_from_ligand; everything else is value-initialised.
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

// Every invariant below is asserted on BOTH assembly paths, because they are
// different code: build_from_ca assembles the TORSIONAL Hessian and
// build_from_ligand the Cartesian ANM. A defect injected into one is invisible
// to a test that only drives the other -- verified, that gap existed in the
// first draft of this file and the mutation matrix caught it.
enum class Path { Torsional, Cartesian };

const char* path_name(Path p) {
    return p == Path::Torsional ? "build_from_ca(torsional)"
                                : "build_from_ligand(cartesian)";
}

// Holds the atom backing store alive for the Cartesian path.
struct Built {
    tencm::TorsionalENM       enm;
    std::vector<atom>         backing;
};

void build_path(Built& b, Path path,
                const std::vector<std::array<float,3>>& xyz,
                float cutoff, float k0) {
    if (path == Path::Torsional) {
        b.enm.build_from_ca(xyz, cutoff, k0);
    } else {
        b.backing = carbon_atoms(xyz);
        b.enm.build_from_ligand(b.backing.data(), 0,
                                static_cast<int>(b.backing.size()), cutoff, k0);
    }
}

std::vector<double> sorted_eigs(const tencm::TorsionalENM& enm) {
    std::vector<double> e;
    for (const auto& m : enm.modes()) e.push_back(m.eigenvalue);
    std::sort(e.begin(), e.end());
    return e;
}

// tencm::NormalMode and encom::NormalMode are DISTINCT types that share a name
// across namespaces; the entropy engine consumes the encom one, so the bridge is
// explicit here rather than implicit anywhere. Only the eigenvalue and its
// sqrt are carried, which is exactly what compute_vibrational_entropy reads.
std::vector<encom::NormalMode> to_encom(const std::vector<tencm::NormalMode>& in) {
    std::vector<encom::NormalMode> out;
    out.reserve(in.size());
    int idx = 1;
    for (const auto& m : in) {
        encom::NormalMode e;
        e.index      = idx++;
        e.eigenvalue = m.eigenvalue;
        e.frequency  = m.eigenvalue > 0.0 ? std::sqrt(m.eigenvalue) : 0.0;
        out.push_back(std::move(e));
    }
    return out;
}

double svib(const tencm::TorsionalENM& enm, double T) {
    return encom::ENCoMEngine::compute_vibrational_entropy(to_encom(enm.modes()), T)
               .S_vib_kcal_mol_K;
}

constexpr float kCut = 9.0f;   // tencm::DEFAULT_RC
constexpr float kK0  = 1.0f;   // tencm::DEFAULT_K0

} // namespace

// ─── Rigid-body invariance ──────────────────────────────────────────────────
// The class of bug this catches is the one that just cost us the docking gate:
// a coordinate frame applied on one side and not the other.

TEST(SvibInvariants, SpectrumIsTranslationInvariant) {
    // TOLERANCE IS CALIBRATED TO A MEASURED LIMIT, NOT CHOSEN TO PASS.
    // Coordinates are stored as float, so absolute precision degrades as |x|
    // grows and a large translation costs significant digits. Measured on this
    // build (torsional path), max relative eigenvalue residual vs shift:
    //     0 A -> 0.0 exactly    1 A -> 3.1e-07   100 A -> 1.5e-06
    //  1000 A -> 3.1e-05    10000 A -> 1.8e-04
    // The residual scales with the shift, which is float storage, not a frame
    // bug: a genuine frame error moves the spectrum by O(1) relative, four
    // decades above this.
    for (Path path : {Path::Torsional, Path::Cartesian}) {
        auto a = helix_ca(20);
        auto b = a;
        for (auto& q : b) { q[0] += 100.0f; q[1] -= 50.0f; q[2] += 33.0f; }

        Built ba, bb;
        build_path(ba, path, a, kCut, kK0);
        build_path(bb, path, b, kCut, kK0);
        ASSERT_TRUE(ba.enm.is_built()) << path_name(path);
        ASSERT_TRUE(bb.enm.is_built()) << path_name(path);

        const auto sa = sorted_eigs(ba.enm), sb = sorted_eigs(bb.enm);
        ASSERT_EQ(sa.size(), sb.size()) << path_name(path);
        for (size_t i = 0; i < sa.size(); ++i)
            EXPECT_NEAR(sa[i], sb[i], 1e-5 * std::max(1.0, std::fabs(sa[i])))
                << path_name(path) << " mode " << i
                << " moved under pure translation";
    }
}

TEST(SvibInvariants, SpectrumIsRotationInvariant) {
    for (Path path : {Path::Torsional, Path::Cartesian}) {
        auto a = helix_ca(20);
        auto b = a;
        const double th = 0.7;                  // rotation about z
        for (auto& q : b) {
            const double x = q[0], y = q[1];
            q[0] = static_cast<float>(x * std::cos(th) - y * std::sin(th));
            q[1] = static_cast<float>(x * std::sin(th) + y * std::cos(th));
        }

        Built ba, bb;
        build_path(ba, path, a, kCut, kK0);
        build_path(bb, path, b, kCut, kK0);

        const auto sa = sorted_eigs(ba.enm), sb = sorted_eigs(bb.enm);
        ASSERT_EQ(sa.size(), sb.size()) << path_name(path);
        for (size_t i = 0; i < sa.size(); ++i)
            EXPECT_NEAR(sa[i], sb[i], 1e-5 * std::max(1.0, std::fabs(sa[i])))
                << path_name(path) << " mode " << i
                << " moved under pure rotation";
    }
}

// Permutation invariance is asserted on the CARTESIAN ligand path, not on
// build_from_ca. For a CA trace the index order IS the chain connectivity, so
// permuting it describes a different molecule and the spectrum SHOULD change --
// asserting invariance there would be asserting something false. The Cartesian
// ANM is purely distance/contact based, so atom order carries no information
// and the spectrum must be invariant.
TEST(SvibInvariants, CartesianSpectrumIsAtomOrderInvariant) {
    auto xyz = helix_ca(18);
    auto atoms_a = carbon_atoms(xyz);

    auto shuffled = xyz;
    std::reverse(shuffled.begin(), shuffled.end());
    auto atoms_b = carbon_atoms(shuffled);

    tencm::TorsionalENM ea, eb;
    ea.build_from_ligand(atoms_a.data(), 0, static_cast<int>(atoms_a.size()), kCut, kK0);
    eb.build_from_ligand(atoms_b.data(), 0, static_cast<int>(atoms_b.size()), kCut, kK0);
    ASSERT_TRUE(ea.is_built()); ASSERT_TRUE(eb.is_built());

    const auto sa = sorted_eigs(ea), sb = sorted_eigs(eb);
    ASSERT_EQ(sa.size(), sb.size());
    for (size_t i = 0; i < sa.size(); ++i)
        EXPECT_NEAR(sa[i], sb[i], 1e-6 * std::max(1.0, std::fabs(sa[i])))
            << "mode " << i << " depends on atom ORDER";
}

// ─── Spectrum well-formedness ───────────────────────────────────────────────

TEST(SvibInvariants, HessianSpectrumIsPositiveSemiDefinite) {
    for (Path path : {Path::Torsional, Path::Cartesian}) {
        Built b;
        build_path(b, path, helix_ca(20), kCut, kK0);
        ASSERT_TRUE(b.enm.is_built()) << path_name(path);

        const auto e = sorted_eigs(b.enm);
        ASSERT_FALSE(e.empty()) << path_name(path);
        const double tol = 1e-8 * std::max(1.0, std::fabs(e.back()));
        for (size_t i = 0; i < e.size(); ++i)
            EXPECT_GE(e[i], -tol)
                << path_name(path) << " negative eigenvalue at " << i
                << " => imaginary frequency, unstable network";
    }
}

// The zero-mode count is a CONVENTION of the basis, taken from the source
// (tencm.h: n_expected receives 6 on the Cartesian path and 0 on the torsional
// one) rather than assumed to be Cartesian.
TEST(SvibInvariants, DroppedModeCountMatchesBasisConvention) {
    auto ca = helix_ca(24);
    tencm::TorsionalENM tors;
    tors.build_from_ca(ca, kCut, kK0);
    int dropped = -1, expected = -1;
    (void)tors.vibrational_eigenvalues(&dropped, &expected);
    EXPECT_EQ(expected, 0) << "torsional basis should declare 0 rigid modes";
    EXPECT_EQ(dropped, expected)
        << "dropped != expected means the contact graph fragmented";

    auto atoms = carbon_atoms(helix_ca(18));
    tencm::TorsionalENM cart;
    cart.build_from_ligand(atoms.data(), 0, static_cast<int>(atoms.size()), kCut, kK0);
    int d2 = -1, e2 = -1;
    (void)cart.vibrational_eigenvalues(&d2, &e2);
    EXPECT_EQ(e2, 6) << "Cartesian ANM should declare 6 rigid-body modes";
    EXPECT_EQ(d2, e2)
        << "dropped != expected means the contact graph fragmented";
}

TEST(SvibInvariants, BuildIsDeterministic) {
    for (Path path : {Path::Torsional, Path::Cartesian}) {
        const auto ca = helix_ca(20);
        Built b1, b2;
        build_path(b1, path, ca, kCut, kK0);
        build_path(b2, path, ca, kCut, kK0);

        const auto s1 = sorted_eigs(b1.enm), s2 = sorted_eigs(b2.enm);
        ASSERT_EQ(s1.size(), s2.size()) << path_name(path);
        for (size_t i = 0; i < s1.size(); ++i)
            EXPECT_DOUBLE_EQ(s1[i], s2[i])
                << path_name(path) << " mode " << i << " is not reproducible";
        EXPECT_DOUBLE_EQ(svib(b1.enm, 300.0), svib(b2.enm, 300.0)) << path_name(path);
    }
}

// ─── Thermodynamic monotonicity ─────────────────────────────────────────────

// C_v = T dS/dT >= 0, so S must be non-decreasing in T. A violation breaks
// thermodynamics, not a convention.
TEST(SvibInvariants, EntropyIncreasesWithTemperature) {
    for (Path path : {Path::Torsional, Path::Cartesian}) {
        Built b;
        build_path(b, path, helix_ca(20), kCut, kK0);
        ASSERT_TRUE(b.enm.is_built()) << path_name(path);

        double prev = -std::numeric_limits<double>::infinity();
        for (double T : {100.0, 200.0, 300.0, 400.0, 600.0}) {
            const double s = svib(b.enm, T);
            EXPECT_GT(s, prev) << path_name(path) << " S_vib not increasing at T="
                               << T << " (C_v = T dS/dT would be negative)";
            prev = s;
        }
    }
}

// Stiffer springs raise frequencies and lower entropy. Direction must hold.
TEST(SvibInvariants, StifferNetworkHasLowerEntropy) {
    for (Path path : {Path::Torsional, Path::Cartesian}) {
        const auto ca = helix_ca(20);
        Built soft, stiff;
        build_path(soft,  path, ca, kCut, 1.0f);
        build_path(stiff, path, ca, kCut, 4.0f);
        ASSERT_TRUE(soft.enm.is_built()) << path_name(path);
        ASSERT_TRUE(stiff.enm.is_built()) << path_name(path);

        const auto ss = sorted_eigs(soft.enm), st = sorted_eigs(stiff.enm);
        ASSERT_EQ(ss.size(), st.size()) << path_name(path);
        for (size_t i = 0; i < ss.size(); ++i)
            EXPECT_GE(st[i], ss[i] - 1e-9)
                << path_name(path) << " stiffening lowered eigenvalue " << i;

        EXPECT_LT(svib(stiff.enm, 300.0), svib(soft.enm, 300.0))
            << path_name(path)
            << " a stiffer network must have LOWER vibrational entropy";
    }
}

// ─── Additivity: the strongest available check ──────────────────────────────
// Two non-interacting copies must give S_total = S_1 + S_2. This probes
// ABSOLUTE correctness with no reference implementation, and catches
// normalisation and double-counting errors that every self-consistency check
// passes. Separation is far beyond the contact cutoff, so the two blocks share
// no springs and the combined spectrum is the union of the two.
TEST(SvibInvariants, NonInteractingSubsystemsAreAdditive) {
    const auto base = helix_ca(12);

    auto a_atoms = carbon_atoms(base);
    auto far = base;
    for (auto& p : far) p[0] += 1000.0f;          // >> kCut
    auto b_atoms = carbon_atoms(far);

    std::vector<std::array<float,3>> both;
    both.insert(both.end(), base.begin(), base.end());
    both.insert(both.end(), far.begin(),  far.end());
    auto ab_atoms = carbon_atoms(both);

    tencm::TorsionalENM ea, eb, eab;
    ea.build_from_ligand(a_atoms.data(), 0, static_cast<int>(a_atoms.size()), kCut, kK0);
    eb.build_from_ligand(b_atoms.data(), 0, static_cast<int>(b_atoms.size()), kCut, kK0);
    eab.build_from_ligand(ab_atoms.data(), 0, static_cast<int>(ab_atoms.size()), kCut, kK0);
    ASSERT_TRUE(ea.is_built()); ASSERT_TRUE(eb.is_built()); ASSERT_TRUE(eab.is_built());

    const double sa  = svib(ea,  300.0);
    const double sb  = svib(eb,  300.0);
    const double sab = svib(eab, 300.0);

    // Tolerance is relative to the magnitude being summed, not absolute: the
    // quantity is a sum of ~N logs and the eigen solver carries rounding.
    const double scale = std::max(1e-12, std::fabs(sa) + std::fabs(sb));
    EXPECT_NEAR(sab, sa + sb, 1e-6 * scale)
        << "S(A+B) != S(A) + S(B) for non-interacting blocks: "
        << "S(AB)=" << sab << " S(A)=" << sa << " S(B)=" << sb
        << " residual=" << (sab - (sa + sb));
}
