// tests/test_vcontacts_calc_areas_degenerate.cpp
//
// Degenerate-face guard for the P2.1/P2.3 float spherical-excess area path
// added in PR#370 (the FlexAID-Fast recoding).  Vcontacts.cpp:981-1000 moved
// the l'Huilier tangent formula for the contact-area triangle from double to
// float:
//
//     U,V,W,X   = sqrt( (1±c0)(1±c1)(1±c2)/8 )                 // float
//     tansqrS   = (1 - U + V + W + X) / (1 + U - V - W - X)    // float ratio
//     ... (three sibling ratios) ...
//     area     += 4 r^2 * atan( sqrt( sqrt(tS*tA*tB*tC) ) )    // accum in double
//
// The accumulator and the cosine inputs stay double (bounds the error on
// generic faces), but the RATIO DENOMINATORS `(1 + U - V - W - X)` etc. go
// toward zero on near-degenerate / near-coplanar faces — precisely the
// geometry where a Voronoi contact sits on the knife-edge of existing.
// float32 loses relative precision there that double absorbed, and the
// downstream atan(sqrt(...)) amplifies it.
//
// So the drift is NOT uniform epsilon: measured on this formula it is ~2e-7
// (relative) on a well-conditioned face and climbs ~4 orders of magnitude,
// toward ~1e-3, as the vertex cosines approach 1.  On a marginal face that
// error is what can flip whether a contact (and its area) is counted, which
// is a pose-RANKING risk — and the Astex-85 A/B can pass on average while
// individual poses reorder.  This test pins the fragile-face behaviour as a
// tripwire: it fires if the near-degenerate relative error ever exceeds a
// documented ceiling, i.e. if the float path is made materially less faithful
// to the double reference than it is today (or than any future double
// fallback would be).
//
// The two formulas are replicated inline (not linked from Vcontacts.cpp) so
// the test is a pure numerical-property check with no engine dependency; the
// float path mirrors Vcontacts.cpp:981-1000 exactly, and the double path
// mirrors the pre-PR reference it replaced.
//
// Origin: PR#370 independent review (Fizz), review ask #3.
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>

namespace {

// Pre-PR reference: spherical-excess triangle area in double precision.
double area_double(double c0, double c1, double c2, double rado) {
    const double U = std::sqrt((1 + c0) * (1 + c1) * (1 + c2) / 8.0);
    const double V = std::sqrt((1 - c0) * (1 - c1) * (1 + c2) / 8.0);
    const double W = std::sqrt((1 - c0) * (1 + c1) * (1 - c2) / 8.0);
    const double X = std::sqrt((1 + c0) * (1 - c1) * (1 - c2) / 8.0);
    const double tS = (1 - U + V + W + X) / (1 + U - V - W - X);
    const double tA = (1 - U - V - W + X) / (1 + U + V + W - X);
    const double tB = (1 - U - V + W - X) / (1 + U + V - W + X);
    const double tC = (1 - U + V - W - X) / (1 + U - V + W + X);
    const double tp = tS * tA * tB * tC;
    if (tp <= 0.0) return 0.0;
    return 4.0 * rado * rado * std::atan(std::sqrt(std::sqrt(tp)));
}

// PR#370 float path: mirrors LIB/Vcontacts.cpp:981-1000 exactly.
double area_float(double c0d, double c1d, double c2d, double rado) {
    const float c0 = static_cast<float>(c0d);
    const float c1 = static_cast<float>(c1d);
    const float c2 = static_cast<float>(c2d);
    float rad4[4] = {
        (1.0f + c0) * (1.0f + c1) * (1.0f + c2) * 0.125f,
        (1.0f - c0) * (1.0f - c1) * (1.0f + c2) * 0.125f,
        (1.0f - c0) * (1.0f + c1) * (1.0f - c2) * 0.125f,
        (1.0f + c0) * (1.0f - c1) * (1.0f - c2) * 0.125f,
    };
    float root4[4];
    for (int q = 0; q < 4; ++q) root4[q] = std::sqrt(rad4[q]);
    const float U = root4[0], V = root4[1], W = root4[2], X = root4[3];
    const float tS = (1.0f - U + V + W + X) / (1.0f + U - V - W - X);
    const float tA = (1.0f - U - V - W + X) / (1.0f + U + V + W - X);
    const float tB = (1.0f - U - V + W - X) / (1.0f + U + V - W + X);
    const float tC = (1.0f - U + V - W - X) / (1.0f + U - V + W + X);
    const float tp = std::sqrt(tS * tA * tB * tC);
    if (tp > 0.0f) return 4.0 * rado * rado * static_cast<double>(std::atan(std::sqrt(tp)));
    return 0.0;
}

double rel_err(double ref, double got) {
    return ref != 0.0 ? std::fabs(ref - got) / std::fabs(ref) : std::fabs(got);
}

constexpr double kRado = 1.9;  // representative contact radius (Å)

}  // namespace

// A well-conditioned face: the float path must track double to near-epsilon.
// If this ever loosens, the float move has cost precision on GENERIC faces,
// not just fragile ones — a different and larger regression.
TEST(CalcAreasFloat, GenericFaceMatchesDoubleToEpsilon) {
    const double d = area_double(0.30, 0.20, 0.25, kRado);
    const double f = area_float(0.30, 0.20, 0.25, kRado);
    EXPECT_GT(d, 0.0);
    EXPECT_LT(rel_err(d, f), 1e-5)
        << "generic-face float area diverged from double by more than 1e-5 relative";
}

// Near-degenerate faces (vertex cosines -> 1, denominators -> 0).  The float
// error concentrates here; today it peaks near ~1e-3.  This ceiling is the
// tripwire: crossing 1e-2 means the fragile-face precision has degraded past
// the point where "rank/pose-equivalent" can be assumed without an Astex A/B.
TEST(CalcAreasFloat, DegenerateFaceStaysWithinCeiling) {
    double worst = 0.0;
    for (double c = 0.99; c <= 0.99999; c = 1.0 - (1.0 - c) * 0.1) {
        const double d = area_double(c, c, c, kRado);
        const double f = area_float(c, c, c, kRado);
        const double re = rel_err(d, f);
        worst = std::max(worst, re);
        EXPECT_LT(re, 1e-2)
            << "near-degenerate face c=" << c
            << " float area diverged from double by " << re
            << " relative (ceiling 1e-2)";
    }
    // Document the concentration: the fragile-face error is orders of
    // magnitude above the generic-face epsilon measured above.  This is not a
    // pass/fail assertion (a future double fallback would legitimately shrink
    // it); it records the measured worst case for reviewers.
    RecordProperty("worst_degenerate_rel_err_x1e6",
                   static_cast<int>(worst * 1e6));
}
