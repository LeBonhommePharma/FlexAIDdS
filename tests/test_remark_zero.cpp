// tests/test_remark_zero.cpp -- regression test for strip_rendered_negzero().
//
// WHAT THIS GUARDS.  cluster.cpp and BindingMode.cpp emit `REMARK entropy` at
// %.8f. That value is (E_avg - F)/T, a near-cancellation whose residue is
// ~1e-17 and whose SIGN depends on summation order, so it flips when unrelated
// code is added to the translation unit. Under %.8f the residue renders as
// "0.00000000" or "-0.00000000", and the flip shows up as a byte difference in
// a pose file whose coordinates are identical -- which voided a byte-identity
// inertness comparison (measured: 91/98 identical, 7 differing, every one a
// lone entropy line, flipping in both directions).
//
// THE INVARIANT UNDER TEST, and the reason this is a test rather than a
// comment: the helper must strip the sign ONLY when every digit the format
// printed is zero. Any nonzero displayed digit means the sign carries
// information and must survive. Two earlier designs failed exactly here --
// a signed-zero bit-pattern test fired on nothing (the values are not -0.0),
// and a magnitude clamp would silently start hiding digits if the format
// precision were ever raised. This test fails loudly in both directions.

#include <gtest/gtest.h>

#include "remark_zero.h"

#include <cstdio>
#include <string>

namespace {

std::string stripped(const char* in)
{
	char buf[256];
	std::snprintf(buf, sizeof buf, "%s", in);
	strip_rendered_negzero(buf);
	return std::string(buf);
}

// Format a double exactly as the engine does, then canonicalise.
std::string emitted(double v, const char* fmt = "REMARK entropy = %.8f\n")
{
	char buf[256];
	std::snprintf(buf, sizeof buf, fmt, v);
	strip_rendered_negzero(buf);
	return std::string(buf);
}

} // namespace

// --- the sign is dropped when every printed digit is zero -------------------

TEST(RemarkZero, StripsRenderedNegativeZero)
{
	EXPECT_EQ(stripped("REMARK entropy = -0.00000000\n"),
	                   "REMARK entropy = 0.00000000\n");
	EXPECT_EQ(stripped("REMARK entropy = -0\n"), "REMARK entropy = 0\n");
	EXPECT_EQ(stripped("REMARK entropy = -0.0\n"), "REMARK entropy = 0.0\n");
}

TEST(RemarkZero, LeavesPositiveZeroUntouched)
{
	EXPECT_EQ(stripped("REMARK entropy = 0.00000000\n"),
	                   "REMARK entropy = 0.00000000\n");
}

// --- the sign SURVIVES whenever a printed digit is nonzero ------------------

TEST(RemarkZero, KeepsSignOnAnyDisplayedDigit)
{
	// The smallest magnitude %.8f can show. Stripping this would be data loss.
	EXPECT_EQ(stripped("REMARK entropy = -0.00000001\n"),
	                   "REMARK entropy = -0.00000001\n");
	EXPECT_EQ(stripped("REMARK entropy = -1.50000000\n"),
	                   "REMARK entropy = -1.50000000\n");
	EXPECT_EQ(stripped("REMARK enthalpy = -15.945917\n"),
	                   "REMARK enthalpy = -15.945917\n");
}

// This is the property a magnitude clamp would have broken: raise the
// precision and the sign must come back, with no change to the helper.
TEST(RemarkZero, PrecisionAgnostic)
{
	const double residue = -5.9211894646675019e-18;   // measured, 1JD0_32.pdb
	EXPECT_EQ(emitted(residue), "REMARK entropy = 0.00000000\n");
	EXPECT_EQ(emitted(residue, "REMARK entropy = %.17g\n"),
	          "REMARK entropy = -5.9211894646675019e-18\n");
}

// Values measured from a real 1JD0 run at %.17g: the residue ladder renders as
// zero under %.8f and must therefore emit no sign, while the genuine value in
// the same run must keep its digits.
TEST(RemarkZero, MeasuredResidueLadder)
{
	for (double v : {-2.3684757858670008e-17, 4.7369515717340015e-17,
	                 -2.960594732333751e-18,  8.8580994391425824e-14,
	                 -4.7369515717340015e-17, 5.5135747819197906e-13}) {
		EXPECT_EQ(emitted(v), "REMARK entropy = 0.00000000\n")
			<< "residue " << v << " should render as unsigned zero at %.8f";
	}
	EXPECT_EQ(emitted(0.0021548977890150432),
	          "REMARK entropy = 0.00215490\n");
	EXPECT_EQ(emitted(-0.0021548977890150432),
	          "REMARK entropy = -0.00215490\n");
}

// --- scoping and robustness -------------------------------------------------

TEST(RemarkZero, ScansOnlyAfterTheLastEquals)
{
	// A minus sign in the field name must not be mistaken for the value's.
	EXPECT_EQ(stripped("REMARK soft-beta = -0.00000000\n"),
	                   "REMARK soft-beta = 0.00000000\n");
	// No '=' at all: the whole string is scanned, still only a rendered zero.
	EXPECT_EQ(stripped("-0.000\n"), "0.000\n");
}

TEST(RemarkZero, HandlesDegenerateInput)
{
	strip_rendered_negzero(nullptr);              // must not crash
	EXPECT_EQ(stripped("REMARK entropy = -\n"), "REMARK entropy = -\n");
	EXPECT_EQ(stripped("REMARK entropy = -.\n"), "REMARK entropy = -.\n");
	EXPECT_EQ(stripped(""), "");
}

TEST(RemarkZero, IsIdempotent)
{
	char buf[64];
	std::snprintf(buf, sizeof buf, "REMARK entropy = -0.00000000\n");
	strip_rendered_negzero(buf);
	strip_rendered_negzero(buf);
	EXPECT_STREQ(buf, "REMARK entropy = 0.00000000\n");
}
