#ifndef FLEXAIDDS_REMARK_ZERO_H
#define FLEXAIDDS_REMARK_ZERO_H

#include <cstring>

// ---------------------------------------------------------------------------
// strip_rendered_negzero -- drop the minus sign from an already-FORMATTED
// numeric field whose every displayed digit is zero.
//
// THE PROBLEM, MEASURED.  `REMARK entropy` is (E_avg - F)/T (statmech.cpp:321,
// 391).  For tightly-clustered poses that difference is a near-cancellation, so
// the quotient is floating-point residue whose SIGN depends on summation order
// -- and summation order shifts when unrelated code is added to the translation
// unit.  Re-emitting one 1JD0 run at %.17g instead of %.8f showed exactly that
// (98 poses, seed 12345):
//
//     26 poses   entropy = 0                     (exactly zero, positive)
//      4 poses   entropy = -2.3684757858670008e-17
//      2 poses   entropy = +4.7369515717340015e-17
//      2 poses   entropy = -2.960594732333751e-18
//      1 pose    entropy = +8.8580994391425824e-14
//      ...       a continuous ladder up to +0.0021548977890150432
//
// Under %.8f every magnitude below 5e-9 renders as "0.00000000" or
// "-0.00000000", so a sign flip in the residue shows up as a byte difference in
// a pose file while the pose itself is untouched.  That is what voided a
// byte-identity inertness comparison: 91 of 98 poses identical, 7 differing,
// every one a lone entropy line, flipping in BOTH directions.
//
// WHY NOT A SIGNED-ZERO TEST.  The first attempt collapsed IEEE-754 negative
// zero via its bit pattern.  It was correct and it fired on nothing, because
// -5.92e-18 is not negative zero -- it is a real negative number.  (An exact
// cancellation cannot produce -0.0 anyway: x - x is +0.0 in round-to-nearest.)
//
// WHY NOT A MAGNITUDE CLAMP.  A threshold like |v| < 1e-9 works only while the
// format stays %.8f; raise the precision and the same threshold starts hiding
// digits the format would have displayed.  The threshold and the format would
// have to be kept in agreement by hand, and nothing would fail if they drifted.
//
// WHAT THIS DOES INSTEAD.  It canonicalises the RENDERED TEXT, so the rule is
// precision-agnostic and self-maintaining: strip the sign only when every digit
// the format actually printed is zero.  By construction it cannot hide
// information -- if any displayed digit is nonzero, the sign is meaningful and
// is kept.  At %.17g the values above render their mantissa, so the sign
// survives; at %.8f they render as zero, so it does not.
//
// Pure string work, so no floating-point relaxation applies: cluster.cpp and
// BindingMode.cpp compile -O3 -ffast-math (confirmed in the build's own
// compile_commands.json), which permits the optimiser to delete a
// signed-zero comparison but has no bearing on memmove and character tests.
//
// Scoped to the value: it scans only after the last '=', because a REMARK field
// name never contains one.  Exponent forms that render a zero mantissa with an
// exponent (e.g. "-0.0e+00") are deliberately left ALONE -- conservative, since
// keeping a sign never hides a digit, and no emitter in this codebase uses %e.
// ---------------------------------------------------------------------------
static inline void strip_rendered_negzero(char* s)
{
	if (!s) return;

	char* eq = std::strrchr(s, '=');
	char* p  = eq ? eq + 1 : s;
	while (*p == ' ' || *p == '\t') ++p;
	if (*p != '-') return;

	bool any_digit = false;
	for (const char* q = p + 1; *q && *q != '\n' && *q != '\r'; ++q) {
		if (*q == '0') { any_digit = true; continue; }
		if (*q == '.') continue;
		return;              // a significant digit or an exponent: sign is real
	}
	if (!any_digit) return;  // "-" or "-." alone: not a rendered number

	std::memmove(p, p + 1, std::strlen(p + 1) + 1);
}

#endif // FLEXAIDDS_REMARK_ZERO_H
