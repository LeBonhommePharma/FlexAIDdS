// test_shell_exec.cpp — unit tests for path-safe shell quoting / argv exec helpers
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>

#include "shell_exec.h"

#include <stdexcept>
#include <string>

using flexaids::shell_exec::PathReject;
using flexaids::shell_exec::is_safe_exec_path;
using flexaids::shell_exec::run_argv_first_token;
using flexaids::shell_exec::shell_quote;
using flexaids::shell_exec::shell_quote_checked;
using flexaids::shell_exec::shell_quote_raw;
using flexaids::shell_exec::validate_exec_path;

// ---------------------------------------------------------------------------
// shell_quote_raw — single-quote rule
// ---------------------------------------------------------------------------

TEST(ShellExec, QuoteSimplePath) {
    EXPECT_EQ(shell_quote_raw("/tmp/data"), "'/tmp/data'");
}

TEST(ShellExec, QuoteEmbeddedSingleQuote) {
    // foo'bar → 'foo'\''bar'
    EXPECT_EQ(shell_quote_raw("foo'bar"), "'foo'\\''bar'");
}

TEST(ShellExec, QuoteSpacesAndMetacharacters) {
    // Single-quoting neutralizes spaces, $, `, ;, |, &, etc.
    EXPECT_EQ(shell_quote_raw("a b;c$(x)"), "'a b;c$(x)'");
    EXPECT_EQ(shell_quote_raw("x`id`y"), "'x`id`y'");
    EXPECT_EQ(shell_quote_raw("p|q&r"), "'p|q&r'");
}

TEST(ShellExec, QuoteEmptyIsEmptyQuoted) {
    EXPECT_EQ(shell_quote_raw(""), "''");
}

// ---------------------------------------------------------------------------
// validate_exec_path — refuse NUL / newlines / controls
// ---------------------------------------------------------------------------

TEST(ShellExec, AcceptNormalPaths) {
    EXPECT_EQ(validate_exec_path("/Users/me/data/MC_st0r5.2_6.dat"), PathReject::Ok);
    EXPECT_EQ(validate_exec_path("relative/path with spaces.pdb"), PathReject::Ok);
    EXPECT_TRUE(is_safe_exec_path("/tmp/ligand.sdf"));
}

TEST(ShellExec, AcceptTabInPath) {
    EXPECT_EQ(validate_exec_path("name\twith\ttab"), PathReject::Ok);
}

TEST(ShellExec, RejectEmpty) {
    EXPECT_EQ(validate_exec_path(""), PathReject::Empty);
    EXPECT_FALSE(is_safe_exec_path(""));
}

TEST(ShellExec, RejectNul) {
    std::string p = "evil";
    p.push_back('\0');
    p += "payload";
    EXPECT_EQ(validate_exec_path(p), PathReject::ContainsNul);
}

TEST(ShellExec, RejectNewlineAndCr) {
    EXPECT_EQ(validate_exec_path("a\nb"), PathReject::ContainsNewline);
    EXPECT_EQ(validate_exec_path("a\rb"), PathReject::ContainsNewline);
}

TEST(ShellExec, RejectOtherControls) {
    EXPECT_EQ(validate_exec_path(std::string("a") + '\x01' + "b"),
              PathReject::ContainsControl);
    EXPECT_EQ(validate_exec_path(std::string("a") + '\x7f' + "b"),
              PathReject::ContainsControl);
}

// ---------------------------------------------------------------------------
// shell_quote / shell_quote_checked — validate then quote
// ---------------------------------------------------------------------------

TEST(ShellExec, ShellQuoteThrowsOnNewline) {
    EXPECT_THROW(
        { auto q = shell_quote("a\nb"); (void)q; },
        std::invalid_argument);
}

TEST(ShellExec, ShellQuoteCheckedNulloptOnNul) {
    std::string p = "x";
    p.push_back('\0');
    p += "y";
    EXPECT_FALSE(shell_quote_checked(p).has_value());
}

TEST(ShellExec, ShellQuoteCheckedOk) {
    auto q = shell_quote_checked("/data/dir");
    ASSERT_TRUE(q.has_value());
    EXPECT_EQ(*q, "'/data/dir'");
}

TEST(ShellExec, ShellQuoteHandlesEmbeddedQuote) {
    EXPECT_EQ(shell_quote("o'reilly"), "'o'\\''reilly'");
}

// ---------------------------------------------------------------------------
// run_argv_first_token — no shell; injection cannot expand
// ---------------------------------------------------------------------------

TEST(ShellExec, RunArgvFirstTokenEcho) {
    // argv form: metacharacters are literal args, not shell syntax.
    auto t = run_argv_first_token({"printf", "%s", "hello_world"});
    EXPECT_EQ(t, "hello_world");
}

TEST(ShellExec, RunArgvRejectsUnsafeArg) {
    auto t = run_argv_first_token({"printf", "%s", "a\nb"});
    EXPECT_TRUE(t.empty());
}
