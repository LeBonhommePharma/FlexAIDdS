#ifndef FLEXAIDS_FS_SAFE_H
#define FLEXAIDS_FS_SAFE_H

#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <initializer_list>
#include <string>
#include <system_error>

// Non-throwing wrappers for the std::filesystem queries whose THROWING
// overload aborts the engine on a path that has merely gone away.
//
// WHY THIS EXISTS
// ---------------
// Companion to temp_dir.h, same root cause, wider blast radius. Most of
// std::filesystem has two overloads: one that reports through a
// std::error_code and one that throws std::filesystem_error. The throwing one
// is the default -- you get it by writing the obvious thing -- and an
// uncaught filesystem_error terminates the process. On 2026-09-20 exactly one
// such call (temp_directory_path) cost nine cells of an 85-target Astex
// campaign, and the driver reported it only as "Incomplete docking", so it
// read like a docking failure on specific targets rather than a startup abort.
//
// A tree-wide count taken while writing this header found 220 std::filesystem
// calls in LIB/, of which 59 were non-predicate throwing calls and 48 of those
// sat outside any try block. Those 48 are the abort surface. This header
// exists to shrink the part of it that can be shrunk WITHOUT converting a loud
// failure into a quiet one.
//
// THE RULE THAT MATTERS MORE THAN THE WRAPPERS
// --------------------------------------------
// `remove(p)` and `remove(p, ec)` are DIFFERENT PROGRAMS. Mechanically
// swapping every throwing call for its error_code sibling does not make the
// engine safer; it makes it silent, which is worse, because a campaign that
// aborts is one you re-run and a campaign that quietly skips is one you
// publish. Every wrapper below therefore has an explicit, documented answer to
// "what happens instead of the abort", and each answer is the CONSERVATIVE
// branch the call site already takes for a missing file:
//
//   file_size_or()  -> an unstattable file counts as SIZE ZERO, i.e. absent.
//                      Every call site in LIB guards file_size with exists()
//                      and asks "is there a non-trivially-sized file here";
//                      "cannot stat" and "not there" already led to the same
//                      branch, so this removes an abort and changes nothing.
//   out_of_date()   -> a product whose mtime cannot be compared is STALE and
//                      gets regenerated. Never the other way round: treating
//                      an unverifiable cache as fresh is how you silently
//                      publish a stale receptor.
//   canonical_or()  -> an unresolvable path falls back to the string as
//                      given, which is still subject to whatever validation
//                      the call site already applies.
//   ensure_dir()    -> returns false and says why. The CALLER decides; this
//                      function deliberately does not, because "carry on
//                      without an output directory" is right in one place and
//                      catastrophic in the next.
//
// WHAT IS DELIBERATELY NOT HERE
// -----------------------------
// No wrapper for remove(): [fs.op.remove] already returns false rather than
// throwing when the path does not exist, so the vanishing-path failure mode
// does not apply to it. It throws only when the filesystem actively refuses a
// delete (EACCES, EBUSY, a non-empty directory), and that is a real fault
// worth surfacing loudly. Wrapping it would suppress a genuine error to buy
// nothing.
//
// No wrapper for exists()/is_directory()/is_regular_file(): these return false
// for a missing path and throw only on a non-ENOENT stat error.
//
// No class, no RAII type, no policy object. The goal is fewer aborts, not an
// abstraction layer.

namespace flexaids {
namespace fs_safe {

namespace detail {
inline void warn(const char* what, const std::filesystem::path& p,
                 const std::error_code& ec, const char* consequence) {
    // One line, names the path, names what the engine does next. The 09-20
    // post-mortem cost six hours because the failure surfaced as a generic
    // driver-level message with no path in it.
    std::fprintf(stderr, "[FS] %s failed on '%s' (%s); %s\n",
                 what, p.string().c_str(), ec.message().c_str(), consequence);
}
}  // namespace detail

// Size in bytes, or `fallback` when the file is missing, unreadable, or not a
// regular file. Silent by design: at every call site the caller is asking a
// yes/no question about usability and immediately branches on the answer, so
// a diagnostic here would fire once per cached file per target.
inline std::uintmax_t file_size_or(const std::filesystem::path& p,
                                   std::uintmax_t fallback = 0) {
    std::error_code ec;
    const std::uintmax_t n = std::filesystem::file_size(p, ec);
    return ec ? fallback : n;
}

// True when `product` must be regenerated from `sources`.
//
// Conservative in one direction on purpose: any stat failure -- on the product
// OR on a source -- returns true. We can never prove a cache current from
// missing information, and regenerating an already-current file costs CPU
// while reusing a stale one costs a wrong number in a results table.
//
// Replaces the idiom
//     !fs::exists(prod) || fs::last_write_time(prod) < fs::last_write_time(src)
// whose SECOND operand had no existence guard anywhere in LIB: a source file
// pruned between runs aborted the engine instead of triggering a rebuild.
inline bool out_of_date(const std::filesystem::path& product,
                        std::initializer_list<std::filesystem::path> sources) {
    std::error_code ec;
    const auto t_prod = std::filesystem::last_write_time(product, ec);
    if (ec) return true;                    // missing or unstattable -> rebuild
    for (const auto& s : sources) {
        std::error_code sec;
        const auto t_src = std::filesystem::last_write_time(s, sec);
        if (sec) {
            // Distinct from the product case and worth a line: the INPUT is
            // unreadable, so the rebuild we are about to trigger will probably
            // fail too, and the reader should know which file to look at.
            detail::warn("last_write_time", s, sec,
                         "treating the cached product as stale and rebuilding");
            return true;
        }
        if (t_prod < t_src) return true;
    }
    return false;
}

// Canonical form, or the input unchanged when it cannot be resolved.
//
// canonical() throws when any component does not exist. Call sites here feed
// the result into a shell command built once per docking target, so an abort
// lands in the middle of a campaign. Falling back to the given string is safe
// only because the call sites re-apply their own validation (is_safe_exec_path)
// to whatever comes back -- this function does NOT sanitise.
inline std::string canonical_or(const std::string& p) {
    std::error_code ec;
    const auto c = std::filesystem::canonical(p, ec);
    if (!ec) return c.string();
    // weakly_canonical tolerates non-existent trailing components and still
    // collapses ".." and symlinks, so it is a strictly better second try.
    std::error_code wec;
    const auto w = std::filesystem::weakly_canonical(p, wec);
    if (!wec) return w.string();
    detail::warn("canonical", p, ec, "using the path as given, unresolved");
    return p;
}

// Absolute form, or the input unchanged.
//
// absolute() is not a pure string operation: it consults current_path(), which
// throws when the process working directory has been DELETED. That is the
// precise object an agent-harness sweep destroys, and it is why a function
// that looks like string concatenation is on this list at all.
inline std::string absolute_or(const std::string& p) {
    std::error_code ec;
    const auto a = std::filesystem::absolute(p, ec);
    if (ec) {
        detail::warn("absolute", p, ec,
                     "using the relative path as given (cwd may be deleted)");
        return p;
    }
    return a.string();
}

// Create `d` and its parents. Returns true when the directory exists and is
// usable afterwards. On failure fills *why (when non-null) and returns false.
//
// Returns a value rather than deciding, because the right reaction is not
// uniform: a CLI that cannot create its output directory should stop, while a
// per-target sidecar directory inside a long campaign may be skippable. Note
// create_directories() returns false WITHOUT error when the directory already
// exists, so the error_code -- not the return value -- is the success test.
inline bool ensure_dir(const std::filesystem::path& d, std::string* why = nullptr) {
    std::error_code ec;
    std::filesystem::create_directories(d, ec);
    if (ec) {
        if (why) *why = ec.message();
        return false;
    }
    std::error_code sec;
    if (!std::filesystem::is_directory(d, sec)) {
        // create_directories can report success while the path is occupied by
        // a regular file left behind by an earlier run.
        if (why) *why = "path exists but is not a directory";
        return false;
    }
    return true;
}

}  // namespace fs_safe
}  // namespace flexaids

#endif  // FLEXAIDS_FS_SAFE_H
