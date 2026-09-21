#include "temp_dir.h"

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <mutex>
#include <string>
#include <system_error>
#include <vector>

#if defined(_WIN32)
#include <process.h>
#define FLEXAIDS_GETPID _getpid
#else
#include <unistd.h>
#define FLEXAIDS_GETPID getpid
#endif

namespace fs = std::filesystem;

namespace flexaids {
namespace {

struct Resolved {
    std::string dir = ".";
    std::string source = "unresolved";
    bool ok = false;
};

// Writability is established by CREATING a file, not by asking. access(2) and
// fs::status() both report on metadata that can be stale or, under a sandbox,
// simply wrong: a directory can be reported writable and still reject every
// write. The probe is the only check that answers the question actually being
// asked, which is "can this process put a scratch file here right now".
bool probe_writable(const fs::path& d, std::string* why) {
    std::error_code ec;

    if (!fs::exists(d, ec) || ec) { *why = "does not exist"; return false; }
    if (!fs::is_directory(d, ec) || ec) { *why = "not a directory"; return false; }

    const fs::path probe =
        d / (".flexaidds_wtest_" + std::to_string(static_cast<long>(FLEXAIDS_GETPID())));
    // Distinct name per process so two concurrent engines cannot race here; a
    // campaign runs several docking processes at once and a shared probe name
    // would make one delete the other's file mid-check.
    std::FILE* f = std::fopen(probe.string().c_str(), "wb");
    if (!f) { *why = "not writable"; return false; }
    const bool wrote = std::fputc('x', f) != EOF;
    std::fclose(f);
    fs::remove(probe, ec);          // best effort; a leftover probe is harmless
    if (!wrote) { *why = "write failed"; return false; }
    return true;
}

Resolved resolve() {
    struct Candidate { fs::path path; std::string source; };
    std::vector<Candidate> cands;

    if (const char* ov = std::getenv("FLEXAIDDS_TMPDIR"); ov && *ov)
        cands.push_back({fs::path(ov), "FLEXAIDDS_TMPDIR"});

    // The error_code overload NEVER throws -- this is the whole point of the
    // file. The throwing sibling is what aborted nine campaign cells.
    std::error_code ec;
    const fs::path sys = fs::temp_directory_path(ec);
    if (!ec && !sys.empty())
        cands.push_back({sys, "system_temp"});  // LABEL, not the function name: it is printed
                                        // on the recovery path, and a monitor grepping
                                        // the function name matched BOTH the abort and
                                        // the healthy recovery (measured 2026-09-20).
    else
        std::fprintf(stderr,
                     "[TEMPDIR] system temp lookup returned no usable path (%s); "
                     "continuing down the fallback chain [this is RECOVERY, not a failure]\n",
                     ec ? ec.message().c_str() : "empty");

#if !defined(_WIN32)
    cands.push_back({fs::path("/tmp"), "/tmp"});
#endif
    cands.push_back({fs::current_path(ec), "cwd"});

    for (const auto& c : cands) {
        if (c.path.empty()) continue;
        std::string why;
        if (probe_writable(c.path, &why)) {
            // Announce ONLY when the answer is surprising -- i.e. when the
            // first choice was not usable. A line on every healthy start is
            // noise that trains readers to ignore the channel.
            if (c.source != cands.front().source)
                std::fprintf(stderr,
                             "[TEMPDIR] using %s (%s) after earlier candidates failed\n",
                             c.path.string().c_str(), c.source.c_str());
            return {c.path.string(), c.source, true};
        }
        std::fprintf(stderr, "[TEMPDIR] rejected %s (%s): %s\n",
                     c.path.string().c_str(), c.source.c_str(), why.c_str());
    }

    std::fprintf(stderr,
                 "[TEMPDIR] FATAL-SOFT: no writable temporary directory found; "
                 "falling back to \".\". Set FLEXAIDDS_TMPDIR to a durable path.\n");
    return {".", "none", false};
}

Resolved& cached() {
    // Resolved once per process, under call_once so concurrent OpenMP regions
    // cannot race the probe. Caching is deliberate: if the chosen directory is
    // deleted mid-run, a long docking job keeps a consistent answer instead of
    // silently migrating its scratch somewhere else between restarts.
    static Resolved r;
    static std::once_flag once;
    std::call_once(once, [] { r = resolve(); });
    return r;
}

}  // namespace

const std::string& temp_dir() { return cached().dir; }
bool temp_dir_ok() { return cached().ok; }
const std::string& temp_dir_source() { return cached().source; }

void temp_dir_reset_for_testing() {
    // Deliberately bypasses the once_flag by replacing the cached value. Only
    // safe single-threaded, which is what the test harness guarantees.
    cached() = resolve();
}

}  // namespace flexaids
