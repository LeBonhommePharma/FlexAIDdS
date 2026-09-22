"""Source contract + behaviour for the non-throwing std::filesystem wrappers.

WHAT THIS PREVENTS
------------------
Companion to test_temp_dir_contract.py, same failure mode one level out. Most
of std::filesystem has a throwing overload and an error_code overload, and the
throwing one is what you get by writing the obvious thing. An uncaught
std::filesystem_error terminates the engine. On 2026-09-20 a single such call
(temp_directory_path) killed nine cells of an 85-target Astex campaign; a
tree-wide count taken afterwards found 48 more non-predicate throwing calls in
LIB/ sitting outside any try block.

29 of those 48 were converted to flexaids::fs_safe wrappers. The conversions
are NOT uniform and that is the point: the dangerous move would have been to
swap all 48 mechanically, because `remove(p)` and `remove(p, ec)` are
different programs and the second one is silent. These tests pin the
conversions that were made and, more importantly, pin the BEHAVIOUR each one
promised -- a wrapper that returns a wrong answer quietly is worse than the
abort it replaced.

WHY A SOURCE TEST AND A BEHAVIOUR TEST
--------------------------------------
The source test proves the throwing overload has not crept back into the
functions that were fixed. It cannot prove the replacement behaves. The
behaviour test compiles the real header and exercises the three paths that
actually bite: a file that vanishes, a cache SOURCE that vanishes (the operand
that had no existence guard anywhere in LIB), and a directory that cannot be
created.

WHY THE BEHAVIOUR TEST SHELLS OUT TO A COMPILER
-----------------------------------------------
tests/ is scanned by scripts/validate_sources.py, so a new .cpp that is not
referenced by a CMake target fails a STRICT configure -- and wiring one in
means editing the root CMakeLists.txt. This file therefore compiles its own
translation unit into a temp dir outside the repo. See MISSING COVERAGE below
for what that costs.

MISSING COVERAGE -- read this before trusting a green run
---------------------------------------------------------
  * No coverage of the 19 throwing calls deliberately LEFT in place (7
    fs::remove, 3 directory_iterator, 5 create_directories at CLI startup, 4
    absolute in benchmark_datasets.cpp which is another lane's file). Those
    are decisions, not oversights; see fs_hardening_inventory.csv.
  * No coverage of DatasetRunner.cpp's create_directories(pb_dir): triaged
    HIGH and knowingly left, because a silent skip there would leave a pose
    unvalidated by mandatory PoseBusters.
  * The behaviour test exercises the WRAPPERS, not the call sites. Nothing
    here proves DatasetRunner still elects the same pose; that needs the
    docking regression suite.
  * No concurrency coverage. The TOCTOU window between exists() and
    file_size() is narrowed by these changes, not closed, and no test here
    races two processes.
  * No coverage of a read-only or full filesystem (EROFS/ENOSPC); the
    ensure_dir failure is provoked with a path occupied by a regular file.
  * Source scanning is lexical. It strips comments and string literals but
    does not parse C++, so a call assembled by a macro would be invisible.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
LIB = REPO / "LIB"
HEADER = LIB / "fs_safe.h"

# Call sites converted in this migration: file -> functions that must no longer
# appear in their THROWING form anywhere in that file.
CONVERTED = {
    "LIB/DatasetRunner.cpp": ("file_size", "last_write_time", "canonical"),
    "LIB/LibrarySplitter.cpp": ("create_directories",),
    "LIB/ParallelCampaign.cpp": ("create_directories",),
    "LIB/PoseBust/BustCli.cpp": ("absolute",),
}


# ---------------------------------------------------------------------------
# lexical helpers
# ---------------------------------------------------------------------------
def _blank_comments_and_strings(text: str) -> str:
    """Replace comment and string-literal bytes with spaces, keeping offsets.

    Offsets are preserved so reported line numbers stay true to the file.
    """
    out = list(text)
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c == "/" and i + 1 < n and text[i + 1] == "/":
            while i < n and text[i] != "\n":
                out[i] = " "
                i += 1
        elif c == "/" and i + 1 < n and text[i + 1] == "*":
            out[i] = out[i + 1] = " "
            i += 2
            while i + 1 < n and not (text[i] == "*" and text[i + 1] == "/"):
                if text[i] != "\n":
                    out[i] = " "
                i += 1
            if i + 1 < n:
                out[i] = out[i + 1] = " "
                i += 2
        elif c in "\"'":
            quote = c
            out[i] = " "
            i += 1
            while i < n and text[i] != quote:
                if text[i] == "\\":
                    if text[i] != "\n":
                        out[i] = " "
                    i += 1
                if i < n and text[i] != "\n":
                    out[i] = " "
                i += 1
            if i < n:
                out[i] = " "
                i += 1
        else:
            i += 1
    return "".join(out)


def _args_of(code: str, open_paren: int) -> str | None:
    depth, i = 0, open_paren
    while i < len(code):
        if code[i] == "(":
            depth += 1
        elif code[i] == ")":
            depth -= 1
            if depth == 0:
                return code[open_paren + 1 : i]
        i += 1
    return None


def _last_arg(args: str) -> str:
    depth, cur, parts = 0, [], []
    for ch in args:
        if ch in "([{<":
            depth += 1
        elif ch in ")]}>":
            depth -= 1
        if ch == "," and depth == 0:
            parts.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    parts.append("".join(cur))
    return parts[-1].strip().lstrip("&*").strip()


def throwing_calls(rel: str, funcs) -> list[str]:
    """Every call to `funcs` in `rel` that does NOT pass a std::error_code."""
    text = (REPO / rel).read_text(encoding="utf-8", errors="replace")
    code = _blank_comments_and_strings(text)
    lines = text.splitlines()
    # names actually declared as error_code in this file -- avoids guessing
    ec_names = set(re.findall(r"(?:std\s*::\s*)?error_code\s+([A-Za-z_]\w*)", code))
    pattern = re.compile(
        r"(?<![\w:])(?:std\s*::\s*filesystem|fs)\s*::\s*(" + "|".join(funcs) + r")\s*\("
    )
    hits = []
    for m in pattern.finditer(code):
        args = _args_of(code, m.end() - 1)
        if args is None:
            continue
        last = _last_arg(args)
        if "error_code" in last or (re.fullmatch(r"[A-Za-z_]\w*", last) and last in ec_names):
            continue
        lineno = code.count("\n", 0, m.start()) + 1
        hits.append(f"{rel}:{lineno}: {lines[lineno - 1].strip()[:90]}")
    return hits


# ---------------------------------------------------------------------------
# source contract
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("rel,funcs", sorted(CONVERTED.items()))
def test_converted_sites_stay_converted(rel, funcs):
    """No throwing overload of a converted function may return to these files."""
    bad = throwing_calls(rel, funcs)
    assert not bad, (
        f"throwing std::filesystem call reintroduced in {rel}; "
        f"use flexaids::fs_safe instead:\n  " + "\n  ".join(bad)
    )


def test_the_source_guard_can_actually_fail():
    """Non-vacuity: the scanner must flag a throwing call that really exists.

    A guard never observed failing is not a guard. fs::remove() was triaged
    LOW and deliberately left throwing in DatasetRunner.cpp, so it is a real,
    checked-in example the scanner must catch. If this returns nothing, the
    regex or the file walk is broken and the tests above are vacuous.
    """
    hits = throwing_calls("LIB/DatasetRunner.cpp", ("remove",))
    assert hits, (
        "scanner found no throwing fs::remove() in DatasetRunner.cpp -- either "
        "someone converted them (update this test and the inventory CSV) or "
        "the scanner is broken and every other test here is vacuous"
    )


def test_scanner_does_not_flag_the_error_code_form():
    """Inverse non-vacuity: a call WITH an error_code must not be reported.

    Guards against the opposite failure -- a scanner that flags everything and
    would make the tests above unfixable rather than informative.
    """
    assert not throwing_calls("LIB/temp_dir.cpp", ("temp_directory_path", "remove")), (
        "scanner flagged the error_code overload in temp_dir.cpp as throwing"
    )


def test_header_documents_the_control_flow_decisions():
    """Rationale travels with the code; the next editor reads the header."""
    hdr = HEADER.read_text(encoding="utf-8")
    for token in (
        "DIFFERENT PROGRAMS",      # the semantic rule
        "WHAT IS DELIBERATELY NOT HERE",
        "remove()",                # why remove is not wrapped
        "current_path()",          # why absolute() is dangerous
        "STALE",                   # the conservative direction for caches
    ):
        assert token in hdr, f"fs_safe.h no longer documents {token!r}"


def test_no_wrapper_for_remove_exists():
    """remove() must stay unwrapped: absence is not an error for it.

    Pins the decision so a later 'completeness' pass cannot quietly add a
    remove wrapper and convert a refused delete into a no-op.
    """
    hdr = HEADER.read_text(encoding="utf-8")
    assert not re.search(r"inline\s+\w+\s+remove_or\b|inline\s+bool\s+safe_remove\b", hdr)


# ---------------------------------------------------------------------------
# behaviour
# ---------------------------------------------------------------------------
PROBE = r"""
#include "fs_safe.h"
#include <cassert>
#include <cstdio>
#include <fstream>
#include <filesystem>
namespace fs = std::filesystem;
using namespace flexaids::fs_safe;

static void write(const fs::path& p, const char* s) {
    std::ofstream o(p); o << s; o.close();
}

int main(int argc, char** argv) {
    const fs::path root = argv[1];
    int checks = 0;

    // 1. file_size_or: the vanishing file. This is the case the old
    //    `exists(p) && file_size(p)` idiom lost to a TOCTOU.
    const fs::path gone = root / "gone.txt";
    write(gone, "1234567890");
    assert(file_size_or(gone) == 10); ++checks;
    fs::remove(gone);
    assert(file_size_or(gone) == 0); ++checks;          // no throw, reads as absent
    assert(file_size_or(gone, 99) == 99); ++checks;     // fallback honoured

    // 1b. file_size_or on a DIRECTORY -- exists() says true, file_size throws.
    //     The old idiom's exists() guard did not cover this at all.
    assert(file_size_or(root) == 0); ++checks;

    // 2. out_of_date: the operand that had NO existence guard anywhere in LIB
    //    was the SOURCE. A pruned source must rebuild, not abort.
    const fs::path src = root / "src.pdb";
    const fs::path prod = root / "prod.pdb";
    write(src, "s");
    write(prod, "p");
    fs::last_write_time(prod, fs::last_write_time(src) + std::chrono::seconds(10));
    assert(out_of_date(prod, {src}) == false); ++checks;   // product newer -> fresh
    fs::last_write_time(prod, fs::last_write_time(src) - std::chrono::seconds(10));
    assert(out_of_date(prod, {src}) == true); ++checks;    // product older -> stale
    fs::remove(src);
    assert(out_of_date(prod, {src}) == true); ++checks;    // SOURCE gone -> stale, no abort
    write(src, "s");
    fs::remove(prod);
    assert(out_of_date(prod, {src}) == true); ++checks;    // product gone -> stale
    // never claims fresh from missing information:
    assert(out_of_date(root / "nope_a", {root / "nope_b"}) == true); ++checks;

    // 3. canonical_or / absolute_or on a path that does not exist.
    const std::string missing = (root / "no" / "such").string();
    assert(!canonical_or(missing).empty()); ++checks;      // no throw
    assert(!absolute_or("rel/path").empty()); ++checks;
    // an existing path still canonicalises to something real
    assert(fs::exists(canonical_or(root.string()))); ++checks;

    // 4. ensure_dir: success, idempotence, and the occupied-path failure.
    std::string why;
    assert(ensure_dir(root / "a" / "b", &why)); ++checks;
    assert(ensure_dir(root / "a" / "b", &why)); ++checks;  // already there -> still true
    const fs::path blocker = root / "blocker";
    write(blocker, "x");
    why.clear();
    assert(!ensure_dir(blocker, &why)); ++checks;          // file in the way -> false
    assert(!why.empty()); ++checks;                        // and says why

    std::printf("OK %d\n", checks);
    return 0;
}
"""


def _compiler():
    for cc in ("clang++", "g++"):
        if shutil.which(cc):
            return cc
    return None


@pytest.mark.skipif(_compiler() is None, reason="no C++ compiler on PATH")
def test_fs_safe_behaviour(tmp_path):
    """Compile the real header and exercise the three highest-risk paths."""
    src = tmp_path / "probe.cpp"
    src.write_text(textwrap.dedent(PROBE), encoding="utf-8")
    exe = tmp_path / "probe"
    compile_ = subprocess.run(
        [_compiler(), "-std=c++20", "-O0", f"-I{LIB}", str(src), "-o", str(exe)],
        capture_output=True, text=True,
    )
    assert compile_.returncode == 0, f"probe failed to compile:\n{compile_.stderr}"

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    run = subprocess.run([str(exe), str(scratch)], capture_output=True, text=True, timeout=60)
    assert run.returncode == 0, (
        f"fs_safe behaviour probe aborted (rc={run.returncode}); "
        f"stdout={run.stdout!r} stderr={run.stderr!r}"
    )
    assert run.stdout.startswith("OK "), run.stdout


@pytest.mark.skipif(_compiler() is None, reason="no C++ compiler on PATH")
def test_behaviour_probe_would_catch_a_broken_wrapper(tmp_path):
    """Non-vacuity for the behaviour test.

    Rebuild the probe against a header whose out_of_date() has been sabotaged
    in the one direction that matters -- reporting a cache FRESH when the
    source cannot be stat'ed. That is the silent-wrong-answer failure the
    wrapper exists to prevent, and the probe must fail on it.
    """
    broken = (tmp_path / "hdr")
    broken.mkdir()
    text = HEADER.read_text(encoding="utf-8")
    needle = "        if (sec) {"
    assert needle in text, "fs_safe.h shape changed; update this sabotage"
    sabotaged = text.replace(needle, "        if (sec) { return false;", 1)
    assert sabotaged != text
    (broken / "fs_safe.h").write_text(sabotaged, encoding="utf-8")

    src = tmp_path / "probe.cpp"
    src.write_text(textwrap.dedent(PROBE), encoding="utf-8")
    exe = tmp_path / "probe_broken"
    compile_ = subprocess.run(
        [_compiler(), "-std=c++20", "-O0", f"-I{broken}", str(src), "-o", str(exe)],
        capture_output=True, text=True,
    )
    assert compile_.returncode == 0, compile_.stderr
    scratch = tmp_path / "scratch2"
    scratch.mkdir()
    run = subprocess.run([str(exe), str(scratch)], capture_output=True, text=True, timeout=60)
    assert run.returncode != 0, (
        "probe PASSED against a header that reports a cache fresh when its "
        "source cannot be stat'ed -- the behaviour test is vacuous"
    )
