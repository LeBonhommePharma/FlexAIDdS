"""Source contract for temporary-directory resolution.

WHAT THIS PREVENTS
------------------
std::filesystem::temp_directory_path() has a throwing overload. When $TMPDIR
names a directory that no longer exists it does not return an error -- it
aborts the process:

    Fatal error: filesystem error: in temp_directory_path: <path>

On 2026-09-20 that destroyed nine cells of an 85-target Astex campaign. $TMPDIR
pointed into an agent session workspace, the harness deleted that workspace
after a few hours idle, and every engine invocation launched afterwards died at
startup having written no poses. The driver surfaced it only as "Incomplete
docking", so it read like a docking failure on nine specific targets. It was
not: the same nine had docked normally in an earlier campaign on the same
engine, and the distribution of casualties was simply whichever cells happened
to start after the sweep.

The repair was a single validated resolver, flexaids::temp_dir(), which never
throws and probes writability by creating a file rather than trusting
metadata. These tests keep it that way: a future edit that reintroduces a bare
temp_directory_path() call fails here instead of in a campaign six hours after
someone walks away.

WHY A SOURCE TEST AND NOT ONLY A UNIT TEST
------------------------------------------
A unit test proves the resolver behaves. It cannot prove the resolver is the
only path to a temp directory -- that is a property of the whole source tree,
so it needs a tree-wide assertion.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
LIB = REPO / "LIB"

# The resolver itself must call the std function -- that is its job. CifReader
# uses the non-throwing error_code overload, which is already correct and is
# deliberately left alone rather than churned.
ALLOWED = {"LIB/temp_dir.cpp", "LIB/temp_dir.h", "LIB/CifReader.cpp"}

# A call is SAFE when the error_code overload is used: temp_directory_path(ec).
# It is UNSAFE when called with no arguments, which is the throwing form.
THROWING_CALL = re.compile(r"temp_directory_path\s*\(\s*\)")


def _cpp_sources():
    return sorted(p for p in LIB.rglob("*") if p.suffix in {".cpp", ".h", ".hpp", ".cc"})


def _offenders(extra_allowed=frozenset()):
    allowed = ALLOWED | set(extra_allowed)
    bad = []
    for p in _cpp_sources():
        rel = p.relative_to(REPO).as_posix()
        if rel in allowed:
            continue
        text = p.read_text(encoding="utf-8", errors="replace")
        for i, line in enumerate(text.splitlines(), 1):
            stripped = line.lstrip()
            if stripped.startswith("//") or stripped.startswith("*"):
                continue          # a comment naming the function is not a call
            if THROWING_CALL.search(line):
                bad.append(f"{rel}:{i}: {line.strip()[:90]}")
    return bad


def test_no_throwing_temp_directory_path_outside_the_resolver():
    """The throwing overload must not appear anywhere but the resolver."""
    bad = _offenders()
    assert not bad, (
        "throwing temp_directory_path() found; use flexaids::temp_dir() instead:\n  "
        + "\n  ".join(bad)
    )


def test_the_guard_can_actually_fail():
    """Non-vacuity: the guard must flag the resolver itself once un-allowlisted.

    A guard that has never been observed failing is not a guard. temp_dir.cpp
    genuinely contains the pattern (inside a comment and as the error_code
    form), so removing it from the allowlist must produce a hit -- otherwise
    the regex or the file walk is broken and the passing test above is vacuous.
    """
    text = (LIB / "temp_dir.cpp").read_text(encoding="utf-8")
    assert "temp_directory_path" in text, "resolver no longer calls the std function at all"
    # the resolver must use the ERROR_CODE form, never the throwing one, in code
    code_lines = [
        ln for ln in text.splitlines()
        if "temp_directory_path" in ln and not ln.lstrip().startswith(("//", "*"))
    ]
    assert code_lines, "no non-comment call found in the resolver"
    assert all("ec" in ln or "error" in ln for ln in code_lines), (
        "resolver itself uses the THROWING overload:\n  " + "\n  ".join(code_lines))


def test_every_former_call_site_now_uses_the_resolver():
    """The eight migrated sites must reference flexaids::temp_dir()."""
    migrated = [
        "LIB/top.cpp", "LIB/direct_input.cpp", "LIB/read_input.cpp",
        "LIB/LibrarySplitter.cpp", "LIB/PoseBust/ChecksChemistry.cpp",
        "LIB/PoseBust/Engine.cpp",
    ]
    missing = []
    for rel in migrated:
        text = (REPO / rel).read_text(encoding="utf-8", errors="replace")
        if "flexaids::temp_dir()" not in text:
            missing.append(rel)
        if "temp_dir.h" not in text:
            missing.append(f"{rel} (no include)")
    assert not missing, f"sites not migrated: {missing}"


def test_resolver_is_in_the_core_library_build():
    """A missing link is better than a silent fallback to the throwing call."""
    lib_cmake = (LIB / "CMakeLists.txt").read_text(encoding="utf-8")
    assert "temp_dir.cpp" in lib_cmake, "temp_dir.cpp not in LIB/CMakeLists.txt"


def test_resolver_documents_why_it_exists():
    """Rationale must travel with the code, not only with the test.

    The next person editing temp_dir.cpp sees the header, not this file.
    """
    hdr = (LIB / "temp_dir.h").read_text(encoding="utf-8")
    for token in ("never throws", "WRITABLE", "FALLBACK CHAIN", "FLEXAIDDS_TMPDIR"):
        assert token in hdr, f"header no longer documents {token!r}"
