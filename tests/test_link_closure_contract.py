"""Every unit-test fixture that compiles a caller of flexaids::temp_dir() must
also compile LIB/temp_dir.cpp.

WHAT THIS PREVENTS
------------------
On 2026-09-22 the PoseBust migration added flexaids::temp_dir() calls to
LIB/PoseBust/ChecksChemistry.cpp and LIB/PoseBust/Engine.cpp. The fixture
test_posebust compiles both of those translation units but did not list
LIB/temp_dir.cpp, so it failed to LINK:

    Undefined symbols for architecture arm64:
      "flexaids::temp_dir()", referenced from:
          flexaids::posebust::check_chemistry_sanity(...)
          flexaids::posebust::validate_elected_pose(...)

That broke three CI jobs -- macOS CPU tests (blocking), tsan-linux and
Coverage Analysis -- ALL at the Build step, because each builds every target.
The Release job passed throughout: it builds named targets only and never
reaches test_posebust. A named-target build is therefore NOT a check on link
closure, which is why this guard reads the build graph instead.

flexaids_add_unit_test links only GTest::gtest_main plus explicit
LINK_LIBRARIES -- there is no core library to inherit definitions from -- so a
fixture's SOURCES list must be link-closed on its own.

fs_safe.h needs no such rule: its helpers are inline in the header. Only
temp_dir has a separate definition TU, so only temp_dir can produce this.

SCOPE
-----
Only flexaids_add_unit_test() fixtures are checked. Bare add_executable()
targets that compile the PoseBust callers (test_dataset_runner,
test_cofactor_blacklist, benchmark_datasets) obtain the definition through
target_link_libraries(... flexaid_core) -- the OBJECT library that compiles
LIB/temp_dir.cpp -- which a static reading of their ${VAR} source lists cannot
see. The full-target build (cmake --build with NO --target) is the check for
those. Calls made from a header (an inline helper) are also out of scope: only
.cpp files under LIB/ and tests/ are scanned for call sites.
"""
from __future__ import annotations

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CMAKE = REPO / "CMakeLists.txt"
SCAN_DIRS = ("LIB", "tests")

# Symbols that are DECLARED in a header but DEFINED in a separate .cpp, so a
# caller alone is not link-closed. Extend this map if another such pair appears.
DEFINED_IN = {"flexaids::temp_dir": "LIB/temp_dir.cpp"}

# The cmake_parse_arguments keyword sets of flexaids_add_unit_test, verbatim
# from cmake/FlexAIDHelpers.cmake (options / oneValueArgs / multiValueArgs).
# Any of these ENDS the SOURCES list. TEST_NAME is a oneValueArgs, so a file
# listed after it becomes an UNPARSED argument and add_executable never sees
# it. The first version of this guard searched the whole block and therefore
# PASSED on exactly that mistake -- LIB/temp_dir.cpp sat after TEST_NAME, was
# silently dropped, and test_posebust still failed to link.
_KEYWORDS = (
    "CONFIGURE_SIMD", "LINK_OPENMP", "GTEST_MAIN_ONLY", "NO_DEFAULT_COMPILE_OPTS",
    "TEST_NAME",
    "SOURCES", "LINK_LIBRARIES", "DEFINES", "INCLUDES", "COMPILE_OPTIONS",
)

_FIXTURE = re.compile(r"flexaids_add_unit_test\(\s*(\w+)(.*?)\n\s*\)", re.S)
_SRC_PREFIXES = ("${CMAKE_CURRENT_SOURCE_DIR}/", "${CMAKE_SOURCE_DIR}/")


def _code_tokens(body: str) -> list[str]:
    """Non-empty lines of a fixture body with `#` comments removed."""
    out = []
    for raw in body.splitlines():
        tok = raw.split("#")[0].strip()
        if tok:
            out.append(tok)
    return out


def _norm(path: str) -> str:
    for pre in _SRC_PREFIXES:
        if path.startswith(pre):
            return path[len(pre):]
    return path


def _sources_section(body: str) -> list[str]:
    """The files add_executable will actually receive from one fixture body.

    Walks the body line by line. Collection starts at SOURCES and STOPS at the
    next cmake_parse_arguments keyword; anything after TEST_NAME is an
    unparsed argument, not a source, and is not returned.

    THIS is the function the placement tests below exercise. _blocks() must
    keep calling it -- a private copy of the loop inside a test would test the
    copy, not the parser the real checks run on.
    """
    toks = _code_tokens(body)
    if not any(t.split()[0] == "SOURCES" for t in toks):
        # bare form: flexaids_add_unit_test(name src1 src2 ...) with no keyword
        return [_norm(w) for t in toks if t.split()[0] not in _KEYWORDS
                for w in t.split()]
    srcs: list[str] = []
    collecting = False
    for tok in toks:
        words = tok.split()
        head = words[0]
        if head == "SOURCES":
            collecting = True
            srcs.extend(_norm(w) for w in words[1:])
            continue
        if head in _KEYWORDS:
            collecting = False
            continue
        if collecting:
            srcs.extend(_norm(w) for w in words)
    return srcs


def _blocks(text: str | None = None) -> list[tuple[str, list[str]]]:
    """(fixture name, sources) for every fixture in CMakeLists.txt.

    `text` lets a test run the PRODUCTION parser on mutated CMake text.
    """
    s = CMAKE.read_text(encoding="utf-8") if text is None else text
    return [(m.group(1), _sources_section(m.group(2))) for m in _FIXTURE.finditer(s)]


def _strip_c_comments(t: str) -> str:
    t = re.sub(r"/\*.*?\*/", "", t, flags=re.S)
    return re.sub(r"//[^\n]*", "", t)


def _call_pattern(symbol: str) -> re.Pattern[str]:
    """Matches `flexaids::temp_dir(` AND an unqualified `temp_dir(` made from
    inside the namespace; rejects member access (`.temp_dir(`, `->temp_dir(`)."""
    last = symbol.rsplit("::", 1)[-1]
    return re.compile(r"(?<![\w.>])" + re.escape(last) + r"\s*\(")


def _callers(symbol: str) -> list[str]:
    """Translation units under SCAN_DIRS that call `symbol` (comments stripped).
    The definition TU itself is not a caller."""
    pat = _call_pattern(symbol)
    out = []
    for d in SCAN_DIRS:
        for p in sorted((REPO / d).rglob("*.cpp")):
            rel = p.relative_to(REPO).as_posix()
            if rel == DEFINED_IN.get(symbol):
                continue
            try:
                t = p.read_text(encoding="utf-8", errors="replace")
            except OSError:
                continue
            if pat.search(_strip_c_comments(t)):
                out.append(rel)
    return out


def _violations(blocks: list[tuple[str, list[str]]]) -> list[tuple[str, list[str], str]]:
    """(fixture, callers it compiles, missing definition TU) for every fixture
    that will fail to link."""
    problems = []
    for symbol, defn in DEFINED_IN.items():
        callers = _callers(symbol)
        assert callers, f"no caller of {symbol} found -- the scan is vacuous"
        for name, srcs in blocks:
            compiled = [c for c in callers if c in srcs]
            if compiled and defn not in srcs:
                problems.append((name, compiled, defn))
    return problems


# ── the rule ────────────────────────────────────────────────────────────────

def test_every_fixture_compiling_a_caller_also_defines_the_symbol():
    problems = _violations(_blocks())
    assert not problems, "fixtures that will fail to link: " + "; ".join(
        f"{n} compiles {c} without {d}" for n, c, d in problems
    )


def test_test_posebust_is_wired_for_it():
    """The specific fixture from the incident, pinned."""
    hit = [srcs for n, srcs in _blocks() if n == "test_posebust"]
    assert hit, "test_posebust fixture not found in CMakeLists.txt"
    assert "LIB/temp_dir.cpp" in hit[0], "test_posebust lost its temp_dir.cpp source"


def test_fs_safe_is_header_only_so_needs_no_rule():
    """If fs_safe ever gains a .cpp, this test fails and DEFINED_IN must grow."""
    assert not (REPO / "LIB" / "fs_safe.cpp").exists(), (
        "LIB/fs_safe.cpp now exists -- add 'fs_safe::' to DEFINED_IN or fixtures "
        "compiling its callers will fail to link exactly as test_posebust did"
    )


# ── non-vacuity: the instruments must detect a present case ─────────────────

def test_the_scan_finds_the_known_callers():
    """The PoseBust TUs that caused the incident, and the fixture's own test
    file, must be detected as callers."""
    callers = _callers("flexaids::temp_dir")
    for expect in ("LIB/PoseBust/ChecksChemistry.cpp", "LIB/PoseBust/Engine.cpp",
                   "tests/test_temp_dir.cpp"):
        assert expect in callers, f"{expect} not detected as a caller (scan broken)"


def test_the_scan_ignores_comment_mentions_and_member_access():
    pat = _call_pattern("flexaids::temp_dir")
    assert pat.search("x = flexaids::temp_dir();")
    assert pat.search("x = temp_dir();")                       # inside namespace
    assert not pat.search("lib.temp_dir = other;")             # member, no call
    assert not pat.search("cfg->temp_dir(1);")                 # member call
    assert not pat.search(_strip_c_comments("// flexaids::temp_dir() is called\n"))
    assert not pat.search(_strip_c_comments("/* temp_dir() */"))


def test_the_parser_sees_every_fixture_in_the_file():
    """Fixture count by regex must equal the count of call lines, so a fixture
    the parser cannot see cannot hide a violation."""
    s = CMAKE.read_text(encoding="utf-8")
    call_lines = [l for l in s.splitlines()
                  if l.split("#")[0].strip().startswith("flexaids_add_unit_test(")]
    assert len(_blocks(text=s)) == len(call_lines) > 0


def test_a_source_after_TEST_NAME_is_not_counted_as_a_source():
    """PLACEMENT, not mere presence -- run on the PRODUCTION parser.

    Simulates the broken wiring (LIB/temp_dir.cpp listed after TEST_NAME) and
    requires _blocks() to NOT report it as a source. A parser that greps the
    whole block passes the broken wiring and is worthless here.
    """
    broken = """
    flexaids_add_unit_test(test_posebust
        SOURCES
            tests/test_posebust.cpp
            LIB/PoseBust/Engine.cpp
        INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/LIB/PoseBust
        TEST_NAME PoseBustTests
        LIB/temp_dir.cpp
    )
"""
    [(name, srcs)] = _blocks(text=broken)
    assert name == "test_posebust"
    assert "LIB/temp_dir.cpp" not in srcs, (
        "the parser counted a file listed after TEST_NAME as a source; it would "
        "pass the exact wiring that broke three CI jobs"
    )
    assert srcs == ["tests/test_posebust.cpp", "LIB/PoseBust/Engine.cpp"], srcs

    # positive control: the same file BEFORE the keywords is a source
    fixed = broken.replace("        LIB/temp_dir.cpp\n", "").replace(
        "            LIB/PoseBust/Engine.cpp\n",
        "            LIB/PoseBust/Engine.cpp\n            LIB/temp_dir.cpp\n")
    [(_, srcs_ok)] = _blocks(text=fixed)
    assert "LIB/temp_dir.cpp" in srcs_ok, "parser dropped a correctly placed source"


def _posebust_block(text: str) -> str:
    m = re.search(r"flexaids_add_unit_test\(\s*test_posebust\b.*?\n\s*\)", text, re.S)
    assert m, "test_posebust fixture not found"
    return m.group(0)


def test_the_broken_placement_in_the_real_file_is_reported():
    """End-to-end injection on the REAL CMakeLists.txt text: move the
    LIB/temp_dir.cpp line to after TEST_NAME and the rule must go red on
    test_posebust. The unmodified file must stay green (positive control)."""
    s = CMAKE.read_text(encoding="utf-8")
    block = _posebust_block(s)
    without = re.sub(r"\n[ \t]*LIB/temp_dir\.cpp[ \t]*(?=\n)", "", block, count=1)
    assert without != block, "could not remove the LIB/temp_dir.cpp source line"
    # Anchored to a CODE line: the fixture's own comment block also says
    # "TEST_NAME", and an unanchored pattern re-inserted the file INSIDE
    # SOURCES after that comment -- a valid placement, so the injection was
    # not injecting anything. Instruments get checked too.
    moved, n = re.subn(r"(?m)^([ \t]*TEST_NAME\b[^\n]*\n)",
                       r"\1            LIB/temp_dir.cpp\n", without, count=1)
    assert n == 1, "TEST_NAME keyword line not found in test_posebust"
    [(_, moved_srcs)] = _blocks(text=moved)
    assert "LIB/temp_dir.cpp" not in moved_srcs, (
        "injection failed: the moved file still parses as a source"
    )
    mutated = s.replace(block, moved, 1)

    red = [n for n, _, _ in _violations(_blocks(text=mutated))]
    assert "test_posebust" in red, (
        "the guard did not report LIB/temp_dir.cpp placed after TEST_NAME -- "
        "it would pass the wiring that broke CI"
    )
    assert not _violations(_blocks(text=s)), "positive control: real file is not clean"


def test_a_missing_definition_in_the_real_file_is_reported():
    """End-to-end injection on the REAL CMakeLists.txt text: delete the
    LIB/temp_dir.cpp line and the rule must go red on test_posebust."""
    s = CMAKE.read_text(encoding="utf-8")
    block = _posebust_block(s)
    without = re.sub(r"\n[ \t]*LIB/temp_dir\.cpp[ \t]*(?=\n)", "", block, count=1)
    assert without != block, "could not remove the LIB/temp_dir.cpp source line"
    red = [n for n, _, _ in _violations(_blocks(text=s.replace(block, without, 1)))]
    assert "test_posebust" in red, "the guard did not report a missing LIB/temp_dir.cpp"
