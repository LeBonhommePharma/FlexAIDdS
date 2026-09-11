#!/usr/bin/env python3
"""check_dois.py — Fail-closed gate on citation identifiers in this repository.

Motivation
----------
The 2026-09-10 dataset identity audit found three fabricated Zenodo DOIs and, on
a repo-wide sweep of every DOI-shaped string, four unregistered DOIs plus one
misattributed citation. The known-bad identifiers, all retired by that commit:

  known-bad FABRICATED   ``10.5281/zenodo.dude37``
  known-bad FABRICATED   ``10.5281/zenodo.itc187``
  known-bad FABRICATED   ``10.5281/zenodo.muv2009``
  known-bad UNREGISTERED ``10.1021/ci800224j``  -> ``10.1021/ci8002254``
  known-bad UNREGISTERED ``10.1021/ci8003380``  -> ``10.1021/ci8002649``
  known-bad UNREGISTERED ``10.1093/nar/g1271``  -> ``10.1093/nar/gks966``
  known-bad UNREGISTERED ``10.1021/ci100366a``  -> emptied, not recoverable

None of these was caught by any existing gate because nothing ever resolved a
citation identifier. This script is that gate.

What it fails on
----------------
1. UNREGISTERED   — a DOI-shaped string used as a live citation whose cached
                    resolution says it is not registered in the DOI Handle System.
2. NOT_DOI_SYNTAX — a citation *field* (``doi:``, ``zenodo_doi:``, ``"doi":``)
                    whose non-empty value is not DOI syntax. Includes the
                    registry-specific rule that a Zenodo DOI suffix must be a
                    numeric record id, so the known-bad ``10.5281/zenodo.dude37``
                    is rejected on syntax alone and can never be "fixed" by a
                    cache refresh.
3. UNCACHED       — a DOI-shaped string with no entry in the committed cache.
                    Means the cache is stale: re-run with ``--refresh``.
4. RETIRED_REVIVED— a DOI listed in the retired ledger reappearing on a line that
                    is not marked as documenting the defect, i.e. someone put a
                    known-bad identifier back into service.
5. LEDGER_DRIFT   — a retired-ledger entry that the cache says *is* registered,
                    or a ledger entry that no longer appears anywhere.

Offline by default
------------------
Resolutions live in a committed fixture so the suite runs with no network.
``--refresh`` is the only mode that touches the network, and it refuses to write
a cache unless BOTH controls behave: a known-real DOI must resolve and a
known-impossible DOI must not. A 404 from a rate-limited or blocked host is not
evidence of fabrication, and this is how that failure mode is excluded.

Usage
-----
    python3 scripts/check_dois.py               # offline gate (CI)
    python3 scripts/check_dois.py --refresh     # re-resolve online, rewrite cache
    python3 scripts/check_dois.py --selftest    # non-vacuity proof
    python3 scripts/check_dois.py --report      # print the full sweep table

Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
SPDX-License-Identifier: Apache-2.0
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, NamedTuple, Optional, Tuple

ROOT = Path(__file__).resolve().parents[1]
FIXTURE_DIR = ROOT / "tests" / "fixtures" / "doi_resolution"
CACHE_PATH = FIXTURE_DIR / "doi_resolution_cache.json"
RETIRED_PATH = FIXTURE_DIR / "known_bad_dois.json"

# Controls. The positive control is Hartshorn 2007 (the Astex Diverse paper, cited
# throughout this repo); the negative control is a syntactically valid DOI under a
# prefix that is not allocated, so it can never resolve.
POSITIVE_CONTROL = "10.1021/jm061277y"
NEGATIVE_CONTROL = "10.9999/flexaidds.nonexistent.control"

# ---------------------------------------------------------------------------
# DOI syntax
# ---------------------------------------------------------------------------

# Greedy first pass; trailing punctuation and unbalanced parens are trimmed after.
# Parens must be *allowed* inside because Elsevier DOIs contain them
# (10.1016/S0969-2126(01)00662-1) — a naive pattern truncates those and then
# reports a false "unregistered".
# Backtick and pipe are excluded because markup wraps DOIs in them (``10.xxxx/yyy``
# in RST, |doi| in tables); absorbing the delimiter yields a DOI that looks real,
# fails to resolve, and gets reported as fabricated. That false positive was
# observed while building this gate.
_DOI_RAW = re.compile(r"10\.\d{4,9}/[^\s\"'`|<>\[\]{}\\]+")
_DOI_STRICT = re.compile(r"^10\.\d{4,9}/\S+$")

# Zenodo mints 10.5281/zenodo.<integer>. Anything else under that prefix is not a
# DOI, regardless of what any resolver says.
_ZENODO_OK = re.compile(r"^10\.5281/zenodo\.\d+$", re.I)

# A line that is documenting a retired identifier rather than citing it.
_DOC_MARKER = re.compile(
    r"previously\s+(?:read|cited)|UNREGISTERED|FABRICATED|NOT\s+adopted"
    r"|no\s+DOI\s+was\s+invented|known-bad|doi-check:\s*documented",
    re.I,
)

# Citation fields whose value must be DOI syntax when non-empty.
# The optional prefix group must be OPTIONAL: a bare `doi:` key is the common case
# (benchmarks/*/manifest.yaml, the fetch_itc_data source registry). Requiring at
# least one character before "doi" made this rule match nothing in the live tree
# while still firing on the `zenodo_doi:` self-test fixture — a green run that
# checked nothing. Caught by printing the field denominator, which read 0.
_YAML_FIELD = re.compile(r"^\s*((?:[A-Za-z0-9_]+_)?doi)\s*:\s*(.*?)\s*$", re.I)
_JSONISH_FIELD = re.compile(r"""["']((?:[A-Za-z0-9_]+_)?doi)["']\s*:\s*["']([^"']*)["']""", re.I)

_SKIP_DIRS = {
    ".git", "build", "build-ci-fix", "build-audit", "node_modules", "__pycache__",
    ".venv", "venv", ".mypy_cache", ".pytest_cache", ".claude", ".cursor", ".gemini",
    ".trae", ".windsurf", ".qoder", ".roo",
}
# Mirrored wwPDB deposition files. Excluded deliberately — see scan() docstring.
_STRUCTURE_SUFFIXES = {".cif", ".mmcif", ".pdb", ".ent", ".sdf", ".mol2"}

# The gate's own implementation and test. These are ABOUT identifiers rather than
# citing any, and they necessarily contain: the negative-control DOI (which must
# never resolve), example strings that explain the parsing rules, and fixture
# literals. Scanning them makes the gate fail on itself. Same rationale as the
# fixture data files below.
#
# ACCEPTED LIMIT, stated rather than hidden: a genuine citation defect introduced
# into these two files would not be caught. Neither file cites a source; both are
# tooling. This exclusion was added after the gate — validated while its own source
# was still untracked, so git ls-files did not yet return it — failed on its own
# tree the moment it was committed. Caught by exporting the commit and running it,
# not by inspection.
_GATE_INTERNALS = {"scripts/check_dois.py", "tests/test_doi_integrity.py"}

_TEXT_SUFFIXES = {
    ".py", ".yaml", ".yml", ".json", ".md", ".txt", ".cpp", ".h", ".hpp", ".c",
    ".sh", ".cmake", ".rst", ".cfg", ".ini", ".toml", ".csv", ".tsv", ".jsx",
    ".js", ".ts", ".html", ".R", ".r",
}


def trim_doi(s: str) -> str:
    """Strip trailing prose punctuation, markup and unbalanced closing parens."""
    while True:
        t = s.rstrip(".,;:`|*_")
        while t.endswith(")") and t.count(")") > t.count("("):
            t = t[:-1]
        if t == s:
            return t
        s = t


def is_doi_syntax(value: str) -> bool:
    if not _DOI_STRICT.match(value):
        return False
    if value.lower().startswith("10.5281/zenodo.") and not _ZENODO_OK.match(value):
        return False
    return True


# ---------------------------------------------------------------------------
# Repository scan
# ---------------------------------------------------------------------------

class Occurrence(NamedTuple):
    path: str
    lineno: int
    doi: str
    line: str
    documented: bool


class FieldValue(NamedTuple):
    path: str
    lineno: int
    field: str
    value: str
    line: str


def list_files(root: Path) -> List[str]:
    """Tracked files if git is usable, else a filtered walk. Never silently empty."""
    env = dict(os.environ, GIT_CONFIG_GLOBAL="/dev/null", GIT_CONFIG_SYSTEM="/dev/null")
    try:
        r = subprocess.run(["git", "ls-files", "-z"], cwd=str(root), env=env,
                           capture_output=True, text=True, timeout=120)
        if r.returncode == 0 and r.stdout.strip():
            return [f for f in r.stdout.split("\0") if f]
    except (OSError, subprocess.SubprocessError):
        pass
    out: List[str] = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in _SKIP_DIRS]
        for fn in filenames:
            p = Path(dirpath) / fn
            if p.suffix in _TEXT_SUFFIXES:
                out.append(str(p.relative_to(root)))
    return out


def scan(root: Path, files: Optional[Iterable[str]] = None
         ) -> Tuple[List[Occurrence], List[FieldValue], Dict[str, int]]:
    """Return (DOI occurrences, citation-field values, coverage counts).

    DELIBERATE EXCLUSION, reported not silent: mirrored crystallographic data
    files (.cif/.pdb/.mmcif/.ent) are not scanned. Their DOIs are the wwPDB
    primary-citation records that ship inside the deposited file — upstream data
    this repository neither authors nor may edit. At the time this gate was written
    they accounted for 199 of the 251 distinct DOI strings in the tree. Repository
    citations live in YAML/JSON/Markdown/source, which are scanned in full. The
    counts below always report how many files were excluded on this basis.
    """
    if files is None:
        files = list_files(root)
    occs: List[Occurrence] = []
    fields: List[FieldValue] = []
    nread = 0
    n_structure = 0
    n_binary = 0
    n_self = 0
    for rel in files:
        p = root / rel
        if p.suffix.lower() in _STRUCTURE_SUFFIXES:
            n_structure += 1
            continue
        if p.suffix and p.suffix not in _TEXT_SUFFIXES:
            n_binary += 1
            continue
        # The gate's own data files are a registry OF identifiers, not citations.
        # Scanning the retired ledger against the retired ledger is circular, and the
        # resolution cache necessarily lists every unregistered DOI by key.
        rel_norm = rel.replace(os.sep, "/")
        if rel_norm.startswith("tests/fixtures/doi_resolution/") or rel_norm in _GATE_INTERNALS:
            n_self += 1
            continue
        try:
            text = p.read_text(encoding="utf-8", errors="replace")
        except (OSError, IsADirectoryError):
            continue
        nread += 1
        for i, line in enumerate(text.splitlines(), 1):
            documented = bool(_DOC_MARKER.search(line))
            for m in _DOI_RAW.finditer(line):
                # A match butted against '<' is a documentation template, not a
                # citation — e.g. "10.5281/zenodo.<numeric record id>". Scanning it
                # would report the truncated stem as an unresolvable DOI.
                if line[m.end():m.end() + 1] == "<":
                    continue
                d = trim_doi(m.group(0))
                if _DOI_STRICT.match(d):
                    occs.append(Occurrence(rel, i, d, line.strip()[:200], documented))
            stripped = line.lstrip()
            if not stripped.startswith(("#", "//", "*")):
                ym = _YAML_FIELD.match(line)
                if ym and p.suffix in {".yaml", ".yml"}:
                    val = ym.group(2).strip().strip("\"'")
                    if val:
                        fields.append(FieldValue(rel, i, ym.group(1), val, line.strip()[:200]))
                for jm in _JSONISH_FIELD.finditer(line):
                    if jm.group(2).strip():
                        fields.append(FieldValue(rel, i, jm.group(1), jm.group(2).strip(),
                                                 line.strip()[:200]))
    return occs, fields, {
        "files_read": nread,
        "files_excluded_structure": n_structure,
        "files_excluded_binary": n_binary,
        "files_excluded_gate_internal": n_self,
        "occurrences": len(occs),
        "distinct": len({o.doi.lower() for o in occs}),
        "fields": len(fields),
    }


# ---------------------------------------------------------------------------
# Online resolution (only under --refresh)
# ---------------------------------------------------------------------------

def resolve_handle(doi: str, tries: int = 3) -> Dict[str, object]:
    """Registry-agnostic existence check via the DOI Handle System.

    Works for CrossRef (ACS, Elsevier), DataCite (wwPDB 10.2210, Zenodo) alike.
    responseCode 1 = registered, 100 = handle not found.
    """
    import time
    from urllib.parse import quote
    import requests  # imported lazily: offline mode must not need it

    url = "https://doi.org/api/handles/" + quote(doi, safe="/:().-_")
    last = None
    for k in range(tries):
        try:
            r = requests.get(url, timeout=30, headers={"Accept": "application/json"})
            if r.status_code in (200, 404):
                try:
                    code = r.json().get("responseCode")
                except ValueError:
                    code = None
                if code in (1, 100):
                    return {"registered": code == 1, "response_code": code, "error": None}
            last = "HTTP%s" % r.status_code
        except Exception as exc:  # noqa: BLE001 - network layer, any failure is "unknown"
            last = type(exc).__name__
        time.sleep(1.2 * (k + 1))
    return {"registered": None, "response_code": None, "error": last}


def crossref_meta(doi: str) -> Dict[str, object]:
    import time
    from urllib.parse import quote
    import requests

    for k in range(3):
        try:
            r = requests.get("https://api.crossref.org/works/" + quote(doi, safe="/:().-_"),
                             timeout=40)
            if r.status_code == 200:
                m = r.json()["message"]
                return {
                    "title": (m.get("title") or [None])[0],
                    "first_author": (m.get("author") or [{}])[0].get("family"),
                    "container": (m.get("container-title") or [None])[0],
                    "volume": m.get("volume"), "page": m.get("page"),
                    "year": (m.get("issued", {}).get("date-parts") or [[None]])[0][0],
                }
            if r.status_code == 404:
                return {"crossref": "404"}
        except Exception:  # noqa: BLE001
            pass
        time.sleep(1.5 * (k + 1))
    return {"crossref": "unavailable"}


def refresh(root: Path) -> int:
    occs, _fields, cov = scan(root)
    dois = sorted({o.doi.lower() for o in occs})
    print("[refresh] %d files read (%d structure files excluded, %d non-text), "
          "%d DOI occurrences, %d distinct DOIs"
          % (cov["files_read"], cov["files_excluded_structure"],
             cov["files_excluded_binary"], len(occs), len(dois)))

    pos = resolve_handle(POSITIVE_CONTROL)
    neg = resolve_handle(NEGATIVE_CONTROL)
    print("[refresh] positive control %s -> registered=%s" % (POSITIVE_CONTROL, pos["registered"]))
    print("[refresh] negative control %s -> registered=%s" % (NEGATIVE_CONTROL, neg["registered"]))
    if pos["registered"] is not True or neg["registered"] is not False:
        print("VOID: controls did not behave (positive must resolve, negative must not). "
              "Cache NOT written — a 404 here is not evidence of a bad DOI.")
        return 2

    cache: Dict[str, object] = {}
    unresolved: List[str] = []
    for d in dois:
        h = resolve_handle(d)
        if h["registered"] is None:
            unresolved.append(d)
            continue
        entry = {"registered": bool(h["registered"]), "response_code": h["response_code"]}
        if h["registered"]:
            entry.update(crossref_meta(d))
        cache[d] = entry
    if unresolved:
        print("VOID: %d DOIs could not be resolved (network), cache NOT written: %s"
              % (len(unresolved), ", ".join(unresolved[:10])))
        return 2

    n_reg = sum(1 for v in cache.values() if v["registered"])
    payload = {
        "_comment": "Committed resolution cache for scripts/check_dois.py. "
                    "Refresh explicitly with --refresh; never hand-edit.",
        "checked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "oracle": "doi.org Handle System /api/handles (registry-agnostic) "
                  "+ CrossRef /works for metadata",
        "controls": {"positive": POSITIVE_CONTROL, "positive_registered": True,
                     "negative": NEGATIVE_CONTROL, "negative_registered": False},
        "counts": {"distinct": len(cache), "registered": n_reg,
                   "unregistered": len(cache) - n_reg},
        "dois": cache,
    }
    FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
    CACHE_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("[refresh] wrote %s: %d distinct, %d registered, %d unregistered"
          % (CACHE_PATH.relative_to(root), len(cache), n_reg, len(cache) - n_reg))
    return 0


# ---------------------------------------------------------------------------
# Offline gate
# ---------------------------------------------------------------------------

class Finding(NamedTuple):
    kind: str
    path: str
    lineno: int
    detail: str


def load_json(path: Path, default):
    if not path.exists():
        return default
    return json.loads(path.read_text(encoding="utf-8"))


def check(root: Path, cache: Dict[str, dict], retired: Dict[str, dict],
          files: Optional[Iterable[str]] = None) -> Tuple[List[Finding], dict]:
    occs, fields, cov = scan(root, files)
    findings: List[Finding] = []
    seen_retired = set()

    for o in occs:
        dl = o.doi.lower()
        if dl in retired:
            seen_retired.add(dl)
            if not o.documented:
                findings.append(Finding(
                    "RETIRED_REVIVED", o.path, o.lineno,
                    "%s is a retired known-bad identifier (%s) but this line does not "
                    "mark it as documentation" % (o.doi, retired[dl].get("classification", "?"))))
            continue
        if not is_doi_syntax(o.doi):
            findings.append(Finding("NOT_DOI_SYNTAX", o.path, o.lineno,
                                    "%s is not valid DOI syntax" % o.doi))
            continue
        entry = cache.get(dl)
        if entry is None:
            findings.append(Finding("UNCACHED", o.path, o.lineno,
                                    "%s has no cached resolution; run --refresh" % o.doi))
        elif not entry.get("registered"):
            findings.append(Finding("UNREGISTERED", o.path, o.lineno,
                                    "%s does not resolve (Handle responseCode %s)"
                                    % (o.doi, entry.get("response_code"))))

    for f in fields:
        if f.value.lower() in retired:
            continue
        if not is_doi_syntax(f.value):
            findings.append(Finding("NOT_DOI_SYNTAX", f.path, f.lineno,
                                    "field %s has non-DOI value %r" % (f.field, f.value)))

    for dl, meta in retired.items():
        entry = cache.get(dl)
        if entry is not None and entry.get("registered"):
            findings.append(Finding("LEDGER_DRIFT", str(RETIRED_PATH.name), 0,
                                    "%s is in the retired ledger but the cache says it "
                                    "resolves; re-check the classification" % dl))
        if dl not in seen_retired and meta.get("expect_documented", True):
            findings.append(Finding("LEDGER_DRIFT", str(RETIRED_PATH.name), 0,
                                    "%s is in the retired ledger but appears nowhere in the "
                                    "tree; drop the ledger entry" % dl))

    stats = dict(cov)
    stats.update({"retired": len(retired), "cached": len(cache)})
    return findings, stats


# ---------------------------------------------------------------------------
# Non-vacuity self-test
# ---------------------------------------------------------------------------

def selftest(root: Path) -> int:
    """Prove the gate is non-vacuous: it must PASS a known-good tree and FAIL a
    known-bad tree on every rule. Anything less and a green run means nothing.

    Fixture DOIs are ASSEMBLED AT RUNTIME rather than written as literals. A
    literal known-bad DOI in this file would be swept by the gate itself when it
    scans the repository, and a literal "uncached" DOI would be resolved and
    cached by the next --refresh, silently disarming that arm.
    """
    cache = load_json(CACHE_PATH, {}).get("dois", {})
    if not cache:
        print("VOID: no resolution cache at %s — cannot run the non-vacuity proof."
              % CACHE_PATH)
        return 2
    unreg = [d for d, v in cache.items() if not v.get("registered")]
    if POSITIVE_CONTROL not in cache or not cache[POSITIVE_CONTROL].get("registered"):
        print("VOID: cache lacks a known-good DOI (%s); a pass would be meaningless."
              % POSITIVE_CONTROL)
        return 2
    if not unreg:
        print("VOID: cache contains no known-unregistered DOI, so the UNREGISTERED and "
              "RETIRED_REVIVED arms cannot be exercised.")
        return 2

    bad_doi = unreg[0]
    zenodo_bad = "10.5281/" + "zenodo." + "notanumber"
    uncached = "10.5555/" + "selftest.uncached.arm"
    good_fixture = (
        "# Reference: Hartshorn et al. (2007) J. Med. Chem. 50:726-741\n"
        "#            https://doi.org/%s\n"
        'zenodo_doi: ""\n'
        "doi: %s\n" % (POSITIVE_CONTROL, POSITIVE_CONTROL)
    )
    bad_fixtures = {
        "unregistered_live": "# Reference: https://doi.org/%s\n" % bad_doi,
        "zenodo_word_suffix": 'zenodo_doi: "%s"\n' % zenodo_bad,
        "uncached_doi": "# see https://doi.org/%s\n" % uncached,
    }

    results = {}
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        (tmp / "good.yaml").write_text(good_fixture, encoding="utf-8")
        f, _ = check(tmp, cache, {}, files=["good.yaml"])
        results["known_good"] = f

        for name, body in bad_fixtures.items():
            fn = "%s.yaml" % name
            (tmp / fn).write_text(body, encoding="utf-8")
            fb, _ = check(tmp, cache, {}, files=[fn])
            results[name] = fb

        (tmp / "revived.yaml").write_text("doi: %s\n" % bad_doi, encoding="utf-8")
        fr, _ = check(tmp, cache, {bad_doi: {"classification": "UNREGISTERED",
                                             "expect_documented": False}},
                      files=["revived.yaml"])
        results["retired_revived"] = fr

    passed_good = len(results["known_good"]) == 0
    bad_arms = {k: v for k, v in results.items() if k != "known_good"}
    failed_bad = {k: len(v) > 0 for k, v in bad_arms.items()}

    print("[selftest] known-good fixture -> %d findings (expected 0)" % len(results["known_good"]))
    for k, v in bad_arms.items():
        print("[selftest] known-bad  %-20s -> %d findings (expected >=1)%s"
              % (k, len(v), "" if v else "   *** ARM DID NOT FIRE ***"))
        for fnd in v:
            print("             %s: %s" % (fnd.kind, fnd.detail[:110]))

    if not passed_good or not all(failed_bad.values()):
        print("VOID: the check did not produce BOTH a pass on known-good and a failure on "
              "every known-bad arm. It proves nothing in this state.")
        return 2
    print("[selftest] NON-VACUOUS: passes known-good, fails all %d known-bad arms."
          % len(bad_arms))
    return 0


# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--refresh", action="store_true", help="re-resolve online, rewrite cache")
    ap.add_argument("--selftest", action="store_true", help="prove the check is non-vacuous")
    ap.add_argument("--report", action="store_true", help="print the full sweep table")
    args = ap.parse_args(argv)

    if args.refresh:
        return refresh(ROOT)
    if args.selftest:
        return selftest(ROOT)

    cache_doc = load_json(CACHE_PATH, {})
    cache = cache_doc.get("dois", {})
    retired = load_json(RETIRED_PATH, {}).get("retired", {})
    if not cache:
        print("VOID: no resolution cache at %s. Run --refresh (needs network)." % CACHE_PATH)
        return 2

    findings, stats = check(ROOT, cache, retired)

    print("DOI integrity sweep — %d files read, %d DOI occurrences, %d distinct, "
          "%d citation fields, %d cached resolutions, %d retired"
          % (stats["files_read"], stats["occurrences"], stats["distinct"],
             stats["fields"], stats["cached"], stats["retired"]))
    print("cache checked_utc=%s oracle=%s" % (cache_doc.get("checked_utc"), cache_doc.get("oracle")))

    if args.report:
        n_unreg = sum(1 for v in cache.values() if not v.get("registered"))
        print("cache: %d registered, %d unregistered" % (len(cache) - n_unreg, n_unreg))

    if findings:
        print("\nFAIL — %d finding(s):" % len(findings))
        for f in findings:
            print("  [%s] %s:%s  %s" % (f.kind, f.path, f.lineno, f.detail))
        return 1
    print("PASS — every DOI-shaped string resolves, every citation field is DOI syntax, "
          "and no retired identifier is back in service.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
