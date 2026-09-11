#!/usr/bin/env python3
"""Fail-closed gate on citation identifiers — pytest face of scripts/check_dois.py.

The 2026-09-10 dataset identity audit found three fabricated Zenodo DOIs and four
unregistered DOIs that had survived indefinitely because no gate in this repository
ever resolved a citation identifier. This module is that gate, wired into the suite
that already runs the claim-firewall checks (tests/test_thermo_claim_firewall.py,
tests/test_check_published_astex_rates.py).

Runs OFFLINE: resolutions come from a committed fixture under
tests/fixtures/doi_resolution/, refreshed explicitly with

    python3 scripts/check_dois.py --refresh

The denominator assertions below exist because the field rule was briefly
LIVE-VACUOUS during development — it matched nothing in the real tree while still
firing on its own fixture, so the suite was green and checking nothing. A gate that
cannot say how many things it inspected is not a gate.

Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
SPDX-License-Identifier: Apache-2.0
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "check_dois.py"
FIXTURE_DIR = ROOT / "tests" / "fixtures" / "doi_resolution"
CACHE = FIXTURE_DIR / "doi_resolution_cache.json"
RETIRED = FIXTURE_DIR / "known_bad_dois.json"

sys.path.insert(0, str(ROOT / "scripts"))
import check_dois  # noqa: E402


def _cache() -> dict:
    return json.loads(CACHE.read_text(encoding="utf-8"))


def _retired() -> dict:
    return json.loads(RETIRED.read_text(encoding="utf-8"))["retired"]


def test_script_and_fixtures_exist() -> None:
    assert SCRIPT.is_file(), "scripts/check_dois.py is missing"
    assert CACHE.is_file(), "resolution cache fixture is missing; run --refresh"
    assert RETIRED.is_file(), "retired-DOI ledger fixture is missing"


def test_cache_controls_recorded() -> None:
    """A cache written without both controls passing is not evidence of anything."""
    doc = _cache()
    ctrl = doc.get("controls", {})
    assert ctrl.get("positive_registered") is True, "cache lacks a passing positive control"
    assert ctrl.get("negative_registered") is False, "cache lacks a failing negative control"


def test_sweep_denominators_are_nonzero() -> None:
    """Guard against a rule that inspects nothing and reports success."""
    occs, fields, stats = check_dois.scan(ROOT)
    assert stats["files_read"] > 500, "scanned implausibly few files: %r" % (stats,)
    assert len(occs) > 0, "no DOI-shaped strings found at all — the scanner is broken"
    assert len(fields) > 0, (
        "zero citation fields inspected: the NOT_DOI_SYNTAX rule is vacuous. "
        "This is the exact failure mode the rule was written to catch."
    )


def test_no_unresolvable_or_fabricated_dois() -> None:
    """The gate itself: every live citation identifier must resolve."""
    findings, stats = check_dois.check(ROOT, _cache()["dois"], _retired())
    assert not findings, "DOI integrity findings:\n" + "\n".join(
        "  [%s] %s:%s  %s" % (f.kind, f.path, f.lineno, f.detail) for f in findings
    )


def test_gate_is_non_vacuous() -> None:
    """Prove the gate produces BOTH a pass and a failure. Otherwise it proves nothing."""
    rc = check_dois.selftest(ROOT)
    assert rc == 0, "non-vacuity proof failed or printed VOID (rc=%d)" % rc


def test_retired_ledger_entries_are_actually_unregistered() -> None:
    cache = _cache()["dois"]
    for doi, meta in _retired().items():
        entry = cache.get(doi)
        if entry is None:
            continue
        assert not entry.get("registered"), (
            "%s sits in the retired ledger as %s but the cache says it resolves"
            % (doi, meta.get("classification"))
        )


# Values are ASSEMBLED from fragments rather than written as DOI literals. A literal
# here would be swept by the gate when it scans the repository, so this test file
# would silently change the very input surface it is meant to be checking — and an
# unregistered literal would fail the gate outright.
_SYNTAX_CASES = [
    ("10." + "1021/jm061277y", True),           # registered CrossRef DOI
    ("10." + "5281/zenodo.16861024", True),     # Zenodo DOI, numeric record id
    ("10." + "2210/pdb1g9v/pdb", True),         # wwPDB DataCite DOI, multi-segment
    ("10." + "5281/zenodo.dude37", False),      # word suffix: not Zenodo DOI syntax
    ("10." + "5281/zenodo.", False),            # empty record id
    ("not-a-doi", False),
    ("", False),
]


@pytest.mark.parametrize("value,ok", _SYNTAX_CASES)
def test_doi_syntax_rule(value: str, ok: bool) -> None:
    assert check_dois.is_doi_syntax(value) is ok
