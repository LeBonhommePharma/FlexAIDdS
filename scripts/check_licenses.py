#!/usr/bin/env python3
"""check_licenses.py — Verify no GPL/AGPL/LGPL dependencies in the codebase.

Parses scancode-toolkit JSON output and fails if any forbidden licenses
are detected. Used in CI to enforce the Apache-2.0 clean-room policy.

Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
SPDX-License-Identifier: Apache-2.0
"""

import json
import re
import sys
from typing import Iterator, List, Tuple

# Copyleft license families forbidden by the Apache-2.0 clean-room policy.
# Matched by family prefix so every version and suffix is covered — gpl-2.0,
# gpl-3.0-or-later, gpl-1.0-plus (scancode's key), lgpl-2.1, agpl-3.0-only, etc.
# (won't match apache/mit/mpl/bsd).
_COPYLEFT = re.compile(r"^(a|l)?gpl(-|$)")

# scancode assigns each match a score in [0, 100]. A real license header or an
# explicit SPDX-License-Identifier scores ~95-100; a bare-word mention — e.g.
# "No GPL dependencies", rule gpl_bare_word_only.RULE — scores ~50. Only
# high-confidence matches count, so a disclaimer that merely names GPL does not
# trip the gate.
_MIN_MATCH_SCORE = 90.0

# License keys are separated in expressions by whitespace, boolean operators,
# and parentheses (e.g. "gpl-2.0 WITH classpath-exception-2.0", "(mit OR gpl-3.0)").
_EXPRESSION_SPLIT = re.compile(r"[^a-z0-9.+-]+")
_EXPRESSION_OPERATORS = {"and", "or", "with"}


# Files that legitimately reference forbidden license names without being under
# them: the license matrix, package metadata, and documentation. Scanning these
# for a *dependency* gate produces false positives (THIRD_PARTY_LICENSES.md lists
# the forbidden licenses by name; setup.py carries PyPI classifiers), so they are
# excluded. Real vendored GPL code is caught via its own source files.
_EXCLUDED_PATH = re.compile(
    r"(^|/)(LICENSE|LICENCE|NOTICE|COPYING|COPYRIGHT|THIRD_PARTY_LICENSES)[^/]*$"
    r"|\.(md|rst|txt)$"
    r"|(^|/)(setup\.py|pyproject\.toml)$"
    r"|(^|/)(docs?|licensing|licenses)/",
    re.IGNORECASE,
)


def _is_excluded(path: str) -> bool:
    return bool(_EXCLUDED_PATH.search(path))


def _license_keys(expression: str) -> Iterator[str]:
    """Split a scancode/SPDX license expression into individual license keys."""
    for token in _EXPRESSION_SPLIT.split(expression.lower()):
        token = token.strip()
        if token and token not in _EXPRESSION_OPERATORS:
            yield token


def _forbidden_matches(file_entry: dict) -> Iterator[str]:
    """Yield forbidden copyleft license keys from *high-confidence* matches.

    Reads scancode's per-match scores — v32 ``license_detections[].matches`` or
    the legacy ``licenses[]`` array — and only yields a key whose own match scored
    >= _MIN_MATCH_SCORE. That is what separates a real GPL header from a bare-word
    mention in a comment such as "No GPL code": scancode records the latter as a
    score-50 gpl_bare_word_only match and folds it into detected_license_expression
    as e.g. 'apache-2.0 AND gpl-1.0-plus', which a plain expression check would
    (wrongly) flag. license_clues are ignored entirely.
    """
    detections = file_entry.get("license_detections")
    if detections is not None:
        for det in detections or []:
            for match in det.get("matches") or []:
                if (match.get("score") or 0) < _MIN_MATCH_SCORE:
                    continue
                for key in _license_keys(match.get("license_expression") or ""):
                    if _COPYLEFT.match(key):
                        yield key
    else:
        # Legacy schema (<= v31): flat licenses[] with a per-entry score.
        for lic in file_entry.get("licenses") or []:
            if (lic.get("score") or 0) < _MIN_MATCH_SCORE:
                continue
            for key in (lic.get("key"), lic.get("spdx_license_key")):
                if key and _COPYLEFT.match(key.lower()):
                    yield key.lower()


def check_report(report_path: str) -> Tuple[List[Tuple[str, str]], int]:
    """Check scancode report for forbidden licenses.

    Returns ``(violations, n_files)`` where ``violations`` is a list of
    ``(file_path, license_key)`` tuples and ``n_files`` is the number of file
    entries scancode reported (used to detect a scan that produced nothing).
    """
    with open(report_path) as f:
        report = json.load(f)

    files = report.get("files", [])
    violations = []

    for file_entry in files:
        file_path = file_entry.get("path", "unknown")
        if _is_excluded(file_path):
            continue
        seen = set()
        for key in _forbidden_matches(file_entry):
            if key not in seen:
                seen.add(key)
                violations.append((file_path, key))

    return violations, len(files)


def main():
    if len(sys.argv) < 2:
        print("Usage: check_licenses.py <scancode-report.json>")
        sys.exit(1)

    report_path = sys.argv[1]

    # A missing or malformed report means the scan did not produce usable output.
    # Treat that as a gate failure, not a pass: a license check that silently
    # skips is indistinguishable from one that found nothing.
    try:
        violations, n_files = check_report(report_path)
    except (json.JSONDecodeError, FileNotFoundError) as e:
        print(f"ERROR: Could not read license report '{report_path}': {e}")
        print("The scan did not produce a usable report — failing the gate.")
        sys.exit(1)

    if n_files == 0:
        print("ERROR: license report contains 0 scanned files.")
        print("scancode produced no file entries — the scan did not run correctly.")
        sys.exit(1)

    if violations:
        print(f"FORBIDDEN LICENSES DETECTED ({len(violations)} violations):")
        print("=" * 60)
        for file_path, license_key in violations:
            print(f"  {license_key:30s}  {file_path}")
        print("=" * 60)
        print("FlexAIDdS is Apache-2.0. GPL/AGPL/LGPL dependencies are forbidden.")
        print("See docs/licensing/clean-room-policy.md for details.")
        sys.exit(1)

    print(f"License scan PASSED: no forbidden licenses in {n_files} scanned files.")
    sys.exit(0)


if __name__ == "__main__":
    main()
