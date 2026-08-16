#!/usr/bin/env python3
"""Fail-closed scanner: primary docs must not re-assert withdrawn Astex rates.

Extends the claim-firewall policy in ``tests/test_thermo_claim_firewall.py``
(provenance for entropy.help) and the JCIM label pin in
``tests/test_comparative_benchmark_methodology.py`` (top-1 45.2% / top-10 66.7%).

This module watches README.md, docs/BENCHMARK.md, REPRODUCIBILITY.md, and
scripts/reproduce_astex85.sh. Hits are allowed only with an adjacent
withdrawal / oracle / not-docking-power qualifier.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

SURFACES = (
    ROOT / "README.md",
    ROOT / "docs" / "BENCHMARK.md",
    ROOT / "REPRODUCIBILITY.md",
    ROOT / "scripts" / "reproduce_astex85.sh",
)

# Tokens that must not appear as current FlexAID(∆S) docking power.
CLAIM_TOKEN_RE = re.compile(
    r"(?<!\d)(?:91\.8|94\.1|88\.2)(?!\d)"
    r"|78/85|80/85"
    r"|r\s*=\s*0\.93"
    r"|SEED_ELITISM=1(?!\d)"
    r"|NATIVE_SEED_FRAC=0\.90"
    r"|run_dataset\.py"
    r"|analyze_affinity\.py",
    re.I,
)

# Same-line or ±N-line window must mark the hit as withdrawn / historical / negative.
QUALIFIER_RE = re.compile(
    r"(?:"
    r"withdrawn|withdrawal|oracle|ceiling|seed-echo|seed echo|"
    r"unverified|misquote|disqualifying|former|formerly|"
    r"previously|historical|"
    r"not docking power|not a published|not the default|not current|"
    r"do not cite|do not export|do not invoke|do not document|do not treat|"
    r"does not exist|do not exist|"
    r"no receipt|pending receipt|forbids reporting|"
    r"live oracle lever|dead knob|not_docking_power"
    r")",
    re.I,
)

WINDOW_RADIUS = 3


def _window(lines: list[str], index: int, radius: int = WINDOW_RADIUS) -> str:
    lo = max(0, index - radius)
    hi = min(len(lines), index + radius + 1)
    return "\n".join(lines[lo:hi])


def unqualified_claim_hits(text: str, *, source: str = "") -> list[str]:
    """Return human-readable hits that lack a withdrawal/oracle qualifier."""
    lines = text.splitlines()
    hits: list[str] = []
    for i, line in enumerate(lines):
        for match in CLAIM_TOKEN_RE.finditer(line):
            ctx = _window(lines, i)
            if QUALIFIER_RE.search(line) or QUALIFIER_RE.search(ctx):
                continue
            loc = f"{source}:{i + 1}" if source else f"line {i + 1}"
            hits.append(f"{loc}: unqualified {match.group(0)!r}: {line.strip()[:160]}")
    return hits


def test_primary_surfaces_do_not_reassert_withdrawn_astex_rates() -> None:
    failures: list[str] = []
    for path in SURFACES:
        assert path.is_file(), f"missing claim surface {path}"
        failures.extend(
            unqualified_claim_hits(path.read_text(encoding="utf-8"), source=path.name)
        )
    assert failures == [], "unqualified Astex docking-power claims:\n" + "\n".join(
        failures
    )


def test_scanner_rejects_unqualified_rate() -> None:
    text = "FlexAID∆S achieves 91.8% S1 on Astex-85.\n"
    hits = unqualified_claim_hits(text, source="fixture")
    assert hits, "scanner must fail an unqualified 91.8% docking-power sentence"


def test_scanner_allows_withdrawn_rate() -> None:
    text = "The figures previously stated here (78/85 = 91.8%) are withdrawn.\n"
    assert unqualified_claim_hits(text, source="fixture") == []


def test_benchmark_md_pins_jcim_table2_not_88_2_as_2015_s1() -> None:
    text = (ROOT / "docs" / "BENCHMARK.md").read_text(encoding="utf-8")
    assert "45.2%" in text and "66.7%" in text
    assert re.search(
        r"top-1[^\n]{0,80}45\.2%|45\.2%[^\n]{0,80}[Tt]op-1",
        text,
    ), "45.2% must be labelled JCIM top-1"
    assert re.search(
        r"top-10[^\n]{0,80}66\.7%|66\.7%[^\n]{0,80}[Tt]op-10",
        text,
    ), "66.7% must be labelled JCIM top-10"
    assert not re.search(
        r"FlexAID \(original\)\s*\|\s*88\.2%",
        text,
    ), "do not list 88.2% as FlexAID 2015 S1 in the comparison table"


def test_reproduce_script_default_is_blind() -> None:
    text = (ROOT / "scripts" / "reproduce_astex85.sh").read_text(encoding="utf-8")
    assert re.search(r"^SEED_ELITISM=0\s*$", text, re.M)
    assert re.search(r"^NATIVE_SEED_FRAC=0\s*$", text, re.M)
    assert re.search(r"^\s+SEED_ELITISM=1\s*$", text, re.M)
    assert not re.search(r"^SEED_ELITISM=1\s*$", text, re.M)
    assert re.search(
        r'^export FLEXAIDDS_SEED_ELITISM="\$\{SEED_ELITISM\}"',
        text,
        re.M,
    )
    assert re.search(
        r'^export FLEXAIDDS_NATIVE_SEED_FRAC="\$\{NATIVE_SEED_FRAC\}"',
        text,
        re.M,
    )
    assert not re.search(r"^export FLEXAIDDS_SEED_ELITISM=1\b", text, re.M)
    assert not re.search(r"^export FLEXAIDDS_NATIVE_SEED_FRAC=0\.90\b", text, re.M)
    assert "--oracle-ceiling" in text
    assert "ORACLE CEILING — not docking power" in text
    assert "N_DENOM = 85" in text
    assert "{'80':>12}" not in text


def test_named_runners_exist_and_phantoms_do_not() -> None:
    assert (ROOT / ".grok" / "skills" / "flexaidds" / "scripts" / "dataset_runner.py").is_file()
    cmake = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
    assert "add_executable(benchmark_datasets" in cmake
    assert not (ROOT / "scripts" / "run_dataset.py").exists()
    assert not (ROOT / "scripts" / "analyze_affinity.py").exists()


def test_methodology_success_operator_is_inclusive_le_2() -> None:
    text = (ROOT / "METHODOLOGY.md").read_text(encoding="utf-8")
    assert "Success ⇔ rank-0 in-place RMSD" in text
    assert "<= 2.0 Å" in text


def test_readme_does_not_publish_unreceipted_session_rates() -> None:
    text = (ROOT / "README.md").read_text(encoding="utf-8")
    for token in ("24.1%", "48.8%", "25.3%", "91.8%", "94.1%", "88.2%"):
        assert token not in text, f"README.md must not publish {token}"
