"""Regression tests: CF scoring proxy must not be labeled as free energy / ΔG.

These are lightweight source-audits (no C++ build required). They guard the
AGENTS.md rule that separates the CF/contact-function scoring proxy from true
thermodynamic free energy claims.

Live campaign CSVs may keep historical column names (`best_score`,
`predicted_dG`); comments and docs must describe them accurately.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

# Repo root: python/tests/ -> python/ -> repo
REPO_ROOT = Path(__file__).resolve().parents[2]


def _read(rel: str) -> str:
    path = REPO_ROOT / rel
    assert path.is_file(), f"missing required file: {rel}"
    return path.read_text(encoding="utf-8", errors="replace")


class TestDatasetRunnerHeaderNaming:
    """LIB/DatasetRunner.h field comments must not call best_score free energy."""

    def test_best_score_not_called_free_energy(self):
        text = _read("LIB/DatasetRunner.h")
        # Historical bug: `float best_score{0.0f}; // FlexAIDdS free energy`
        offenders = []
        for i, line in enumerate(text.splitlines(), start=1):
            if "best_score" not in line:
                continue
            lower = line.lower()
            # Same-line free-energy claim on the field declaration is banned
            if re.search(r"best_score\s*\{", line) or "float best_score" in line:
                if "free energy" in lower or "Δg" in lower or "delta g" in lower:
                    offenders.append(f"{i}: {line.strip()}")
            # Inline comment on the declaration line only (after //)
            if "float best_score" in line and "//" in line:
                comment = line.split("//", 1)[1].lower()
                if "free energy" in comment or "Δg" in comment:
                    if "not" not in comment and "≠" not in comment and "!=" not in comment:
                        offenders.append(f"{i}: {line.strip()}")
        assert not offenders, (
            "best_score must be documented as CF/contact-function scoring proxy, "
            f"not free energy:\n" + "\n".join(offenders)
        )

    def test_best_score_documents_cf_proxy(self):
        text = _read("LIB/DatasetRunner.h")
        assert "best_score" in text
        assert re.search(
            r"best_score[^\n]{0,200}CF|CF[^\n]{0,200}best_score|"
            r"contact-function scoring proxy[^\n]{0,400}best_score|"
            r"best_score[^\n]{0,400}contact-function scoring proxy|"
            r"best_score\s*≡\s*CF",
            text,
            re.IGNORECASE | re.DOTALL,
        ), "DockingResult docs must state that best_score is the CF scoring proxy"

    def test_predicted_dg_has_ledger_caveat(self):
        text = _read("LIB/DatasetRunner.h")
        assert "predicted_dG" in text
        banned = re.compile(
            r"predicted_dG[^\n]{0,80}\bis\b[^\n]{0,40}experimental|"
            r"experimental binding free energy[^\n]{0,40}predicted_dG|"
            r"float predicted_dG[^\n]*//[^\n]*experimental",
            re.IGNORECASE,
        )
        hits = []
        for i, line in enumerate(text.splitlines(), start=1):
            if banned.search(line):
                hits.append(f"{i}: {line.strip()}")
        assert not hits, "predicted_dG must not be labeled experimental ΔG:\n" + "\n".join(hits)
        assert re.search(
            r"predicted_dG.{0,500}(ensemble|estimate|NOT experimental|not experimental|fallback)",
            text,
            re.IGNORECASE | re.DOTALL,
        ), "predicted_dG must be documented as ensemble estimate / not experimental ΔG"


class TestDatasetRunnerCppAssignmentComments:
    """Assignment site comments in DatasetRunner.cpp must stay accurate."""

    def test_best_score_assignment_mentions_cf(self):
        text = _read("LIB/DatasetRunner.cpp")
        m = re.search(
            r"result\.best_score\s*=\s*best_cf\s*;[^\n]*",
            text,
        )
        assert m, "expected result.best_score = best_cf assignment"
        start = max(0, m.start() - 1200)
        window = text[start : m.end() + 200]
        assert re.search(r"CF|contact-function|scoring proxy", window, re.I), (
            "best_score assignment context must describe CF scoring proxy"
        )
        bad = re.findall(
            r".{0,40}best_score.{0,40}free energy.{0,40}",
            window,
            flags=re.I,
        )
        bad = [b for b in bad if "not" not in b.lower()]
        assert not bad, f"best_score must not be called free energy: {bad}"

    def test_predicted_dg_assignment_has_fallback_caveat(self):
        text = _read("LIB/DatasetRunner.cpp")
        m = re.search(r"result\.predicted_dG\s*=\s*have_free_energy", text)
        assert m, "expected predicted_dG assignment with free_energy / CF fallback"
        start = max(0, m.start() - 800)
        window = text[start : m.end() + 120]
        assert re.search(
            r"NOT experimental|not experimental|fallback|ensemble",
            window,
            re.I,
        ), "predicted_dG assignment must document ensemble estimate / CF fallback"


class TestDocsAndGuidance:
    """Docs and skill guidance must encode the CF vs ΔG contract."""

    def test_guidance_has_csv_column_contract(self):
        text = _read(".grok/skills/flexaidds/references/flexaidds-guidance.md")
        assert "best_score" in text
        assert "predicted_dG" in text
        assert re.search(r"CF/contact-function scoring proxy", text)
        assert "experimental" in text.lower() or "NOT experimental" in text

    def test_thermodynamics_md_csv_contract(self):
        text = _read("docs/thermodynamics.md")
        assert "best_score" in text
        assert re.search(r"contact-function scoring proxy", text, re.I)
        assert "predicted_dG" in text

    def test_known_limitations_mentions_best_score(self):
        text = _read("docs/KNOWN_LIMITATIONS.md")
        assert "best_score" in text
        assert "predicted_dG" in text
        assert "not" in text.lower()


class TestPythonModelDocstrings:
    """Python pose CF fields must not be described as free energy."""

    def test_pose_result_cf_docstring(self):
        text = _read("python/flexaidds/models.py")
        m = re.search(
            r"class PoseResult:.*?cf:.*?(?:\n\s{8}\w|\n    [a-z_]+:)",
            text,
            re.DOTALL,
        )
        assert m, "PoseResult.cf docstring not found"
        block = m.group(0)
        assert re.search(r"CF|contact-function|scoring proxy", block, re.I)
        for line in block.splitlines():
            if re.search(r"\bcf\b", line, re.I) and "free energy" in line.lower():
                if "not" not in line.lower():
                    pytest.fail(f"cf field must not be free energy: {line.strip()}")
