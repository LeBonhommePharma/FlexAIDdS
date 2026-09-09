"""GitHub workflow YAML must parse the way GitHub parses it.

Duplicate mapping keys (especially duplicated job ids) make GitHub reject
the file at startup: 0s, no jobs, workflow name = the file path.
yaml.safe_load cannot catch that class — it keeps the last duplicate.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import check_github_workflows as gw  # noqa: E402


def test_repo_workflows_and_actions_parse_with_unique_keys():
    errors = gw.collect_errors(REPO_ROOT)
    assert errors == [], "\n".join(errors)


def test_packaging_job_ids_are_unique():
    text = (REPO_ROOT / ".github/workflows/packaging.yml").read_text(encoding="utf-8")
    data = gw.load_unique(text)
    jobs = data["jobs"]
    assert "setup-action-smoke" in jobs
    assert list(jobs).count("setup-action-smoke") == 1


def test_duplicate_job_id_is_rejected(tmp_path: Path):
    workflows = tmp_path / ".github" / "workflows"
    workflows.mkdir(parents=True)
    (workflows / "dup.yml").write_text(
        "name: Dup\n"
        "on: push\n"
        "jobs:\n"
        "  lint:\n"
        "    runs-on: ubuntu-latest\n"
        "    steps:\n"
        "      - run: true\n"
        "  lint:\n"
        "    runs-on: ubuntu-latest\n"
        "    steps:\n"
        "      - run: true\n",
        encoding="utf-8",
    )
    errors = gw.collect_errors(tmp_path)
    assert errors, "duplicate job id must fail the gate"
    assert any("duplicate YAML key" in e and "lint" in e for e in errors)


def test_safe_load_misses_the_duplicate_that_github_rejects():
    text = (
        "jobs:\n"
        "  a:\n"
        "    runs-on: ubuntu-latest\n"
        "  a:\n"
        "    runs-on: windows-latest\n"
    )
    silent = yaml.safe_load(text)
    assert list(silent["jobs"]) == ["a"]
    with pytest.raises(yaml.YAMLError, match="duplicate YAML key"):
        gw.load_unique(text)


def test_cli_ok_on_this_repo():
    assert gw.main(["--root", str(REPO_ROOT)]) == 0
