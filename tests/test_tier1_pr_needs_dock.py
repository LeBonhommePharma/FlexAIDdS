"""Tier-1 skip gate: PoseBust/docs cannot move docking_power_top1."""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import tier1_pr_needs_dock as gate  # noqa: E402


def test_posebust_audit_readme_skips_dock():
    files = [
        "LIB/PoseBust/Engine.cpp",
        "LIB/PoseBust/BustCli.cpp",
        "docs/audit/2026-08-18_posebust_science_and_code_audit.md",
        "README.md",
        "tests/test_posebust.cpp",
        "tests/test_posebust_upstream_parity.py",
    ]
    assert all(gate.is_post_election_path(p) for p in files)
    assert gate.files_need_dock(files) is False


def test_dataset_runner_still_docks():
    files = [
        "LIB/PoseBust/Engine.cpp",
        "LIB/DatasetRunner.cpp",
        "docs/audit/note.md",
    ]
    assert gate.files_need_dock(files) is True


def test_gaboom_still_docks():
    assert gate.files_need_dock(["LIB/gaboom.cpp"]) is True


def test_python_package_still_docks():
    assert gate.files_need_dock(["python/flexaidds/results.py"]) is True


def test_benchmark_yaml_still_docks():
    assert gate.files_need_dock(["benchmarks/datasets/astex_diverse.yaml"]) is True


def test_empty_diff_skips():
    assert gate.files_need_dock([]) is False


def test_skip_gate_ci_files_do_not_dock():
    files = [
        "scripts/tier1_pr_needs_dock.py",
        "tests/test_tier1_pr_needs_dock.py",
        "tests/test_classify_diff.py",
        ".github/workflows/benchmark-tier1.yml",
        ".github/workflows/claude-code-review.yml",
        ".github/workflows/ci.yml",
    ]
    assert gate.files_need_dock(files) is False


def test_cli_files_flag(capsys, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    github_output = tmp_path / "github_output"
    monkeypatch.setenv("GITHUB_OUTPUT", str(github_output))
    rc = gate.main(["--files", "LIB/PoseBust/Engine.cpp", "README.md"])
    assert rc == 0
    captured = capsys.readouterr()
    assert "dock=false" in captured.out
    assert github_output.read_text(encoding="utf-8").strip() == "dock=false"


def test_cli_fail_closed_on_bad_git(capsys, monkeypatch):
    monkeypatch.setattr(
        gate,
        "_name_only",
        lambda base, head: ([], "cannot compute merge-base"),
    )
    rc = gate.main(["origin/main", "HEAD"])
    assert rc == 0
    assert "dock=true" in capsys.readouterr().out
