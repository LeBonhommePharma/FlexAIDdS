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


def test_install_only_paths_skip_dock():
    files = [
        "Formula/flexaidds.rb",
        "packaging/spack/package.py",
        "conda/meta.yaml",
        "containers/Dockerfile.locked",
        ".devcontainer/devcontainer.json",
        "python/flexaidds/updater.py",
        "python/README.md",
        "python/pyproject.toml",
        "python/setup.py",
        "scripts/install.sh",
        "scripts/update.sh",
        "site/FlexAIDdS/index.html",
    ]
    assert all(gate.is_post_election_path(p) for p in files)
    assert gate.files_need_dock(files) is False


def test_updater_plus_results_still_docks():
    assert gate.files_need_dock(
        ["python/flexaidds/updater.py", "python/flexaidds/results.py"]
    ) is True


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


def test_docs_skills_receipt_paths_skip_dock():
    files = [
        "METHODOLOGY.md",
        "AGENTS.md",
        ".grok/skills/flexaidds/SKILL.md",
        ".grok/skills/flexaidds-dataset-runner/SKILL.md",
        ".agents/skills/flexaidds-benchmarking/SKILL.md",
        "benchmarks/protocols/three_engine_entropy_comparison.md",
        "benchmarks/protocols/science_exclusions.md",
        "scripts/blind_astex85_receipt_protocol.py",
        "scripts/check_run_receipt.py",
        "scripts/check_run_receipt_codes.py",
        "python/flexaidds/__init__.py",
    ]
    assert all(gate.is_post_election_path(p) for p in files)
    assert gate.files_need_dock(files) is False


def test_pr497_matrix_pin_docs_skips_dock():
    """#497 (9dc9 pin) is skills + METHODOLOGY + protocols + receipt printer."""
    files = [
        ".grok/skills/flexaidds-dataset-runner/SKILL.md",
        ".grok/skills/flexaidds/SKILL.md",
        "METHODOLOGY.md",
        "benchmarks/protocols/three_engine_entropy_comparison.md",
        "docs/implementation/3dsig_red_pair_protocol.md",
        "docs/implementation/BLIND_ASTEX85_RECEIPT_PROTOCOL.md",
        "scripts/blind_astex85_receipt_protocol.py",
        "tests/test_blind_astex85_receipt_protocol.py",
    ]
    assert gate.files_need_dock(files) is False


def test_pr501_init_docstring_skips_dock():
    assert gate.files_need_dock(["python/flexaidds/__init__.py"]) is False


def test_pr499_softbeta_ghost_flags_skips_dock():
    files = [
        "docs/implementation/FORWARD_SUCCESS_RATE_PLAN.md",
        "docs/implementation/SCORING_PROVENANCE.md",
        "docs/implementation/softbeta_election_policy.md",
        "scripts/acf_strict_offline_reelect.py",
        "scripts/check_run_receipt.py",
        "scripts/e10_election_vs_scoring.py",
    ]
    assert gate.files_need_dock(files) is False


def test_engine_paths_still_dock():
    for path in (
        "LIB/gaboom.cpp",
        "LIB/DatasetRunner.cpp",
        "LIB/Vcontacts.cpp",
        "src/backends/webgpu/webgpu_eval.cpp",
        "python/flexaidds/results.py",
        "python/flexaidds/docking.py",
        "python/flexaidds/dataset_runner/runner.py",
        "benchmarks/datasets/astex_diverse.yaml",
        "CMakeLists.txt",
        "scripts/generate_flexaid_inp.py",
    ):
        assert gate.files_need_dock([path]) is True, path


def test_skills_plus_gaboom_still_docks():
    assert gate.files_need_dock(
        [".grok/skills/flexaidds/SKILL.md", "LIB/gaboom.cpp"]
    ) is True


def test_init_plus_results_still_docks():
    assert gate.files_need_dock(
        ["python/flexaidds/__init__.py", "python/flexaidds/results.py"]
    ) is True


def test_protocol_plus_dataset_yaml_still_docks():
    assert gate.files_need_dock(
        [
            "benchmarks/protocols/three_engine_entropy_comparison.md",
            "benchmarks/datasets/astex_diverse.yaml",
        ]
    ) is True
