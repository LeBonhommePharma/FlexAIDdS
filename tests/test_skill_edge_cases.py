"""Edge-case tests for flexaidds skill agent autonomy."""
from __future__ import annotations

import json
import os
import stat
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SKILL_ROOT = REPO_ROOT / ".grok/skills/flexaidds"
sys.path.insert(0, str(SKILL_ROOT))
sys.path.insert(0, str(REPO_ROOT / "python"))
sys.path.insert(0, str(SKILL_ROOT / "scripts"))

import resolve_build as rb  # noqa: E402
from flexaidds_skill import (  # noqa: E402
    figure_parameters_from_mapping,
    generate_flexaids_figure,
    invoke_manifest_action,
)


def _touch_exec(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"\xcf\xfa\xed\xfe" + b"\x00" * 4)  # mach-o-ish header stub
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


def test_stale_env_path_rejected_when_binaries_missing(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    stale = tmp_path / "empty_build"
    stale.mkdir()
    good = tmp_path / "good_build"
    good.mkdir()
    _touch_exec(good / rb.ENGINE_NAME)
    _touch_exec(good / rb.RUNNER_NAME)
    for name in rb.RUNTIME_DATA_FILES:
        (good / name).write_text("ok")

    monkeypatch.setenv("FLEXAIDDS_BUILD", str(stale))
    monkeypatch.setenv("FLEXAIDDS_BUILD_DIR", str(good))
    monkeypatch.delenv("FLEXAIDDS_ENGINE_SHA256", raising=False)
    monkeypatch.setenv("HOME", str(tmp_path / "home"))

    resolution = rb.resolve_build(repo_root=tmp_path)
    assert Path(resolution.build_dir) == good.resolve()
    assert any("empty_build" in path for path, _ in resolution.rejected)


def test_malformed_active_manifest_ignored(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    home = tmp_path / "home"
    manifest = home / ".flexaidds"
    manifest.mkdir(parents=True)
    (manifest / "active_build.json").write_text("{not json", encoding="utf-8")

    build = tmp_path / "build"
    build.mkdir()
    _touch_exec(build / rb.ENGINE_NAME)
    _touch_exec(build / rb.RUNNER_NAME)
    for name in rb.RUNTIME_DATA_FILES:
        (build / name).write_text("x")

    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("FLEXAIDDS_BUILD", str(build))
    monkeypatch.delenv("FLEXAIDDS_ENGINE_SHA256", raising=False)

    resolution = rb.resolve_build(repo_root=tmp_path)
    assert resolution.engine_sha256


def test_require_fresh_rejects_old_binaries(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    build = tmp_path / "build"
    build.mkdir()
    _touch_exec(build / rb.ENGINE_NAME)
    _touch_exec(build / rb.RUNNER_NAME)
    for name in rb.RUNTIME_DATA_FILES:
        (build / name).write_text("x")

    source = tmp_path / "LIB/benchmark_datasets.cpp"
    source.parent.mkdir(parents=True)
    source.write_text("// newer source")
    os.utime(source, (2_000_000_000, 2_000_000_000))

    monkeypatch.setenv("FLEXAIDDS_BUILD", str(build))
    monkeypatch.delenv("FLEXAIDDS_ENGINE_SHA256", raising=False)

    with pytest.raises(SystemExit):
        rb.resolve_build(repo_root=tmp_path, require_fresh=True)


def test_figure_unknown_manifest_action():
    with pytest.raises(ValueError, match="Unknown manifest action"):
        invoke_manifest_action("not_a_real_action", {})


def test_figure_both_tds_and_enthalpy_alias():
    params = figure_parameters_from_mapping(
        {"entropy_value": 0.5, "enthalpy_value": 9.9, "tds_value": 1.1, "index_value": 0.5}
    )
    assert params.tds_value == 1.1


def test_figure_empty_output_dir_rejected():
    def fake_gen(_prompt: str):
        tmp = Path("/tmp/flexaidds_figure_edge.png")
        tmp.write_bytes(b"png")
        return {"path": str(tmp)}

    with pytest.raises(ValueError, match="output_dir"):
        generate_flexaids_figure(
            params={"index_value": 0.5},
            image_generator=fake_gen,
            output_dir="   ",
        )
    # Generator must not run when output_dir is invalid
    calls = {"n": 0}

    def counting_gen(_prompt: str):
        calls["n"] += 1
        return {"path": "/tmp/should-not-run.png"}

    with pytest.raises(ValueError, match="output_dir"):
        generate_flexaids_figure(image_generator=counting_gen, output_dir="")
    assert calls["n"] == 0


def test_dataset_runner_dry_run_without_build(monkeypatch: pytest.MonkeyPatch):
    """Dry-run must not require a C++ build tree."""
    monkeypatch.setenv("FLEXAIDDS_BUILD", "/nonexistent/stale/path")
    proc = subprocess.run(
        [
            sys.executable,
            str(SKILL_ROOT / "scripts/dataset_runner.py"),
            "--dataset",
            "astex_diverse",
            "--tier",
            "1",
            "--dry-run",
            "--no-ensure-data",
            "--results-dir",
            "results/edge_case_dry_run",
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        timeout=120,
        env={**os.environ, "PYTHONPATH": str(REPO_ROOT / "python")},
    )
    assert proc.returncode == 0, proc.stderr[-2000:] if proc.stderr else proc.stdout[-2000:]


def test_launch_script_usage_not_fatal_env():
    proc = subprocess.run(
        ["bash", str(SKILL_ROOT / "scripts/launch_full_benchmark.sh")],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert proc.returncode == 1
    assert "Usage:" in proc.stdout
    assert "FATAL: Could not source" not in proc.stdout


def test_bin_wrappers_not_symlinks():
    for name in (
        "validate-skill",
        "resolve-build",
        "ensure-docking-data",
        "dataset-runner",
    ):
        path = SKILL_ROOT / "bin" / name
        assert path.is_file(), name
        assert not path.is_symlink(), f"{name} must not be a symlink"
        assert os.access(path, os.X_OK), f"{name} must be executable"


def test_resolve_cli_json_roundtrip():
    proc = subprocess.run(
        [sys.executable, str(SKILL_ROOT / "scripts/resolve_build.py"), "--json", "--repo-root", str(REPO_ROOT)],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=REPO_ROOT,
    )
    if proc.returncode != 0:
        pytest.skip("no production build on this runner")
    data = json.loads(proc.stdout)
    assert len(data["engine_sha256"]) == 64
    assert data["runner_path"].endswith("benchmark_datasets")