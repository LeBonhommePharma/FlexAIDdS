"""Tests for resolve_build.py stale-path rejection and SHA pinning."""
from __future__ import annotations

import json
import os
import stat
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
RESOLVE = REPO_ROOT / ".grok/skills/flexaidds/scripts/resolve_build.py"
sys.path.insert(0, str(RESOLVE.parent))

import resolve_build as rb  # noqa: E402


def _touch_exec(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"\x00")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


def test_rejects_stale_tmp_without_pin(monkeypatch: pytest.MonkeyPatch):
    build = Path("/tmp/flexaidds_resolve_test_stale")
    if build.exists():
        import shutil
        shutil.rmtree(build)
    build.mkdir()
    _touch_exec(build / rb.ENGINE_NAME)
    _touch_exec(build / rb.RUNNER_NAME)
    for name in rb.RUNTIME_DATA_FILES:
        (build / name).write_text("x")

    monkeypatch.delenv("FLEXAIDDS_BUILD", raising=False)
    monkeypatch.delenv("FLEXAIDDS_BUILD_DIR", raising=False)
    monkeypatch.delenv("FLEXAIDDS_ENGINE_SHA256", raising=False)
    home = Path("/tmp/flexaidds_resolve_test_home")
    home.mkdir(exist_ok=True)
    monkeypatch.setenv("HOME", str(home))

    with pytest.raises(SystemExit):
        rb.resolve_build(repo_root=build.parent, pin_sha=None)


def test_pin_requires_exact_sha(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    build = tmp_path / "build"
    build.mkdir()
    _touch_exec(build / rb.ENGINE_NAME)
    _touch_exec(build / rb.RUNNER_NAME)
    for name in rb.RUNTIME_DATA_FILES:
        (build / name).write_text("data")

    engine_sha = rb._sha256(build / rb.ENGINE_NAME)
    monkeypatch.setenv("FLEXAIDDS_BUILD", str(build))

    resolution = rb.resolve_build(repo_root=tmp_path, pin_sha=engine_sha)
    assert resolution.engine_sha256 == engine_sha
    assert resolution.pinned is True

    with pytest.raises(SystemExit):
        rb.resolve_build(repo_root=tmp_path, pin_sha="0" * 64)


@pytest.mark.skipif(
    not (Path.home() / "Projects/FlexAIDdS/build_lto/benchmark_datasets").is_file(),
    reason="production build not present on this machine",
)
def test_resolves_main_checkout_build():
    proc = subprocess.run(
        [sys.executable, str(RESOLVE), "--check", "--repo-root", str(REPO_ROOT)],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=REPO_ROOT,
    )
    assert proc.returncode == 0, proc.stderr or proc.stdout


def test_export_shell_is_safe(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    build = tmp_path / "build"
    build.mkdir()
    _touch_exec(build / rb.ENGINE_NAME)
    _touch_exec(build / rb.RUNNER_NAME)
    for name in rb.RUNTIME_DATA_FILES:
        (build / name).write_text("data")
    monkeypatch.setenv("FLEXAIDDS_BUILD", str(build))

    resolution = rb.resolve_build(repo_root=tmp_path)
    shell = rb.export_shell(resolution)
    assert "export FLEXAIDDS_ENGINE_SHA256=" in shell
    assert "';" not in shell