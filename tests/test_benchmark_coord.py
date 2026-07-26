"""Drive shipped scripts/benchmark_coord.py — hold, lock, disk, workers, stamp."""
from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
MOD_PATH = REPO / "scripts" / "benchmark_coord.py"


def _load():
    import sys

    name = "benchmark_coord_under_test"
    spec = importlib.util.spec_from_file_location(name, MOD_PATH)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    # dataclasses on 3.14 requires module registered before @dataclass runs
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture()
def coord(tmp_path, monkeypatch):
    mod = _load()
    root = tmp_path / "results"
    root.mkdir()
    monkeypatch.setenv("FLEXAIDDS_LOCAL_ROOT", str(root))
    # ensure disk floor does not trip on tmp (usually huge free)
    monkeypatch.delenv("FLEXAIDDS_DISK_FLOOR_OVERRIDE", raising=False)
    return mod, root


def test_hold_refuses_preflight(coord, tmp_path):
    mod, root = coord
    mod.write_hold(root, owner="ops", reason="test hold")
    binary = tmp_path / "bin" / "engine"
    binary.parent.mkdir(parents=True)
    binary.write_bytes(b"\x7fELFfake")
    out = tmp_path / "run1"
    r = mod.preflight(
        root=root,
        out=out,
        workers=2,
        binary=binary,
        owner="tester",
    )
    assert r.ok is False
    assert "HOLD" in r.reason or "hold" in r.reason.lower()
    assert not (out / "bin").exists()


def test_mkdir_lock_exclusive(coord, tmp_path):
    mod, root = coord
    ok1, msg1, tok1 = mod.acquire_lock(root, owner="a", purpose="dock")
    assert ok1 and tok1
    ok2, msg2, tok2 = mod.acquire_lock(root, owner="b", purpose="dock")
    assert ok2 is False
    assert tok2 is None
    assert "lock" in msg2.lower()
    ok_rel, _ = mod.release_lock(root, tok1)
    assert ok_rel
    ok3, _, tok3 = mod.acquire_lock(root, owner="c", purpose="dock")
    assert ok3 and tok3
    mod.release_lock(root, tok3)


def test_release_token_mismatch(coord):
    mod, root = coord
    ok, _, tok = mod.acquire_lock(root, owner="a")
    assert ok
    bad, msg = mod.release_lock(root, "wrong-token")
    assert bad is False
    assert "mismatch" in msg.lower()
    good, _ = mod.release_lock(root, tok)
    assert good


def test_workers_cap(coord):
    mod, root = coord
    ok, msg = mod.check_workers(4)
    assert ok
    ok2, msg2 = mod.check_workers(5)
    assert ok2 is False
    assert "MAX_WORKERS" in msg2 or "exceeds" in msg2


def test_disk_floor_refuse(coord, tmp_path, monkeypatch):
    mod, root = coord
    monkeypatch.setenv("FLEXAIDDS_DISK_FLOOR_GB", "999999")
    monkeypatch.delenv("FLEXAIDDS_DISK_FLOOR_OVERRIDE", raising=False)
    ok, msg, free = mod.check_disk(root)
    assert ok is False
    assert free < 999999


def test_stamp_binary(coord, tmp_path):
    mod, root = coord
    binary = tmp_path / "FlexAIDdS"
    binary.write_bytes(b"mach-o-fake-binary-content")
    out = tmp_path / "run_out"
    stamped = mod.stamp_binary(binary, out)
    assert stamped.is_file()
    assert stamped.read_bytes() == binary.read_bytes()
    assert (out / "bin" / "BINARY_STAMP.json").is_file()


def test_preflight_success_stamps_and_locks(coord, tmp_path):
    mod, root = coord
    binary = tmp_path / "engine"
    binary.write_bytes(b"\x00binary")
    out = tmp_path / "run_ok"
    r = mod.preflight(
        root=root,
        out=out,
        workers=2,
        binary=binary,
        owner="unit",
    )
    assert r.ok, r.reason
    assert r.lock_token
    assert r.stamped_binary and Path(r.stamped_binary).is_file()
    assert mod.lock_dir(root).is_dir()
    # second preflight fails while locked
    r2 = mod.preflight(
        root=root,
        out=tmp_path / "run2",
        workers=2,
        binary=binary,
        owner="other",
    )
    assert r2.ok is False
    mod.release_lock(root, r.lock_token)


def test_cli_status_and_hold(tmp_path, monkeypatch):
    monkeypatch.setenv("FLEXAIDDS_LOCAL_ROOT", str(tmp_path / "res"))
    (tmp_path / "res").mkdir()
    env = os.environ.copy()
    env["FLEXAIDDS_LOCAL_ROOT"] = str(tmp_path / "res")
    r = subprocess.run(
        [sys.executable, str(MOD_PATH), "status"],
        capture_output=True,
        text=True,
        env=env,
        timeout=30,
        check=False,
    )
    assert r.returncode == 0, r.stderr
    body = json.loads(r.stdout)
    assert "may_dock" in body
    r2 = subprocess.run(
        [
            sys.executable,
            str(MOD_PATH),
            "hold",
            "--owner",
            "cli",
            "--reason",
            "unit",
        ],
        capture_output=True,
        text=True,
        env=env,
        timeout=30,
        check=False,
    )
    assert r2.returncode == 0
    r3 = subprocess.run(
        [sys.executable, str(MOD_PATH), "status"],
        capture_output=True,
        text=True,
        env=env,
        timeout=30,
        check=False,
    )
    st = json.loads(r3.stdout)
    assert st["hold"] is not None
    assert st["may_dock"] is False
