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


def test_cli_preflight_hold_then_lock_stamp_then_second_abort(tmp_path, monkeypatch):
    """Verification plan step 2: clean-temp hold→abort, lock+stamp, second abort."""
    root = tmp_path / "res"
    root.mkdir()
    binary = tmp_path / "engine"
    binary.write_bytes(b"\x7fELF-fake-engine")
    binary.chmod(0o755)
    env = os.environ.copy()
    env["FLEXAIDDS_LOCAL_ROOT"] = str(root)
    env.pop("FLEXAIDDS_DISK_FLOOR_OVERRIDE", None)

    def run_coord(*args):
        return subprocess.run(
            [sys.executable, str(MOD_PATH), *args],
            capture_output=True,
            text=True,
            env=env,
            timeout=30,
            check=False,
        )

    # hold present → preflight abort (exit 2), no stamp dir
    rh = run_coord("hold", "--owner", "ops", "--reason", "cli-demo-hold")
    assert rh.returncode == 0, rh.stderr
    out1 = tmp_path / "run_hold"
    r_hold = run_coord(
        "preflight",
        "--out",
        str(out1),
        "--workers",
        "2",
        "--binary",
        str(binary),
        "--owner",
        "cli-demo",
    )
    assert r_hold.returncode != 0
    body_hold = json.loads(r_hold.stdout)
    assert body_hold["ok"] is False
    assert "hold" in body_hold["reason"].lower() or "HOLD" in body_hold["reason"]
    assert not (out1 / "bin").exists()

    # clear hold → first preflight acquires lock + stamps
    ru = run_coord("unhold")
    assert ru.returncode == 0
    out2 = tmp_path / "run_ok"
    r_ok = run_coord(
        "preflight",
        "--out",
        str(out2),
        "--workers",
        "2",
        "--binary",
        str(binary),
        "--owner",
        "cli-demo",
    )
    assert r_ok.returncode == 0, r_ok.stdout + r_ok.stderr
    body_ok = json.loads(r_ok.stdout)
    assert body_ok["ok"] is True
    assert body_ok["lock_token"]
    assert body_ok["stamped_binary"]
    assert Path(body_ok["stamped_binary"]).is_file()
    assert (root / "BENCHMARK_DOCK_LOCK").is_dir()

    # second preflight while lock held → abort
    out3 = tmp_path / "run_second"
    r2 = run_coord(
        "preflight",
        "--out",
        str(out3),
        "--workers",
        "2",
        "--binary",
        str(binary),
        "--owner",
        "other-session",
    )
    assert r2.returncode != 0
    body2 = json.loads(r2.stdout)
    assert body2["ok"] is False
    assert "lock" in body2["reason"].lower()

    # release
    rr = run_coord("release", "--token", body_ok["lock_token"])
    assert rr.returncode == 0


def _load_lib_launch():
    name = "lib_launch_under_test"
    path = REPO / "scripts" / "lib_launch.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    # Ensure scripts/ is on path for benchmark_coord hard import
    scripts = str(REPO / "scripts")
    if scripts not in sys.path:
        sys.path.insert(0, scripts)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_lib_launch_import_requires_coord():
    """Sol #9: lib_launch must hard-import benchmark_coord (no soft None)."""
    mod = _load_lib_launch()
    assert mod.dock_preflight is not None
    assert callable(mod.dock_preflight)
    assert mod.MAX_WORKERS == 4


def test_launch_session_isolated_refuses_on_hold(tmp_path, monkeypatch):
    """Shipped launch path: hold present → RuntimeError before any dock child."""
    root = tmp_path / "results"
    root.mkdir()
    monkeypatch.setenv("FLEXAIDDS_LOCAL_ROOT", str(root))
    mod = _load_lib_launch()
    mod.write_hold = None  # not on lib_launch
    # write hold via coord
    coord = _load()
    coord.write_hold(root, owner="ops", reason="block launch_session")
    binary = tmp_path / "engine"
    binary.write_bytes(b"fake")
    binary.chmod(0o755)
    out = tmp_path / "out_launch"
    env = {"FLEXAIDDS_BINARY": str(binary), "PATH": os.environ.get("PATH", "")}
    with pytest.raises(RuntimeError, match="[Ss]ol #9|preflight|HOLD|hold"):
        mod.launch_session_isolated(
            cmd=[sys.executable, "-c", "print('should-not-run')"],
            env=env,
            output_dir=str(out),
            workers=2,
            binary=str(binary),
            owner="unit-test",
        )
    # no dock child pid file from a successful launch
    assert not (out / "benchmark.pid").exists()


def test_launch_session_isolated_preflight_success_path(tmp_path, monkeypatch):
    """Preflight passes → token + stamp written; use skip after preflight by
    checking DOCK_PREFLIGHT without completing the double-fork daemon.

    We only assert the preflight side-effects by calling dock_preflight through
    the same helpers launch_session_isolated uses, then verify a second call
    to launch_session_isolated fails on the held lock (no real dock needed).
    """
    root = tmp_path / "results"
    root.mkdir()
    monkeypatch.setenv("FLEXAIDDS_LOCAL_ROOT", str(root))
    mod = _load_lib_launch()
    binary = tmp_path / "engine"
    binary.write_bytes(b"mach-o-fake")
    binary.chmod(0o755)
    out = tmp_path / "out_ok"
    # First: direct preflight via exported dock_preflight (same as launch path)
    r = mod.dock_preflight(
        out=out,
        workers=2,
        binary=binary,
        owner="unit",
    )
    assert r.ok, r.reason
    assert (out / "bin").exists() or Path(r.stamped_binary).is_file()
    # Second launch_session_isolated must refuse while lock held (no fork success)
    env = {"FLEXAIDDS_BINARY": str(binary), "PATH": os.environ.get("PATH", "")}
    with pytest.raises(RuntimeError, match="[Ss]ol #9|preflight|lock"):
        mod.launch_session_isolated(
            cmd=[sys.executable, "-c", "print('no')"],
            env=env,
            output_dir=str(tmp_path / "out2"),
            workers=2,
            binary=str(binary),
            owner="other",
        )
    mod.dock_release_lock(root, r.lock_token)
