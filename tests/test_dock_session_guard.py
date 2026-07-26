#!/usr/bin/env python3
"""Unit tests for scripts/dock_session_guard.py (Sol #9 multi-session preflight).

Drives the shipped module paths with temporary hold/lock/binary dirs.
"""
from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
GUARD_PATH = REPO / "scripts" / "dock_session_guard.py"


def _load_guard():
    name = "dock_session_guard"
    spec = importlib.util.spec_from_file_location(name, GUARD_PATH)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    # Python 3.14 dataclasses require the module to be registered before exec.
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


class TestDockSessionGuard(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.g = _load_guard()

    def test_hold_blocks_when_file_present(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hold = root / "BENCHMARK_HOLD.json"
            hold.write_text('{"active": true, "reason": "unit-test"}\n', encoding="utf-8")
            blocked, msg = self.g.hold_blocks_launch(root)
            self.assertTrue(blocked)
            self.assertIn(str(hold), msg)
            self.assertIn("BENCHMARK_HOLD", msg)

    def test_hold_absent_allows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            blocked, msg = self.g.hold_blocks_launch(Path(tmp))
            self.assertFalse(blocked)
            self.assertEqual(msg, "")

    def test_hold_env_override_path(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            hold = Path(tmp) / "custom_hold.json"
            hold.write_text('{"active": true}\n', encoding="utf-8")
            old = os.environ.get("FLEXAIDDS_BENCHMARK_HOLD_PATH")
            os.environ["FLEXAIDDS_BENCHMARK_HOLD_PATH"] = str(hold)
            try:
                blocked, msg = self.g.hold_blocks_launch(Path(tmp) / "empty_repo")
                self.assertTrue(blocked)
                self.assertIn(str(hold), msg)
            finally:
                if old is None:
                    os.environ.pop("FLEXAIDDS_BENCHMARK_HOLD_PATH", None)
                else:
                    os.environ["FLEXAIDDS_BENCHMARK_HOLD_PATH"] = old

    def test_mkdir_lock_second_acquire_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            lock = Path(tmp) / "dock_session.lock"
            ok1, path1, msg1 = self.g.try_acquire_lock(lock, owner="first", pid=111)
            self.assertTrue(ok1, msg1)
            self.assertEqual(path1, lock)
            self.assertTrue((lock / "owner.json").is_file())
            ok2, path2, msg2 = self.g.try_acquire_lock(lock, owner="second", pid=222)
            self.assertFalse(ok2)
            self.assertEqual(path2, lock)
            self.assertIn("already held", msg2.lower())
            # first owner can release
            ok_r, msg_r = self.g.release_lock(lock)
            # release may refuse if pid != current; force if needed
            if not ok_r:
                ok_r, msg_r = self.g.release_lock(lock, force=True)
            self.assertTrue(ok_r, msg_r)
            ok3, _, msg3 = self.g.try_acquire_lock(lock, owner="third", pid=os.getpid())
            self.assertTrue(ok3, msg3)
            self.g.release_lock(lock, force=True)

    def test_disk_ok_injectable_free_space(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            ok, msg, avail = self.g.disk_ok(path, min_free_gb=20, free_gib=25.0)
            self.assertTrue(ok, msg)
            self.assertEqual(avail, 25.0)
            ok2, msg2, avail2 = self.g.disk_ok(path, min_free_gb=20, free_gib=12.0)
            self.assertFalse(ok2)
            self.assertIn("Disk floor fail", msg2)
            self.assertEqual(avail2, 12.0)

    def test_workers_cap(self) -> None:
        ok, msg = self.g.workers_ok(4, max_workers=4)
        self.assertTrue(ok, msg)
        bad, bmsg = self.g.workers_ok(6, max_workers=4)
        self.assertFalse(bad)
        self.assertIn("Workers cap fail", bmsg)

    def test_copy_binary_isolates_from_source_mutation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            src = tmp_path / "FlexAIDdS"
            src.write_bytes(b"BINARY-V1")
            run = tmp_path / "run_ns"
            pin = self.g.copy_binary_to_run_namespace(src, run)
            dest = Path(pin["dest"])
            self.assertTrue(dest.is_file())
            self.assertEqual(dest.read_bytes(), b"BINARY-V1")
            digest1 = pin["sha256"]
            self.assertEqual(self.g.sha256_file(dest), digest1)
            # Mutate shared source; namespace copy must stay V1
            src.write_bytes(b"BINARY-V2-REBUILD")
            self.assertEqual(dest.read_bytes(), b"BINARY-V1")
            self.assertEqual(self.g.sha256_file(dest), digest1)
            self.assertNotEqual(self.g.sha256_file(src), digest1)
            pin_file = run / "bin" / "FlexAIDdS.SHA256.json"
            self.assertTrue(pin_file.is_file())
            data = json.loads(pin_file.read_text(encoding="utf-8"))
            self.assertEqual(data["sha256"], digest1)

    def test_preflight_hold_refuses(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "BENCHMARK_HOLD.json").write_text('{"active":true}\n', encoding="utf-8")
            out = root / "out"
            out.mkdir()
            result = self.g.preflight_dock(
                out_dir=out,
                binary=None,
                workers=1,
                acquire_lock=False,
                copy_binary=False,
                repo_root=root,
                free_gib=50.0,
            )
            self.assertFalse(result.ok)
            self.assertEqual(result.exit_code, self.g.EXIT_HOLD)

    def test_preflight_dual_lock_refuses_second(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "out"
            out.mkdir()
            lock = root / "dock_session.lock"
            r1 = self.g.preflight_dock(
                out_dir=out,
                binary=None,
                workers=2,
                acquire_lock=True,
                copy_binary=False,
                repo_root=root,
                lock_dir=lock,
                free_gib=50.0,
                max_workers=4,
                owner="owner-a",
            )
            self.assertTrue(r1.ok, r1.messages)
            r2 = self.g.preflight_dock(
                out_dir=out,
                binary=None,
                workers=1,
                acquire_lock=True,
                copy_binary=False,
                repo_root=root,
                lock_dir=lock,
                free_gib=50.0,
                owner="owner-b",
            )
            self.assertFalse(r2.ok)
            self.assertEqual(r2.exit_code, self.g.EXIT_LOCK)
            self.g.release_lock(lock, force=True)

    def test_cli_check_hold_exit_code(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "BENCHMARK_HOLD.json").write_text('{"active":true}\n', encoding="utf-8")
            proc = subprocess.run(
                [sys.executable, str(GUARD_PATH), "check-hold", "--repo-root", str(root)],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc.returncode, 78)
            self.assertIn("BENCHMARK_HOLD", proc.stderr)

    def test_cli_check_hold_absent_ok(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            proc = subprocess.run(
                [sys.executable, str(GUARD_PATH), "check-hold", "--repo-root", tmp],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc.returncode, 0)


if __name__ == "__main__":
    unittest.main()
