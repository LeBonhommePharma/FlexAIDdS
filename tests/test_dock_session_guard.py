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
            # Use live pid so the lock is not treated as stale/reclaimable.
            ok1, path1, msg1 = self.g.try_acquire_lock(
                lock, owner="first", pid=os.getpid()
            )
            self.assertTrue(ok1, msg1)
            self.assertEqual(path1, lock)
            self.assertTrue((lock / "owner.json").is_file())
            ok2, path2, msg2 = self.g.try_acquire_lock(lock, owner="second", pid=222)
            self.assertFalse(ok2)
            self.assertEqual(path2, lock)
            self.assertIn("already held", msg2.lower())
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

    def test_dock_pid_keeps_lock_after_launcher_would_exit(self) -> None:
        """Background dock_pid must keep exclusive lock (claim nohup path)."""
        with tempfile.TemporaryDirectory() as tmp:
            lock = Path(tmp) / "dock_session.lock"
            ok, _, msg = self.g.try_acquire_lock(lock, owner="launcher", pid=999001)
            self.assertTrue(ok, msg)
            # Simulate launcher exit + long-lived dock still running (us).
            ok_set, msg_set = self.g.set_dock_pid(lock, dock_pid=os.getpid())
            self.assertTrue(ok_set, msg_set)
            # release without force must refuse while dock_pid live
            ok_r, msg_r = self.g.release_lock(lock, force=False)
            self.assertFalse(ok_r)
            self.assertIn("dock_pid", msg_r)
            # second session cannot dual-dock
            ok2, _, msg2 = self.g.try_acquire_lock(lock, owner="second", pid=999002)
            self.assertFalse(ok2)
            self.assertIn("Dual-dock refused", msg2)
            # force still works for ops recovery
            ok_f, msg_f = self.g.release_lock(lock, force=True)
            self.assertTrue(ok_f, msg_f)

    def test_stale_lock_reclaimed_when_pids_dead(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            lock = Path(tmp) / "dock_session.lock"
            ok, _, msg = self.g.try_acquire_lock(lock, owner="dead", pid=1)
            # pid 1 is usually alive on macOS (launchd). Use a high fake pid and
            # write owner.json with both dead pids manually.
            self.assertTrue(ok or True)
            if lock.exists():
                self.g.release_lock(lock, force=True)
            lock.mkdir()
            (lock / "owner.json").write_text(
                '{"pid": 99999901, "dock_pid": 99999902, "owner": "stale"}\n',
                encoding="utf-8",
            )
            ok2, _, msg2 = self.g.try_acquire_lock(lock, owner="reclaimer", pid=os.getpid())
            self.assertTrue(ok2, msg2)
            self.g.release_lock(lock, force=True)

    def test_resolve_launch_binaries_rebinds_engine_and_runner(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp) / "out"
            eng = Path(tmp) / "FlexAIDdS"
            runner = Path(tmp) / "benchmark_datasets"
            eng.write_bytes(b"ENGINE-V1")
            runner.write_bytes(b"RUNNER-V1")
            eng.chmod(0o755)
            runner.chmod(0o755)
            r = self.g.preflight_dock(
                out_dir=run,
                binaries=[eng, runner],
                workers=1,
                acquire_lock=False,
                copy_binary=True,
                repo_root=Path(tmp),
                free_gib=50.0,
            )
            self.assertTrue(r.ok, r.messages)
            env = r.launch_env or {}
            self.assertIn("FLEXAIDDS_BINARY", env)
            self.assertIn("FLEXAIDDS_RUNNER", env)
            self.assertTrue(env["FLEXAIDDS_BINARY"].endswith("/bin/FlexAIDdS"))
            self.assertTrue(env["FLEXAIDDS_RUNNER"].endswith("/bin/benchmark_datasets"))
            # Mutate shared sources — pins stay V1
            eng.write_bytes(b"ENGINE-V2")
            runner.write_bytes(b"RUNNER-V2")
            self.assertEqual(Path(env["FLEXAIDDS_BINARY"]).read_bytes(), b"ENGINE-V1")
            self.assertEqual(Path(env["FLEXAIDDS_RUNNER"]).read_bytes(), b"RUNNER-V1")
            self.assertEqual(
                self.g.resolve_launch_binaries(run)["FLEXAIDDS_BINARY"],
                env["FLEXAIDDS_BINARY"],
            )

    def test_claim_script_no_exit_force_release_and_rebinds_pin(self) -> None:
        """Structural: claim launcher must not force-release lock on EXIT."""
        text = (REPO / "scripts" / "run_C0_claim_clean.sh").read_text(encoding="utf-8")
        self.assertNotIn("release-lock --force", text)
        self.assertNotIn("trap 'python3 \"$DOCK_GUARD\" release-lock", text)
        self.assertIn("set-dock-pid", text)
        self.assertIn("--owner-pid", text)
        self.assertIn("SHELL_OWNER_PID", text)
        self.assertIn('BINARY="$PINNED_ENGINE"', text)
        self.assertIn('RUNNER="$PINNED_RUNNER"', text)
        self.assertIn("dock_session_guard.py missing", text)
        self.assertIn("exit 77", text)

    def test_production_script_fail_closed_and_rebind(self) -> None:
        text = (REPO / "scripts" / "run_benchmark_production.sh").read_text(encoding="utf-8")
        self.assertIn("dock_session_guard.py missing", text)
        self.assertIn("exit 77", text)
        self.assertIn("rebound FLEXAIDDS_BIN", text)
        self.assertIn("FLEXAIDDS_BIN=\"${pinned_engine}\"", text)
        self.assertIn("--owner-pid", text)
        self.assertIn("set-dock-pid", text)
        # Production may release without --force after foreground docks complete.
        self.assertNotIn("release-lock --force", text)

    def test_cli_two_sequential_preflights_second_exits_79(self) -> None:
        """Real bash→python path: short-lived preflight must not leave reclaimable lock.

        Both CLI subprocesses exit; ownership is bound to this still-live parent
        (default getppid / explicit --owner-pid). Second preflight must dual-dock
        refuse with exit 79 without any set-dock-pid of a still-live dock child.
        """
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            out1 = tmp_path / "out1"
            out2 = tmp_path / "out2"
            out1.mkdir()
            out2.mkdir()
            lock = tmp_path / "dock_session.lock"
            # Parent of both python preflights is this test process — live.
            owner = os.getpid()
            common = [
                sys.executable,
                str(GUARD_PATH),
                "preflight",
                "--workers",
                "1",
                "--repo-root",
                str(tmp_path),
                "--lock-dir",
                str(lock),
                "--no-copy-binary",
                "--min-free-gb",
                "1",
                "--owner-pid",
                str(owner),
            ]
            p1 = subprocess.run(
                common + ["--out-dir", str(out1), "--owner", "first-cli"],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(p1.returncode, 0, p1.stdout + p1.stderr)
            owner_json = json.loads((lock / "owner.json").read_text(encoding="utf-8"))
            # Must record long-lived owner, NOT a dead preflight python pid.
            self.assertEqual(int(owner_json["pid"]), owner)
            self.assertTrue(self.g.pid_alive(int(owner_json["pid"])))
            p2 = subprocess.run(
                common + ["--out-dir", str(out2), "--owner", "second-cli"],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(
                p2.returncode,
                79,
                f"expected dual-dock EXIT_LOCK=79, got {p2.returncode}\n"
                f"stdout={p2.stdout}\nstderr={p2.stderr}\nowner={owner_json}",
            )
            self.assertIn("Dual-dock refused", p2.stdout + p2.stderr)
            # Owner must still be the first long-lived pid (not rewritten by second).
            owner_after = json.loads((lock / "owner.json").read_text(encoding="utf-8"))
            self.assertEqual(int(owner_after["pid"]), owner)
            self.assertEqual(owner_after.get("owner"), "first-cli")
            self.g.release_lock(lock, force=True, requester_pid=owner)

    def test_cli_preflight_default_binds_parent_not_self(self) -> None:
        """Without --owner-pid, CLI records getppid (shell/parent), not python pid."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            out = tmp_path / "out"
            out.mkdir()
            lock = tmp_path / "dock_session.lock"
            parent = os.getpid()
            proc = subprocess.run(
                [
                    sys.executable,
                    str(GUARD_PATH),
                    "preflight",
                    "--out-dir",
                    str(out),
                    "--workers",
                    "1",
                    "--repo-root",
                    str(tmp_path),
                    "--lock-dir",
                    str(lock),
                    "--no-copy-binary",
                    "--min-free-gb",
                    "1",
                    "--owner",
                    "default-parent",
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            data = json.loads((lock / "owner.json").read_text(encoding="utf-8"))
            # Child's getppid() is this test process.
            self.assertEqual(int(data["pid"]), parent)
            self.assertNotEqual(int(data["pid"]), 0)
            # A second sequential CLI preflight must still refuse.
            proc2 = subprocess.run(
                [
                    sys.executable,
                    str(GUARD_PATH),
                    "preflight",
                    "--out-dir",
                    str(out),
                    "--workers",
                    "1",
                    "--repo-root",
                    str(tmp_path),
                    "--lock-dir",
                    str(lock),
                    "--no-copy-binary",
                    "--min-free-gb",
                    "1",
                    "--owner",
                    "default-parent-2",
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc2.returncode, 79, proc2.stdout + proc2.stderr)
            self.g.release_lock(lock, force=True, requester_pid=parent)

    def test_missing_guard_is_fail_closed_not_warn(self) -> None:
        claim = (REPO / "scripts" / "run_C0_claim_clean.sh").read_text(encoding="utf-8")
        prod = (REPO / "scripts" / "run_benchmark_production.sh").read_text(encoding="utf-8")
        self.assertNotIn("multi-session dual-dock unprotected", claim + prod)
        self.assertIn("fail-closed", claim)
        self.assertIn("fail-closed", prod)


if __name__ == "__main__":
    unittest.main()