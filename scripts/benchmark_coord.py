#!/usr/bin/env python3
"""Fail-closed multi-session dock coordination (Sol #9 / Codex).

Gates new dock launches so two independent agents cannot OOM or stomp binaries:

  1. BENCHMARK_HOLD.json  — operator hold; any presence refuses launch
  2. Atomic mkdir lock    — portable exclusive owner (APFS-safe; not flock)
  3. Disk floor           — refuse if free space < floor (default 20 GiB)
  4. WORKERS ≤ 4          — hard refuse above cap
  5. Binary stamp         — copy engine into run OUT; use stamped path

Does NOT kill or re-parent an in-flight dock. Live G4.1 (or any owner) may
already hold the box; this module only gates *new* launches.

Usage:
  python3 scripts/benchmark_coord.py status
  python3 scripts/benchmark_coord.py hold --reason "…" --owner "session-name"
  python3 scripts/benchmark_coord.py unhold
  python3 scripts/benchmark_coord.py preflight --out DIR --workers N --binary PATH
  python3 scripts/benchmark_coord.py release --token TOKEN

Environment:
  FLEXAIDDS_LOCAL_ROOT   results root (default ~/flexaidds_results)
  FLEXAIDDS_DISK_FLOOR_GB  free-space floor (default 20)
  FLEXAIDDS_DISK_FLOOR_OVERRIDE=1  emergency skip disk refuse (logged)

Apache-2.0 · Le Bonhomme Pharma
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import socket
import sys
import time
import uuid
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

DEFAULT_FLOOR_GB = 20.0
MAX_WORKERS = 4
HOLD_NAME = "BENCHMARK_HOLD.json"
LOCK_DIRNAME = "BENCHMARK_DOCK_LOCK"
LOCK_META = "owner.json"


def default_results_root() -> Path:
    """Local dock coordination root — never iCloud thin mirror.

    Prefer FLEXAIDDS_LOCAL_ROOT; fall back to ~/flexaidds_results.
    Do **not** use FLEXAIDDS_RESULTS (often a campaigns/results subpath).
    """
    env = os.environ.get("FLEXAIDDS_LOCAL_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    return (Path.home() / "flexaidds_results").resolve()


def hold_path(root: Path) -> Path:
    return root / HOLD_NAME


def lock_dir(root: Path) -> Path:
    return root / LOCK_DIRNAME


def disk_floor_gb() -> float:
    raw = os.environ.get("FLEXAIDDS_DISK_FLOOR_GB")
    if raw is None or raw.strip() == "":
        return DEFAULT_FLOOR_GB
    try:
        return float(raw)
    except ValueError:
        return DEFAULT_FLOOR_GB


def free_gib(path: Path) -> float:
    """Free space in GiB for the filesystem containing path (must exist)."""
    p = path if path.exists() else path.parent
    usage = shutil.disk_usage(str(p))
    return usage.free / (1024.0**3)


@dataclass
class PreflightResult:
    ok: bool
    reason: str = ""
    lock_token: Optional[str] = None
    stamped_binary: Optional[str] = None
    free_gib: Optional[float] = None
    hold_path: Optional[str] = None
    lock_path: Optional[str] = None
    details: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def check_hold(root: Path) -> Optional[Dict[str, Any]]:
    """Return hold metadata if present, else None."""
    p = hold_path(root)
    if not p.is_file():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {"raw": True, "path": str(p), "error": "unreadable hold file"}


def write_hold(
    root: Path,
    *,
    owner: str,
    reason: str,
    extra: Optional[Dict[str, Any]] = None,
) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    p = hold_path(root)
    body = {
        "schema": "benchmark_hold/v1",
        "owner": owner,
        "reason": reason,
        "written_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "host": socket.gethostname(),
        "pid": os.getpid(),
    }
    if extra:
        body.update(extra)
    p.write_text(json.dumps(body, indent=2) + "\n", encoding="utf-8")
    return p


def clear_hold(root: Path) -> bool:
    p = hold_path(root)
    if p.is_file():
        p.unlink()
        return True
    return False


def acquire_lock(
    root: Path,
    *,
    owner: str,
    purpose: str = "dock",
    out: Optional[str] = None,
) -> Tuple[bool, str, Optional[str]]:
    """Atomic exclusive lock via mkdir. Returns (ok, message, token)."""
    root.mkdir(parents=True, exist_ok=True)
    d = lock_dir(root)
    token = str(uuid.uuid4())
    try:
        d.mkdir(exist_ok=False)
    except FileExistsError:
        meta_path = d / LOCK_META
        meta: Dict[str, Any] = {}
        if meta_path.is_file():
            try:
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                meta = {"error": "unreadable"}
        return (
            False,
            f"dock lock held at {d} (owner={meta.get('owner')!r} purpose={meta.get('purpose')!r})",
            None,
        )
    meta = {
        "schema": "benchmark_dock_lock/v1",
        "token": token,
        "owner": owner,
        "purpose": purpose,
        "out": out,
        "pid": os.getpid(),
        "host": socket.gethostname(),
        "acquired_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    (d / LOCK_META).write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    return True, f"lock acquired {d}", token


def release_lock(root: Path, token: Optional[str] = None, *, force: bool = False) -> Tuple[bool, str]:
    """Release lock if token matches (or force=True for operator recovery)."""
    d = lock_dir(root)
    if not d.is_dir():
        return True, "no lock present"
    meta_path = d / LOCK_META
    if not force and token is not None and meta_path.is_file():
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return False, "lock meta unreadable; use force=True"
        if meta.get("token") != token:
            return False, "token mismatch; refusing release"
    # remove contents then dir
    for child in d.iterdir():
        if child.is_file():
            child.unlink()
        elif child.is_dir():
            shutil.rmtree(child)
    d.rmdir()
    return True, f"lock released {d}"


def check_workers(workers: int, max_workers: int = MAX_WORKERS) -> Tuple[bool, str]:
    if workers < 1:
        return False, f"workers={workers} invalid (need >=1)"
    if workers > max_workers:
        return False, f"workers={workers} exceeds MAX_WORKERS={max_workers} (refuse)"
    return True, f"workers={workers} ok"


def check_disk(path: Path, floor_gb: Optional[float] = None) -> Tuple[bool, str, float]:
    floor = disk_floor_gb() if floor_gb is None else floor_gb
    free = free_gib(path)
    if os.environ.get("FLEXAIDDS_DISK_FLOOR_OVERRIDE", "").strip() in ("1", "true", "yes"):
        return True, f"disk override active free={free:.2f}GiB floor={floor}", free
    if free < floor:
        return False, f"free={free:.2f}GiB < floor={floor}GiB", free
    return True, f"free={free:.2f}GiB >= floor={floor}", free


def stamp_binary(binary: Path, run_out: Path, name: str = "FlexAIDdS.stamped") -> Path:
    """Copy engine binary into run OUT; return stamped path."""
    binary = binary.expanduser().resolve()
    if not binary.is_file():
        raise FileNotFoundError(f"binary not found: {binary}")
    run_out.mkdir(parents=True, exist_ok=True)
    dest = run_out / "bin" / name
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(binary, dest)
    try:
        dest.chmod(dest.stat().st_mode | 0o111)
    except OSError:
        pass
    # receipt
    receipt = {
        "source": str(binary),
        "stamped": str(dest),
        "size": dest.stat().st_size,
        "mtime_src": binary.stat().st_mtime,
        "stamped_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    (run_out / "bin" / "BINARY_STAMP.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )
    return dest


def preflight(
    *,
    root: Optional[Path] = None,
    out: Path,
    workers: int,
    binary: Path,
    owner: str = "agent",
    purpose: str = "dock",
    max_workers: int = MAX_WORKERS,
    floor_gb: Optional[float] = None,
    acquire: bool = True,
    stamp: bool = True,
) -> PreflightResult:
    """Fail-closed preflight for a new dock launch.

    On success (and acquire=True): holds mkdir lock and optionally stamps binary.
    Caller must release_lock(root, token) when the dock finishes.
    """
    root = (root or default_results_root()).expanduser().resolve()
    out = out.expanduser().resolve()

    hold = check_hold(root)
    if hold is not None:
        return PreflightResult(
            ok=False,
            reason=f"BENCHMARK_HOLD present: {hold_path(root)} owner={hold.get('owner')!r} reason={hold.get('reason')!r}",
            hold_path=str(hold_path(root)),
            details={"hold": hold},
        )

    ok_w, msg_w = check_workers(workers, max_workers=max_workers)
    if not ok_w:
        return PreflightResult(ok=False, reason=msg_w)

    # disk against out parent / results root
    probe = out if out.parent.exists() else root
    if not probe.exists():
        probe.mkdir(parents=True, exist_ok=True)
    ok_d, msg_d, free = check_disk(probe, floor_gb=floor_gb)
    if not ok_d:
        return PreflightResult(
            ok=False,
            reason=msg_d,
            free_gib=free,
            details={"disk": msg_d},
        )

    token = None
    lock_path = str(lock_dir(root))
    if acquire:
        ok_l, msg_l, token = acquire_lock(
            root, owner=owner, purpose=purpose, out=str(out)
        )
        if not ok_l:
            return PreflightResult(
                ok=False,
                reason=msg_l,
                free_gib=free,
                lock_path=lock_path,
            )

    stamped = None
    if stamp:
        try:
            stamped = stamp_binary(binary, out)
        except (OSError, FileNotFoundError) as e:
            if token:
                release_lock(root, token)
            return PreflightResult(
                ok=False,
                reason=f"binary stamp failed: {e}",
                free_gib=free,
                lock_token=None,
                lock_path=lock_path,
            )

    return PreflightResult(
        ok=True,
        reason="preflight ok",
        lock_token=token,
        stamped_binary=str(stamped) if stamped else None,
        free_gib=free,
        lock_path=lock_path,
        details={"workers": workers, "disk": msg_d, "owner": owner},
    )


def status(root: Optional[Path] = None) -> Dict[str, Any]:
    root = (root or default_results_root()).expanduser().resolve()
    hold = check_hold(root)
    ld = lock_dir(root)
    lock_meta = None
    if (ld / LOCK_META).is_file():
        try:
            lock_meta = json.loads((ld / LOCK_META).read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            lock_meta = {"error": "unreadable"}
    free = free_gib(root if root.exists() else Path.home())
    return {
        "results_root": str(root),
        "hold": hold,
        "hold_path": str(hold_path(root)),
        "lock_held": ld.is_dir(),
        "lock_meta": lock_meta,
        "lock_path": str(ld),
        "free_gib": free,
        "disk_floor_gb": disk_floor_gb(),
        "max_workers": MAX_WORKERS,
        "may_dock": hold is None and not ld.is_dir() and free >= disk_floor_gb(),
    }


def offline_queue_doc() -> str:
    return """# Offline Benchmarks queue (no docking)

While another session holds the dock lock / BENCHMARK_HOLD, this role is **offline only**.

## Rules
1. **Do not launch docks** while `~/flexaidds_results/BENCHMARK_HOLD.json` exists
   or `~/flexaidds_results/BENCHMARK_DOCK_LOCK/` is present.
2. **Read the other session's commits first** (`git log -5 --oneline`) — do not re-derive
   HEM / population / baseline debates already on main.
3. Owner (dock session) only: WORKERS≤4; stamp binary into OUT; never rebuild while the
   other session's run is live.

## Offline queue (highest value first)
1. **Full-population ceiling** on SEARCH-MISS (1J3J 1K3U 1L7F 1N1M 1M2Z) — frozen poses only.
2. **Per-term CF decomposition** of SCORING-LOCKED gaps (+17.9 / +28.8 / +70.2) — frozen poses.
3. **Matrix pin resolve** — `72d7` vs `9dc9` canonical (`md5` / git only).
4. **Crystal-reference PoseBusters ceiling** — no docking; bounds strict claims.

## Preflight for the dock owner
```bash
python3 scripts/benchmark_coord.py status
python3 scripts/benchmark_coord.py preflight \\
  --out ~/flexaidds_results/my_run --workers 2 \\
  --binary build/FlexAIDdS --owner "session-name"
# on success: export stamped path from JSON; on exit: release --token …
```
"""


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--root",
        type=Path,
        default=None,
        help="results root (default FLEXAIDDS_LOCAL_ROOT or ~/flexaidds_results)",
    )
    sub = ap.add_subparsers(dest="cmd", required=True)

    sub.add_parser("status", help="print hold/lock/disk status JSON")

    p_hold = sub.add_parser("hold", help="create BENCHMARK_HOLD.json")
    p_hold.add_argument("--owner", required=True)
    p_hold.add_argument("--reason", required=True)

    sub.add_parser("unhold", help="remove BENCHMARK_HOLD.json")

    p_pre = sub.add_parser("preflight", help="fail-closed dock preflight")
    p_pre.add_argument("--out", type=Path, required=True)
    p_pre.add_argument("--workers", type=int, required=True)
    p_pre.add_argument("--binary", type=Path, required=True)
    p_pre.add_argument("--owner", default="agent")
    p_pre.add_argument("--purpose", default="dock")
    p_pre.add_argument("--no-acquire", action="store_true")
    p_pre.add_argument("--no-stamp", action="store_true")

    p_rel = sub.add_parser("release", help="release dock lock")
    p_rel.add_argument("--token", default=None)
    p_rel.add_argument("--force", action="store_true")

    sub.add_parser("offline-queue", help="print offline role queue markdown")

    args = ap.parse_args(argv)
    root = args.root.expanduser().resolve() if args.root else default_results_root()

    if args.cmd == "status":
        print(json.dumps(status(root), indent=2))
        return 0
    if args.cmd == "hold":
        p = write_hold(root, owner=args.owner, reason=args.reason)
        print(json.dumps({"ok": True, "hold": str(p)}, indent=2))
        return 0
    if args.cmd == "unhold":
        cleared = clear_hold(root)
        print(json.dumps({"ok": True, "cleared": cleared}, indent=2))
        return 0
    if args.cmd == "preflight":
        r = preflight(
            root=root,
            out=args.out,
            workers=args.workers,
            binary=args.binary,
            owner=args.owner,
            purpose=args.purpose,
            acquire=not args.no_acquire,
            stamp=not args.no_stamp,
        )
        print(json.dumps(r.to_dict(), indent=2))
        return 0 if r.ok else 2
    if args.cmd == "release":
        ok, msg = release_lock(root, args.token, force=args.force)
        print(json.dumps({"ok": ok, "message": msg}, indent=2))
        return 0 if ok else 2
    if args.cmd == "offline-queue":
        print(offline_queue_doc())
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
