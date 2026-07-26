#!/usr/bin/env python3
"""Fail-closed multi-session docking preflight (Sol #9 / Codex recommendation).

Prevents two independent agent sessions from dual-docking an 18 GB MacBook:

1. **BENCHMARK_HOLD.json** — any presence refuses launch (exit 78).
2. **mkdir lock** — APFS-atomic exclusive ownership (not flock).
3. **Disk floor** — default 20 GiB free on the OUT filesystem.
4. **Binary isolation** — copy engine/runner into the run namespace + SHA-256 pin
   so a rebuild of a shared build tree cannot rewrite a live run's binary.
5. **Workers cap** — default WORKERS ≤ 4 on this box.

Env overrides (tests use temp paths):
  FLEXAIDDS_BENCHMARK_HOLD_PATH   absolute path to hold file (checked first)
  FLEXAIDDS_DOCK_LOCK_DIR         directory whose *creation* is the lock
  FLEXAIDDS_MIN_FREE_GB           default 20
  FLEXAIDDS_MAX_WORKERS           default 4

Copyright 2026 Le Bonhomme Pharma
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, Optional

EXIT_OK = 0
EXIT_HOLD = 78
EXIT_LOCK = 79
EXIT_DISK = 80
EXIT_WORKERS = 81
EXIT_USAGE = 2

DEFAULT_MIN_FREE_GB = 20
DEFAULT_MAX_WORKERS = 4
HOLD_FILENAME = "BENCHMARK_HOLD.json"
LOCK_DIRNAME = "dock_session.lock"


def repo_root_from_here() -> Path:
    return Path(__file__).resolve().parents[1]


def hold_candidates(repo_root: Optional[Path] = None) -> list[Path]:
    """Paths checked for an active hold (first existing wins for messaging)."""
    paths: list[Path] = []
    env = os.environ.get("FLEXAIDDS_BENCHMARK_HOLD_PATH", "").strip()
    if env:
        paths.append(Path(os.path.expandvars(env)).expanduser())
    root = repo_root or repo_root_from_here()
    paths.append(root / HOLD_FILENAME)
    local = os.environ.get("FLEXAIDDS_LOCAL_ROOT", "").strip()
    if local:
        paths.append(Path(os.path.expandvars(local)).expanduser() / HOLD_FILENAME)
    # Common durable local results root (no CloudDocs walk).
    paths.append(Path.home() / "flexaidds_results" / HOLD_FILENAME)
    return paths


def find_hold_file(repo_root: Optional[Path] = None) -> Optional[Path]:
    for path in hold_candidates(repo_root):
        try:
            if path.is_file():
                return path
        except OSError:
            continue
    return None


def hold_blocks_launch(repo_root: Optional[Path] = None) -> tuple[bool, str]:
    found = find_hold_file(repo_root)
    if found is None:
        return False, ""
    try:
        text = found.read_text(encoding="utf-8", errors="replace")[:500]
    except OSError as exc:
        text = f"(unreadable: {exc})"
    msg = (
        f"BENCHMARK_HOLD active: refuse dock launch while {found} exists.\n"
        f"Remove/rename the hold only after LP releases the multi-session hold.\n"
        f"Preview:\n{text}"
    )
    return True, msg


def default_lock_dir() -> Path:
    env = os.environ.get("FLEXAIDDS_DOCK_LOCK_DIR", "").strip()
    if env:
        return Path(os.path.expandvars(env)).expanduser()
    local = os.environ.get("FLEXAIDDS_LOCAL_ROOT", "").strip()
    if local:
        return Path(os.path.expandvars(local)).expanduser() / LOCK_DIRNAME
    return Path.home() / "flexaidds_results" / LOCK_DIRNAME


@dataclass
class LockInfo:
    lock_dir: str
    pid: int
    owner: str
    created_utc: str
    out_dir: str = ""
    note: str = ""


def try_acquire_lock(
    lock_dir: Optional[Path] = None,
    *,
    owner: str = "",
    out_dir: str = "",
    note: str = "",
    pid: Optional[int] = None,
) -> tuple[bool, Path, str]:
    """Acquire exclusive dock ownership via mkdir (atomic on APFS/POSIX).

    Returns (ok, lock_dir, message). On success the caller owns the lock until
    release_lock(). On failure another session holds the box — go offline.
    """
    path = Path(lock_dir) if lock_dir is not None else default_lock_dir()
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        path.mkdir(exist_ok=False)
    except FileExistsError:
        holder = ""
        meta = path / "owner.json"
        if meta.is_file():
            try:
                holder = meta.read_text(encoding="utf-8", errors="replace")[:400]
            except OSError:
                holder = "(owner.json unreadable)"
        return (
            False,
            path,
            f"Dock lock already held at {path}. Dual-dock refused.\n{holder}",
        )
    except OSError as exc:
        return False, path, f"Dock lock mkdir failed at {path}: {exc}"

    info = LockInfo(
        lock_dir=str(path),
        pid=int(pid if pid is not None else os.getpid()),
        owner=owner or os.environ.get("USER", "unknown"),
        created_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        out_dir=out_dir,
        note=note,
    )
    try:
        (path / "owner.json").write_text(
            json.dumps(asdict(info), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    except OSError as exc:
        # Best-effort cleanup if we cannot record ownership.
        try:
            shutil.rmtree(path, ignore_errors=True)
        except OSError:
            pass
        return False, path, f"Dock lock acquired but owner.json write failed: {exc}"
    return True, path, f"Dock lock acquired: {path}"


def release_lock(lock_dir: Optional[Path] = None, *, force: bool = False) -> tuple[bool, str]:
    path = Path(lock_dir) if lock_dir is not None else default_lock_dir()
    if not path.exists():
        return True, f"No lock present at {path}"
    meta = path / "owner.json"
    if meta.is_file() and not force:
        try:
            data = json.loads(meta.read_text(encoding="utf-8"))
            owner_pid = int(data.get("pid", -1))
            if owner_pid not in (-1, os.getpid()):
                # Stale-PID cleanup: if process is gone, allow release.
                try:
                    os.kill(owner_pid, 0)
                    alive = True
                except OSError:
                    alive = False
                if alive:
                    return (
                        False,
                        f"Refuse release of lock owned by live pid={owner_pid} "
                        f"(use force=True only for explicit ops recovery)",
                    )
        except (OSError, json.JSONDecodeError, TypeError, ValueError):
            pass
    try:
        shutil.rmtree(path)
    except OSError as exc:
        return False, f"Failed to remove lock {path}: {exc}"
    return True, f"Dock lock released: {path}"


def free_disk_gib(path: Path, *, statvfs: Optional[Callable[[str], os.statvfs_result]] = None) -> float:
    """Return free space in GiB for the filesystem containing path."""
    probe = Path(path)
    # Walk up to an existing ancestor so pre-create OUT dirs still work.
    while not probe.exists() and probe != probe.parent:
        probe = probe.parent
    if not probe.exists():
        probe = Path.cwd()
    fn = statvfs or os.statvfs
    st = fn(str(probe))
    return (st.f_bavail * st.f_frsize) / (1024.0**3)


def disk_ok(
    path: Path,
    *,
    min_free_gb: Optional[float] = None,
    free_gib: Optional[float] = None,
) -> tuple[bool, str, float]:
    floor = float(
        min_free_gb
        if min_free_gb is not None
        else os.environ.get("FLEXAIDDS_MIN_FREE_GB", DEFAULT_MIN_FREE_GB)
    )
    available = float(free_gib) if free_gib is not None else free_disk_gib(path)
    if available + 1e-9 < floor:
        return (
            False,
            f"Disk floor fail: {available:.2f} GiB free at {path} "
            f"(need ≥ {floor:.0f} GiB). Refuse dock launch.",
            available,
        )
    return True, f"Disk OK: {available:.2f} GiB free (≥ {floor:.0f} GiB)", available


def workers_ok(workers: int, *, max_workers: Optional[int] = None) -> tuple[bool, str]:
    cap = int(
        max_workers
        if max_workers is not None
        else os.environ.get("FLEXAIDDS_MAX_WORKERS", DEFAULT_MAX_WORKERS)
    )
    if workers < 1:
        return False, f"Invalid workers={workers} (must be ≥ 1)"
    if workers > cap:
        return (
            False,
            f"Workers cap fail: requested {workers} > max {cap} on 18 GB box. "
            f"Use WORKERS≤{cap}.",
        )
    return True, f"Workers OK: {workers} ≤ {cap}"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def copy_binary_to_run_namespace(
    source: Path,
    run_dir: Path,
    *,
    dest_name: Optional[str] = None,
) -> dict[str, str]:
    """Copy binary into run_dir/bin/ and write SHA256 pin. Source may later change."""
    src = Path(source).resolve()
    if not src.is_file():
        raise FileNotFoundError(f"Binary not found: {src}")
    run = Path(run_dir)
    bin_dir = run / "bin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    name = dest_name or src.name
    dest = bin_dir / name
    shutil.copy2(src, dest)
    dest.chmod(dest.stat().st_mode | 0o111)
    digest = sha256_file(dest)
    pin = {
        "source": str(src),
        "dest": str(dest),
        "sha256": digest,
        "bytes": str(dest.stat().st_size),
        "copied_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    pin_path = bin_dir / f"{name}.SHA256.json"
    pin_path.write_text(json.dumps(pin, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (bin_dir / f"{name}.SHA256").write_text(digest + "\n", encoding="utf-8")
    return pin


@dataclass
class PreflightResult:
    ok: bool
    exit_code: int
    messages: list[str]
    hold_path: str = ""
    lock_dir: str = ""
    free_gib: float = -1.0
    binary_pin: Optional[dict[str, str]] = None


def preflight_dock(
    *,
    out_dir: Path,
    binary: Optional[Path] = None,
    workers: int = 1,
    acquire_lock: bool = True,
    copy_binary: bool = True,
    repo_root: Optional[Path] = None,
    lock_dir: Optional[Path] = None,
    min_free_gb: Optional[float] = None,
    max_workers: Optional[int] = None,
    free_gib: Optional[float] = None,
    owner: str = "",
    note: str = "",
) -> PreflightResult:
    """Single entry used by launch scripts and unit tests."""
    messages: list[str] = []
    blocked, hold_msg = hold_blocks_launch(repo_root)
    if blocked:
        return PreflightResult(
            ok=False,
            exit_code=EXIT_HOLD,
            messages=[hold_msg],
            hold_path=str(find_hold_file(repo_root) or ""),
        )

    ok_w, msg_w = workers_ok(workers, max_workers=max_workers)
    messages.append(msg_w)
    if not ok_w:
        return PreflightResult(ok=False, exit_code=EXIT_WORKERS, messages=messages)

    ok_d, msg_d, avail = disk_ok(Path(out_dir), min_free_gb=min_free_gb, free_gib=free_gib)
    messages.append(msg_d)
    if not ok_d:
        return PreflightResult(
            ok=False, exit_code=EXIT_DISK, messages=messages, free_gib=avail
        )

    lock_path = Path(lock_dir) if lock_dir is not None else default_lock_dir()
    if acquire_lock:
        ok_l, lock_path, msg_l = try_acquire_lock(
            lock_path,
            owner=owner,
            out_dir=str(out_dir),
            note=note,
        )
        messages.append(msg_l)
        if not ok_l:
            return PreflightResult(
                ok=False,
                exit_code=EXIT_LOCK,
                messages=messages,
                lock_dir=str(lock_path),
                free_gib=avail,
            )

    pin: Optional[dict[str, str]] = None
    if copy_binary and binary is not None:
        try:
            pin = copy_binary_to_run_namespace(Path(binary), Path(out_dir))
            messages.append(
                f"Binary isolated: {pin['dest']} sha256={pin['sha256'][:16]}…"
            )
        except (OSError, FileNotFoundError) as exc:
            if acquire_lock:
                release_lock(lock_path, force=True)
            return PreflightResult(
                ok=False,
                exit_code=EXIT_USAGE,
                messages=messages + [f"Binary copy failed: {exc}"],
                lock_dir=str(lock_path),
                free_gib=avail,
            )

    return PreflightResult(
        ok=True,
        exit_code=EXIT_OK,
        messages=messages + ["Preflight OK — exclusive dock ownership granted."],
        lock_dir=str(lock_path) if acquire_lock else "",
        free_gib=avail,
        binary_pin=pin,
    )


def _cmd_check_hold(args: argparse.Namespace) -> int:
    blocked, msg = hold_blocks_launch(
        Path(args.repo_root) if args.repo_root else None
    )
    if blocked:
        print(msg, file=sys.stderr)
        return EXIT_HOLD
    print("No BENCHMARK_HOLD present.")
    return EXIT_OK


def _cmd_preflight(args: argparse.Namespace) -> int:
    result = preflight_dock(
        out_dir=Path(args.out_dir),
        binary=Path(args.binary) if args.binary else None,
        workers=int(args.workers),
        acquire_lock=not args.no_lock,
        copy_binary=bool(args.binary) and not args.no_copy_binary,
        repo_root=Path(args.repo_root) if args.repo_root else None,
        lock_dir=Path(args.lock_dir) if args.lock_dir else None,
        min_free_gb=float(args.min_free_gb) if args.min_free_gb is not None else None,
        max_workers=int(args.max_workers) if args.max_workers is not None else None,
        owner=args.owner or "",
        note=args.note or "",
    )
    for line in result.messages:
        print(line, file=sys.stderr if not result.ok else sys.stdout)
    if result.ok and result.binary_pin:
        print(json.dumps({"binary_pin": result.binary_pin}, indent=2))
    return result.exit_code


def _cmd_release(args: argparse.Namespace) -> int:
    ok, msg = release_lock(
        Path(args.lock_dir) if args.lock_dir else None,
        force=bool(args.force),
    )
    print(msg, file=sys.stdout if ok else sys.stderr)
    return EXIT_OK if ok else EXIT_LOCK


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    h = sub.add_parser("check-hold", help="Exit 78 if BENCHMARK_HOLD is present")
    h.add_argument("--repo-root", default="")
    h.set_defaults(func=_cmd_check_hold)

    pf = sub.add_parser("preflight", help="Hold + lock + disk + workers + binary copy")
    pf.add_argument("--out-dir", required=True)
    pf.add_argument("--binary", default="")
    pf.add_argument("--workers", type=int, default=1)
    pf.add_argument("--repo-root", default="")
    pf.add_argument("--lock-dir", default="")
    pf.add_argument("--min-free-gb", type=float, default=None)
    pf.add_argument("--max-workers", type=int, default=None)
    pf.add_argument("--no-lock", action="store_true")
    pf.add_argument("--no-copy-binary", action="store_true")
    pf.add_argument("--owner", default="")
    pf.add_argument("--note", default="")
    pf.set_defaults(func=_cmd_preflight)

    rl = sub.add_parser("release-lock", help="Release mkdir dock lock")
    rl.add_argument("--lock-dir", default="")
    rl.add_argument("--force", action="store_true")
    rl.set_defaults(func=_cmd_release)
    return p


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    sys.exit(main())
