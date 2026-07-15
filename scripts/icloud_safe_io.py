#!/usr/bin/env python3
"""icloud_safe_io — production I/O helpers that never hang agents on CloudDocs.

Problem
-------
macOS iCloud Drive (``Mobile Documents/com~apple~CloudDocs``) is a FileProvider
filesystem. ``open()``, ``stat()``, ``Path.rglob()``, and ``find`` can block
indefinitely waiting for download / coordination. Agent sessions then pile up
hung processes and starve real docking work.

Policy
------
1. **Live work** lives under ``~/flexaidds_results`` (local APFS).
2. **iCloud** is a durable *mirror* only (thin ``result.csv`` + receipts).
3. Any CloudDocs access is **timeout-bounded** and preferably **materialized**
   into a local pin-cache before hashing/reading.
4. Never ``rglob`` / ``find`` trees under CloudDocs from ops/monitor code.

Usage
-----
::

    from icloud_safe_io import is_clouddocs, safe_md5, materialize, safe_glob_result_csvs

    p = materialize(Path("…/MC_st0r5.2_6.dat"))  # local pin path
    print(safe_md5(p))

Copyright 2026 Le Bonhomme Pharma
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import hashlib
import os
import shutil
import signal
import sys
from concurrent.futures import ProcessPoolExecutor, TimeoutError as FuturesTimeout
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Union

PathLike = Union[str, Path]

CLOUDDOCS_MARKERS = (
    "Mobile Documents/com~apple~CloudDocs",
    "Mobile Documents/com~apple~CloudDocs/",
)


def local_root() -> Path:
    env = os.environ.get("FLEXAIDDS_LOCAL_ROOT", "").strip()
    if env:
        return Path(env).expanduser()
    return Path.home() / "flexaidds_results"


def pin_cache_dir() -> Path:
    env = os.environ.get("FLEXAIDDS_PIN_CACHE", "").strip()
    d = Path(env).expanduser() if env else (local_root() / "pins" / "materialize")
    d.mkdir(parents=True, exist_ok=True)
    return d


def is_clouddocs(path: PathLike) -> bool:
    """True if path is under iCloud Drive CloudDocs (FileProvider)."""
    try:
        s = str(Path(path).expanduser())
    except Exception:
        s = str(path)
    # Do not resolve() — resolve can itself hang on CloudDocs.
    return any(m in s for m in CLOUDDOCS_MARKERS)


def is_local_apfs(path: PathLike) -> bool:
    return not is_clouddocs(path)


def _worker_read(path_s: str, max_bytes: int) -> bytes:
    with open(path_s, "rb") as f:
        if max_bytes <= 0:
            return f.read()
        return f.read(max_bytes)


def _worker_md5(path_s: str) -> str:
    h = hashlib.md5()
    with open(path_s, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _worker_sha256(path_s: str) -> str:
    h = hashlib.sha256()
    with open(path_s, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _worker_exists(path_s: str) -> bool:
    return os.path.isfile(path_s)


def _run_isolated(fn, arg: str, timeout_s: float, default=None):
    """Run *fn(arg)* in a child process; kill on timeout (CloudDocs safe)."""
    # ProcessPoolExecutor is more reliable than threads for stuck syscalls.
    try:
        with ProcessPoolExecutor(max_workers=1) as ex:
            fut = ex.submit(fn, arg)
            return fut.result(timeout=timeout_s)
    except FuturesTimeout:
        return default
    except Exception:
        return default


def safe_exists(path: PathLike, timeout_s: float = 5.0) -> bool:
    p = Path(path).expanduser()
    if not is_clouddocs(p):
        try:
            return p.is_file() or p.is_dir()
        except OSError:
            return False
    return bool(_run_isolated(_worker_exists, str(p), timeout_s, default=False))


def safe_read_bytes(
    path: PathLike,
    *,
    timeout_s: float = 15.0,
    max_bytes: int = 0,
) -> Optional[bytes]:
    """Read file bytes with a hard timeout. Returns None on hang/error."""
    p = Path(path).expanduser()
    if not is_clouddocs(p):
        try:
            with p.open("rb") as f:
                return f.read() if max_bytes <= 0 else f.read(max_bytes)
        except OSError:
            return None
    return _run_isolated(
        lambda s: _worker_read(s, max_bytes),
        str(p),
        timeout_s,
        default=None,
    )


def safe_md5(path: PathLike, timeout_s: float = 20.0) -> Optional[str]:
    """MD5 hex digest with hard timeout. Prefer materialize() first for CloudDocs."""
    p = Path(path).expanduser()
    if not is_clouddocs(p):
        try:
            return _worker_md5(str(p))
        except OSError:
            return None
    return _run_isolated(_worker_md5, str(p), timeout_s, default=None)


def safe_sha256(path: PathLike, timeout_s: float = 30.0) -> Optional[str]:
    p = Path(path).expanduser()
    if not is_clouddocs(p):
        try:
            return _worker_sha256(str(p))
        except OSError:
            return None
    return _run_isolated(_worker_sha256, str(p), timeout_s, default=None)


def materialize(
    path: PathLike,
    *,
    timeout_s: float = 30.0,
    force: bool = False,
) -> Optional[Path]:
    """Copy CloudDocs file into local pin-cache; return local path.

    Local paths are returned unchanged. On timeout/failure returns None.
    """
    p = Path(path).expanduser()
    if not is_clouddocs(p):
        return p if p.exists() else None

    # Stable cache name from absolute string hash (avoid re-reading for name).
    key = hashlib.sha256(str(p).encode()).hexdigest()[:16]
    name = p.name
    dest = pin_cache_dir() / f"{key}_{name}"
    if dest.is_file() and not force:
        return dest

    data = safe_read_bytes(p, timeout_s=timeout_s, max_bytes=0)
    if data is None:
        return None
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    try:
        tmp.write_bytes(data)
        tmp.replace(dest)
    except OSError:
        try:
            tmp.unlink(missing_ok=True)  # type: ignore[arg-type]
        except TypeError:
            if tmp.exists():
                tmp.unlink()
        return None
    # sidecar provenance
    try:
        (dest.with_suffix(dest.suffix + ".src.txt")).write_text(str(p) + "\n")
    except OSError:
        pass
    return dest


def safe_glob_result_csvs(
    campaign_dir: PathLike,
    *,
    timeout_s: float = 20.0,
) -> List[Path]:
    """List ``*/result.csv`` under a campaign dir without recursive rglob.

    Uses a one-level glob only. For CloudDocs, runs glob in a child process
    with timeout so a wedged FileProvider cannot stick the parent forever.
    """
    root = Path(campaign_dir).expanduser()

    def _glob(path_s: str) -> List[str]:
        r = Path(path_s)
        out: List[str] = []
        try:
            for rc in sorted(r.glob("*/result.csv")):
                # Skip incomplete archive dirs
                if "incomplete" in rc.parent.name:
                    continue
                out.append(str(rc))
        except OSError:
            pass
        return out

    if not is_clouddocs(root):
        return [Path(s) for s in _glob(str(root))]

    res = _run_isolated(_glob, str(root), timeout_s, default=None)
    if res is None:
        return []
    return [Path(s) for s in res]


def prefer_local_campaign(
    rel: str,
    *,
    local_campaigns: Optional[Path] = None,
    icloud_results: Optional[Path] = None,
) -> Path:
    """Prefer ``~/flexaidds_results/campaigns/<rel>`` over iCloud results.

    *rel* is like ``campaigns/C0_full85_…`` or just the campaign id basename.
    """
    rel_p = Path(rel)
    name = rel_p.name if rel_p.parts[-1] else str(rel_p)
    # allow full relative campaigns/foo
    local_base = local_campaigns or (local_root() / "campaigns")
    candidates = [
        local_base / name,
        local_base / rel_p,
    ]
    if rel_p.parts and rel_p.parts[0] == "campaigns":
        candidates.insert(0, local_root() / rel_p)

    for c in candidates:
        try:
            if c.is_dir() and any(c.glob("*/result.csv")):
                return c
        except OSError:
            continue
        try:
            if c.is_dir():
                return c  # empty but local is still preferred for live runs
        except OSError:
            continue

    if icloud_results is not None:
        ic = icloud_results / rel if not str(rel).startswith("/") else Path(rel)
        return ic
    # last resort: local path even if missing
    return local_base / name


def local_log_dir() -> Path:
    d = local_root() / "logs" / "ops"
    d.mkdir(parents=True, exist_ok=True)
    return d


# ─── CLI ──────────────────────────────────────────────────────────────────────


def _cli(argv: Sequence[str]) -> int:
    if len(argv) < 2:
        print(
            "Usage: icloud_safe_io.py md5|sha256|materialize|is-cloud <path>",
            file=sys.stderr,
        )
        return 2
    cmd, path = argv[0], argv[1]
    if cmd == "is-cloud":
        print("yes" if is_clouddocs(path) else "no")
        return 0
    if cmd == "md5":
        # materialize first if cloud
        p = materialize(path) if is_clouddocs(path) else Path(path)
        if p is None:
            print("TIMEOUT_OR_ERROR", file=sys.stderr)
            return 1
        dig = safe_md5(p)
        print(dig or "ERROR")
        return 0 if dig else 1
    if cmd == "sha256":
        p = materialize(path) if is_clouddocs(path) else Path(path)
        if p is None:
            print("TIMEOUT_OR_ERROR", file=sys.stderr)
            return 1
        dig = safe_sha256(p)
        print(dig or "ERROR")
        return 0 if dig else 1
    if cmd == "materialize":
        p = materialize(path)
        if p is None:
            print("TIMEOUT_OR_ERROR", file=sys.stderr)
            return 1
        print(p)
        return 0
    print(f"unknown cmd {cmd}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    # Avoid SIGINT hanging children on Ctrl-C in agent sessions
    try:
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    except Exception:
        pass
    sys.exit(_cli(sys.argv[1:]))
