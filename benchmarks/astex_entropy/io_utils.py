from __future__ import annotations

import concurrent.futures
import csv
import hashlib
import io
import json
import os
import shutil
import subprocess
from pathlib import Path
from typing import Callable, Iterable, Sequence, TypeVar

from .models import PoseRecord, TargetRecord

T = TypeVar("T")

_LOG_CACHE_ROOT = Path(os.environ.get("ASTEX_ENTROPY_LOG_CACHE", "/tmp/astex_entropy_logs"))
_SCRATCH_ROOT = Path(os.environ.get("ASTEX_ENTROPY_LOCAL_SCRATCH", "/tmp/astex_entropy_scratch"))
_FS_TIMEOUT = float(os.environ.get("ASTEX_ENTROPY_FS_TIMEOUT_SEC", "30"))


def _run_fs(func: Callable[..., T], /, *args, timeout: float | None = None) -> T:
    limit = _FS_TIMEOUT if timeout is None else timeout
    with concurrent.futures.ThreadPoolExecutor(max_workers=1) as pool:
        future = pool.submit(func, *args)
        try:
            return future.result(timeout=limit)
        except concurrent.futures.TimeoutError as exc:
            raise TimeoutError(f"filesystem operation timed out after {limit}s") from exc


def _is_icloud_path(path: Path) -> bool:
    text = str(path)
    return "Mobile Documents" in text or "com~apple~CloudDocs" in text


def _local_scratch_root() -> Path:
    return _SCRATCH_ROOT


def _materialize_local(path: Path, *, timeout: float = 120) -> Path:
    """Stage slow iCloud inputs on local disk for rescoring hot paths."""
    if not _is_icloud_path(path):
        return path
    digest = hashlib.sha256(str(path).encode()).hexdigest()[:16]
    dest = _local_scratch_root() / "inputs" / digest / path.name
    if dest.exists() and dest.stat().st_size > 0:
        try:
            src_st = path.stat()
            dst_st = dest.stat()
            if dst_st.st_size == src_st.st_size and dst_st.st_mtime >= src_st.st_mtime:
                return dest
        except OSError:
            pass
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        _run_fs(shutil.copy2, path, dest, timeout=timeout)
        return dest
    except (TimeoutError, OSError):
        pass
    try:
        content = _run_fs(lambda: path.read_text(errors="ignore"), timeout=timeout)
        if content:
            dest.write_text(content)
            return dest
    except (TimeoutError, OSError):
        pass
    return dest if dest.exists() and dest.stat().st_size > 0 else path


def _mirror_path(path: Path) -> Path:
    digest = hashlib.sha256(str(path).encode()).hexdigest()[:16]
    return _LOG_CACHE_ROOT / "mirrors" / digest / path.name


def _safe_exists(path: Path) -> bool:
    try:
        return bool(_run_fs(path.exists))
    except (TimeoutError, OSError):
        return False


def _safe_stat_size(path: Path) -> int:
    try:
        return int(_run_fs(lambda: path.stat().st_size))
    except (TimeoutError, OSError):
        return -1


def _safe_mtime(path: Path) -> float:
    try:
        return float(_run_fs(lambda: path.stat().st_mtime))
    except (TimeoutError, OSError):
        return 0.0


def _safe_unlink(path: Path) -> None:
    try:
        _run_fs(path.unlink)
    except (TimeoutError, OSError, FileNotFoundError):
        pass


def _safe_read_text(path: Path, *, default: str = "") -> str:
    for candidate in (path, _mirror_path(path)):
        try:
            if _safe_exists(candidate):
                return _run_fs(lambda: candidate.read_text(errors="ignore"))
        except (TimeoutError, OSError):
            continue
    return default


def _safe_write_text(path: Path, content: str) -> Path:
    """Write log text locally first; mirror to target path when iCloud allows."""
    mirror = _mirror_path(path)
    try:
        mirror.parent.mkdir(parents=True, exist_ok=True)
        mirror.write_text(content)
    except (TimeoutError, OSError):
        pass
    if _is_icloud_path(path):
        return mirror if mirror.exists() else path
    try:
        _run_fs(lambda: (path.parent.mkdir(parents=True, exist_ok=True), path.write_text(content)))
        return path
    except (TimeoutError, OSError):
        return mirror if mirror.exists() else path


def _safe_read_csv(path: Path):
    import pandas as pd

    for candidate in (path, _mirror_path(path)):
        try:
            if _safe_exists(candidate) and _safe_stat_size(candidate) > 0:
                return _run_fs(pd.read_csv, candidate)
        except (TimeoutError, OSError):
            continue
    return None


def _safe_write_csv(df, path: Path) -> Path:
    mirror = _mirror_path(path)
    try:
        mirror.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(mirror, index=False, quoting=csv.QUOTE_MINIMAL)
    except (TimeoutError, OSError):
        pass
    if _is_icloud_path(path):
        return mirror if mirror.exists() else path
    try:
        _run_fs(lambda: (path.parent.mkdir(parents=True, exist_ok=True), df.to_csv(path, index=False, quoting=csv.QUOTE_MINIMAL)))
        return path
    except (TimeoutError, OSError):
        return mirror if mirror.exists() else path


def read_targets(path: str | Path) -> list[TargetRecord]:
    with Path(path).open(newline="") as fh:
        return [TargetRecord.from_row(row) for row in csv.DictReader(fh)]


def write_targets(records: Iterable[TargetRecord], csv_path: str | Path) -> Path:
    rows = [record.to_dict() for record in records]
    path = Path(csv_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(TargetRecord.__dataclass_fields__.keys())
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    path.with_suffix(".json").write_text(json.dumps(rows, indent=2) + "\n")
    return path


def read_poses(path: str | Path) -> list[PoseRecord]:
    local_path = _materialize_local(Path(path))
    text = _safe_read_text(local_path)
    if not text.strip() and local_path != Path(path):
        text = _safe_read_text(Path(path))
    if not text.strip():
        return []
    return [PoseRecord.from_row(row) for row in csv.DictReader(io.StringIO(text))]


def write_poses(records: Iterable[PoseRecord], csv_path: str | Path) -> Path:
    rows = [record.to_dict() for record in records]
    path = Path(csv_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(PoseRecord.__dataclass_fields__.keys())
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return path


def append_poses(records: Iterable[PoseRecord], csv_path: str | Path) -> Path:
    existing = read_poses(csv_path)
    existing.extend(records)
    return write_poses(existing, csv_path)


def run_command(
    args: Sequence[str],
    *,
    cwd: str | Path | None = None,
    log_path: str | Path | None = None,
    check: bool = True,
    timeout: int | None = 7200,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    run_env = os.environ.copy()
    if env:
        run_env.update(env)
    cmd_args = [str(a) for a in args]
    cmd = " ".join(cmd_args)
    try:
        result = subprocess.run(
            cmd_args,
            cwd=str(cwd) if cwd else None,
            env=run_env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
            timeout=timeout,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(f"Command executable not found: {cmd_args[0]}") from exc
    except subprocess.TimeoutExpired as exc:
        stdout = exc.stdout or ""
        if isinstance(stdout, bytes):
            stdout = stdout.decode(errors="replace")
        if log_path:
            _safe_write_text(Path(log_path), stdout)
        raise RuntimeError(f"Command timed out after {timeout}s: {cmd}\n{stdout[-4000:]}") from exc
    if log_path:
        _safe_write_text(Path(log_path), result.stdout)
    if check and result.returncode != 0:
        raise RuntimeError(f"Command failed ({result.returncode}): {cmd}\n{result.stdout[-4000:]}")
    return result


def first(items: Iterable[T]) -> T | None:
    for item in items:
        return item
    return None
