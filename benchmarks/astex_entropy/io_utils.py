from __future__ import annotations

import csv
import json
import os
import subprocess
from pathlib import Path
from typing import Iterable, Sequence, TypeVar

from .models import PoseRecord, TargetRecord

T = TypeVar("T")


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
    if not Path(path).exists():
        return []
    with Path(path).open(newline="") as fh:
        return [PoseRecord.from_row(row) for row in csv.DictReader(fh)]


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
            Path(log_path).parent.mkdir(parents=True, exist_ok=True)
            Path(log_path).write_text(stdout)
        raise RuntimeError(f"Command timed out after {timeout}s: {cmd}\n{stdout[-4000:]}") from exc
    if log_path:
        Path(log_path).parent.mkdir(parents=True, exist_ok=True)
        Path(log_path).write_text(result.stdout)
    if check and result.returncode != 0:
        raise RuntimeError(f"Command failed ({result.returncode}): {cmd}\n{result.stdout[-4000:]}")
    return result


def first(items: Iterable[T]) -> T | None:
    for item in items:
        return item
    return None
