# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
"""Phase P0: local-first layout + JCIM matrix pin check."""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .gates import MATRIX_MD5_PIN

REQUIRED_SUBDIRS = [
    "campaigns",
    "campaigns/three_engine",
    "campaigns/three_engine/A",
    "campaigns/three_engine/B0",
    "campaigns/three_engine/B",
    "campaigns/three_engine/C",
    "campaigns/three_engine/analysis",
    "campaigns/three_engine/receipts",
    "logs/ops",
    "logs/ops_monitor",
    "pins/materialize",
    "three_engine_entropy_q1/bin/A",
    "three_engine_entropy_q1/bin/B",
    "three_engine_entropy_q1/bin/C",
    "three_engine_entropy_q1/data",
    "three_engine_entropy_q1/inputs",
]


def repo_root() -> Path:
    env = os.environ.get("FLEXAIDDS_ROOT")
    if env:
        return Path(env).resolve()
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
        return Path(out)
    except (subprocess.CalledProcessError, FileNotFoundError):
        # package at python/flexaidds/comparative_phases/ → repo root parents[3]
        return Path(__file__).resolve().parents[3]


def local_root(override: Optional[str] = None) -> Path:
    if override:
        return Path(override).expanduser().resolve()
    return Path(
        os.environ.get("FLEXAIDDS_LOCAL_ROOT", str(Path.home() / "flexaidds_results"))
    ).expanduser().resolve()


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def ensure_matrix(local: Path, root: Path) -> Tuple[Path, str]:
    """Ensure matrix on live path; return (path, md5). Raises FileNotFoundError / ValueError."""
    dest = local / "three_engine_entropy_q1" / "data" / "MC_st0r5.2_6.dat"
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.is_file():
        candidates = [
            root / "MC_st0r5.2_6.dat",
            root / "WRK" / "MC_st0r5.2_6.dat",
            root / ".grok" / "skills" / "flexaidds" / "data" / "MC_st0r5.2_6.dat",
        ]
        src = next((c for c in candidates if c.is_file()), None)
        if src is None:
            raise FileNotFoundError(
                f"matrix missing and no source under {root}; need MC_st0r5.2_6.dat"
            )
        shutil.copy2(src, dest)
    digest = md5_file(dest)
    if digest != MATRIX_MD5_PIN:
        raise ValueError(
            f"matrix MD5 {digest} != pin {MATRIX_MD5_PIN} path={dest}"
        )
    return dest, digest


def run_p0(
    local_root_path: Optional[str] = None,
    *,
    call_shell_layout: bool = True,
) -> Dict[str, Any]:
    """Execute P0: layout dirs + matrix pin. Returns status dict."""
    root = repo_root()
    local = local_root(local_root_path)
    messages: List[str] = []

    if call_shell_layout:
        script = root / "scripts" / "ensure_local_first_layout.sh"
        if script.is_file():
            env = os.environ.copy()
            env["FLEXAIDDS_ROOT"] = str(root)
            env["FLEXAIDDS_LOCAL_ROOT"] = str(local)
            try:
                subprocess.run(
                    ["bash", str(script)],
                    check=True,
                    env=env,
                    capture_output=True,
                    text=True,
                )
                messages.append("ensure_local_first_layout.sh OK")
            except subprocess.CalledProcessError as exc:
                messages.append(
                    f"ensure_local_first_layout.sh warn: {exc.stderr or exc.stdout}"
                )

    created: List[str] = []
    for rel in REQUIRED_SUBDIRS:
        d = local / rel
        d.mkdir(parents=True, exist_ok=True)
        if d.is_dir():
            created.append(rel)

    try:
        mat_path, digest = ensure_matrix(local, root)
        status = "pass"
        reason = f"matrix md5={digest}"
    except (FileNotFoundError, ValueError) as exc:
        status = "fail"
        reason = str(exc)
        mat_path, digest = local / "three_engine_entropy_q1/data/MC_st0r5.2_6.dat", ""

    return {
        "phase": "P0",
        "status": status,
        "reason": reason,
        "local_root": str(local),
        "repo_root": str(root),
        "matrix_path": str(mat_path),
        "matrix_md5": digest,
        "matrix_md5_pin": MATRIX_MD5_PIN,
        "dirs_ok": created,
        "messages": messages,
    }
