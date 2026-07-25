#!/usr/bin/env python3
"""Validate docking campaign RUN_RECEIPT / provenance for scoring auditability.

Usage:
  python3 scripts/check_run_receipt.py ~/flexaidds_results/v_autonomous_...
  python3 scripts/check_run_receipt.py path/to/RUN_RECEIPT.json --require-matrix-9dc9
  python3 scripts/check_run_receipt.py path/ --json-out report.json

Exit: 0 OK, 1 missing/invalid fields, 2 usage.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Optional

MATRIX_9DC9 = "9dc93717dfed0698006d88dd6a9627bc"

REQUIRED = ("matrix_md5", "git_commit")
BINARY_KEYS = ("binary_sha256", "binary_path")
SCORING_ENV_KEYS = (
    "FLEXAIDDS_ACF_STRICT",
    "FLEXAIDDS_COM_BURIAL_CAP",
    "FLEXAIDDS_COM_FLOOR",
    "FLEXAIDDS_VCT_NORM",
    "FLEXAIDDS_SOFTBETA_ELECTION",
    "FLEXAIDDS_ELECTION_ENTROPY",
)


def load_receipt(path: Path) -> dict[str, Any]:
    if path.is_file() and path.suffix == ".json":
        return json.loads(path.read_text(encoding="utf-8", errors="replace"))
    if path.is_dir():
        for name in ("RUN_RECEIPT.json", "provenance.json", "out/RUN_RECEIPT.json"):
            p = path / name if not name.startswith("out/") else path / "out" / "RUN_RECEIPT.json"
            if name == "out/RUN_RECEIPT.json":
                p = path / "out" / "RUN_RECEIPT.json"
            else:
                p = path / name
            if p.is_file():
                return json.loads(p.read_text(encoding="utf-8", errors="replace"))
    raise FileNotFoundError(f"no RUN_RECEIPT.json / provenance.json under {path}")


def check_receipt(
    data: dict[str, Any],
    *,
    require_matrix_9dc9: bool = False,
    require_scoring_env: bool = False,
) -> list[str]:
    errs: list[str] = []
    for k in REQUIRED:
        v = data.get(k)
        if v is None or str(v).strip() == "":
            errs.append(f"missing_or_empty:{k}")
    if not any(str(data.get(k, "")).strip() for k in BINARY_KEYS):
        errs.append("missing_or_empty:binary_sha256_or_binary_path")
    md = str(data.get("matrix_md5", "")).strip().lower()
    if require_matrix_9dc9:
        if md != MATRIX_9DC9:
            errs.append(f"matrix_md5_not_9dc9:got={md or 'empty'}")
    elif md and len(md) not in (0, 32) and not md.startswith("9dc9") and not md.startswith("72d7"):
        # soft warn only for weird lengths — not hard error
        pass

    env = data.get("scoring_env")
    if require_scoring_env:
        if not isinstance(env, dict) or not env:
            errs.append("missing_or_empty:scoring_env")
        else:
            # at least one known scoring knob recorded (even if value is "unset")
            if not any(k in env for k in SCORING_ENV_KEYS):
                errs.append("scoring_env_missing_known_keys")
    return errs


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("path", type=Path, help="Campaign dir or RUN_RECEIPT.json")
    ap.add_argument("--require-matrix-9dc9", action="store_true")
    ap.add_argument(
        "--require-scoring-env",
        action="store_true",
        help="Require scoring_env object with known FLEXAIDDS_* keys",
    )
    ap.add_argument("--json-out", type=Path, default=None)
    args = ap.parse_args(argv)

    try:
        data = load_receipt(args.path.expanduser())
    except (OSError, ValueError, FileNotFoundError, json.JSONDecodeError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    errs = check_receipt(
        data,
        require_matrix_9dc9=args.require_matrix_9dc9,
        require_scoring_env=args.require_scoring_env,
    )
    report = {
        "path": str(args.path),
        "ok": len(errs) == 0,
        "errors": errs,
        "matrix_md5": data.get("matrix_md5"),
        "git_commit": data.get("git_commit"),
        "has_scoring_env": isinstance(data.get("scoring_env"), dict),
    }
    if args.json_out:
        args.json_out.expanduser().parent.mkdir(parents=True, exist_ok=True)
        args.json_out.expanduser().write_text(json.dumps(report, indent=2) + "\n")

    if errs:
        print("FAIL: " + ", ".join(errs))
        return 1
    print(
        f"OK: matrix_md5={data.get('matrix_md5')} git_commit={str(data.get('git_commit'))[:12]}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
