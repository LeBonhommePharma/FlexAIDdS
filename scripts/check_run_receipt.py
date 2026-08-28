#!/usr/bin/env python3
"""Validate docking campaign RUN_RECEIPT / provenance for scoring auditability.

TWO DIALECTS
  engine    Written by DatasetRunner via LIB/RunReceipt.cpp. Has schema_version
            and protocol_config. Carries git_commit. Documented in
            docs/run-uniformity/RUN_RECEIPT_CONTRACT.md.
  campaign  Written by hand by a campaign driver. Has a string `schema` (e.g.
            "1jd0_ga_wal400_v1") and frozen_utc, and deliberately has no
            git_commit or protocol_config — engine identity is carried by
            binary_sha256 / engine_id instead.

Both are legitimate. Validating one against the other's contract produces a pile
of missing-key errors that send the reader looking for the wrong problem, so the
dialect is detected structurally and each is checked against its own
expectations. Use --require-engine-dialect when the caller specifically needs an
engine receipt; a campaign file then fails with one clear message rather than
several misleading ones.

Usage:
  python3 scripts/check_run_receipt.py ~/flexaidds_results/v_autonomous_...
  python3 scripts/check_run_receipt.py path/to/RUN_RECEIPT.json --require-matrix-9dc9
  python3 scripts/check_run_receipt.py path/ --require-engine-dialect
  python3 scripts/check_run_receipt.py path/ --json-out report.json

Exit: 0 OK, 1 missing/invalid fields or wrong dialect, 2 usage/load failure.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Optional

MATRIX_9DC9 = "9dc93717dfed0698006d88dd6a9627bc"

# Resolution order inside a campaign directory. RUN_RECEIPT.json stays first so
# an engine receipt still wins where both are present; PREREGISTRATION.json is
# the name a campaign pre-registration takes once it is no longer overloading
# the engine's filename.
RECEIPT_NAMES = (
    "RUN_RECEIPT.json",
    "PREREGISTRATION.json",
    "provenance.json",
    "out/RUN_RECEIPT.json",
)

DIALECT_ENGINE = "engine"
DIALECT_CAMPAIGN = "campaign"
DIALECT_UNKNOWN = "unknown"

REQUIRED = ("matrix_md5", "git_commit")
CAMPAIGN_REQUIRED = ("schema", "frozen_utc", "matrix_md5")
BINARY_KEYS = ("binary_sha256", "binary_path")
SCORING_ENV_KEYS = (
    "FLEXAIDDS_ACF_STRICT",
    "FLEXAIDDS_COM_BURIAL_CAP",
    "FLEXAIDDS_COM_FLOOR",
    "FLEXAIDDS_VCT_NORM",
    "FLEXAIDDS_SOFTBETA_ELECTION",
    "FLEXAIDDS_ELECTION_ENTROPY",
)


def detect_dialect(data: Any) -> str:
    """Classify a parsed receipt structurally, not by filename.

    Structural rather than by the literal schema string: the next campaign will
    invent a new one, and a checker that has to be edited per campaign is a
    checker that will be wrong at exactly the wrong moment.

    UNKNOWN routes to the engine checks, which is the historical default and
    keeps every pre-existing caller behaving as before.
    """
    if not isinstance(data, dict):
        return DIALECT_UNKNOWN
    if "schema_version" in data:
        return DIALECT_ENGINE
    if isinstance(data.get("schema"), str) and data["schema"].strip():
        return DIALECT_CAMPAIGN
    if "protocol_config" in data:
        # legacy slim provenance.json — engine family, no schema_version
        return DIALECT_ENGINE
    return DIALECT_UNKNOWN


def resolve_receipt_path(path: Path) -> Path:
    """Return the file a directory argument resolves to. Raises if none exists."""
    if path.is_file():
        return path
    if path.is_dir():
        for name in RECEIPT_NAMES:
            p = path / name
            if p.is_file():
                return p
    raise FileNotFoundError(
        f"no {' / '.join(RECEIPT_NAMES)} under {path}"
    )


def load_receipt_with_source(path: Path) -> tuple[dict[str, Any], Path]:
    src = resolve_receipt_path(path)
    return json.loads(src.read_text(encoding="utf-8", errors="replace")), src


def load_receipt(path: Path) -> dict[str, Any]:
    return load_receipt_with_source(path)[0]


def _check_engine(
    data: dict[str, Any],
    *,
    require_matrix_9dc9: bool,
    require_scoring_env: bool,
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


def _check_campaign(data: dict[str, Any], *, require_matrix_9dc9: bool) -> list[str]:
    """Campaign dialect checked against its own expectations.

    Deliberately does NOT require git_commit: the campaign dialect has no such
    key by design. Requiring it here is the exact miscategorisation this split
    exists to prevent.
    """
    errs: list[str] = []
    for k in CAMPAIGN_REQUIRED:
        v = data.get(k)
        if v is None or str(v).strip() == "":
            errs.append(f"missing_or_empty:{k}")
    if not any(str(data.get(k, "")).strip() for k in BINARY_KEYS):
        errs.append("missing_or_empty:binary_sha256_or_binary_path")
    md = str(data.get("matrix_md5", "")).strip().lower()
    if require_matrix_9dc9 and md != MATRIX_9DC9:
        errs.append(f"matrix_md5_not_9dc9:got={md or 'empty'}")
    return errs


def check_receipt(
    data: dict[str, Any],
    *,
    require_matrix_9dc9: bool = False,
    require_scoring_env: bool = False,
    dialect: Optional[str] = None,
) -> list[str]:
    resolved = dialect or detect_dialect(data)
    if resolved == DIALECT_CAMPAIGN:
        return _check_campaign(data, require_matrix_9dc9=require_matrix_9dc9)
    return _check_engine(
        data,
        require_matrix_9dc9=require_matrix_9dc9,
        require_scoring_env=require_scoring_env,
    )


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("path", type=Path, help="Campaign dir or receipt .json")
    ap.add_argument("--require-matrix-9dc9", action="store_true")
    ap.add_argument(
        "--require-scoring-env",
        action="store_true",
        help="Require scoring_env object with known FLEXAIDDS_* keys",
    )
    ap.add_argument(
        "--require-engine-dialect",
        action="store_true",
        help="Fail if the file is a campaign pre-registration rather than an engine receipt",
    )
    ap.add_argument("--json-out", type=Path, default=None)
    args = ap.parse_args(argv)

    try:
        data, source = load_receipt_with_source(args.path.expanduser())
    except (OSError, ValueError, FileNotFoundError, json.JSONDecodeError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    dialect = detect_dialect(data)
    schema_label = (
        str(data.get("schema", "")) if dialect == DIALECT_CAMPAIGN
        else str(data.get("schema_version", ""))
    )

    wrong_dialect = args.require_engine_dialect and dialect == DIALECT_CAMPAIGN
    if wrong_dialect:
        errs = [f"campaign_dialect_not_engine_dialect:schema={schema_label or 'unnamed'}"]
    else:
        errs = check_receipt(
            data,
            require_matrix_9dc9=args.require_matrix_9dc9,
            require_scoring_env=args.require_scoring_env,
            dialect=dialect,
        )

    report = {
        "path": str(args.path),
        "source": str(source),
        "dialect": dialect,
        "schema": schema_label,
        "ok": len(errs) == 0,
        "errors": errs,
        "matrix_md5": data.get("matrix_md5"),
        "git_commit": data.get("git_commit"),
        "has_scoring_env": isinstance(data.get("scoring_env"), dict),
    }
    if args.json_out:
        args.json_out.expanduser().parent.mkdir(parents=True, exist_ok=True)
        args.json_out.expanduser().write_text(json.dumps(report, indent=2) + "\n")

    if wrong_dialect:
        print(
            f"FAIL: campaign dialect, not engine dialect "
            f"(schema={schema_label or 'unnamed'}, source={source.name}). "
            f"This file is a campaign pre-registration; it is valid, but it is "
            f"not an engine receipt. Drop --require-engine-dialect to validate "
            f"it against the campaign contract."
        )
        return 1

    if errs:
        print(f"FAIL [{dialect} dialect]: " + ", ".join(errs))
        return 1

    if dialect == DIALECT_CAMPAIGN:
        print(
            f"OK [campaign dialect]: schema={schema_label} "
            f"matrix_md5={data.get('matrix_md5')} "
            f"binary_sha256={str(data.get('binary_sha256'))[:12]}"
        )
        return 0

    print(
        f"OK: matrix_md5={data.get('matrix_md5')} git_commit={str(data.get('git_commit'))[:12]}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
