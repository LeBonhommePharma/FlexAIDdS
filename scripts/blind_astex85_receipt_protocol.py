#!/usr/bin/env python3
"""Receipt-gated blind Astex-85 republish protocol (Wave 4).

This is the protocol, not the 85-target run. It can emit a RUN_RECEIPT for a
blind campaign and must refuse to print a success percentage without a receipt.

Pins (METHODOLOGY.md §0 / §3; Wave 4 handoff):
  - N = 85
  - native_pose_seeded = 0
  - seed_echo = 0
  - matrix MD5 = 72d7c7396702331d96ff12d18f831796 (MC_st0r5.2_6.dat)
  - default SEED_ELITISM = 0, NATIVE_SEED_FRAC = 0
  - claim success = rank-0 in-place RMSD <= 2.0 Å

Do not launch docking from this script. Use --dry-run / emit / claim only.

Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Optional

ROOT = Path(__file__).resolve().parents[1]

N_TARGETS = 85
MATRIX_MD5_PIN = "72d7c7396702331d96ff12d18f831796"
MATRIX_NAME = "MC_st0r5.2_6.dat"
CLAIM_CUTOFF_A = 2.0
SCHEMA_VERSION = 1

# Tokens that must never appear as a published docking-power default.
ORACLE_SEED_ELITISM = "1"
ORACLE_NATIVE_SEED_FRAC = "0.90"


class ProtocolError(RuntimeError):
    pass


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def default_blind_receipt(**overrides: Any) -> dict[str, Any]:
    rec: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "protocol": "blind_astex85_receipt",
        "methodology": "METHODOLOGY.md §0 / §3",
        "dataset": "astex_diverse",
        "n_targets": N_TARGETS,
        "native_pose_seeded": 0,
        "seed_echo": 0,
        "seed_elitism": 0,
        "native_seed_frac": 0,
        "matrix_name": MATRIX_NAME,
        "matrix_md5": MATRIX_MD5_PIN,
        "arm": "blind",
        "claim_validity": "blind",
        "rmsd_instrument": "rank-0 in-place RMSD <= 2.0 Å",
        "dry_run": False,
        "started_utc": utc_now(),
        "git_commit": "",
        "binary_sha256": "",
        "binary_path": "",
    }
    rec.update(overrides)
    rec["inert_note"] = (
        "Receipt does not authorise a success % until result.csv is present "
        "and validate_blind_receipt() returns no errors."
    )
    return rec


def _as_int(value: Any, default: int = -1) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _as_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def validate_blind_receipt(
    data: Mapping[str, Any],
    *,
    allow_oracle: bool = False,
) -> list[str]:
    errs: list[str] = []
    if _as_int(data.get("n_targets")) != N_TARGETS:
        errs.append(f"n_targets_not_85:got={data.get('n_targets')}")
    if _as_int(data.get("native_pose_seeded"), 1) != 0:
        errs.append("native_pose_seeded_not_0")
    if _as_int(data.get("seed_echo"), 1) != 0:
        errs.append("seed_echo_not_0")
    md = str(data.get("matrix_md5", "")).strip().lower()
    if md != MATRIX_MD5_PIN:
        errs.append(f"matrix_md5_not_pin:got={md or 'empty'}")
    seed_elitism = _as_int(data.get("seed_elitism"), 1)
    native_frac = str(data.get("native_seed_frac", "")).strip()
    if seed_elitism == 1 or native_frac in {ORACLE_NATIVE_SEED_FRAC, "0.9"}:
        if not allow_oracle:
            errs.append("oracle_seed_not_allowed_on_blind_claim")
    git = str(data.get("git_commit", "")).strip()
    if not git:
        errs.append("missing_or_empty:git_commit")
    binary_ok = any(
        str(data.get(k, "")).strip() for k in ("binary_sha256", "binary_path")
    )
    if not binary_ok:
        errs.append("missing_or_empty:binary_sha256_or_binary_path")
    return errs


def load_receipt(path: Path) -> dict[str, Any]:
    if path.is_file():
        return json.loads(path.read_text(encoding="utf-8"))
    for name in ("RUN_RECEIPT.json", "provenance.json"):
        cand = path / name
        if cand.is_file():
            return json.loads(cand.read_text(encoding="utf-8"))
    raise FileNotFoundError(f"no RUN_RECEIPT.json under {path}")


def find_results_csv(campaign_dir: Path) -> Optional[Path]:
    for name in (
        "astex_crossdock_85_results.csv",
        "result.csv",
        "astex_crossdock_85_summary.csv",
    ):
        cand = campaign_dir / name
        if cand.is_file():
            return cand
    return None


def pick_rank0_rmsd(row: Mapping[str, str]) -> tuple[float, str]:
    for key in ("rmsd_hungarian", "rmsd_top1", "rmsd_to_crystal"):
        raw = (row.get(key) or "").strip()
        if raw:
            try:
                return float(raw), key
            except ValueError:
                continue
    return 9999.0, "missing"


def s1_from_csv(csv_path: Path) -> dict[str, Any]:
    rows = list(csv.DictReader(csv_path.open(encoding="utf-8")))
    picked = [pick_rank0_rmsd(r) for r in rows]
    n_ok = sum(1 for value, name in picked if name != "missing" and value <= CLAIM_CUTOFF_A)
    return {
        "n_rows": len(rows),
        "n_ok": n_ok,
        "n_denominator": N_TARGETS,
        "s1_percent": 100.0 * n_ok / N_TARGETS,
        "instrument": ",".join(sorted({name for _, name in picked})),
    }


def emit_receipt(out_dir: Path, **overrides: Any) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    rec = default_blind_receipt(**overrides)
    path = out_dir / "RUN_RECEIPT.json"
    path.write_text(json.dumps(rec, indent=2) + "\n", encoding="utf-8")
    return path


def format_s1_line(stats: Mapping[str, Any]) -> str:
    return (
        f"Observed S1 (rank-0 in-place RMSD <= {CLAIM_CUTOFF_A:.1f} Å) = "
        f"{stats['n_ok']}/{stats['n_denominator']} "
        f"({stats['s1_percent']:.1f}%)  instrument={stats['instrument']}"
    )


def claim_from_dir(
    campaign_dir: Path,
    *,
    allow_oracle: bool = False,
) -> tuple[int, str]:
    """Return (exit_code, message). Never prints a % when refusing."""
    try:
        receipt = load_receipt(campaign_dir)
    except (OSError, ValueError, FileNotFoundError, json.JSONDecodeError) as exc:
        return 1, f"REFUSE: no RUN_RECEIPT ({exc}). Not printing a success %."

    errs = validate_blind_receipt(receipt, allow_oracle=allow_oracle)
    if errs:
        return 1, "REFUSE: receipt failed blind Astex-85 gates: " + ", ".join(errs)

    csv_path = find_results_csv(campaign_dir)
    if csv_path is None:
        return 2, (
            "Receipt OK for a blind Astex-85 protocol, but no result.csv. "
            "Not printing a success %."
        )
    stats = s1_from_csv(csv_path)
    if stats["n_rows"] != N_TARGETS:
        return 1, (
            f"REFUSE: CSV has {stats['n_rows']} rows, need {N_TARGETS}. "
            "Not printing a success %."
        )
    return 0, format_s1_line(stats)


def refuse_oracle_defaults(seed_elitism: Any, native_seed_frac: Any) -> None:
    if str(seed_elitism).strip() == ORACLE_SEED_ELITISM and str(native_seed_frac).strip() == ORACLE_NATIVE_SEED_FRAC:
        raise ProtocolError(
            "REFUSE: SEED_ELITISM=1 / NATIVE_SEED_FRAC=0.90 is the oracle ceiling, "
            "not the blind default. Pass --oracle-ceiling explicitly; do not cite as S1."
        )


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    emit = sub.add_parser("emit", help="Write RUN_RECEIPT.json (no docking)")
    emit.add_argument("--out", type=Path, required=True)
    emit.add_argument("--git-commit", default="unknown")
    emit.add_argument("--binary-path", default="unspecified")
    emit.add_argument("--binary-sha256", default="")
    emit.add_argument("--dry-run", action="store_true")
    emit.add_argument("--oracle-ceiling", action="store_true")

    claim = sub.add_parser("claim", help="Print S1 only if receipt+CSV pass")
    claim.add_argument("--dir", type=Path, required=True)
    claim.add_argument("--allow-oracle", action="store_true")

    sub.add_parser("validate-defaults", help="Assert the blind default is not oracle")
    return ap


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.cmd == "validate-defaults":
        rec = default_blind_receipt()
        refuse_oracle_defaults(rec["seed_elitism"], rec["native_seed_frac"])
        errs = validate_blind_receipt(
            {**rec, "git_commit": "deadbeef", "binary_path": "/bin/FlexAIDdS"}
        )
        if errs:
            print("REFUSE: " + ", ".join(errs), file=sys.stderr)
            return 1
        print(
            "OK: blind default n=85 seed_elitism=0 native_seed_frac=0 "
            f"native_pose_seeded=0 seed_echo=0 matrix_md5={MATRIX_MD5_PIN}"
        )
        return 0

    if args.cmd == "emit":
        overrides: dict[str, Any] = {
            "git_commit": args.git_commit,
            "binary_path": args.binary_path,
            "binary_sha256": args.binary_sha256,
            "dry_run": bool(args.dry_run),
        }
        if args.oracle_ceiling:
            overrides.update(
                {
                    "seed_elitism": 1,
                    "native_seed_frac": 0.90,
                    "arm": "oracle_ceiling",
                    "claim_validity": "oracle_ceiling_not_docking_power",
                    "native_pose_seeded": 1,
                }
            )
            print("ORACLE CEILING — not docking power", file=sys.stderr)
        else:
            refuse_oracle_defaults(0, 0)
        path = emit_receipt(args.out.expanduser(), **overrides)
        print(f"wrote {path}")
        if args.dry_run:
            print("dry-run: receipt only; not launching an 85-target dock")
        return 0

    if args.cmd == "claim":
        code, msg = claim_from_dir(args.dir.expanduser(), allow_oracle=args.allow_oracle)
        stream = sys.stdout if code == 0 else sys.stderr
        print(msg, file=stream)
        return 0 if code == 0 else 1

    return 2


if __name__ == "__main__":
    raise SystemExit(main())
