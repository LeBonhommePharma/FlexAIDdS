#!/usr/bin/env python3
"""Aggregate claim-table metrics under the admission + S1/S2/S3 contract.

Normative contract:
  benchmarks/protocols/admission_metrics_contract.md
  benchmarks/protocols/three_engine_entropy_comparison.md §1.4–§5

Claim admission (all required):
  seed_echo == 0
  native_pose_seeded == 0
  matrix_md5 == PIN  (default 72d7c7396702331d96ff12d18f831796,
                      or RUN_RECEIPT / provenance / --matrix-md5)

Metrics (always reported separately on claim rows):
  S1  elected RMSD ≤ 2.0 Å          — primary / headline KPI
  S2  S1 ∧ PoseBusters pass         — modern secondary
  S3  best_cluster_rmsd ≤ 2.0 Å     — diagnostic only (BCR / any-pose)

Never treat S3 as abstract success. --headline s3 requires --diagnostic-only
or the process exits nonzero.

Usage:
  python3 scripts/aggregate_claim_metrics.py <campaign_dir> [--json out.json]
  python3 scripts/aggregate_claim_metrics.py --c0-full85
  python3 scripts/aggregate_claim_metrics.py <dir> --headline s1
  python3 scripts/aggregate_claim_metrics.py <dir> --headline s3 --diagnostic-only

Exit codes:
  0  OK (N_claim > 0)
  1  no claim-eligible rows / empty campaign
  2  usage or contract violation
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from pathlib import Path
from typing import Any, Iterable

# Production matrix pin used by three-engine / C0 full85 (see RUN_RECEIPT).
DEFAULT_MATRIX_MD5 = "72d7c7396702331d96ff12d18f831796"
RMSD_SUCCESS_A = 2.0
C0_FULL85_REL = "campaigns/C0_full85_defined_cleft_nativeseed_forbidden"

SUMMARY_CSV_NAMES = (
    "astex_diverse_results.csv",
    "astex_crossdock_85_results.csv",
    "results.csv",
    "summary.csv",
    "claim_summary.csv",
)


# ── helpers ──────────────────────────────────────────────────────────────────


def _f(row: dict[str, str], *keys: str) -> float:
    for k in keys:
        v = row.get(k)
        if v is None or v == "" or str(v).upper() == "NA":
            continue
        try:
            x = float(v)
            if math.isfinite(x):
                return x
        except (TypeError, ValueError):
            continue
    return float("nan")


def _truth(row: dict[str, str], key: str) -> bool:
    return str(row.get(key, "")).strip() in ("1", "True", "true", "YES", "yes")


def _flag0(row: dict[str, str], key: str) -> bool:
    """True when the admission flag is explicitly zero / false / missing-as-zero-ok.

    Missing key is treated as 0 (claim-pass) so older CSVs without the column
    can still be filtered by the other gates. Explicit non-zero fails.
    """
    if key not in row or row.get(key) is None or str(row.get(key, "")).strip() == "":
        return True  # missing → treat as 0 for admission
    return str(row.get(key, "")).strip() in ("0", "False", "false", "NO", "no")


def _pdb_id(row: dict[str, str]) -> str:
    return str(row.get("pdb_id") or row.get("pdb") or row.get("target") or "?").strip()


def resolve_c0_full85_dir() -> Path:
    """Resolve C0 full85 campaign dir from FLEXAIDDS_RESULTS (or FLEXAIDDS_ICLOUD)."""
    results = os.environ.get("FLEXAIDDS_RESULTS", "").strip()
    if results:
        return Path(results) / "campaigns" / "C0_full85_defined_cleft_nativeseed_forbidden"
    icloud = os.environ.get("FLEXAIDDS_ICLOUD", "").strip()
    if icloud:
        # FLEXAIDDS_ICLOUD may be the benchmarks root or the nested astex_entropy path
        base = Path(icloud)
        candidates = [
            base / "results" / C0_FULL85_REL,
            base.parent / "results" / C0_FULL85_REL if base.name else None,
        ]
        for c in candidates:
            if c is not None and c.is_dir():
                return c
        return base / "results" / C0_FULL85_REL
    # Documented default relative to home iCloud Drive (not hardcoded user)
    home = Path.home()
    return (
        home
        / "Library"
        / "Mobile Documents"
        / "com~apple~CloudDocs"
        / "FlexAIDdS_benchmarks"
        / "results"
        / "campaigns"
        / "C0_full85_defined_cleft_nativeseed_forbidden"
    )


def load_matrix_pin(campaign_dir: Path, cli_pin: str | None) -> tuple[str, str]:
    """Return (md5, source_label)."""
    if cli_pin:
        return cli_pin.strip().lower(), "cli"
    for name in ("RUN_RECEIPT.json", "provenance.json"):
        p = campaign_dir / name
        if not p.is_file():
            continue
        try:
            data = json.loads(p.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        md = data.get("matrix_md5")
        if md:
            return str(md).strip().lower(), name
    return DEFAULT_MATRIX_MD5, "default_pin"


def load_campaign_rows(out_dir: Path) -> list[dict[str, str]]:
    """Load per-target result.csv trees first, then flat summary CSVs."""
    rows: list[dict[str, str]] = []
    for rc in sorted(out_dir.glob("*/result.csv")):
        try:
            batch = list(csv.DictReader(rc.open(newline="")))
            if batch:
                # One authoritative row per target dir (first row)
                rows.append(dict(batch[0]))
        except OSError:
            continue
    if rows:
        return rows
    for name in SUMMARY_CSV_NAMES:
        p = out_dir / name
        if p.is_file():
            try:
                return [dict(r) for r in csv.DictReader(p.open(newline=""))]
            except OSError:
                continue
    # Also accept a single CSV path masquerading as "dir" (handled by caller)
    return rows


def load_rows_from_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh)]


def elected_rmsd(row: dict[str, str]) -> float:
    """Preferred elected RMSD: Hungarian, then crystal."""
    return _f(row, "rmsd_hungarian", "rmsd_to_crystal")


def is_s1(row: dict[str, str]) -> bool:
    """S1: elected pose RMSD ≤ 2.0 Å (not seed echo)."""
    if _truth(row, "seed_echo"):
        return False
    if "success_s1" in row and str(row.get("success_s1", "")).strip() != "":
        return _truth(row, "success_s1")
    # Prefer recomputation from Hungarian (protocol) over engine success_rmsd,
    # which historically used min(crystal, hungarian). When success_rmsd is the
    # only signal and RMSD columns are absent, honour it.
    rh = elected_rmsd(row)
    if math.isfinite(rh):
        return 0.0 <= rh <= RMSD_SUCCESS_A
    if "success_rmsd" in row and str(row.get("success_rmsd", "")).strip() != "":
        return _truth(row, "success_rmsd")
    if "success" in row and str(row.get("success", "")).strip() != "":
        return _truth(row, "success")
    return False


def is_s2(row: dict[str, str], s1: bool) -> bool:
    """S2: S1 ∧ PoseBusters pass."""
    if not s1:
        return False
    if "success_pb" in row and str(row.get("success_pb", "")).strip() != "":
        # success_pb is defined as success_rmsd && pb_pass in DatasetRunner;
        # still require s1 so S2 never exceeds S1 under our S1 definition.
        return _truth(row, "success_pb") and s1
    if "pb_pass" in row and str(row.get("pb_pass", "")).strip() != "":
        return _truth(row, "pb_pass")
    return False


def is_s3(row: dict[str, str]) -> bool:
    """S3: any-pose / BCR ceiling (diagnostic only)."""
    if "success_s3" in row and str(row.get("success_s3", "")).strip() != "":
        return _truth(row, "success_s3")
    bc = _f(row, "best_cluster_rmsd", "rmsd_bcr")
    return math.isfinite(bc) and 0.0 <= bc <= RMSD_SUCCESS_A


def row_matrix_ok(row: dict[str, str], pin: str) -> bool:
    md = str(row.get("matrix_md5", "")).strip().lower()
    if not md:
        return True  # campaign-level pin applies
    return md == pin


def is_claim_eligible(row: dict[str, str], matrix_pin: str) -> tuple[bool, list[str]]:
    """Apply admission gates. Returns (ok, list of fail reasons)."""
    reasons: list[str] = []
    if not _flag0(row, "seed_echo"):
        reasons.append("seed_echo!=0")
    if not _flag0(row, "native_pose_seeded"):
        reasons.append("native_pose_seeded!=0")
    if not row_matrix_ok(row, matrix_pin):
        reasons.append(
            f"matrix_md5_mismatch(got={row.get('matrix_md5')!r}, pin={matrix_pin})"
        )
    # Honour engine claim flag when present and false
    if "protocol_claim_eligible" in row and str(
        row.get("protocol_claim_eligible", "")
    ).strip() != "":
        if not _truth(row, "protocol_claim_eligible"):
            reasons.append("protocol_claim_eligible=0")
    return (len(reasons) == 0, reasons)


def aggregate_rows(
    rows: Iterable[dict[str, str]],
    matrix_pin: str,
    matrix_pin_source: str,
    campaign_dir: str | None = None,
) -> dict[str, Any]:
    all_rows = list(rows)
    claim: list[dict[str, str]] = []
    dropped: list[dict[str, Any]] = []

    for r in all_rows:
        ok, reasons = is_claim_eligible(r, matrix_pin)
        if ok:
            claim.append(r)
        else:
            dropped.append({"pdb_id": _pdb_id(r), "reasons": reasons})

    n = len(claim)
    s1_ids: list[str] = []
    s2_ids: list[str] = []
    s3_ids: list[str] = []
    election_gap_ids: list[str] = []
    s1_fail_ids: list[str] = []

    for r in claim:
        pid = _pdb_id(r)
        s1 = is_s1(r)
        s2 = is_s2(r, s1)
        s3 = is_s3(r)
        if s1:
            s1_ids.append(pid)
        else:
            s1_fail_ids.append(pid)
        if s2:
            s2_ids.append(pid)
        if s3:
            s3_ids.append(pid)
        if s3 and not s1:
            election_gap_ids.append(pid)

    def rate(k: int) -> float:
        return (k / n) if n else 0.0

    report: dict[str, Any] = {
        "contract": "admission_metrics_contract",
        "contract_doc": "benchmarks/protocols/admission_metrics_contract.md",
        "campaign_dir": campaign_dir,
        "matrix_md5_pin": matrix_pin,
        "matrix_md5_pin_source": matrix_pin_source,
        "N_raw": len(all_rows),
        "N_claim": n,
        "N_dropped": len(dropped),
        "dropped_rows": dropped,
        "metrics": {
            "S1": {
                "definition": "elected RMSD <= 2.0 A (Hungarian preferred)",
                "role": "primary_headline",
                "n": len(s1_ids),
                "rate": rate(len(s1_ids)),
                "ids": s1_ids,
            },
            "S2": {
                "definition": "S1 AND PoseBusters pass",
                "role": "secondary",
                "n": len(s2_ids),
                "rate": rate(len(s2_ids)),
                "ids": s2_ids,
            },
            "S3": {
                "definition": "best_cluster_rmsd <= 2.0 A (BCR / any-pose ceiling)",
                "role": "diagnostic_only",
                "n": len(s3_ids),
                "rate": rate(len(s3_ids)),
                "ids": s3_ids,
                "warning": "Do not report S3 as abstract / headline success.",
            },
        },
        "election_gap": {
            "definition": "S3=1 and S1=0 (sampling found near-native; elector missed)",
            "n": len(election_gap_ids),
            "rate": rate(len(election_gap_ids)),
            "ids": election_gap_ids,
        },
        "S1_fail_ids": s1_fail_ids,
        "headline": {
            "metric": "S1",
            "n": len(s1_ids),
            "N": n,
            "rate": rate(len(s1_ids)),
            "label": "S1 elected RMSD <= 2.0 A (claim-eligible only)",
        },
        "admission": {
            "seed_echo": 0,
            "native_pose_seeded": 0,
            "matrix_md5": matrix_pin,
        },
    }
    return report


def format_text_report(report: dict[str, Any]) -> str:
    m = report["metrics"]
    n = report["N_claim"]
    lines = [
        f"campaign_dir: {report.get('campaign_dir')}",
        f"matrix_md5_pin: {report['matrix_md5_pin']} (source={report['matrix_md5_pin_source']})",
        f"N_raw={report['N_raw']}  N_claim={n}  N_dropped={report['N_dropped']}",
        "",
        f"S1 (primary):     {m['S1']['n']}/{n} = {100.0 * m['S1']['rate']:.2f}%",
        f"S2 (S1∧PB):       {m['S2']['n']}/{n} = {100.0 * m['S2']['rate']:.2f}%",
        f"S3 (diagnostic):  {m['S3']['n']}/{n} = {100.0 * m['S3']['rate']:.2f}%  "
        f"[NOT abstract success]",
        f"election_gap (S3∧¬S1): {report['election_gap']['n']}/{n} = "
        f"{100.0 * report['election_gap']['rate']:.2f}%",
    ]
    if report["N_dropped"]:
        lines.append("")
        lines.append("dropped (non-claim):")
        for d in report["dropped_rows"][:20]:
            lines.append(f"  {_pdb_id_from_drop(d)}: {', '.join(d['reasons'])}")
        if report["N_dropped"] > 20:
            lines.append(f"  ... +{report['N_dropped'] - 20} more")
    return "\n".join(lines) + "\n"


def _pdb_id_from_drop(d: dict[str, Any]) -> str:
    return str(d.get("pdb_id") or "?")


def apply_headline(
    report: dict[str, Any],
    headline: str,
    diagnostic_only: bool,
) -> tuple[dict[str, Any], int | None]:
    """Mutate report headline; return (report, error_exit_code_or_None)."""
    h = headline.strip().lower()
    if h not in ("s1", "s2", "s3"):
        return report, 2
    if h == "s3" and not diagnostic_only:
        print(
            "CONTRACT VIOLATION: --headline s3 is not allowed as primary without "
            "--diagnostic-only. S3 is BCR/any-pose diagnostic only; use S1 for "
            "abstract / claim success.",
            file=sys.stderr,
        )
        return report, 2

    m = report["metrics"][h.upper()]
    report["headline"] = {
        "metric": h.upper(),
        "n": m["n"],
        "N": report["N_claim"],
        "rate": m["rate"],
        "label": m["definition"],
        "role": m["role"],
        "diagnostic_only": h == "s3" or diagnostic_only and h != "s1",
    }
    if h == "s3":
        report["headline"]["warning"] = (
            "S3 is diagnostic only — do not quote as abstract success rate."
        )
    return report, None


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "campaign_dir",
        type=Path,
        nargs="?",
        default=None,
        help="Campaign directory with */result.csv and/or summary CSV",
    )
    ap.add_argument(
        "--c0-full85",
        action="store_true",
        help=f"Use $FLEXAIDDS_RESULTS/{C0_FULL85_REL}",
    )
    ap.add_argument(
        "--csv",
        type=Path,
        default=None,
        help="Aggregate a single summary CSV instead of a campaign tree",
    )
    ap.add_argument(
        "--matrix-md5",
        type=str,
        default=None,
        help=f"Matrix pin (default: receipt/provenance or {DEFAULT_MATRIX_MD5})",
    )
    ap.add_argument(
        "--headline",
        type=str,
        default="s1",
        choices=("s1", "s2", "s3", "S1", "S2", "S3"),
        help="Primary headline metric (default s1). s3 requires --diagnostic-only.",
    )
    ap.add_argument(
        "--diagnostic-only",
        action="store_true",
        help="Allow --headline s3 (still labelled diagnostic, never abstract success).",
    )
    ap.add_argument("--json", type=Path, default=None, help="Write full JSON report")
    ap.add_argument(
        "--quiet",
        action="store_true",
        help="JSON only to stdout (no human text block)",
    )
    args = ap.parse_args(argv)

    campaign: Path | None = None
    rows: list[dict[str, str]] = []

    if args.csv is not None:
        if not args.csv.is_file():
            print(f"not a file: {args.csv}", file=sys.stderr)
            return 2
        rows = load_rows_from_csv(args.csv)
        campaign = args.csv.parent
    elif args.c0_full85:
        campaign = resolve_c0_full85_dir()
        if not campaign.is_dir():
            print(
                f"C0 full85 campaign not found: {campaign}\n"
                "Source ~/.flexaidds_env or set FLEXAIDDS_RESULTS.",
                file=sys.stderr,
            )
            return 2
        rows = load_campaign_rows(campaign)
    elif args.campaign_dir is not None:
        campaign = args.campaign_dir
        if campaign.is_file() and campaign.suffix.lower() == ".csv":
            rows = load_rows_from_csv(campaign)
            campaign = campaign.parent
        elif campaign.is_dir():
            rows = load_campaign_rows(campaign)
        else:
            print(f"not a directory or CSV: {campaign}", file=sys.stderr)
            return 2
    else:
        ap.print_help()
        print(
            "\nerror: provide campaign_dir, --csv, or --c0-full85",
            file=sys.stderr,
        )
        return 2

    assert campaign is not None
    pin, pin_src = load_matrix_pin(campaign, args.matrix_md5)
    report = aggregate_rows(
        rows,
        matrix_pin=pin,
        matrix_pin_source=pin_src,
        campaign_dir=str(campaign.resolve()),
    )
    report, err = apply_headline(report, args.headline, args.diagnostic_only)
    if err is not None:
        return err

    if not args.quiet:
        print(format_text_report(report), end="")
        print(
            f"headline: {report['headline']['metric']} "
            f"{report['headline']['n']}/{report['headline']['N']} "
            f"= {100.0 * report['headline']['rate']:.2f}%"
        )

    text = json.dumps(report, indent=2)
    if args.json:
        args.json.write_text(text + "\n")
    if args.quiet:
        print(text)

    return 0 if report["N_claim"] > 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
