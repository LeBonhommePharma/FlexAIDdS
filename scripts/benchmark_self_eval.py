#!/usr/bin/env python3
"""Enforce a priori / a posteriori benchmark self-eval contract.

Defers floors to PHASE4 / METHODOLOGY narrative; validates structure and
on-disk metrics without reimplementing docking.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

REQUIRED_APRIORI = [
    "one_variable",
    "panel_class",
    "codes",
    "matrix_pin",
    "no_sec",
    "sol9",
    "matched_control",
    "magnitude_floor",
    "report_tiers_separately",
]

VALID_STATUS = {
    "PASS",
    "FAIL",
    "PASS_LIVENESS",
    "VOID",
    "INVALID",
    "MISSING_OUT",
    "IN_PROGRESS",
    "NOT_RUN",
}

MATRIX_9DC9 = "9dc93717dfed0698006d88dd6a9627bc"


def load_apriori(path: Path) -> dict:
    return json.loads(path.read_text())


def validate_apriori(ap: dict) -> list[str]:
    errs = []
    for k in REQUIRED_APRIORI:
        if k not in ap:
            errs.append(f"missing a priori field: {k}")
    if ap.get("matrix_pin") not in (MATRIX_9DC9, "9dc9", "9dc93717dfed0698006d88dd6a9627bc"):
        errs.append(f"matrix_pin must be 9dc9, got {ap.get('matrix_pin')}")
    if ap.get("panel_class") in ("SEARCH_MISS", "NEAR_MISS", "GROSS_MISS"):
        if not ap.get("no_sec"):
            errs.append("Phase-4 SEARCH sampling docks require no_sec=true")
    codes = ap.get("codes") or []
    if isinstance(codes, str):
        codes = [c.strip() for c in codes.split(",") if c.strip()]
    if ap.get("panel_class") == "NEAR_MISS":
        bad = set(codes) - {"1N1M", "1L7F"}
        if bad:
            errs.append(f"NEAR_MISS codes must be subset of 1N1M,1L7F; extra {bad}")
    if ap.get("report_tiers_separately") is not True:
        errs.append("report_tiers_separately must be true")
    return errs


def result_metrics(out: Path, codes: list[str]) -> dict[str, dict]:
    m = {}
    for c in codes:
        p = out / c / "result.csv"
        if not p.is_file():
            m[c] = {"missing": True}
            continue
        row = next(csv.DictReader(p.open()))
        m[c] = {
            "bcr": float(row["best_cluster_rmsd"]),
            "s3": float(row["conditional_scanned_pool_ceiling"]),
            "elect": float(row["rmsd_to_crystal"]),
            "missing": False,
        }
    return m


def count_marker(out: Path, marker: str) -> int:
    n = 0
    for log in out.rglob("stdout.log"):
        n += log.read_text(errors="replace").count(marker)
    return n


def wipeout_signature(out: Path) -> bool:
    for log in out.rglob("stdout.log"):
        t = log.read_text(errors="replace")
        if re.search(r"wipeout|CF\s*[=~]\s*0\.0+", t, re.I):
            return True
    return False


def posteriori_g4_1_style(
    control: Path,
    treatments: dict[str, Path],
    codes: list[str],
    *,
    marker: str = "[BOOM]",
    magnitude_floor: float = -0.5,
) -> dict:
    ctrl = result_metrics(control, codes)
    if any(ctrl[c].get("missing") for c in codes):
        return {"status": "IN_PROGRESS", "reason": "control result.csv missing"}
    tx_out = {}
    best_name = None
    best_mean = 999.0
    for name, path in treatments.items():
        m = result_metrics(path, codes)
        if any(m[c].get("missing") for c in codes):
            tx_out[name] = {"status": "IN_PROGRESS", "metrics": m}
            continue
        deltas = [m[c]["bcr"] - ctrl[c]["bcr"] for c in codes]
        mean_d = sum(deltas) / len(deltas)
        n_sub2 = sum(1 for c in codes if m[c]["bcr"] <= 2.0)
        n_boom = count_marker(path, marker)
        wo = wipeout_signature(path)
        mag = mean_d <= magnitude_floor or n_sub2 >= 1
        l4_ok = n_boom > 0
        st = "PASS" if mag and l4_ok and not wo else (
            "PASS_LIVENESS" if l4_ok and not mag and not wo else "FAIL"
        )
        if not l4_ok:
            st = "FAIL"
        tx_out[name] = {
            "status": st,
            "mean_dBCR": mean_d,
            "n_bcr_sub2": n_sub2,
            "n_markers": n_boom,
            "wipeout": wo,
            "metrics": m,
            "per_target_dBCR": {c: m[c]["bcr"] - ctrl[c]["bcr"] for c in codes},
        }
        if mean_d < best_mean and st != "IN_PROGRESS":
            best_mean = mean_d
            best_name = name
    ctrl_boom = count_marker(control, marker)
    accept = False
    if best_name and tx_out[best_name]["status"] == "PASS" and ctrl_boom == 0:
        accept = True
    overall = "PASS" if accept else "FAIL"
    if any(tx_out[n].get("status") == "IN_PROGRESS" for n in tx_out):
        overall = "IN_PROGRESS"
    return {
        "status": overall,
        "accept_magnitude": accept,
        "best_treatment": best_name,
        "control_markers": ctrl_boom,
        "control": ctrl,
        "treatments": tx_out,
        "codes": codes,
        "magnitude_floor": magnitude_floor,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("preflight")
    p.add_argument("--apriori", type=Path, required=True)
    p.add_argument("--write-out", type=Path, help="copy apriori into OUT/APRIORI.json")

    q = sub.add_parser("posteriori")
    q.add_argument("--control", type=Path, required=True)
    q.add_argument("--treatment", action="append", nargs=2, metavar=("NAME", "PATH"), required=True)
    q.add_argument("--codes", default="1N1M,1L7F")
    q.add_argument("--marker", default="[BOOM]")
    q.add_argument("--json", action="store_true")

    v = sub.add_parser("validate-contract-doc")
    v.add_argument("--path", type=Path, default=Path("workorders/BENCHMARK_SELF_EVAL_CONTRACT.md"))

    args = ap.parse_args(argv)

    if args.cmd == "preflight":
        ap_data = load_apriori(args.apriori)
        errs = validate_apriori(ap_data)
        if errs:
            print("APRIORI_FAIL")
            for e in errs:
                print(" ", e)
            return 2
        print("APRIORI_OK")
        if args.write_out:
            args.write_out.mkdir(parents=True, exist_ok=True)
            (args.write_out / "APRIORI.json").write_text(json.dumps(ap_data, indent=2) + "\n")
        return 0

    if args.cmd == "posteriori":
        codes = [c.strip() for c in args.codes.split(",") if c.strip()]
        treatments = {n: Path(p) for n, p in args.treatment}
        rec = posteriori_g4_1_style(
            args.control, treatments, codes, marker=args.marker
        )
        print(
            f"status={rec['status']} accept={rec['accept_magnitude']} "
            f"best={rec['best_treatment']} ctrl_markers={rec['control_markers']}"
        )
        if args.json:
            print(json.dumps(rec, indent=2))
        return 0 if rec["status"] == "PASS" else 1

    if args.cmd == "validate-contract-doc":
        text = args.path.read_text()
        need = [
            "a priori",
            "a posteriori",
            "PHASE4_GATES_ACTUALIZED",
            "METHODOLOGY.md",
            "magnitude",
            "NEAR_MISS",
            "Sol #9",
            "9dc9",
        ]
        missing = [n for n in need if n.lower() not in text.lower() and n not in text]
        # METHODOLOGY.md might be found; case for a priori
        missing = []
        for n in need:
            if n not in text and n.lower() not in text.lower():
                missing.append(n)
        if missing:
            print("CONTRACT_DOC_FAIL", missing)
            return 2
        print("CONTRACT_DOC_OK", args.path)
        return 0

    return 2


if __name__ == "__main__":
    raise SystemExit(main())
