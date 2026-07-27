#!/usr/bin/env python3
"""Apply campaign 'what would flip the order' rules to G4.1 / election / scoring results.

Does not reimplement docking — reads result.csv BCR/elect columns and optional
BOOM log markers. Encodes ACCEPT floors from PHASE4 + forward_plan_confidence.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path


MAGNITUDE_FLOOR = -0.5  # mean ΔBCR treatment−control
SUB2 = 2.0


def load_result(path: Path) -> dict:
    row = next(csv.DictReader(path.open()))
    return {
        "bcr": float(row["best_cluster_rmsd"]),
        "s3": float(row["conditional_scanned_pool_ceiling"]),
        "elect": float(row["rmsd_to_crystal"]),
        "success_rmsd": int(float(row.get("success_rmsd") or 0)),
    }


def collect_arm(out: Path, codes: list[str]) -> dict[str, dict]:
    arm: dict[str, dict] = {}
    for code in codes:
        p = out / code / "result.csv"
        if not p.is_file():
            raise FileNotFoundError(p)
        arm[code] = load_result(p)
    return arm


def boom_l4(out: Path) -> dict:
    """Count [BOOM] markers; detect wipeout signature CF≈0 early stop if present."""
    n_boom = 0
    wipeout = False
    for log in out.rglob("stdout.log"):
        text = log.read_text(errors="replace")
        n_boom += len(re.findall(r"\[BOOM\]", text))
        if re.search(r"wipeout|CF\s*[=~]\s*0\.0+\b", text, re.I):
            wipeout = True
    return {"n_boom_markers": n_boom, "wipeout_flag": wipeout}


def evaluate_g4_1(
    control: dict[str, dict],
    treatments: dict[str, dict[str, dict]],
    *,
    codes: list[str],
) -> dict:
    """treatments: name -> code -> metrics."""
    per_tx = {}
    for name, arm in treatments.items():
        deltas = []
        elect_reg = []
        sub2 = []
        for c in codes:
            db = arm[c]["bcr"] - control[c]["bcr"]
            de = arm[c]["elect"] - control[c]["elect"]
            deltas.append(db)
            elect_reg.append(de > 0.25)  # elected worsened by >0.25 Å
            sub2.append(arm[c]["bcr"] <= SUB2)
        mean_db = sum(deltas) / len(deltas)
        per_tx[name] = {
            "mean_dBCR_vs_control": mean_db,
            "per_target_dBCR": {c: arm[c]["bcr"] - control[c]["bcr"] for c in codes},
            "n_bcr_sub2": sum(1 for x in sub2 if x),
            "n_elect_regression": sum(1 for x in elect_reg if x),
            "magnitude_pass": mean_db <= MAGNITUDE_FLOOR or any(sub2),
            "no_elect_regression": not any(elect_reg),
        }
    # best treatment by mean dBCR
    best_name = min(per_tx, key=lambda k: per_tx[k]["mean_dBCR_vs_control"])
    best = per_tx[best_name]
    accept = (
        best["magnitude_pass"]
        and best["no_elect_regression"]
    )
    # flip-order rule 1
    if accept:
        flip = {
            "rule": "G4.1_BOOM_hits_magnitude",
            "action": "PROMOTE_BOOM_to_claim_recipe; deprioritize_G4.3",
            "priority_order": ["BOOM_in_claim", "1N1M_election_audit", "G4.3_later"],
        }
    else:
        flip = {
            "rule": "G4.1_null_on_near_miss",
            "action": "Sampling_inject_exhausted_for_now; prioritize_election_fix_then_G4.3",
            "priority_order": ["1N1M_election_fix", "G4.3_mutation", "new_search_arch"],
        }
    return {
        "accept_g4_1": accept,
        "best_treatment": best_name,
        "treatments": per_tx,
        "flip_order": flip,
        "codes": codes,
        "magnitude_floor": MAGNITUDE_FLOOR,
    }


def evaluate_election_offline(
    pop_tsv: Path,
    result_csv: Path,
    *,
    elect_target: float = 2.5,
) -> dict:
    """If pool has near-native, compare elected vs best-available head/population."""
    rows = list(csv.DictReader(pop_tsv.open(), delimiter="\t"))
    pop_best = min(float(r["rmsd_sym"]) for r in rows)
    elects = [r for r in rows if str(r.get("is_elected", "0")) in ("1", "true", "True")]
    emit_best = min(float(r["rmsd_sym"]) for r in elects) if elects else float("nan")
    # best CF among rows with rmsd_sym <= elect_target + 1
    near = [r for r in rows if float(r["rmsd_sym"]) <= elect_target + 1.0]
    if near:
        best_near = min(near, key=lambda r: float(r["cf_total"]))
        oracle_elect_rmsd = float(best_near["rmsd_sym"])
        oracle_elect_cf = float(best_near["cf_total"])
    else:
        oracle_elect_rmsd = pop_best
        oracle_elect_cf = None
    res = next(csv.DictReader(result_csv.open()))
    actual_elect = float(res["rmsd_to_crystal"])
    # Flip if rank-0 elect is far from a pool geometry that is already ≤ elect_target
    # (emitted head best or any population member).
    pool_near = min(pop_best, emit_best, oracle_elect_rmsd)
    would_flip = pool_near <= elect_target and actual_elect - pool_near > 1.0
    flip = None
    if would_flip:
        flip = {
            "rule": "1N1M_elect_fixable_offline",
            "action": "Selection_lever_becomes_P0_for_bellwether",
            "priority_order": [
                "election_fix_P0",
                "G4.1_BOOM_secondary_or_parallel",
                "G4.3_later",
            ],
            "pool_near_rmsd": pool_near,
            "emitted_best": emit_best,
            "oracle_cf_near_rmsd": oracle_elect_rmsd,
            "actual_elect": actual_elect,
            "pop_best": pop_best,
        }
    return {
        "pop_best_rmsd_sym": pop_best,
        "emitted_best_rmsd_sym": emit_best,
        "actual_elect_rmsd": actual_elect,
        "oracle_cf_rank_among_near_rmsd": oracle_elect_rmsd,
        "oracle_cf": oracle_elect_cf,
        "would_flip_to_election_P0": would_flip,
        "flip_order": flip,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    g = sub.add_parser("g4_1", help="Score G4.1 multi-arm OUT tree")
    g.add_argument("--control", type=Path, required=True)
    g.add_argument("--treatment", action="append", nargs=2, metavar=("NAME", "OUT"), required=True)
    g.add_argument("--codes", default="1N1M,1L7F")
    g.add_argument("--json", action="store_true")

    e = sub.add_parser("election_offline", help="1N1M-style election gap offline")
    e.add_argument("--pop-tsv", type=Path, required=True)
    e.add_argument("--result-csv", type=Path, required=True)
    e.add_argument("--json", action="store_true")

    f = sub.add_parser("flip_summary", help="Merge g4_1 + election JSON into final order")
    f.add_argument("--g4-1-json", type=Path)
    f.add_argument("--election-json", type=Path)
    f.add_argument("--json", action="store_true")

    args = ap.parse_args(argv)
    if args.cmd == "g4_1":
        codes = [c.strip() for c in args.codes.split(",") if c.strip()]
        control = collect_arm(args.control, codes)
        treatments = {name: collect_arm(Path(out), codes) for name, out in args.treatment}
        # L4 from first treatment
        l4 = {}
        for name, outp in args.treatment:
            l4[name] = boom_l4(Path(outp))
        ctrl_l4 = boom_l4(args.control)
        rec = evaluate_g4_1(control, treatments, codes=codes)
        rec["l4_treatments"] = l4
        rec["l4_control"] = ctrl_l4
        # fail accept if control has boom markers or treatment has zero
        for name, m in l4.items():
            if m["n_boom_markers"] == 0:
                rec["accept_g4_1"] = False
                rec["treatments"][name]["l4_fail"] = "no [BOOM] markers"
            if m["wipeout_flag"]:
                rec["accept_g4_1"] = False
                rec["treatments"][name]["l4_fail"] = "wipeout signature"
        if ctrl_l4["n_boom_markers"] > 0:
            rec["accept_g4_1"] = False
            rec["l4_control_fail"] = "control must have zero [BOOM]"
        # recompute flip if accept cleared
        if not rec["accept_g4_1"]:
            rec["flip_order"] = {
                "rule": "G4.1_null_or_l4_fail",
                "action": "prioritize_election_fix_then_G4.3",
                "priority_order": ["1N1M_election_fix", "G4.3_mutation", "new_search_arch"],
            }
        print(
            f"accept={rec['accept_g4_1']} best={rec['best_treatment']} "
            f"mean_dBCR={rec['treatments'][rec['best_treatment']]['mean_dBCR_vs_control']:+.4f} "
            f"flip={rec['flip_order']['rule']}"
        )
        if args.json:
            print(json.dumps(rec, indent=2))
        return 0 if rec["accept_g4_1"] else 1

    if args.cmd == "election_offline":
        rec = evaluate_election_offline(args.pop_tsv, args.result_csv)
        print(
            f"pop_best={rec['pop_best_rmsd_sym']:.4f} actual_elect={rec['actual_elect_rmsd']:.4f} "
            f"oracle_near={rec['oracle_cf_rank_among_near_rmsd']:.4f} "
            f"flip_election_P0={rec['would_flip_to_election_P0']}"
        )
        if args.json:
            print(json.dumps(rec, indent=2))
        return 0

    if args.cmd == "flip_summary":
        g4 = json.loads(args.g4_1_json.read_text()) if args.g4_1_json and args.g4_1_json.is_file() else None
        el = json.loads(args.election_json.read_text()) if args.election_json and args.election_json.is_file() else None
        # precedence: election P0 if offline says so and G4.1 null; else G4.1 promote
        if g4 and g4.get("accept_g4_1"):
            order = g4["flip_order"]
        elif el and el.get("would_flip_to_election_P0"):
            order = el.get("flip_order") or {
                "rule": "1N1M_elect_fixable_offline",
                "action": "Selection_lever_P0",
                "priority_order": ["election_fix_P0", "G4.1_secondary"],
            }
        elif g4:
            order = g4["flip_order"]
        else:
            order = {
                "rule": "insufficient_data",
                "action": "run_G4.1_then_reevaluate",
                "priority_order": ["G4.1_BOOM", "election_offline"],
            }
        print(json.dumps({"final_flip_order": order, "g4_1_accept": bool(g4 and g4.get("accept_g4_1"))}, indent=2))
        if args.json:
            pass
        return 0

    return 2


if __name__ == "__main__":
    raise SystemExit(main())
