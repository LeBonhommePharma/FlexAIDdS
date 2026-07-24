#!/usr/bin/env python3
"""Paired offline election: E_CF (FlexAID-2015 style) vs E_S (first-entropy style).

Both rules are applied to the SAME GA population from ONE search, so the only
variable is the election rule -- the entropy contribution is isolated with zero
GA stochastic noise (see docs/COMPARATIVE_BENCHMARK_METHODOLOGY.md Sec. 3b).

Reference semantics (LIB/BindingMode.cpp of LeBonhommePharma/FlexAID @1a6ae0b):

    w_i = exp(-(1/T) * CF_i)          # beta = 1/T, NOT 1/(k_B*T)
    Z   = sum over the WHOLE population of w_i
    per binding mode (cluster) m:
        H_m = sum_{i in m} (w_i/Z) * CF_i
        S_m = -sum_{i in m} (w_i/Z) * ln(w_i/Z)     # Shannon, no reference state
        G_m = H_m - T * S_m

    E_S  : modes sorted ASCENDING by G_m; mode[0] wins.
    E_CF : modes ranked by raw CF, lowest wins (== the .cad TCF order).
    In both cases the elected structure is the mode's lowest-CF member
    (elect_Representative(false)).

NOTE the normalization: w_i/Z is normalized over the whole population, so the
within-mode probabilities do NOT sum to 1. That is faithful to the reference and
is reproduced exactly here.

Numerical stability: clash poses carry CF ~ 1e4, so exp(-CF/T) underflows in
naive float. All weights are computed in log space with a log-sum-exp shift,
which is exact for the ratios w_i/Z that the formulas actually need.

Usage:
    elect_paired.py <run_dir_or_rrd> [--temper 21] [--cutoff 2.0] [--json out.json]

Copyright 2026 Le Bonhomme Pharma
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

# CF at or above this is a clash sentinel, not a physical score.
CLASH_SENTINEL = 10000.0


def parse_rrd(path: Path) -> list[dict]:
    """Parse a FlexAID .rrd, auto-detecting the V1 vs V2 column layout.

    The two reference binaries do NOT write the same columns (verified against
    the pinned sources and real output):

      V1 `b555e0e` : idx cluster clusRMSD rmsd rmsd_hungarian CF [genes]   (6 numeric)
      V2 `1a6ae0b` : idx cluster clusRMSD rmsd                CF [genes]   (5 numeric)
                     (the Hungarian line is commented out on the entropy branch)

    Hardcoding cf=parts[4] silently reads V1's symmetry-corrected RMSD as the
    CF, which corrupts every downstream weight. Detect instead: the numeric
    prefix ends at the '[' that opens the gene vector.
    """
    poses: list[dict] = []
    for line in path.read_text(errors="ignore").splitlines():
        head = line.split("[", 1)[0]  # drop the gene vector
        parts = head.split()
        if len(parts) < 5:
            continue
        try:
            n = len(parts)
            if n >= 6:  # V1 layout: CF is last, Hungarian RMSD precedes it
                cf, rmsd_h = float(parts[5]), float(parts[4])
            else:  # V2 layout: CF is col 4, no Hungarian
                cf, rmsd_h = float(parts[4]), None
            poses.append(
                {
                    "idx": int(parts[0]),
                    "cluster": int(parts[1]),
                    "cluster_rmsd": float(parts[2]),
                    "rmsd_ref": float(parts[3]),
                    "rmsd_hungarian": rmsd_h,
                    "cf": cf,
                    "layout": "V1" if n >= 6 else "V2",
                }
            )
        except ValueError:
            continue  # header/comment line
    return poses


def elect(poses: list[dict], temper: float) -> dict:
    """Apply both election rules to one population. Returns both winners."""
    if not poses:
        return {"error": "empty population"}

    # Poses beyond the MAXRES ceiling are written with cluster == -1 (unassigned).
    # Grouping them would fabricate a single spurious mega-mode (observed: 376 of
    # 500 poses on 1SJ0) that can win an election. Drop them, and count them.
    n_unassigned = sum(1 for p in poses if p["cluster"] < 0)
    poses = [p for p in poses if p["cluster"] >= 0]
    if not poses:
        return {"error": "all poses unassigned (raise MAXRES to NUMCHROM)"}

    by_cluster: dict[int, list[dict]] = defaultdict(list)
    for p in poses:
        by_cluster[p["cluster"]].append(p)

    # --- log-space Boltzmann weights over the WHOLE population -------------
    # log w_i = -CF_i / T ; log Z = logsumexp(log w_i). Ratios w_i/Z are exact
    # under the shift, so clash poses (CF~1e4) simply carry ~zero weight
    # instead of underflowing to a NaN.
    log_w = {p["idx"]: -p["cf"] / temper for p in poses}
    m = max(log_w.values())
    log_Z = m + math.log(sum(math.exp(lw - m) for lw in log_w.values()))

    modes = []
    for cid, members in by_cluster.items():
        H = 0.0
        S = 0.0
        for p in members:
            log_p = log_w[p["idx"]] - log_Z  # ln(w_i/Z)
            prob = math.exp(log_p)
            if prob <= 0.0:
                continue  # zero-weight (clash) pose contributes nothing
            H += prob * p["cf"]
            S -= prob * log_p
        rep = min(members, key=lambda p: p["cf"])  # elect_Representative(false)
        modes.append(
            {
                "cluster": cid,
                "freq": len(members),
                "tcf": rep["cf"],  # lowest CF in the mode (== .cad TCF)
                "H": H,
                "S": S,
                "G": H - temper * S,
                "rep_idx": rep["idx"],
                "rep_rmsd": rep["rmsd_ref"],
            }
        )

    e_cf = min(modes, key=lambda x: x["tcf"])  # V1-style: lowest raw CF
    e_s = min(modes, key=lambda x: x["G"])  # V2-style: lowest G = H - T*S
    return {
        "n_poses": len(poses),
        "n_modes": len(modes),
        "n_unassigned": n_unassigned,
        "layout": poses[0].get("layout"),
        "n_clash": sum(1 for p in poses if p["cf"] >= CLASH_SENTINEL),
        "temper": temper,
        "E_CF": e_cf,
        "E_S": e_s,
        "agree": e_cf["cluster"] == e_s["cluster"],
        "modes": sorted(modes, key=lambda x: x["G"]),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("target", type=Path, help=".rrd file, or a dir searched recursively for *.rrd")
    ap.add_argument("--temper", type=float, default=21.0, help="T (TEMPER); reference ISMB calibration = 21")
    ap.add_argument("--cutoff", type=float, default=2.0, help="sub-X A success cutoff")
    ap.add_argument("--json", type=Path, help="write the full per-target report here")
    args = ap.parse_args()

    rrds = (
        [args.target]
        if args.target.is_file()
        else sorted(args.target.rglob("*.rrd"))
    )
    if not rrds:
        print(f"no .rrd found under {args.target}", file=sys.stderr)
        return 2

    rows, n_cf, n_s, n_both = [], 0, 0, 0
    for rrd in rrds:
        poses = parse_rrd(rrd)
        if not poses:
            continue
        r = elect(poses, args.temper)
        if "error" in r:
            continue
        name = rrd.stem
        ok_cf = r["E_CF"]["rep_rmsd"] < args.cutoff
        ok_s = r["E_S"]["rep_rmsd"] < args.cutoff
        n_cf += ok_cf
        n_s += ok_s
        n_both += 1
        rows.append({"target": name, "ok_cf": ok_cf, "ok_s": ok_s, **r})
        print(
            f"{name:<28} poses={r['n_poses']:<5} modes={r['n_modes']:<4} "
            f"E_CF: rmsd={r['E_CF']['rep_rmsd']:6.2f} {'OK ' if ok_cf else '   '} "
            f"E_S: rmsd={r['E_S']['rep_rmsd']:6.2f} {'OK ' if ok_s else '   '} "
            f"{'same' if r['agree'] else 'DIFFER'}"
        )

    if n_both:
        # McNemar discordant pairs: the only cells that carry information.
        b = sum(1 for r in rows if r["ok_cf"] and not r["ok_s"])  # CF wins
        c = sum(1 for r in rows if r["ok_s"] and not r["ok_cf"])  # S wins
        print(f"\n{'='*70}\nPAIRED ELECTION RESULT  (n={n_both}, T={args.temper}, cutoff={args.cutoff} A)")
        print(f"  E_CF (FlexAID-2015 style) : {n_cf}/{n_both} = {100*n_cf/n_both:.1f}%")
        print(f"  E_S  (first-entropy style): {n_s}/{n_both} = {100*n_s/n_both:.1f}%")
        print(f"  discordant: E_CF-only={b}  E_S-only={c}  (McNemar uses only these)")
        if b + c:
            # exact binomial two-sided p under H0: P(win)=0.5
            k, n = min(b, c), b + c
            p = min(1.0, 2 * sum(math.comb(n, i) for i in range(k + 1)) / (2**n))
            print(f"  McNemar exact p = {p:.4f}")
        else:
            print("  McNemar: no discordant pairs -> the rules are indistinguishable here")

    if args.json:
        args.json.write_text(json.dumps(rows, indent=2))
        print(f"\nwrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
