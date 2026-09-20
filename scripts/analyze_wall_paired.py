#!/usr/bin/env python3
"""Paired c1_on vs c1_off analysis over the wall_paired_85 campaign.

WHAT THIS MEASURES AND WHAT IT CANNOT
-------------------------------------
Primary endpoint is ORACLE-IN-POOL: is a sub-2 A pose present anywhere in the
target's retained pool. That is the endpoint the FLEXAIDDS_WAL_C1 default flip
was ratified on (satisfied then on 12 admissible targets; this is 82).

Secondary endpoint is MIN-CF ELECTION, and it is an INFERENCE, not a reading.
The campaign's result.csv left `elected_pose_path` empty, so which pose the
engine elected was never recorded. Pooled min-CF was measured byte-identical to
the engine's own election on 85/85 Arm A targets once the re-election composite
was removed, so min-CF is a well-grounded stand-in -- but it is labelled as
min-CF throughout and must never be quoted as "the engine elected".

A SEEDING CHECK RUNS FIRST. Every campaign row carries
native_pose_seeded=1 / native_pose_seed_fraction=0.9000, which reads like a 90%
crystal-pose injection. It is not: DatasetRunner.cpp:8751-8755 derives that flag
purely from `mode not in {AUTONOMOUS, DEFINED_CLEFT_REDOCK, ORACLE_CEILING}` and
then stamps 0.90 as a LABEL, while :5997 states UNSET's seed fraction is always
0.0 and :6003 forces seed_elitism off for UNSET. The empirical test below is
independent of that source reading: a real 90% seeding would place a ~0 A pose
in essentially every pool, so a median oracle far above 0 refutes it.

POOL DEPTH IS NOT EQUAL ACROSS ARMS. The arms docked independently, so pools
differ in size (one target differs 5-fold). Oracle-in-pool is monotone in pool
size, so the comparison is reported BOTH over all paired targets AND restricted
to targets whose pools are within 10% of each other -- otherwise a wall-form
effect and a pool-depth effect are quoted as one number.
"""

from __future__ import annotations

import csv
import glob
import json
import os
import statistics
import sys

import numpy as np
from scipy import stats

import argparse

# Repo scripts must not carry a developer's absolute paths. The results root is
# a CLI argument, falling back to $FLEXAIDDS_RESULTS so a campaign driver can set
# it once for a whole run.
R = os.environ.get("FLEXAIDDS_RESULTS", "")
FORMS = ("c1_on", "c1_off")


def bank(arm):
    """(pose_file, restart) -> symmetry-corrected RMSD, per target."""
    out = {}
    for f in glob.glob(f"{R}/pose_bank_{arm}/*.csv"):
        b = os.path.basename(f)
        if "MANIFEST" in b or "REFERENCE" in b:
            continue
        rows = {}
        for r in csv.DictReader(open(f)):
            try:
                rows[(r["pose_file"], r.get("restart", "0"))] = float(r["rmsd_symmcorr"])
            except (KeyError, TypeError, ValueError):
                continue
        if rows:
            out[b[:-4]] = rows
    return out


def mcnemar(a_only, b_only):
    d = a_only + b_only
    if d < 5:
        return f"unreachable (d={d}, exact floor is 5)"
    p = stats.binomtest(max(a_only, b_only), d, 0.5, alternative="greater").pvalue
    return f"exact p={p:.4f} (d={d})"


def main():
    global R
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=R or None, required=not R,
                    help="results root holding pose_bank_<arm>/ (or $FLEXAIDDS_RESULTS)")
    a = ap.parse_args()
    R = a.root
    B = {x: bank(x) for x in FORMS}
    cf = json.load(open("handoff/campaign_cf.json"))
    paired = sorted(set(B["c1_on"]) & set(B["c1_off"]))
    if not paired:
        print("VACUOUS: no paired targets", file=sys.stderr)
        return 3

    rec = []
    for t in paired:
        row = {"target": t}
        for a in FORMS:
            d = B[a][t]
            row[f"n_{a}"] = len(d)
            row[f"oracle_{a}"] = min(d.values())
            best_rmsd, best_cf = None, None
            for k, v in cf.get(a, {}).get(t, {}).items():
                nm, rs = k.rsplit("|", 1)
                if (nm, rs) in d and (best_cf is None or v < best_cf):
                    best_cf, best_rmsd = v, d[(nm, rs)]
            row[f"mincf_{a}"] = best_rmsd
        rec.append(row)
    N = len(rec)

    # ---- seeding check, independent of the source reading --------------------
    print("=== SEEDING CHECK (flag says seeded=1/0.9; source says UNSET => 0.0) ===")
    med_on = statistics.median([r["oracle_c1_on"] for r in rec])
    for a in FORMS:
        o = [r[f"oracle_{a}"] for r in rec]
        print(f"  {a}: oracle min {min(o):.3f}  median {statistics.median(o):.3f}  "
              f"<0.5A {sum(1 for x in o if x < 0.5)}/{N}  <2A {sum(1 for x in o if x < 2)}/{N}")
    print("  A true 90% crystal seeding puts a ~0 A pose in essentially every pool.")
    print("  VERDICT: " + ("NOT SEEDED -- the CSV flag is a derived label, as the source states"
                           if med_on > 0.5 else "*** POOLS LOOK SEEDED -- STOP ***"))

    # ---- per-arm rates -------------------------------------------------------
    print(f"\n=== PER-ARM, n={N} PAIRED TARGETS ===")
    print(f"{'form':<9}{'oracle<2A':>13}{'median':>10}{'minCF<2A':>13}{'median':>10}{'med pool':>10}")
    for a in FORMS:
        o = [r[f"oracle_{a}"] for r in rec]
        m = [r[f"mincf_{a}"] for r in rec if r[f"mincf_{a}"] is not None]
        npool = [r[f"n_{a}"] for r in rec]
        print(f"{a:<9}{sum(1 for x in o if x < 2):>6}/{N:<6}{statistics.median(o):>10.3f}"
              f"{sum(1 for x in m if x < 2):>6}/{len(m):<6}"
              f"{(statistics.median(m) if m else float('nan')):>10.3f}{statistics.median(npool):>10.0f}")

    # ---- paired tests, full and pool-matched ---------------------------------
    matched = [r for r in rec
               if min(r["n_c1_on"], r["n_c1_off"]) / max(r["n_c1_on"], r["n_c1_off"]) >= 0.90]
    for subset, lbl in ((rec, f"ALL {N} PAIRED"),
                        (matched, f"POOL-MATCHED (within 10%), n={len(matched)}")):
        print(f"\n=== PAIRED TESTS -- {lbl} (negative delta = c1_on better) ===")
        for name, key in (("oracle-in-pool", "oracle"), ("min-CF elected", "mincf")):
            pp = [(r[f"{key}_c1_on"], r[f"{key}_c1_off"]) for r in subset
                  if r.get(f"{key}_c1_on") is not None and r.get(f"{key}_c1_off") is not None]
            if len(pp) < 6:
                print(f"  {name:<16} n={len(pp)} -- too few to test")
                continue
            d = np.array([x - y for x, y in pp])
            nz = d[d != 0]
            p = stats.wilcoxon(nz).pvalue if len(nz) >= 6 else float("nan")
            a_only = sum(1 for x, y in pp if x < 2 <= y)
            b_only = sum(1 for x, y in pp if y < 2 <= x)
            print(f"  {name:<16} n={len(pp)}  median d={np.median(d):+.3f}  "
                  f"better {int((d < 0).sum())} worse {int((d > 0).sum())} tied {int((d == 0).sum())}  "
                  f"Wilcoxon p={p:.4f}")
            print(f"  {'':<16} crossings: c1_on-only {a_only}  c1_off-only {b_only}  "
                  f"McNemar {mcnemar(a_only, b_only)}")

    with open("wall_paired_82.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rec[0]))
        w.writeheader()
        w.writerows(rec)
    json.dump({"n": N, "n_pool_matched": len(matched),
               "oracle_sub2": {a: sum(1 for r in rec if r[f"oracle_{a}"] < 2) for a in FORMS},
               "oracle_median": {a: statistics.median([r[f"oracle_{a}"] for r in rec]) for a in FORMS}},
              open("handoff/wall_paired_82_facts.json", "w"), indent=1)
    print(f"\nwrote wall_paired_82.csv ({N} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
