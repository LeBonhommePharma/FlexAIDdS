#!/usr/bin/env python3
"""Decompose the paired wall effect into SEARCH and ELECTION components.

WHY THIS EXISTS
---------------
Two experiments on FLEXAIDDS_WAL_C1 looked contradictory:

  (b) 84-target rescore over Arm A's FROZEN pool: c1_off 24/84 vs c1_on 25/84,
      exact p=0.5000 -- c1_off nominally WORST of three forms.
  (c) 82-target paired campaign, each arm DOCKED SEPARATELY: c1_off better on
      the elected pose, Wilcoxon p=0.0107, sign test p=0.00064.

They are not contradictory: they measure different causal paths. A rescore holds
the pool fixed, so only ELECTION can move. Separate docking lets the wall form
change the CF that drives the GA, so the POOL ITSELF differs -- SEARCH -- and the
elected-pose difference is search+election combined.

The decomposition is exact, per target:

    total    = mincf_on - mincf_off
    search   = oracle_on - oracle_off
    election = (mincf_on - oracle_on) - (mincf_off - oracle_off)
             = total - search

`election` is the difference in SELECTION GAP: how far above its own pool's best
available pose each arm actually elects. That isolates the elector from the pool
it is electing out of, which is the only way to compare against a rescore.

REPORTING RULE ENCODED HERE. The election term's MEAN and MEDIAN disagree by
roughly an order of magnitude (outlier-driven), so this script prints both plus
a sign test, and the figure shows the distribution rather than a bar. A mean
alone would make election look like the driver when the median and sign test
say it is not.
"""

from __future__ import annotations

import csv
import json
import sys

import numpy as np
from scipy import stats


def load(path="wall_paired_82.csv"):
    rows = []
    for r in csv.DictReader(open(path)):
        d = {"target": r["target"]}
        try:
            for k in ("oracle_c1_on", "oracle_c1_off", "mincf_c1_on",
                      "mincf_c1_off", "n_c1_on", "n_c1_off"):
                d[k] = float(r[k])
        except (TypeError, ValueError):
            continue          # a target without a min-CF election in both arms
        rows.append(d)
    return rows


def describe(d, label):
    nz = d[d != 0]
    w = stats.wilcoxon(nz).pvalue if len(nz) >= 6 else float("nan")
    s = stats.binomtest(int((nz > 0).sum()), len(nz), 0.5).pvalue if len(nz) else float("nan")
    return {"label": label, "n": int(len(d)), "median": float(np.median(d)),
            "mean": float(d.mean()), "off_better": int((d > 0).sum()),
            "on_better": int((d < 0).sum()), "tied": int((d == 0).sum()),
            "wilcoxon_p": float(w), "sign_p": float(s)}


def main():
    rows = load()
    if len(rows) < 6:
        print("VACUOUS: too few paired targets", file=sys.stderr)
        return 3

    total = np.array([r["mincf_c1_on"] - r["mincf_c1_off"] for r in rows])
    search = np.array([r["oracle_c1_on"] - r["oracle_c1_off"] for r in rows])
    gap_on = np.array([r["mincf_c1_on"] - r["oracle_c1_on"] for r in rows])
    gap_off = np.array([r["mincf_c1_off"] - r["oracle_c1_off"] for r in rows])
    election = gap_on - gap_off

    # the decomposition must be exact; assert rather than trust
    assert np.allclose(total, search + election, atol=1e-9), "decomposition is not exact"

    res = [describe(total, "total (elected)"),
           describe(search, "search (oracle)"),
           describe(election, "election (selection gap)")]
    print(f"n = {len(rows)} paired targets   (positive delta = c1_on WORSE)\n")
    print(f"{'component':<28}{'median':>9}{'mean':>9}{'off>':>6}{'on>':>6}{'Wilcoxon':>11}{'sign':>10}")
    for x in res:
        print(f"{x['label']:<28}{x['median']:>+9.4f}{x['mean']:>+9.4f}"
              f"{x['off_better']:>6}{x['on_better']:>6}{x['wilcoxon_p']:>11.4f}{x['sign_p']:>10.5f}")

    print(f"\nselection gap   c1_on  mean {gap_on.mean():.3f}  median {np.median(gap_on):.3f} A")
    print(f"                c1_off mean {gap_off.mean():.3f}  median {np.median(gap_off):.3f} A")
    print(f"mean additivity: {search.mean():+.4f} + {election.mean():+.4f} "
          f"= {search.mean()+election.mean():+.4f}  vs total {total.mean():+.4f}")

    with open("wall_effect_decomposition.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["target", "total", "search", "election",
                    "gap_c1_on", "gap_c1_off", "n_c1_on", "n_c1_off"])
        for i, r in enumerate(rows):
            w.writerow([r["target"], f"{total[i]:.4f}", f"{search[i]:.4f}",
                        f"{election[i]:.4f}", f"{gap_on[i]:.4f}", f"{gap_off[i]:.4f}",
                        int(r["n_c1_on"]), int(r["n_c1_off"])])
    json.dump({"n": len(rows), "components": res,
               "gap_mean": {"c1_on": float(gap_on.mean()), "c1_off": float(gap_off.mean())},
               "gap_median": {"c1_on": float(np.median(gap_on)),
                              "c1_off": float(np.median(gap_off))}},
              open("handoff/wall_decomposition_facts.json", "w"), indent=1)
    print("\nwrote wall_effect_decomposition.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
