#!/usr/bin/env python3
"""Detect and classify targets claimed by more than one shard worker.

WHY THIS EXISTS
---------------
The replicate arm was launched twice with INCOMPATIBLE shard splits: shards 0,1
of a 2-way split (running since 00:23) and shards 2,3 of a 4-way split (added
01:33). The 4-way subsets are strict subsets of the 2-way ones -- shard 2's
targets are all in shard 0, shard 3's all in shard 1 -- so the later pair adds
no coverage, only the chance that two workers dock the same target into the same
directory at the same time.

The launcher's only guard is a CHECK-THEN-ACT at line 94:

    if [ -f "$d/result.csv" ]; then ... continue; fi
    mkdir -p "$d"

There is no atomic claim, so two workers that both see no result.csv both
proceed. `kill` is EPERM in this sandbox, so the extra workers cannot be
retracted; the only sound response is to make every collision VISIBLE and
repairable instead of silent. This project has already lost one cell that way
(c1_on/1G9V in the main campaign), and it was only caught because a pose count
disagreed.

WHAT COUNTS AS EVIDENCE
-----------------------
A target is CONTAMINATED if two different shard logs record a start for it with
no intervening completion. A target is CLEAN-DUPLICATE if the second worker's
log shows "skip <T> (done)", i.e. the resume gate fired as intended.

The ledger is advisory: it names cells to re-run, and never edits or deletes
campaign data.
"""

from __future__ import annotations

import argparse
import collections
import csv
import glob
import json
import os
import re
import sys

# The completion line carries the wall time, and WALL TIME is what separates a
# real docking run from a no-op pass. MEASURED on 1HNN: shard 1 logged 1063 s
# (a genuine 5-restart run, 7 write bursts on disk) and shard 3 logged 1 s for
# the same target with identical poses and evals. The engine refuses to re-dock
# into a populated output directory, so a duplicate worker CANNOT duplicate pose
# data -- it can only overwrite the wrapper's provenance fields (wall_s,
# evals_actual, finished_utc). Classifying on the log line alone called that
# "contaminated" and demanded a 20-minute re-dock; classifying on wall time
# calls it what it is, a receipt overwrite needing a receipt repair.
NOOP_SECONDS = 30

START = re.compile(r"^\s*\[(\d+)/(\d+)\]\s+([0-9A-Z]{4})\b.*?(\d+)s\s*$")
SKIP = re.compile(r"^\s*skip\s+([0-9A-Z]{4})\s+\(done\)")


def parse_logs(campaign):
    """shard -> ordered list of (kind, target) as the worker reported them."""
    per = {}
    for p in sorted(glob.glob(f"{campaign}/log_s*.log")):
        sh = os.path.basename(p)[5:-4]
        ev = []
        for ln in open(p, errors="replace"):
            m = START.match(ln.rstrip())
            if m:
                # printed AFTER the cell finishes; group(4) is its wall time
                ev.append(("done", m.group(3), int(m.group(4))))
                continue
            m = SKIP.match(ln)
            if m:
                ev.append(("skip", m.group(1), None))
        per[sh] = ev
    return per


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--campaign", required=True)
    ap.add_argument("--out", default="shard_collision_ledger.csv")
    a = ap.parse_args()

    per = parse_logs(a.campaign)
    if not per:
        print("VACUOUS: no shard logs found", file=sys.stderr)
        return 3

    real = collections.defaultdict(list)    # target -> [(shard, wall_s), ...] genuine runs
    noop = collections.defaultdict(list)    # target -> [(shard, wall_s), ...] engine-guarded passes
    skipped = collections.defaultdict(list)
    for sh, ev in per.items():
        for kind, t, wall in ev:
            if kind == "skip":
                skipped[t].append(sh)
            elif wall is not None and wall <= NOOP_SECONDS:
                noop[t].append((sh, wall))
            else:
                real[t].append((sh, wall))
    completed = {t: [s for s, _ in v] for t, v in real.items()}

    # a target whose directory exists but has no receipt is IN FLIGHT right now
    inflight = []
    for d in sorted(glob.glob(f"{a.campaign}/c1_on/*/")):
        t = os.path.basename(d.rstrip("/"))
        if len(t) != 4:
            continue
        if not os.path.exists(os.path.join(d, "CELL_PROVENANCE.json")):
            inflight.append(t)

    rows = []
    for t in sorted(set(real) | set(noop) | set(skipped) | set(inflight)):
        rs, ns, sk = real.get(t, []), noop.get(t, []), skipped.get(t, [])
        if len(rs) > 1:
            verdict = "CONTAMINATED: two genuine runs -- re-dock"
        elif rs and ns:
            verdict = "receipt-overwrite: poses intact, repair wall_s/evals"
        elif rs and sk:
            verdict = "clean: resume gate fired on the duplicate"
        elif t in inflight:
            verdict = "in-flight"
        elif rs:
            verdict = "clean"
        else:
            verdict = "no genuine run recorded"
        rows.append({"target": t,
                     "run_by": ",".join(f"{s}:{w}s" for s, w in rs),
                     "noop_by": ",".join(f"{s}:{w}s" for s, w in ns),
                     "skipped_by": ",".join(sk),
                     "in_flight": t in inflight, "verdict": verdict})

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    tally = collections.Counter(r["verdict"].split(":")[0] for r in rows)
    print(f"targets seen: {len(rows)}")
    for k, v in sorted(tally.items()):
        print(f"  {k:<40} {v}")
    bad = [r["target"] for r in rows if r["verdict"].startswith("CONTAMINATED")]
    ovw = [r["target"] for r in rows if r["verdict"].startswith("receipt-overwrite")]
    print(f"\nCONTAMINATED -- two genuine runs, re-dock: {len(bad)} {bad}")
    print(f"RECEIPT-OVERWRITE -- poses intact, repair fields: {len(ovw)} {ovw}")
    print(f"\nper-shard completed counts: "
          f"{ {sh: sum(1 for k, _, _ in ev if k == 'done') for sh, ev in per.items()} }")
    print(f"per-shard skip counts:      "
          f"{ {sh: sum(1 for k, _, _ in ev if k == 'skip') for sh, ev in per.items()} }")
    json.dump({"contaminated": bad, "receipt_overwrite": ovw,
               "tally": dict(tally), "n_targets": len(rows)},
              open("handoff/shard_collisions.json", "w"), indent=1)
    print(f"\nwrote {a.out}  (advisory only -- nothing was edited or deleted)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
