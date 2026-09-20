#!/usr/bin/env python3
"""Recover evals_actual into campaign cell provenance from the engine's own logs.

WHY THIS IS NEEDED
------------------
The paired wall campaign's streams were launched from a launcher version that
grepped the WRAPPER's run.log for the engine's `[EVALS] evals_actual=N` line.
The engine does not write it there: it prints it on the normal exit path of each
RESTART process, to that process's own stdout.log. So every cell recorded
`evals_actual: "unreported"` while the numbers sat one directory down.

The per-target total is the SUM across restart processes, because each restart is
a separate process with its own counter. Measured on 1G9V at pop 1000 / gen 2000
/ restarts 5: 3,524,539 + 3,515,610 + 3,502,372 + 3,522,625 + 3,522,621.

WHY RUN IT BEFORE THE CAMPAIGN ENDS
-----------------------------------
The stdout.log files live inside the per-target directories. Once those are
archived and removed, a provenance file saying "unreported" is the only record,
and the witness is gone. This is idempotent and read-only apart from the
provenance JSON it repairs, so it can be run repeatedly as cells land.

WHAT IT DOES NOT DO
-------------------
It never derives a value. If no restart process reported a count, the field is
set to null and `evals_procs_reporting` to 0 — never to a pop x gen x restarts
product. That derivation is unsound on this engine: the generation loop can exit
early on an operator STOP file without marking the restart short, and the
population is scaled by ligand flexibility (pop_effective), so the nominal grid
is not the search actually performed.

Usage:
  python3 backfill_evals_actual.py --campaign <dir> [--apply]
Dry-run by default.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys

EVALS = re.compile(rb"evals_actual=(\d+)")
SCALE = re.compile(rb"pop_base=(\d+) pop_effective=(\d+) n_gen=(\d+)")
NFLEX = re.compile(rb"n_flex_bonds=(\d+)")


def harvest(cell_dir: str) -> dict:
    """Sum [EVALS] across every restart process's stdout.log under cell_dir."""
    vals, files = [], 0
    for p in glob.glob(f"{cell_dir}/**/stdout.log", recursive=True):
        try:
            blob = open(p, "rb").read()
        except OSError:
            continue
        m = EVALS.findall(blob)
        if m:
            files += 1
            vals.append(int(m[-1]))          # last line per process
    out = {"evals_actual": sum(vals) if vals else None,
           "evals_procs_reporting": files}
    # pop_effective is the correct denominator for any ratio; it appears in the
    # wrapper log and in NO receipt field, which is why the counter exists.
    for p in (f"{cell_dir}/run.log",):
        if not os.path.exists(p):
            continue
        blob = open(p, "rb").read()
        s = SCALE.search(blob)
        if s:
            out["pop_base"] = int(s.group(1))
            out["pop_effective"] = int(s.group(2))
            out["n_gen"] = int(s.group(3))
        f = NFLEX.search(blob)
        if f:
            out["n_flex_bonds"] = int(f.group(1))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--campaign", required=True)
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args(argv)

    cells = sorted(glob.glob(f"{a.campaign}/*/*/CELL_PROVENANCE.json"))
    if not cells:
        print("VACUOUS: no CELL_PROVENANCE.json under --campaign", file=sys.stderr)
        return 3

    fixed = already = nodata = 0
    ratios = []
    for p in cells:
        d = json.load(open(p))
        cur = d.get("evals_actual")
        if isinstance(cur, int) and cur > 0:
            already += 1
            continue
        h = harvest(os.path.dirname(p))
        if not h["evals_actual"]:
            nodata += 1
            print(f"  {d.get('arm','?')}/{d.get('target','?'):<6} NO DATA "
                  f"({h['evals_procs_reporting']} procs reporting)")
            continue
        fixed += 1
        if h.get("pop_effective") and h.get("n_gen"):
            nominal = h["pop_effective"] * h["n_gen"] * int(d.get("restarts") or 5)
            ratios.append(h["evals_actual"] / nominal)
        if a.apply:
            d.update(h)
            d["evals_backfilled"] = True
            d["evals_backfill_note"] = (
                "summed [EVALS] across per-restart stdout.log; the launcher that produced "
                "this cell grepped the wrapper run.log, which never carries the line")
            json.dump(d, open(p, "w"), indent=1)

    print(f"\ncells {len(cells)}   already populated {already}   "
          f"{'backfilled' if a.apply else 'would backfill'} {fixed}   no data {nodata}")
    if ratios:
        rs = sorted(ratios)
        print(f"  measured / (pop_effective x n_gen x restarts): "
              f"median {rs[len(rs)//2]:.3f}  min {rs[0]:.3f}  max {rs[-1]:.3f}  n={len(rs)}")
        print("  (a median near 1.00 confirms pop_effective is the right denominator;")
        print("   pop_base would give ~1.25-1.75 depending on ligand flexibility)")
    if not a.apply:
        print("\nDRY RUN -- pass --apply to write.")
    if fixed == 0 and already == 0:
        print("VACUOUS: nothing populated and nothing backfillable", file=sys.stderr)
        return 3
    return 0


if __name__ == "__main__":
    sys.exit(main())
