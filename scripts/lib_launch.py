#!/usr/bin/env python3
# lib_launch.py — shared helpers for FlexAIDdS benchmark launch scripts
#
# Import from any launch_vN.py:
#   from lib_launch import priority_from_prev_run
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

import csv
import glob
import os


def priority_from_prev_run(prev_result_dir, lo=1.8, hi=2.5):
    """Return comma-separated PDB IDs with best_cluster_rmsd in (lo, hi) from a
    previous run — i.e., near-misses that almost made it under the 2 Å cutoff.

    Pass the result to FLEXAIDDS_PRIORITY_TARGETS so those targets run first in
    the next campaign and get fresh eyes early, before the worker pool disperses
    onto easier wins.

    Args:
        prev_result_dir: root output dir of a prior run (e.g.
            ~/flexaidds_results/v41_20260613_zshannonselect).
        lo: RMSD lower bound (exclusive). Default 1.8 Å.
        hi: RMSD upper bound (exclusive). Default 2.5 Å.

    Returns:
        Comma-separated string of PDB IDs, e.g. "1HP0,1Q4G,1Q41".
        Empty string if no near-misses found or dir doesn't exist.

    Example::
        prio = priority_from_prev_run("~/flexaidds_results/v41_...", lo=1.8, hi=2.5)
        if prio:
            env["FLEXAIDDS_PRIORITY_TARGETS"] = prio
    """
    prev_result_dir = os.path.expanduser(prev_result_dir)
    if not os.path.isdir(prev_result_dir):
        return ""

    near_misses = []
    for csv_path in glob.glob(f"{prev_result_dir}/*/result.csv"):
        pdb_id = os.path.basename(os.path.dirname(csv_path))
        try:
            with open(csv_path, newline="") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    try:
                        rmsd = float(
                            row.get("best_cluster_rmsd")
                            or row.get("rmsd_hungarian")
                            or 999
                        )
                        if lo < rmsd < hi:
                            near_misses.append(pdb_id)
                            break  # one hit per target is enough
                    except (ValueError, TypeError):
                        pass
        except OSError:
            pass

    return ",".join(near_misses)
