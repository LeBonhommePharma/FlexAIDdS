#!/usr/bin/env python3
"""Recover benchmark targets voided by the shared-_dockin.sdf restart race.

A target is *voided* when one restart dies reading the shared
``run/<PDB>/<PDB>_dockin.sdf`` (torn read -> exit code 2). ``ret`` aggregates
any-failure-wins (DatasetRunner.cpp:6440/6467), so ``docking_completed`` goes
false (:6830) and the whole crystal-RMSD + elected-pose-persist block at :6866
is skipped -- even though the surviving restarts produced poses and the pooled
election already chose a winner and logged it as ``[3DSIG-RANK] rank=0``.

The docking is not lost. This script re-derives, read-only, what the run would
have recorded: it reads the elected pose back out of claim.log and scores it
against the crystal ligand.

Why -1 must never be left in the table: ``-1.0 < 2.0``, so a naive
``rmsd < 2`` success filter counts a voided target as a *success*.

Read-only. Safe to run against a live campaign.

Usage:
    python3 scripts/recover_voided_targets.py <campaign_dir> [--cache DIR] [--csv OUT]
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import re
import sys

RANK_RE = re.compile(
    r"\[3DSIG-RANK\] rank=(\d+) .*? cf=(\S+) freq=\d+ nmembers=(\d+) path=(\S+)"
)
SENTINEL = -1.0
SUCCESS_A = 2.0


def read_sdf_heavy(path):
    """Heavy atoms from a V2000 SDF, in file order: [(x, y, z, element)]."""
    lines = open(path).read().splitlines()
    n = int(lines[3][0:3])
    out = []
    for ln in lines[4 : 4 + n]:
        p = ln.split()
        el = p[3]
        if el.upper() != "H":
            out.append((float(p[0]), float(p[1]), float(p[2]), el))
    return out


def read_pose_ligand(path, resname):
    """Heavy ligand atoms from a pose PDB, in file order."""
    out = []
    for ln in open(path, errors="replace"):
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        if ln[17:20].strip() != resname:
            continue
        el = (ln[76:78].strip() or ln[12:16].strip()[:1]).upper()
        if el == "H":
            continue
        out.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54]), el))
    return out


def ordered_rmsd(pose, crystal):
    """Ordered positional RMSD. Returns (rmsd, note); rmsd is None if unusable.

    Pose and crystal encode the same molecule in the same atom order, so no
    alignment or matching is applied -- same convention as
    compute_pose_ligand_rmsd() in the engine. Element agreement is reported as
    the check that the orders really did correspond.
    """
    if len(pose) != len(crystal):
        return None, f"atom count {len(pose)} != crystal {len(crystal)}"
    mismatch = sum(1 for a, b in zip(pose, crystal) if a[3] != b[3].upper())
    ssd = sum(
        (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2
        for a, b in zip(pose, crystal)
    )
    note = "ok" if mismatch == 0 else f"{mismatch} element mismatches"
    return math.sqrt(ssd / len(pose)), note


def ligand_resname(sdf_path):
    """SDF title line carries the residue name the engine writes into the pose."""
    with open(sdf_path) as fh:
        return fh.readline().strip()


def elected_from_log(claim_log, pdb_id):
    """The rank-0 pose the pooled election chose, from the claim.log trace.

    Restricted to paths under run/<pdb_id>/ so interleaved targets in a
    parallel log cannot bleed across.
    """
    best = None
    needle = f"/run/{pdb_id}/"
    with open(claim_log, errors="replace") as fh:
        for line in fh:
            m = RANK_RE.search(line)
            if not m or needle not in m.group(4):
                continue
            rank, cf, nmem, path = int(m.group(1)), m.group(2), int(m.group(3)), m.group(4)
            if rank == 0:
                best = (path, float(cf), nmem)
    return best


def failed_restarts(target_dir):
    """Restarts that died on the torn ligand read."""
    out = []
    for err in sorted(glob.glob(os.path.join(target_dir, "r*", "stderr.log"))):
        try:
            body = open(err, errors="replace").read()
        except OSError:
            continue
        if "truncated at line" in body or "Failed to read ligand file" in body:
            out.append(os.path.basename(os.path.dirname(err)))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("campaign", help="campaign dir containing run/ and claim.log")
    ap.add_argument(
        "--cache",
        default="/Users/lp.more/flexaidds_results/cache_v2/astex_diverse",
        help="dataset cache holding <PDB>/<PDB>_ligand.sdf",
    )
    ap.add_argument("--csv", help="write the recovered table here")
    args = ap.parse_args()

    run_dir = os.path.join(args.campaign, "run")
    claim_log = os.path.join(args.campaign, "claim.log")
    if not os.path.isdir(run_dir):
        print(f"no run/ under {args.campaign}", file=sys.stderr)
        return 2
    if not os.path.exists(claim_log):
        print(f"no claim.log under {args.campaign}", file=sys.stderr)
        return 2

    rows = []
    for path in sorted(glob.glob(os.path.join(run_dir, "*", "result.csv"))):
        with open(path) as fh:
            for r in csv.DictReader(fh):
                r["_dir"] = os.path.dirname(path)
                rows.append(r)

    recovered, voided_unrecovered = [], []
    table = []
    for r in rows:
        pdb = r["pdb_id"]
        try:
            rmsd = float(r.get("rmsd_to_crystal") or "nan")
        except ValueError:
            rmsd = float("nan")
        voided = rmsd == SENTINEL or not (r.get("elected_pose_path") or "")

        rec_rmsd, source, note = None, "", ""
        if voided:
            elected = elected_from_log(claim_log, pdb)
            sdf = os.path.join(args.cache, pdb, f"{pdb}_ligand.sdf")
            if elected and os.path.exists(elected[0]) and os.path.exists(sdf):
                crystal = read_sdf_heavy(sdf)
                pose = read_pose_ligand(elected[0], ligand_resname(sdf))
                rec_rmsd, note = ordered_rmsd(pose, crystal)
                source = elected[0]
            if rec_rmsd is None:
                voided_unrecovered.append(pdb)
                note = note or "no rank-0 pose in claim.log"
            else:
                recovered.append(pdb)

        eff = rec_rmsd if rec_rmsd is not None else (rmsd if rmsd >= 0 else None)
        table.append(
            {
                "pdb_id": pdb,
                "recorded_rmsd": f"{rmsd:.4f}" if rmsd == rmsd else "NA",
                "voided": int(voided),
                "failed_restarts": ",".join(failed_restarts(r["_dir"])) if voided else "",
                "recovered_rmsd": f"{rec_rmsd:.4f}" if rec_rmsd is not None else "",
                "effective_rmsd": f"{eff:.4f}" if eff is not None else "NA",
                "success_lt2": "" if eff is None else int(eff < SUCCESS_A),
                "recovered_from": source,
                "note": note,
            }
        )

    print(f"campaign      : {args.campaign}")
    print(f"targets       : {len(rows)}")
    print(f"voided by race: {len(recovered) + len(voided_unrecovered)}"
          f"  (recovered {len(recovered)}, unrecoverable {len(voided_unrecovered)})")
    if recovered:
        print()
        print(f"{'pdb':6} {'recorded':>9} {'recovered':>10} {'restarts lost':>14}  note")
        for t in table:
            if t["recovered_rmsd"]:
                print(f"{t['pdb_id']:6} {t['recorded_rmsd']:>9} {t['recovered_rmsd']:>10} "
                      f"{t['failed_restarts']:>14}  {t['note']}")
    if voided_unrecovered:
        print(f"\nUNRECOVERABLE (exclude and say so): {', '.join(voided_unrecovered)}")

    scored = [t for t in table if t["success_lt2"] != ""]
    hits = sum(int(t["success_lt2"]) for t in scored)
    n = len(scored)
    naive = sum(
        1 for t in table if t["recorded_rmsd"] != "NA" and float(t["recorded_rmsd"]) < SUCCESS_A
    )
    print()
    if n:
        p = hits / n
        se = math.sqrt(p * (1 - p) / n)
        print(f"top-1 <2A after recovery : {hits}/{n} = {p:.3f} +/- {se:.3f} (1 SE)")
    print(f"top-1 <2A naive 'rmsd<2.0': {naive}/{len(table)} "
          f"  <-- counts every -1.0 sentinel as a SUCCESS")
    print("\nThis is top-1. The 2017 3Dsig baseline is top-10 with a 10k bootstrap "
          "(scripts/bootstrap_3dsig_s_top10.py) -- not comparable as printed.")

    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(table[0].keys()))
            w.writeheader()
            w.writerows(table)
        print(f"\nwrote {args.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
