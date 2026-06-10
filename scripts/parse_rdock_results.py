#!/usr/bin/env python3
# =============================================================================
# parse_rdock_results.py — score rDock Astex output & build the comparison CSV
# =============================================================================
# Reads the per-complex rDock SD output produced by run_rdock_astex.sh, extracts
# the TOP-1 emitted pose (lowest rDock SCORE), computes the Hungarian RMSD vs the
# crystal ligand, applies the same sub-2 A success criterion used by FlexAIDdS,
# and writes:
#
#   rdock_astex_results.csv     rDock-only results (FlexAIDdS-shaped columns)
#   rdock_vs_flexaidds.csv      side-by-side head-to-head comparison
#
# Hungarian RMSD: optimal element-matched atom assignment (scipy
# linear_sum_assignment) — the same family of metric FlexAIDdS uses for its
# success flag, so the two engines are scored identically. No fitting/alignment
# is applied (both poses live in the crystal receptor frame), matching the
# benchmark's "is the pose in the right place" question.
#
# Pure-Python SDF parsing — no RDKit dependency (RDKit isn't installed here).
# =============================================================================
import argparse
import csv
import math
import os
import sys

import numpy as np
from scipy.optimize import linear_sum_assignment

DEFAULT_RDOCK_DIR = os.path.expanduser(
    "~/flexaidds_benchmark_results/rdock_astex")
DEFAULT_FLEXAIDDS_CSV = os.path.expanduser(
    "~/flexaidds_benchmark_results/astex_diverse_results.csv")
DEFAULT_ASTEX_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "benchmarks", "astex_diverse", "astex_diverse")
SUCCESS_THRESHOLD = 2.0  # Angstrom


# --------------------------------------------------------------------------- #
# Minimal V2000 SDF reader
# --------------------------------------------------------------------------- #
def _read_sdf_records(path):
    """Yield (atoms, props) for each record in a (possibly multi-) SDF file.

    atoms: list of (element:str, x, y, z)
    props: dict of SD tag -> first value line (e.g. 'SCORE' -> '-21.4')
    Hydrogens are dropped (heavy-atom RMSD, matching docking convention).
    """
    with open(path, "r", errors="replace") as fh:
        lines = fh.read().splitlines()

    i = 0
    n = len(lines)
    while i < n:
        # A record = header(3) + counts line + atom block + ... + '$$$$'
        if n - i < 4:
            break
        header = lines[i:i + 3]            # noqa: F841 (kept for clarity)
        counts = lines[i + 3]
        try:
            natoms = int(counts[0:3])
            # nbonds = int(counts[3:6])
        except ValueError:
            # Not a valid counts line — advance to next $$$$ and retry
            i = _skip_to_delim(lines, i) + 1
            continue

        atoms = []
        a0 = i + 4
        for a in range(natoms):
            ln = lines[a0 + a]
            try:
                x = float(ln[0:10]); y = float(ln[10:20]); z = float(ln[20:30])
            except (ValueError, IndexError):
                parts = ln.split()
                x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                elem = parts[3]
            else:
                elem = ln[31:34].strip()
            if elem.upper() != "H":
                atoms.append((elem, x, y, z))

        # Parse SD property tags ( > <TAG> \n value )
        props = {}
        j = a0 + natoms
        while j < n and lines[j].strip() != "$$$$":
            s = lines[j].strip()
            if s.startswith(">"):
                lb = s.find("<"); rb = s.find(">", lb + 1)
                if lb != -1 and rb != -1:
                    tag = s[lb + 1:rb]
                    val = lines[j + 1].strip() if j + 1 < n else ""
                    props.setdefault(tag, val)
            j += 1

        yield atoms, props
        i = j + 1  # skip the $$$$ delimiter


def _skip_to_delim(lines, i):
    n = len(lines)
    while i < n and lines[i].strip() != "$$$$":
        i += 1
    return i


# --------------------------------------------------------------------------- #
# Hungarian (element-matched) RMSD
# --------------------------------------------------------------------------- #
def hungarian_rmsd(atoms_a, atoms_b):
    """Optimal element-matched RMSD between two atom lists (no superposition).

    Atoms are matched only within the same element; the assignment minimising
    total squared distance is found per element via the Hungarian algorithm.
    Returns RMSD in Angstrom, or None if the heavy-atom composition differs.
    """
    if not atoms_a or not atoms_b or len(atoms_a) != len(atoms_b):
        return None

    by_elem_a, by_elem_b = {}, {}
    for e, x, y, z in atoms_a:
        by_elem_a.setdefault(e.upper(), []).append((x, y, z))
    for e, x, y, z in atoms_b:
        by_elem_b.setdefault(e.upper(), []).append((x, y, z))

    if set(by_elem_a) != set(by_elem_b):
        return None
    for e in by_elem_a:
        if len(by_elem_a[e]) != len(by_elem_b[e]):
            return None

    total_sq = 0.0
    total_n = 0
    for e in by_elem_a:
        A = np.asarray(by_elem_a[e], dtype=float)
        B = np.asarray(by_elem_b[e], dtype=float)
        # cost[i, j] = squared distance between A_i and B_j
        diff = A[:, None, :] - B[None, :, :]
        cost = np.einsum("ijk,ijk->ij", diff, diff)
        ri, ci = linear_sum_assignment(cost)
        total_sq += cost[ri, ci].sum()
        total_n += len(ri)

    if total_n == 0:
        return None
    return math.sqrt(total_sq / total_n)


# --------------------------------------------------------------------------- #
# Top-1 pose selection from rDock output
# --------------------------------------------------------------------------- #
def top1_pose(docked_sd):
    """Return (atoms, score, num_poses) for the best (lowest SCORE) rDock pose."""
    best_atoms, best_score = None, None
    n = 0
    for atoms, props in _read_sdf_records(docked_sd):
        n += 1
        # rDock total score tag is 'SCORE'; fall back to other common tags.
        raw = props.get("SCORE") or props.get("rDock.Score") or props.get("Score")
        try:
            score = float(raw) if raw is not None else None
        except ValueError:
            score = None
        if score is None:
            continue
        if best_score is None or score < best_score:
            best_score, best_atoms = score, atoms
    return best_atoms, best_score, n


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rdock-dir", default=DEFAULT_RDOCK_DIR,
                    help="Dir with <CODE>/docked.sd (default: %(default)s)")
    ap.add_argument("--astex-dir", default=DEFAULT_ASTEX_DIR,
                    help="Astex dataset dir with <CODE>/<CODE>_ligand.sdf")
    ap.add_argument("--flexaidds-csv", default=DEFAULT_FLEXAIDDS_CSV,
                    help="FlexAIDdS astex_diverse_results.csv for comparison")
    ap.add_argument("--out-dir", default=None,
                    help="Where to write output CSVs (default: --rdock-dir)")
    ap.add_argument("--threshold", type=float, default=SUCCESS_THRESHOLD,
                    help="Success RMSD cutoff in Angstrom (default: %(default)s)")
    args = ap.parse_args()

    out_dir = args.out_dir or args.rdock_dir
    os.makedirs(out_dir, exist_ok=True)

    if not os.path.isdir(args.rdock_dir):
        print(f"ERROR: rdock dir not found: {args.rdock_dir}", file=sys.stderr)
        print("Run scripts/run_rdock_astex.sh first.", file=sys.stderr)
        sys.exit(1)

    codes = sorted(d for d in os.listdir(args.astex_dir)
                   if os.path.isdir(os.path.join(args.astex_dir, d)))

    rows = []
    for code in codes:
        docked = os.path.join(args.rdock_dir, code, "docked.sd")
        crystal = os.path.join(args.astex_dir, code, f"{code}_ligand.sdf")

        score = rmsd = None
        num_poses = 0
        status = "ok"

        if not os.path.isfile(crystal):
            status = "no_crystal_ref"   # e.g. 1TW6 peptide
        elif not os.path.isfile(docked):
            status = "no_rdock_output"
        else:
            try:
                ref_atoms, _ = next(_read_sdf_records(crystal))
            except StopIteration:
                ref_atoms = None
            pose_atoms, score, num_poses = top1_pose(docked)
            if pose_atoms is None:
                status = "no_scored_pose"
            elif not ref_atoms:
                status = "bad_crystal_ref"
            else:
                rmsd = hungarian_rmsd(ref_atoms, pose_atoms)
                if rmsd is None:
                    status = "atom_mismatch"

        success = int(rmsd is not None and rmsd < args.threshold)
        rows.append({
            "pdb_id": code,
            "rdock_score": "" if score is None else f"{score:.4f}",
            "rmsd_to_crystal": "" if rmsd is None else f"{rmsd:.4f}",
            "num_poses": num_poses,
            "success": success,
            "status": status,
        })

    # ---- rDock-only CSV ----------------------------------------------------
    rdock_csv = os.path.join(out_dir, "rdock_astex_results.csv")
    fields = ["pdb_id", "rdock_score", "rmsd_to_crystal",
              "num_poses", "success", "status"]
    with open(rdock_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    n_eval = sum(1 for r in rows if r["rmsd_to_crystal"] != "")
    n_succ = sum(r["success"] for r in rows)
    print(f"rDock: {n_succ}/{n_eval} sub-{args.threshold}A "
          f"(top-1, Hungarian RMSD)  -> {rdock_csv}")

    # ---- Side-by-side comparison ------------------------------------------
    fa = {}
    if os.path.isfile(args.flexaidds_csv):
        with open(args.flexaidds_csv, newline="") as fh:
            for r in csv.DictReader(fh):
                fa[r["pdb_id"]] = r
    else:
        print(f"WARN: FlexAIDdS CSV not found ({args.flexaidds_csv}); "
              "comparison will have blank FlexAIDdS columns.", file=sys.stderr)

    cmp_csv = os.path.join(out_dir, "rdock_vs_flexaidds.csv")
    cmp_fields = ["pdb_id",
                  "flexaidds_rmsd", "flexaidds_success",
                  "rdock_rmsd", "rdock_success",
                  "winner", "rdock_status"]
    fa_succ_total = 0
    with open(cmp_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cmp_fields)
        w.writeheader()
        for r in rows:
            far = fa.get(r["pdb_id"], {})
            fa_rmsd = far.get("rmsd_to_crystal", "")
            # Recompute FlexAIDdS success from RMSD (its 'success' col means
            # "docking ran", not RMSD<2 — see project memory).
            try:
                fa_succ = int(float(fa_rmsd) < args.threshold) if fa_rmsd != "" else 0
            except ValueError:
                fa_succ = 0
            fa_succ_total += fa_succ
            rd_succ = r["success"]

            if rd_succ and not fa_succ:
                winner = "rDock"
            elif fa_succ and not rd_succ:
                winner = "FlexAIDdS"
            elif fa_succ and rd_succ:
                winner = "both"
            else:
                winner = "neither"

            w.writerow({
                "pdb_id": r["pdb_id"],
                "flexaidds_rmsd": fa_rmsd,
                "flexaidds_success": fa_succ,
                "rdock_rmsd": r["rmsd_to_crystal"],
                "rdock_success": rd_succ,
                "winner": winner,
                "rdock_status": r["status"],
            })

    print(f"FlexAIDdS: {fa_succ_total}/{len(rows)} sub-{args.threshold}A "
          f"(recomputed from RMSD)")
    print(f"Side-by-side -> {cmp_csv}")


if __name__ == "__main__":
    main()
