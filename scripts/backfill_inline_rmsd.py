#!/usr/bin/env python3
"""backfill_inline_rmsd.py — recompute -1 RMSD sentinels from pose archives.

Bug 2026-08-22 context: campaign arms 8/9/10 wrote astex_diverse_results.csv
with rmsd_to_crystal = -1 on most/all rows while their pose libraries were
complete and valid (crystal-reference resolution failed at runtime; see
workorders/BUGREPORT_inline_rmsd_minus1_sentinel.md). Their derived summaries
certified 0% success for runs whose offline pooled ceilings were 37-42/85.

This tool writes a SIDE-BY-SIDE backfilled copy (never modifies originals):

    <run-dir>/<name>_backfilled.csv      original columns, -1 cells repaired
    <run-dir>/<name>_backfill_receipt.json   method + provenance + verify stats

Repaired columns (measurement only; success flags are NEVER fabricated):
    rmsd_to_crystal   rank-0/elected pose, ordered direct RMSD (fail-closed)
    rmsd_hungarian    symmetry-corrected diagnostic RMSD, element-grouped
    best_cluster_rmsd oracle best-of-pooled-heads (min ordered RMSD over heads)
    cf_top1_rmsd      RMSD of the CF-best head (when its path column resolves)

Semantics mirror LIB/DatasetRunner.cpp compute_pose_ligand_rmsd exactly:
  - crystal reference parsed from the SDF atom block (heavy atoms, elements
    normalised upper-first; 'M'/'A' pseudo-elements skipped, as in production);
  - docked ligand selected by CONECT serial membership (fallback: HETATM with
    resSeq==1 when no CONECT records exist), hydrogens dropped, serial-sorted;
  - atom-count mismatch => fail-closed -1 (count_mismatch); ordered RMSD
    additionally requires per-index element identity (elem_order_mismatch);
  - Hungarian groups by element with optimal assignment (scipy).

Self-verification is mandatory before any output is written: rows whose
original RMSD was valid are recomputed and must agree within tolerance
(default 0.02 A) at >=95%, otherwise the tool aborts unless --force.

Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

WHOLESALE_REASONS = {"ref_empty", "input_missing", "pose_block_empty"}
DEFAULT_TOL = 0.02
MIN_AGREEMENT = 0.95


# ────────────────────────── crystal reference (SDF) ──────────────────────────

def _norm_elem(tok: str) -> str:
    """Upper-first, rest-lower — mirrors the C++ element normalisation."""
    if not tok:
        return "X"
    return tok[0].upper() + tok[1:].lower()


def parse_crystal_sdf(path: Path):
    """Mirror the DatasetRunner crystal-parse loop.

    Scans lines >39 chars for '<x> <y> <z> <elem>' with an alphabetic element
    that is not 'M'/'A' (pseudo-element guard) and not 'H' (heavy atoms only);
    stops at the first 'M  END' once atoms were seen.
    """
    xyz: list[list[float]] = []
    elem: list[str] = []
    with open(path, errors="replace") as fh:
        for line in fh:
            if len(line.rstrip("\n")) <= 39:
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
            except ValueError:
                continue
            e = parts[3]
            e0 = e[0] if e else ""
            if not e0.isalpha() or e0 in ("M", "A") or e0 == "H":
                continue
            xyz.append([x, y, z])
            elem.append(_norm_elem(e))
            if "M  END" in line:
                break
    return np.asarray(xyz, dtype=float), elem


# ──────────────────────────── pose ligand (PDB) ──────────────────────────────

def _parse_xyz_span(line: str):
    """Strict finite PDB coordinate decode (cols 31-38 / 39-46 / 47-54)."""
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except ValueError:
        return None
    if not (math.isfinite(x) and math.isfinite(y) and math.isfinite(z)):
        return None
    return x, y, z


def parse_pose_ligand(path: Path):
    """Mirror CONECT-membership selection; returns (serial-sorted) arrays."""
    conect: set[int] = set()
    with open(path, errors="replace") as fh:
        raw = fh.readlines()
    for line in raw:
        if line.startswith("CONECT"):
            tail = line[6:]
            for i in range(0, len(tail), 5):
                tok = tail[i:i + 5].strip()
                try:
                    conect.add(int(tok))
                except ValueError:
                    pass
    picked = []  # (serial, xyz, elem)
    for line in raw:
        rec = line[:6]
        if rec not in ("HETATM", "ATOM  "):
            continue
        if len(line.rstrip("\n")) < 54:
            continue
        try:
            serial = int(line[6:11].strip())
        except ValueError:
            continue
        if conect:
            selected = serial in conect
        else:
            resname = line[17:20].strip()
            try:
                resseq = int(line[22:26].strip())
            except ValueError:
                resseq = 0
            selected = rec == "HETATM" and resseq == 1 and len(resname) >= 2
        if not selected:
            continue
        tok = _norm_elem(line.rstrip().split()[-1] if line.split() else "X")
        if tok.lower() in ("h", "d", "du"):
            continue
        xyz = _parse_xyz_span(line)
        if xyz is None:
            continue
        picked.append((serial, xyz, tok))
    picked.sort(key=lambda t: t[0])
    if not picked:
        return None, None, "pose_block_empty"
    xyz = np.asarray([p[1] for p in picked], dtype=float)
    elem = [p[2] for p in picked]
    return xyz, elem, "none"


# ─────────────────────────────── RMSD kernels ────────────────────────────────

def _hungarian(cr_xyz, cr_elem, po_xyz, po_elem) -> float:
    """Element-grouped optimal-assignment RMSD (mirrors hungarian_rmsd)."""
    if cr_xyz.size == 0 or po_xyz.size == 0:
        return -1.0
    from scipy.optimize import linear_sum_assignment

    total_sq = 0.0
    total_n = 0
    for el in sorted(set(cr_elem)):
        ci = [i for i, e in enumerate(cr_elem) if e == el]
        di = [i for i, e in enumerate(po_elem) if e == el]
        n = min(len(ci), len(di))
        if n == 0:
            continue
        # C++ takes the FIRST n indices of each list (prefix truncation).
        c = np.zeros((n, n))
        for a in range(n):
            dxyz = po_xyz[di[a]]
            for b in range(n):
                cxyz = cr_xyz[ci[b]]
                c[a, b] = ((dxyz[0] - cxyz[0]) ** 2 +
                           (dxyz[1] - cxyz[1]) ** 2 +
                           (dxyz[2] - cxyz[2]) ** 2)
        rows, cols = linear_sum_assignment(c)
        total_sq += float(c[rows, cols].sum())
        total_n += n
    return math.sqrt(total_sq / total_n) if total_n > 0 else -1.0


def pose_ligand_rmsd(pose_path: Path, cr_xyz, cr_elem):
    """Mirror compute_pose_ligand_rmsd: {serial, hungarian, fail_reason}."""
    if cr_xyz.size == 0:
        return -1.0, -1.0, "ref_empty"
    po_xyz, po_elem, reason = parse_pose_ligand(pose_path)
    if reason != "none":
        return -1.0, -1.0, reason
    if cr_xyz.shape[0] != po_xyz.shape[0]:
        return -1.0, -1.0, "count_mismatch"
    if len(cr_elem) != cr_xyz.shape[0] or len(po_elem) != po_xyz.shape[0]:
        return -1.0, -1.0, "elem_mismatch"
    rc = -1.0
    reason = "none"
    if po_elem == cr_elem:
        d = po_xyz - cr_xyz
        rc = float(math.sqrt(float((d * d).sum()) / cr_xyz.shape[0]))
    else:
        reason = "elem_order_mismatch"
    rh = _hungarian(cr_xyz, cr_elem, po_xyz, po_elem)
    return rc, rh, reason


# ────────────────────────────────── driver ───────────────────────────────────

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve_cache(cache: Path, dataset: str) -> Path:
    ds = cache / dataset
    return ds if ds.is_dir() else cache


def crystal_for(cache_ds: Path, pdb_id: str):
    for name in (f"{pdb_id}_ligand.sdf",):
        p = cache_ds / pdb_id / name
        if p.is_file():
            return p
    return None


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def _relocate(p: str, run_dir: Path, pdb_id: str):
    """Recorded paths may be absolute paths from a since-moved tree (e.g.
    archive extractions); try the recorded location, then inside this run."""
    if not p:
        return None
    q = Path(p)
    if q.is_file():
        return q
    for alt in (run_dir / pdb_id / q.name,
                run_dir.parent / q.name,
                run_dir / q.name):
        if alt.is_file():
            return alt
    return None


def rank0_pose(run_dir: Path, pdb_id: str, row: dict):
    """Resolve the pose production used for authoritative rmsd_to_crystal.

    Resolution order:
      1. per-complex <ID>/result.csv 'elected_pose_path' (writer #1 schema)
      2. score-matched election among recorded top-1 heads: the elected pose is
         whichever of cf_top1 / entropy_top1 reproduces best_score (SoftBeta
         arms elect by entropy, so cf_top1 alone is NOT authoritative there)
      3. otherwise unresolved — an honest gap beats a wrong-pose fill
    Returns (path | None, source_label).
    """
    per = run_dir / pdb_id / "result.csv"
    if per.is_file():
        try:
            with open(per, newline="") as fh:
                r = csv.DictReader(fh).fieldnames or []
            if "elected_pose_path" in r:
                with open(per, newline="") as fh:
                    prow = next(iter(csv.DictReader(fh)), None)
                cand = _relocate((prow or {}).get("elected_pose_path", ""),
                                 run_dir, pdb_id)
                if cand is not None:
                    return cand, "elected_result_csv"
        except OSError:
            pass

    best = _f(row.get("best_score"))
    cands = []
    for score_col, path_col in (("cf_top1_score", "cf_top1_pose_path"),
                                ("entropy_top1_score", "entropy_top1_pose_path")):
        p = _relocate(row.get(path_col, ""), run_dir, pdb_id)
        s = _f(row.get(score_col))
        if p is not None and math.isfinite(best) and math.isfinite(s):
            cands.append((abs(best - s), p))
    if cands:
        delta, path = min(cands, key=lambda t: t[0])
        if delta <= 1e-3 * max(1.0, abs(best)):
            return path, "elected_via_scores"
        return None, "election_unresolvable"
    return None, "missing"


def backfill_run(run_dir: Path, cache_ds: Path, tol: float,
                 force: bool, verify_only: bool) -> int:
    src = run_dir / "astex_diverse_results.csv"
    if not src.is_file():
        print(f"[skip] no results.csv in {run_dir}")
        return 0
    with open(src, newline="") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])

    out_name = src.stem + "_backfilled.csv"
    out_path = run_dir / out_name
    receipt_path = run_dir / (src.stem + "_backfill_receipt.json")

    neg_cols = ("rmsd_to_crystal", "rmsd_hungarian",
                "best_cluster_rmsd", "cf_top1_rmsd")
    for col in ("backfill_note",):
        if col not in fieldnames:
            fieldnames.append(col)

    crystals: dict[str, tuple] = {}
    verify_ok = verify_bad = verify_skipped = 0
    max_delta = 0.0
    fixed = 0
    still_neg: dict[str, int] = {}

    n_valid_expected = sum(
        1 for row in rows
        if _f(row.get("rmsd_to_crystal")) >= 0 and row.get("pdb_id", "").strip())
    n_neg_expected = sum(
        1 for row in rows
        if _f(row.get("rmsd_to_crystal")) < 0 and row.get("pdb_id", "").strip())

    for row in rows:
        pdb_id = row.get("pdb_id", "").strip()
        note_parts = []

        def get(col):
            v = row.get(col, "")
            try:
                return float(v)
            except (TypeError, ValueError):
                return float("nan")

        orig_rank0 = get("rmsd_to_crystal")

        # Mandatory self-verification on originally-valid rows.
        if orig_rank0 >= 0 and pdb_id:
            if pdb_id not in crystals:
                sdf = crystal_for(cache_ds, pdb_id)
                crystals[pdb_id] = parse_crystal_sdf(sdf) if sdf else (np.zeros((0,)), [])
            cr_xyz, cr_elem = crystals[pdb_id]
            pose, how = rank0_pose(run_dir, pdb_id, row)
            if pose is None:
                verify_skipped += 1
                continue
            if cr_xyz.size:
                rc, _, reason = pose_ligand_rmsd(pose, cr_xyz, cr_elem)
                if reason == "none":
                    d = abs(rc - orig_rank0)
                    max_delta = max(max_delta, d)
                    if d <= tol:
                        verify_ok += 1
                    else:
                        verify_bad += 1
            continue  # valid rows pass through untouched below

        if not (orig_rank0 < 0) or not pdb_id:
            row["backfill_note"] = row.get("backfill_note", "")
            continue

        if pdb_id not in crystals:
            sdf = crystal_for(cache_ds, pdb_id)
            if sdf is None:
                still_neg["ref_sdf_missing"] = still_neg.get("ref_sdf_missing", 0) + 1
                row["backfill_note"] = "unfixable:ref_sdf_missing"
                continue
            crystals[pdb_id] = parse_crystal_sdf(sdf)
        cr_xyz, cr_elem = crystals[pdb_id]

        pose, how = rank0_pose(run_dir, pdb_id, row)
        if pose is not None:
            rc, rh, reason = pose_ligand_rmsd(pose, cr_xyz, cr_elem)
            if reason == "none" or reason == "elem_order_mismatch":
                row["rmsd_to_crystal"] = f"{rc:.4f}"
                row["rmsd_hungarian"] = f"{rh:.4f}"
                fixed += 1
                note_parts.append(f"rank0:{how}")
            else:
                still_neg[reason] = still_neg.get(reason, 0) + 1
                note_parts.append(f"unfixable:{reason}")
        else:
            still_neg[how] = still_neg.get(how, 0) + 1
            note_parts.append(f"rank0:{how}")

        # Oracle best-of-pooled-heads ceiling (ordered-direct, fail-closed).
        # Election-independent: fillable even when rank-0 identity is lost.
        heads = sorted((run_dir / pdb_id).glob(f"{pdb_id}_*.pdb"))
        if heads:
            best = float("inf")
            for head in heads:
                hx, he, hr = parse_pose_ligand(head)
                if hr != "none":
                    continue
                if hx.shape[0] != cr_xyz.shape[0]:
                    continue
                if he != cr_elem:
                    continue
                d = hx - cr_xyz
                r = math.sqrt(float((d * d).sum()) / cr_xyz.shape[0])
                best = min(best, r)
            if math.isfinite(best):
                row["best_cluster_rmsd"] = f"{best:.4f}"
                note_parts.append(f"ceiling_heads={len(heads)}")

        cfp = row.get("cf_top1_pose_path", "")
        if get("cf_top1_rmsd") < 0 and cfp and Path(cfp).is_file():
            crc, _, crr = pose_ligand_rmsd(Path(cfp), cr_xyz, cr_elem)
            if crr == "none":
                row["cf_top1_rmsd"] = f"{crc:.4f}"
                note_parts.append("cftop1")

        row["backfill_note"] = "|".join(note_parts)

    total_valid = verify_ok + verify_bad
    agreement = (verify_ok / total_valid) if total_valid else 0.0
    print(f"  self-check: {verify_ok}/{total_valid} valid rows agree "
          f"(max |delta| = {max_delta:.4f} A, skipped-unresolvable: "
          f"{verify_skipped})")
    # A silent 0/N self-check on a file that HAS valid rows is itself a
    # failure mode: never let "nothing checked" read as "everything agrees".
    if n_valid_expected and total_valid == 0:
        print("  ABORT: could not re-verify any originally-valid row "
              "(pose artifacts unresolvable) — use --force to override")
        return 2
    if total_valid and agreement < MIN_AGREEMENT and not force:
        print(f"  ABORT: agreement {agreement:.3f} < {MIN_AGREEMENT} "
              f"(use --force to override)")
        return 2

    if verify_only:
        print(f"  [verify-only] would repair {fixed} rows; "
              f"still-negative: {still_neg or '{}'}")
        return 0

    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    receipt = {
        "tool": "scripts/backfill_inline_rmsd.py",
        "bug_ref": "workorders/BUGREPORT_inline_rmsd_minus1_sentinel.md",
        "source_csv_sha256": sha256_file(src),
        "output_csv": out_path.name,
        "rows_total": len(rows),
        "rows_repaired": fixed,
        "still_negative_by_reason": still_neg,
        "self_check": {
            "valid_rows_in_csv": n_valid_expected,
            "negative_rows_in_csv": n_neg_expected,
            "valid_rows_checked": total_valid,
            "verify_skipped_unresolvable": verify_skipped,
            "agreement_rate": round(agreement, 5),
            "max_abs_delta_A": round(max_delta, 5),
            "tolerance_A": tol,
        },
        "semantics": "mirrors LIB/DatasetRunner.cpp compute_pose_ligand_rmsd",
        "note": "side-by-side artifact; success flags intentionally untouched",
    }
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"  wrote {out_path.name} ({fixed} rows repaired) "
          f"+ {receipt_path.name}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", action="append", required=True,
                    help="run directory holding astex_diverse_results.csv "
                         "(repeatable)")
    ap.add_argument("--cache", required=True,
                    help="dataset cache root (parent containing the dataset "
                         "dir, or the dataset dir itself)")
    ap.add_argument("--dataset", default="astex_diverse")
    ap.add_argument("--tol", type=float, default=DEFAULT_TOL,
                    help=f"self-check tolerance in A (default {DEFAULT_TOL})")
    ap.add_argument("--verify-only", action="store_true",
                    help="run the self-check, write nothing")
    ap.add_argument("--force", action="store_true",
                    help="write outputs even if the self-check underperforms")
    args = ap.parse_args(argv)

    cache_ds = resolve_cache(Path(args.cache).expanduser(), args.dataset)
    if not cache_ds.is_dir():
        print(f"error: cache dir not found: {cache_ds}", file=sys.stderr)
        return 1

    rc = 0
    for rd in args.run_dir:
        rd = Path(rd).expanduser()
        print(f"[{rd}]")
        rc = max(rc, backfill_run(rd, cache_ds, args.tol, args.force,
                                  args.verify_only))
    return rc


if __name__ == "__main__":
    sys.exit(main())
