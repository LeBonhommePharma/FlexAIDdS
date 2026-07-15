#!/usr/bin/env python3
"""Parse classic FlexAID pilot outputs into a normalized result.csv row.

Reads ``*_r*_0.pdb`` (cluster rank 0) REMARK lines for CF and RMSD when
RMSDST was set. Elects best CF.app among restart rank-0 poses.

Does not claim PoseBusters or thermodynamic ΔG.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

RMSD_RE = re.compile(
    r"REMARK\s+([-+]?\d+\.?\d*)\s+RMSD to ref\. structure\s+\(symmetry corrected\)",
    re.I,
)
RMSD_NS_RE = re.compile(
    r"REMARK\s+([-+]?\d+\.?\d*)\s+RMSD to ref\. structure\s+\(no symmetry",
    re.I,
)
CF_RE = re.compile(r"REMARK\s+CF=([-+]?\d+\.?\d*)", re.I)
CF_APP_RE = re.compile(r"REMARK\s+CF\.app=([-+]?\d+\.?\d*)", re.I)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_pose_pdb(path: Path) -> Dict[str, Optional[float]]:
    text = path.read_text(errors="replace")
    out: Dict[str, Optional[float]] = {
        "rmsd_sym": None,
        "rmsd_nosym": None,
        "cf": None,
        "cf_app": None,
    }
    m = RMSD_RE.search(text)
    if m:
        out["rmsd_sym"] = float(m.group(1))
    m = RMSD_NS_RE.search(text)
    if m:
        out["rmsd_nosym"] = float(m.group(1))
    m = CF_RE.search(text)
    if m:
        out["cf"] = float(m.group(1))
    m = CF_APP_RE.search(text)
    if m:
        out["cf_app"] = float(m.group(1))
    return out


def collect_restart_poses(out_dir: Path, pdb: str) -> List[Tuple[int, Path, Dict]]:
    rows = []
    for rdir_pdb in sorted(out_dir.glob(f"{pdb}_r*_0.pdb")):
        name = rdir_pdb.name
        try:
            r = int(name.split("_r")[1].split("_")[0])
        except (IndexError, ValueError):
            r = -1
        rows.append((r, rdir_pdb, parse_pose_pdb(rdir_pdb)))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arm", required=True)
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--work-dir", type=Path, default=None)
    ap.add_argument("--matrix-md5", default="")
    ap.add_argument("--binary", type=Path, default=None)
    args = ap.parse_args()

    pdb = args.pdb.upper()
    poses = collect_restart_poses(args.out_dir, pdb)
    if not poses:
        for p in sorted(args.out_dir.glob("*_0.pdb")):
            poses.append((0, p, parse_pose_pdb(p)))

    best = None
    for r, path, meta in poses:
        cf = meta.get("cf_app") if meta.get("cf_app") is not None else meta.get("cf")
        if cf is None:
            continue
        if best is None or cf < best[0]:
            best = (cf, r, path, meta)

    rmsd_top1 = None
    score_top1 = None
    elected_path = ""
    restarts_finished = len(poses)
    rmsd_bcr = None
    if best:
        score_top1 = best[0]
        elected_path = str(best[2])
        m = best[3]
        rmsd_top1 = m.get("rmsd_sym") if m.get("rmsd_sym") is not None else m.get("rmsd_nosym")

    all_rmsds = []
    for p in args.out_dir.glob(f"{pdb}_r*_*.pdb"):
        if p.name.endswith("_INI.pdb"):
            continue
        meta = parse_pose_pdb(p)
        v = meta.get("rmsd_sym") if meta.get("rmsd_sym") is not None else meta.get("rmsd_nosym")
        if v is not None:
            all_rmsds.append(v)
    if all_rmsds:
        rmsd_bcr = min(all_rmsds)

    success_s1 = int(rmsd_top1 is not None and rmsd_top1 <= 2.0)
    success_s3 = int(rmsd_bcr is not None and rmsd_bcr <= 2.0)

    bin_sha = ""
    if args.binary and args.binary.is_file():
        bin_sha = sha256_file(args.binary.resolve())

    row = {
        "arm": args.arm,
        "engine_sha": bin_sha,
        "matrix_md5": args.matrix_md5,
        "pdb_id": pdb,
        "rmsd_top1": "" if rmsd_top1 is None else f"{rmsd_top1:.4f}",
        "rmsd_bcr": "" if rmsd_bcr is None else f"{rmsd_bcr:.4f}",
        "success_s1": success_s1,
        "success_s2": "",
        "success_s3": success_s3,
        "rank_native_mode": "",
        "n_poses": restarts_finished,
        "n_modes": restarts_finished,
        "score_top1": "" if score_top1 is None else f"{score_top1:.5f}",
        "H": "",
        "TS": "",
        "F": "",
        "pb_pass": "",
        "tencom_status": "NA",
        "seed_echo": 0,
        "native_pose_seeded": 0,
        "protocol_claim_eligible": 1 if args.matrix_md5 else 0,
        "wall_s": "",
        "restarts_finished": restarts_finished,
        "evals_actual": "",
        "budget_class": "full",
        "elected_path": elected_path,
        "parsed_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }

    out_csv = args.out_dir / "result.csv"
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(row.keys()))
        w.writeheader()
        w.writerow(row)
    print(f"wrote {out_csv} s1={success_s1} rmsd_top1={row['rmsd_top1']} bcr={row['rmsd_bcr']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
