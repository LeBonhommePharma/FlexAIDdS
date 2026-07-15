#!/usr/bin/env python3
"""Parse classic FlexAID pilot outputs into a normalized result.csv row.

Reads ranked pose PDBs::

  {pdb}_r{restart}_{rank}.pdb
  {pdb}_r{restart}_{minPts}_{rank}.pdb   # FO dual-suffix (arm B)

REMARK lines for CF and RMSD when RMSDST was set. Elects best CF.app among
restart rank-0 poses for S1; min RMSD among all poses for BCR (success_s3);
min RMSD among top-10 ranked modes for S_top10 (3Dsig primary).

Emits ``mode_rmsd_0`` … ``mode_rmsd_9`` for ``bootstrap_3dsig_s_top10.py``.

Does not claim PoseBusters or thermodynamic ΔG.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
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
RMSD_KEY_RE = re.compile(r"REMARK\s+rmsd_(?:sym|raw)\s*=\s*([-+]?\d+\.?\d*)", re.I)
CF_RE = re.compile(r"REMARK\s+CF=([-+]?\d+\.?\d*)", re.I)
CF_APP_RE = re.compile(r"REMARK\s+CF\.app=([-+]?\d+\.?\d*)", re.I)
ACF_RE = re.compile(r"REMARK\s+(?:ACF|acf|soft_G|tilde_G)\s*=\s*([-+]?\d+\.?\d*)", re.I)

POSE_NAME_RE = re.compile(
    r"^(?P<pdb>[0-9A-Za-z]{4})_r(?P<restart>\d+)"
    r"(?:_(?P<minpts>\d+))?_(?P<rank>\d+)\.pdb$",
    re.I,
)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _f(m: Optional[re.Match]) -> Optional[float]:
    if not m:
        return None
    try:
        v = float(m.group(1))
    except (TypeError, ValueError):
        return None
    if not math.isfinite(v):
        return None
    return v


def parse_pose_pdb(path: Path) -> Dict[str, Optional[float]]:
    try:
        text = path.read_text(errors="replace")
    except OSError:
        return {
            "rmsd_sym": None,
            "rmsd_nosym": None,
            "cf": None,
            "cf_app": None,
            "acf": None,
        }
    head = "\n".join(text.splitlines()[:220])
    out: Dict[str, Optional[float]] = {
        "rmsd_sym": _f(RMSD_RE.search(head)),
        "rmsd_nosym": _f(RMSD_NS_RE.search(head)),
        "cf": _f(CF_RE.search(head)),
        "cf_app": _f(CF_APP_RE.search(head)),
        "acf": _f(ACF_RE.search(head)),
    }
    if out["rmsd_sym"] is None:
        out["rmsd_sym"] = _f(RMSD_KEY_RE.search(head))
    return out


def rmsd_of(meta: Dict[str, Optional[float]]) -> Optional[float]:
    v = meta.get("rmsd_sym")
    if v is not None and v >= 0:
        return v
    v = meta.get("rmsd_nosym")
    if v is not None and v >= 0:
        return v
    return None


def score_of(meta: Dict[str, Optional[float]]) -> Optional[float]:
    if meta.get("cf_app") is not None:
        return meta["cf_app"]
    if meta.get("cf") is not None:
        return meta["cf"]
    if meta.get("acf") is not None:
        return meta["acf"]
    return None


def collect_all_poses(out_dir: Path, pdb: str) -> List[Tuple[int, int, Path, Dict]]:
    """Return list of (restart, rank, path, meta)."""
    rows: List[Tuple[int, int, Path, Dict]] = []
    for path in sorted(out_dir.glob(f"{pdb}_r*.pdb")):
        if path.name.endswith("_INI.pdb") or path.name.endswith("_prepped.pdb"):
            continue
        m = POSE_NAME_RE.match(path.name)
        if not m:
            continue
        restart = int(m.group("restart"))
        rank = int(m.group("rank"))
        rows.append((restart, rank, path, parse_pose_pdb(path)))
    return rows


def collect_restart_poses(out_dir: Path, pdb: str) -> List[Tuple[int, Path, Dict]]:
    """Rank-0 (elected) pose per restart — legacy helper."""
    all_p = collect_all_poses(out_dir, pdb)
    by_r: Dict[int, Tuple[int, Path, Dict]] = {}
    for restart, rank, path, meta in all_p:
        if rank != 0:
            continue
        by_r[restart] = (restart, path, meta)
    return [by_r[k] for k in sorted(by_r)]


def top10_mode_rmsds(
    poses: List[Tuple[int, int, Path, Dict]], top_n: int = 10
) -> List[float]:
    """Global rank by CF among emission ranks < top_n; fallback restart heads."""
    if not poses:
        return []
    pool = [(r, rank, path, meta) for r, rank, path, meta in poses if rank < top_n]
    if not pool:
        pool = list(poses)
    # prefer multi-rank global; if only rank-0s, still sort by CF
    scored: List[Tuple[float, float, str]] = []
    for _r, rank, path, meta in pool:
        rmsd = rmsd_of(meta)
        sc = score_of(meta)
        if rmsd is None:
            continue
        if sc is None:
            sc = 1.0e30 + rank
        scored.append((sc, rmsd, path.name))
    if not scored:
        return []
    scored.sort(key=lambda t: t[0])
    return [rmsd for _sc, rmsd, _n in scored[:top_n]]


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
    all_poses = collect_all_poses(args.out_dir, pdb)
    poses = collect_restart_poses(args.out_dir, pdb)
    if not poses and not all_poses:
        for p in sorted(args.out_dir.glob("*_0.pdb")):
            poses.append((0, p, parse_pose_pdb(p)))
            all_poses.append((0, 0, p, poses[-1][2]))

    best = None
    for r, path, meta in poses:
        cf = score_of(meta)
        if cf is None:
            continue
        if best is None or cf < best[0]:
            best = (cf, r, path, meta)

    rmsd_top1 = None
    score_top1 = None
    elected_path = ""
    restarts_finished = len(poses) if poses else len({r for r, _, _, _ in all_poses})
    rmsd_bcr = None
    if best:
        score_top1 = best[0]
        elected_path = str(best[2])
        rmsd_top1 = rmsd_of(best[3])

    all_rmsds: List[float] = []
    for _r, _rank, _path, meta in all_poses:
        v = rmsd_of(meta)
        if v is not None:
            all_rmsds.append(v)
    # also scan loose globs (legacy non-matching names)
    if not all_rmsds:
        for p in args.out_dir.glob(f"{pdb}_r*_*.pdb"):
            if p.name.endswith("_INI.pdb"):
                continue
            meta = parse_pose_pdb(p)
            v = rmsd_of(meta)
            if v is not None:
                all_rmsds.append(v)
    if all_rmsds:
        rmsd_bcr = min(all_rmsds)

    mode_rmsds = top10_mode_rmsds(all_poses, top_n=10)
    success_s1 = int(rmsd_top1 is not None and rmsd_top1 <= 2.0)
    success_s3 = int(rmsd_bcr is not None and rmsd_bcr <= 2.0)
    success_s_top10 = int(bool(mode_rmsds) and min(mode_rmsds) < 2.0)

    bin_sha = ""
    if args.binary and args.binary.is_file():
        bin_sha = sha256_file(args.binary.resolve())

    n_cad = len(list(args.out_dir.glob(f"{pdb}_r*.cad")))
    gap = ""
    if not mode_rmsds:
        if n_cad and not all_poses:
            gap = "cad_only_no_pose_pdb"
        elif all_poses and not any(rmsd_of(m) is not None for *_, m in all_poses):
            gap = "poses_without_rmsd_remark"
        elif not all_poses:
            gap = "no_ranked_poses"

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
        "success_s_top10": success_s_top10 if mode_rmsds else "",
        "rank_native_mode": "",
        "n_poses": len(all_poses),
        "n_modes": restarts_finished,
        "n_top10_rmsds": len(mode_rmsds),
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
        "n_cad": n_cad,
        "rmsd_extract_gap": gap,
        "evals_actual": "",
        "budget_class": "full",
        "elected_path": elected_path,
        "parsed_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    for i in range(10):
        row[f"mode_rmsd_{i}"] = (
            f"{mode_rmsds[i]:.4f}" if i < len(mode_rmsds) else ""
        )

    out_csv = args.out_dir / "result.csv"
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(row.keys()))
        w.writeheader()
        w.writerow(row)
    s10 = (
        f"s_top10={success_s_top10} n_top10={len(mode_rmsds)}"
        if mode_rmsds
        else f"s_top10=NA gap={gap or 'no_rmsds'}"
    )
    print(
        f"wrote {out_csv} s1={success_s1} rmsd_top1={row['rmsd_top1']} "
        f"bcr={row['rmsd_bcr']} {s10}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
