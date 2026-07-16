#!/usr/bin/env python3
"""Parse classic FlexAID pilot outputs into a normalized result.csv row.

Reads emitted cluster-head PDBs for CF and RMSD when RMSDST was set.

**S1 / top-1:** elects best CF.app (else CF) among restart rank-0 poses
(``*_r*_0.pdb`` or ``*_0.pdb``).

**S_top10 / mode_rmsd_0..9:** RMSDs of the top-10 modes in **emitted rank
order** (crystal-blind). For multi-restart outputs, uses the elected restart's
rank 0..9 heads. Missing ranks → empty cells (not NaN strings).

Does not claim PoseBusters or thermodynamic ΔG.

Copyright 2026 Le Bonhomme Pharma
SPDX-License-Identifier: Apache-2.0
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

# Success threshold for S1 / S_top10 / BCR (claim contract: ≤ 2.0 Å)
RMSD_SUCCESS_THRESH = 2.0
TOP_N_MODES = 10

# Filename patterns (stem without .pdb):
#   {PDB}_{rank}                  CF/DP emission
#   {PDB}_{minPts}_{rank}         FO dual-suffix
#   {PDB}_r{restart}_{rank}       multi-restart CF/DP
#   {PDB}_r{restart}_{minPts}_{rank}  multi-restart FO
_RESTART_FO = re.compile(
    r"^(?P<pdb>[A-Za-z0-9]+)_r(?P<restart>\d+)_(?P<minpts>\d+)_(?P<rank>\d+)$",
    re.I,
)
_RESTART_CF = re.compile(
    r"^(?P<pdb>[A-Za-z0-9]+)_r(?P<restart>\d+)_(?P<rank>\d+)$",
    re.I,
)
_FO_DUAL = re.compile(
    r"^(?P<pdb>[A-Za-z0-9]+)_(?P<minpts>\d+)_(?P<rank>\d+)$",
    re.I,
)
_CF_RANK = re.compile(
    r"^(?P<pdb>[A-Za-z0-9]+)_(?P<rank>\d+)$",
    re.I,
)


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


def pose_rmsd(meta: Dict[str, Optional[float]]) -> Optional[float]:
    if meta.get("rmsd_sym") is not None:
        return meta["rmsd_sym"]
    return meta.get("rmsd_nosym")


def pose_score(meta: Dict[str, Optional[float]]) -> Optional[float]:
    """Crystal-blind score for election: CF.app preferred, else CF (lower better)."""
    if meta.get("cf_app") is not None:
        return meta["cf_app"]
    return meta.get("cf")


def parse_pose_filename(
    path: Path, pdb: str
) -> Optional[Tuple[Optional[int], int]]:
    """Return (restart_or_None, emitted_rank) or None if not a cluster head.

    Rank is the engine-emitted rank index (0 = top), never derived from RMSD.
    """
    stem = path.stem
    if stem.upper().endswith("_INI") or stem.upper() == f"{pdb.upper()}_INI":
        return None
    if not stem.upper().startswith(pdb.upper()):
        return None

    m = _RESTART_FO.match(stem)
    if m and m.group("pdb").upper() == pdb.upper():
        return int(m.group("restart")), int(m.group("rank"))

    m = _RESTART_CF.match(stem)
    if m and m.group("pdb").upper() == pdb.upper():
        return int(m.group("restart")), int(m.group("rank"))

    m = _FO_DUAL.match(stem)
    if m and m.group("pdb").upper() == pdb.upper():
        # FO dual-suffix: middle token is minPts (not a restart index)
        return None, int(m.group("rank"))

    m = _CF_RANK.match(stem)
    if m and m.group("pdb").upper() == pdb.upper():
        return None, int(m.group("rank"))

    return None


def collect_all_heads(
    out_dir: Path, pdb: str
) -> List[Tuple[Optional[int], int, Path, Dict[str, Optional[float]]]]:
    """All cluster-head PDBs: (restart, rank, path, meta)."""
    rows: List[Tuple[Optional[int], int, Path, Dict[str, Optional[float]]]] = []
    for path in sorted(out_dir.glob(f"{pdb}*.pdb")):
        parsed = parse_pose_filename(path, pdb)
        if parsed is None:
            continue
        restart, rank = parsed
        rows.append((restart, rank, path, parse_pose_pdb(path)))
    return rows


def collect_restart_poses(
    out_dir: Path, pdb: str
) -> List[Tuple[int, Path, Dict[str, Optional[float]]]]:
    """Rank-0 heads per restart (legacy helper; restart index from filename)."""
    rows: List[Tuple[int, Path, Dict[str, Optional[float]]]] = []
    for restart, rank, path, meta in collect_all_heads(out_dir, pdb):
        if rank != 0:
            continue
        r = 0 if restart is None else restart
        rows.append((r, path, meta))
    # Prefer restart-style globs if nothing matched via collect_all_heads
    if not rows:
        for rdir_pdb in sorted(out_dir.glob(f"{pdb}_r*_0.pdb")):
            name = rdir_pdb.name
            try:
                r = int(name.split("_r")[1].split("_")[0])
            except (IndexError, ValueError):
                r = -1
            rows.append((r, rdir_pdb, parse_pose_pdb(rdir_pdb)))
    return rows


def elect_best_rank0(
    heads: List[Tuple[Optional[int], int, Path, Dict[str, Optional[float]]]]
) -> Optional[Tuple[float, Optional[int], Path, Dict[str, Optional[float]]]]:
    """Elect best CF.app/CF among rank-0 heads. Returns (score, restart, path, meta)."""
    best: Optional[Tuple[float, Optional[int], Path, Dict[str, Optional[float]]]] = None
    for restart, rank, path, meta in heads:
        if rank != 0:
            continue
        cf = pose_score(meta)
        if cf is None:
            continue
        if best is None or cf < best[0]:
            best = (cf, restart, path, meta)
    return best


def mode_rmsds_emitted_order(
    heads: List[Tuple[Optional[int], int, Path, Dict[str, Optional[float]]]],
    elected_restart: Optional[int],
    n: int = TOP_N_MODES,
) -> List[Optional[float]]:
    """mode_rmsd_0..n-1 in **emitted rank order** (not RMSD-sorted).

    - Single-emission (no restart token): ranks 0..n-1 from that emission.
    - Multi-restart: ranks 0..n-1 from the **elected** restart only.
    - If a rank is missing → None (empty CSV cell).
    - If multiple heads share a rank (e.g. FO dual minPts), keep lowest score.
    """
    # Filter to the emission group we rank within
    restarts_present = {r for r, _, _, _ in heads if r is not None}
    if restarts_present:
        # Multi-restart: stick to elected restart; if None, fall back to first
        # available restart that has rank-0
        target = elected_restart
        if target is None:
            # Prefer restart with best rank-0 score already encoded as elected;
            # if still None, take the minimum restart id that has any head
            target = min(restarts_present)
        group = [(rank, path, meta) for r, rank, path, meta in heads if r == target]
    else:
        group = [(rank, path, meta) for r, rank, path, meta in heads if r is None]

    # Per-rank: keep best score when duplicates (FO multi minPts should not
    # happen under single-MinPts protocol, but be safe)
    by_rank: Dict[int, Tuple[Optional[float], Optional[float]]] = {}
    # value: (score_for_tiebreak, rmsd)
    for rank, _path, meta in group:
        if rank < 0 or rank >= n:
            continue
        sc = pose_score(meta)
        rmsd = pose_rmsd(meta)
        if rank not in by_rank:
            by_rank[rank] = (sc, rmsd)
        else:
            prev_sc, _ = by_rank[rank]
            if sc is not None and (prev_sc is None or sc < prev_sc):
                by_rank[rank] = (sc, rmsd)

    return [by_rank[i][1] if i in by_rank else None for i in range(n)]


def success_s_top10(
    mode_rmsds: List[Optional[float]], thresh: float = RMSD_SUCCESS_THRESH
) -> int:
    """1 if any finite mode_rmsd_i ≤ thresh."""
    for v in mode_rmsds:
        if v is not None and v == v and v >= 0.0 and v <= thresh:  # v==v → not NaN
            return 1
    return 0


def estimate_wall_s(out_dir: Path, pdb: str) -> Optional[float]:
    """Estimate compute wall when launcher did not pass --wall-s.

    Priority:
      1. ``wall_s.txt`` / ``wall_timing.json`` written by arm launcher
      2. max(mtime) − min(mtime) over ``{pdb}_r*_*.pdb`` (or ``{pdb}_*.pdb``)

    Returns None if unavailable.
    """
    for name in ("wall_s.txt", "wall_timing.json"):
        wp = out_dir / name
        if not wp.is_file():
            continue
        try:
            text = wp.read_text().strip()
            if name.endswith(".json"):
                import json

                d = json.loads(text)
                v = d.get("wall_s")
                if v is not None:
                    return float(v)
            else:
                return float(text.split()[0])
        except (OSError, ValueError, TypeError, KeyError):
            continue

    pdb_u = pdb.upper()
    mt: List[float] = []
    for pat in (f"{pdb_u}_r*_*.pdb", f"{pdb_u}_*.pdb", f"{pdb.lower()}_r*_*.pdb"):
        for p in out_dir.glob(pat):
            try:
                if p.is_file() and p.suffix.lower() == ".pdb":
                    mt.append(p.stat().st_mtime)
            except OSError:
                continue
        if len(mt) >= 2:
            break
    if len(mt) < 2:
        return None
    span = max(mt) - min(mt)
    if span <= 0:
        return None
    return float(span)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arm", required=True)
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--work-dir", type=Path, default=None)
    ap.add_argument("--matrix-md5", default="")
    ap.add_argument("--binary", type=Path, default=None)
    ap.add_argument(
        "--wall-s",
        type=float,
        default=None,
        help="Measured FlexAID wall seconds for this target (from launcher).",
    )
    args = ap.parse_args()

    pdb = args.pdb.upper()
    heads = collect_all_heads(args.out_dir, pdb)
    if not heads:
        # Fallback: any *_0.pdb rank-0 style
        for p in sorted(args.out_dir.glob("*_0.pdb")):
            meta = parse_pose_pdb(p)
            heads.append((0, 0, p, meta))

    best = elect_best_rank0(heads)

    rmsd_top1: Optional[float] = None
    score_top1: Optional[float] = None
    elected_path = ""
    elected_restart: Optional[int] = None
    if best:
        score_top1 = best[0]
        elected_restart = best[1]
        elected_path = str(best[2])
        rmsd_top1 = pose_rmsd(best[3])

    mode_rmsds = mode_rmsds_emitted_order(heads, elected_restart, TOP_N_MODES)
    # Keep S1 aligned with elected rank-0 (mode_rmsd_0 when elected emission used)
    if mode_rmsds and mode_rmsds[0] is not None and rmsd_top1 is None:
        rmsd_top1 = mode_rmsds[0]
    # Prefer elected rmsd_top1 as mode_rmsd_0 when both exist (same pose)
    if rmsd_top1 is not None:
        mode_rmsds[0] = rmsd_top1

    all_rmsds: List[float] = []
    for _r, _rank, _path, meta in heads:
        v = pose_rmsd(meta)
        if v is not None:
            all_rmsds.append(v)
    rmsd_bcr = min(all_rmsds) if all_rmsds else None

    restarts_finished = len({r for r, rank, _, _ in heads if rank == 0})
    if restarts_finished == 0:
        restarts_finished = len(heads)

    success_s1 = int(rmsd_top1 is not None and rmsd_top1 <= RMSD_SUCCESS_THRESH)
    success_s3 = int(rmsd_bcr is not None and rmsd_bcr <= RMSD_SUCCESS_THRESH)
    s_top10 = success_s_top10(mode_rmsds, RMSD_SUCCESS_THRESH)

    bin_sha = ""
    if args.binary and args.binary.is_file():
        bin_sha = sha256_file(args.binary.resolve())

    row: Dict[str, object] = {
        "arm": args.arm,
        "engine_sha": bin_sha,
        "matrix_md5": args.matrix_md5,
        "pdb_id": pdb,
        "rmsd_top1": "" if rmsd_top1 is None else f"{rmsd_top1:.4f}",
        "rmsd_bcr": "" if rmsd_bcr is None else f"{rmsd_bcr:.4f}",
        "success_s1": success_s1,
        "success_s2": "",
        "success_s3": success_s3,
        "success_s_top10": s_top10,
        "rank_native_mode": "",
        "n_poses": len(heads),
        "n_modes": sum(1 for v in mode_rmsds if v is not None),
        "score_top1": "" if score_top1 is None else f"{score_top1:.5f}",
        "H": "",
        "TS": "",
        "F": "",
        "pb_pass": "",
        "tencom_status": "NA",
        "seed_echo": 0,
        "native_pose_seeded": 0,
        "protocol_claim_eligible": 1 if args.matrix_md5 else 0,
        "wall_s": "",  # filled below
        "restarts_finished": restarts_finished,
        "evals_actual": "",
        "budget_class": "full",
        "elected_path": elected_path,
        "parsed_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    for i, v in enumerate(mode_rmsds):
        row[f"mode_rmsd_{i}"] = "" if v is None else f"{v:.4f}"

    # Computational walltime: launcher measurement preferred, else pose-mtime proxy
    wall: Optional[float] = None
    if args.wall_s is not None and args.wall_s >= 0:
        wall = float(args.wall_s)
    else:
        wall = estimate_wall_s(args.out_dir, pdb)
    if wall is not None and wall >= 0:
        row["wall_s"] = f"{wall:.1f}"

    # Stable column order: identity → top1/BCR → mode_rmsd_* → flags
    fieldnames = [
        "arm",
        "engine_sha",
        "matrix_md5",
        "pdb_id",
        "rmsd_top1",
        "rmsd_bcr",
        *[f"mode_rmsd_{i}" for i in range(TOP_N_MODES)],
        "success_s1",
        "success_s2",
        "success_s3",
        "success_s_top10",
        "rank_native_mode",
        "n_poses",
        "n_modes",
        "score_top1",
        "H",
        "TS",
        "F",
        "pb_pass",
        "tencom_status",
        "seed_echo",
        "native_pose_seeded",
        "protocol_claim_eligible",
        "wall_s",
        "restarts_finished",
        "evals_actual",
        "budget_class",
        "elected_path",
        "parsed_utc",
    ]

    out_csv = args.out_dir / "result.csv"
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerow(row)
    print(
        f"wrote {out_csv} s1={success_s1} s_top10={s_top10} "
        f"rmsd_top1={row['rmsd_top1']} bcr={row['rmsd_bcr']} "
        f"wall_s={row['wall_s'] or 'NA'} "
        f"modes={[row[f'mode_rmsd_{i}'] for i in range(TOP_N_MODES)]}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
