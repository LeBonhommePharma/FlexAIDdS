#!/usr/bin/env python3
"""Extract per-case top-10 RMSDs for 3Dsig S_top10 bootstrap.

Reads FlexAID arm OUT trees::

  <arm-dir>/<PDB>/{PDB}_r{restart}_{rank}.pdb
  <arm-dir>/<PDB>/{PDB}_r{restart}_{minPts}_{rank}.pdb   # FO dual-suffix

and optional ``.rrd`` sidecars (crystal RMSDs per chromosome).

Emits cases JSON for ``scripts/bootstrap_3dsig_s_top10.py``.

S_top10 definition (3Dsig 2017): success if min(RMSD among top-10 ranked
modes) < 2.0 Å. Ranking default is CF.app ascending (arm A / CF); use
``--score acf`` when REMARK ACF / soft-β free energy is present (arm B).

Strategies
----------
global
    Pool all emitted modes with emission rank < 10 across restarts, sort by
    score, keep top 10 RMSDs. Preferred when multi-rank pose PDBs exist.
restart_heads
    Rank-0 head from each restart, sorted by score (≤10). Natural when each
    of 10 sims emits only a single elected mode.
auto
    Use global if ≥1 multi-rank pose exists, else restart_heads.

Gap report
----------
Prints why a PDB lacks RMSDs (no poses, segfault after .cad, empty REMARK,
in-progress). Does **not** invent RMSDs from inter-cluster .cad pairs.

Usage::

  python3 scripts/extract_3dsig_s_top10_from_arm.py \\
    --arm-dir ~/flexaidds_results/campaigns/three_engine/A/3dsig_r10 \\
    --json-out /tmp/A_cases.json

  python3 scripts/bootstrap_3dsig_s_top10.py --cases /tmp/A_cases.json

Copyright 2026 Le Bonhomme Pharma
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

RMSD_SYM_RE = re.compile(
    r"REMARK\s+([-+]?\d+\.?\d*)\s+RMSD to ref\. structure\s+\(symmetry corrected\)",
    re.I,
)
RMSD_NS_RE = re.compile(
    r"REMARK\s+([-+]?\d+\.?\d*)\s+RMSD to ref\. structure\s+\(no symmetry",
    re.I,
)
RMSD_KEY_RE = re.compile(r"REMARK\s+rmsd_(?:sym|raw)\s*=\s*([-+]?\d+\.?\d*)", re.I)
CF_APP_RE = re.compile(r"REMARK\s+CF\.app=([-+]?\d+\.?\d*)", re.I)
CF_RE = re.compile(r"REMARK\s+CF=([-+]?\d+\.?\d*)", re.I)
ACF_RE = re.compile(r"REMARK\s+(?:ACF|acf|soft_G|tilde_G)\s*=\s*([-+]?\d+\.?\d*)", re.I)

# {pdb}_r{restart}_{rank}.pdb  OR  {pdb}_r{restart}_{minPts}_{rank}.pdb
POSE_NAME_RE = re.compile(
    r"^(?P<pdb>[0-9A-Za-z]{4})_r(?P<restart>\d+)"
    r"(?:_(?P<minpts>\d+))?_(?P<rank>\d+)\.pdb$",
    re.I,
)


@dataclass
class ModePose:
    path: Path
    restart: int
    rank: int
    minpts: Optional[int]
    rmsd: Optional[float]
    cf: Optional[float]
    acf: Optional[float]


@dataclass
class CaseExtract:
    pdb_id: str
    rmsds: List[float] = field(default_factory=list)
    scores: List[float] = field(default_factory=list)
    sources: List[str] = field(default_factory=list)
    n_pose_pdbs: int = 0
    n_with_rmsd: int = 0
    n_restarts_with_cad: int = 0
    n_restarts_with_pose: int = 0
    gap: str = ""
    strategy_used: str = ""


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


def parse_pose_pdb(path: Path) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    """Return (rmsd_prefer_sym, cf_app_or_cf, acf_if_any)."""
    try:
        text = path.read_text(errors="replace")
    except OSError:
        return None, None, None
    # only scan REMARK block (first ~200 lines)
    head = "\n".join(text.splitlines()[:220])
    rmsd = _f(RMSD_SYM_RE.search(head))
    if rmsd is None:
        rmsd = _f(RMSD_NS_RE.search(head))
    if rmsd is None:
        rmsd = _f(RMSD_KEY_RE.search(head))
    if rmsd is not None and rmsd < 0:
        rmsd = None  # sentinel
    cf = _f(CF_APP_RE.search(head))
    if cf is None:
        cf = _f(CF_RE.search(head))
    acf = _f(ACF_RE.search(head))
    return rmsd, cf, acf


def parse_rrd_cluster_heads(rrd: Path, max_rank: int = 10) -> List[Tuple[int, float, float]]:
    """Parse .rrd → list of (cluster_id, rmsd_sym, cf) for unique cluster heads.

    .rrd columns (write_rrd.cpp): chrom_idx, clus_gapop, clus_rmsd, rmsd,
    rmsd_corrected, evalue, [genes...]
    """
    heads: Dict[int, Tuple[float, float]] = {}
    try:
        lines = rrd.read_text(errors="replace").splitlines()
    except OSError:
        return []
    for line in lines:
        parts = line.split()
        if len(parts) < 6:
            continue
        try:
            clus = int(float(parts[1]))
            rmsd_sym = float(parts[4])
            cf = float(parts[5])
        except ValueError:
            continue
        if not math.isfinite(rmsd_sym) or rmsd_sym < 0:
            continue
        # keep best CF per cluster head id
        prev = heads.get(clus)
        if prev is None or cf < prev[1]:
            heads[clus] = (rmsd_sym, cf)
    ranked = sorted(heads.items(), key=lambda kv: kv[1][1])  # by CF
    out: List[Tuple[int, float, float]] = []
    for clus, (rmsd, cf) in ranked[:max_rank]:
        out.append((clus, rmsd, cf))
    return out


def collect_pose_files(case_dir: Path, pdb: str) -> List[ModePose]:
    poses: List[ModePose] = []
    for path in sorted(case_dir.glob(f"{pdb}_r*.pdb")):
        name = path.name
        if name.endswith("_INI.pdb") or name.endswith("_prepped.pdb"):
            continue
        m = POSE_NAME_RE.match(name)
        if not m:
            continue
        restart = int(m.group("restart"))
        rank = int(m.group("rank"))
        minpts = int(m.group("minpts")) if m.group("minpts") else None
        rmsd, cf, acf = parse_pose_pdb(path)
        poses.append(
            ModePose(
                path=path,
                restart=restart,
                rank=rank,
                minpts=minpts,
                rmsd=rmsd,
                cf=cf,
                acf=acf,
            )
        )
    return poses


def score_of(mode: ModePose, score: str) -> float:
    if score == "acf" and mode.acf is not None:
        return mode.acf
    if mode.cf is not None:
        return mode.cf
    if mode.acf is not None:
        return mode.acf
    # unknown score → push to end but keep deterministic by rank
    return 1.0e30 + mode.rank


def select_top10(
    poses: Sequence[ModePose],
    strategy: str,
    score: str,
    top_n: int = 10,
) -> Tuple[List[ModePose], str]:
    if not poses:
        return [], strategy

    has_multi = any(p.rank > 0 for p in poses)
    used = strategy
    if strategy == "auto":
        used = "global" if has_multi else "restart_heads"

    if used == "restart_heads":
        # best (lowest rank, then score) per restart — prefer rank 0
        by_r: Dict[int, ModePose] = {}
        for p in poses:
            cur = by_r.get(p.restart)
            if cur is None or p.rank < cur.rank or (
                p.rank == cur.rank and score_of(p, score) < score_of(cur, score)
            ):
                by_r[p.restart] = p
        selected = sorted(by_r.values(), key=lambda p: score_of(p, score))[:top_n]
        return selected, used

    # global: emission rank < top_n, then sort by score
    pool = [p for p in poses if p.rank < top_n]
    if not pool:
        pool = list(poses)
    selected = sorted(pool, key=lambda p: score_of(p, score))[:top_n]
    return selected, used


def extract_case(
    case_dir: Path,
    strategy: str = "auto",
    score: str = "cf",
    top_n: int = 10,
) -> CaseExtract:
    pdb = case_dir.name.upper()[:4]
    out = CaseExtract(pdb_id=pdb)
    poses = collect_pose_files(case_dir, pdb)
    out.n_pose_pdbs = len(poses)
    out.n_with_rmsd = sum(1 for p in poses if p.rmsd is not None)
    out.n_restarts_with_pose = len({p.restart for p in poses})
    out.n_restarts_with_cad = len(list(case_dir.glob(f"{pdb}_r*.cad")))

    # Prefer pose PDBs
    if poses:
        selected, used = select_top10(poses, strategy, score, top_n=top_n)
        out.strategy_used = used
        for p in selected:
            if p.rmsd is None:
                continue
            out.rmsds.append(p.rmsd)
            out.scores.append(score_of(p, score))
            out.sources.append(p.path.name)
        if out.rmsds:
            out.gap = ""
            return out
        out.gap = (
            f"found {len(poses)} pose PDB(s) but none have REMARK RMSD "
            f"(RMSDST missing or write failed before RMSD remark)"
        )
        return out

    # Fallback: .rrd (crystal RMSDs) if clustering finished write_rrd
    rrds = sorted(case_dir.glob(f"{pdb}_r*.rrd"))
    if rrds:
        # pool cluster heads from all restarts by CF
        pooled: List[Tuple[float, float, str]] = []  # cf, rmsd, src
        for rrd in rrds:
            for clus, rmsd, cf in parse_rrd_cluster_heads(rrd, max_rank=top_n):
                pooled.append((cf, rmsd, f"{rrd.name}:clus{clus}"))
        pooled.sort(key=lambda t: t[0])
        for cf, rmsd, src in pooled[:top_n]:
            out.rmsds.append(rmsd)
            out.scores.append(cf)
            out.sources.append(src)
        out.strategy_used = "rrd_heads"
        if out.rmsds:
            out.gap = ""
            return out

    # Diagnose empty
    inis = list(case_dir.glob(f"{pdb}_r*_INI.pdb"))
    cads = list(case_dir.glob(f"{pdb}_r*.cad"))
    if cads and not poses:
        out.gap = (
            f"cad_only: {len(cads)} .cad and {len(inis)} _INI.pdb but 0 ranked "
            f"pose PDBs — clustering likely crashed after writing .cad "
            f"(see 'Segmentation fault' near 'clustering all individuals')"
        )
    elif inis and not cads:
        out.gap = f"ini_only: {len(inis)} _INI.pdb, dock in progress or failed before cluster"
    elif not inis and not cads:
        out.gap = "empty: no FlexAID restart outputs yet"
    else:
        out.gap = "no_rmsd_source"
    out.strategy_used = "none"
    return out


def scan_arm(
    arm_dir: Path,
    strategy: str,
    score: str,
    top_n: int,
) -> Dict[str, CaseExtract]:
    cases: Dict[str, CaseExtract] = {}
    if not arm_dir.is_dir():
        return cases
    for child in sorted(arm_dir.iterdir()):
        if not child.is_dir():
            continue
        name = child.name
        if len(name) < 4 or not name[:4].isalnum():
            continue
        # skip obvious non-case dirs
        if name.startswith(".") or name in ("logs", "meta", "tmp"):
            continue
        ex = extract_case(child, strategy=strategy, score=score, top_n=top_n)
        cases[ex.pdb_id] = ex
    return cases


def cases_to_json_payload(
    cases: Dict[str, CaseExtract],
    arm_dir: Path,
    strategy: str,
    score: str,
) -> dict:
    ready = {k: v for k, v in cases.items() if v.rmsds}
    payload = {
        "metric": "S_top10",
        "thresh_A": 2.0,
        "arm_dir": str(arm_dir),
        "strategy": strategy,
        "score": score,
        "n_case_dirs": len(cases),
        "n_cases_with_rmsds": len(ready),
        "n_cases_missing": len(cases) - len(ready),
        "cases": {
            k: {
                "rmsds": v.rmsds,
                "scores": v.scores,
                "sources": v.sources,
                "strategy_used": v.strategy_used,
                "n_pose_pdbs": v.n_pose_pdbs,
                "n_with_rmsd": v.n_with_rmsd,
                "n_restarts_with_cad": v.n_restarts_with_cad,
                "n_restarts_with_pose": v.n_restarts_with_pose,
                "gap": v.gap,
            }
            for k, v in sorted(cases.items())
        },
        "ready_pdb_ids": sorted(ready.keys()),
        "missing": {
            k: v.gap
            for k, v in sorted(cases.items())
            if not v.rmsds
        },
    }
    return payload


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arm-dir", type=Path, required=True, help="e.g. .../A/3dsig_r10")
    ap.add_argument("--json-out", type=Path, default=None, help="cases JSON for bootstrap")
    ap.add_argument(
        "--strategy",
        choices=("auto", "global", "restart_heads"),
        default="auto",
    )
    ap.add_argument(
        "--score",
        choices=("cf", "acf"),
        default="cf",
        help="cf=CF.app (arm A); acf=soft-β free energy when present (arm B)",
    )
    ap.add_argument("--top-n", type=int, default=10)
    ap.add_argument(
        "--require-full-top10",
        action="store_true",
        help="Only include cases with ≥ top-n RMSDs (strict); default keeps partial lists",
    )
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    arm_dir = args.arm_dir.expanduser().resolve()
    cases = scan_arm(arm_dir, args.strategy, args.score, args.top_n)

    if args.require_full_top10:
        for v in cases.values():
            if 0 < len(v.rmsds) < args.top_n:
                v.gap = (
                    f"partial_top{args.top_n}: only {len(v.rmsds)} RMSDs "
                    f"(dropped by --require-full-top10)"
                )
                v.rmsds = []
                v.scores = []
                v.sources = []

    payload = cases_to_json_payload(cases, arm_dir, args.strategy, args.score)

    if not args.quiet:
        print(f"arm_dir={arm_dir}")
        print(
            f"case_dirs={payload['n_case_dirs']} "
            f"with_rmsds={payload['n_cases_with_rmsds']} "
            f"missing={payload['n_cases_missing']}"
        )
        for pid, ex in sorted(cases.items()):
            if ex.rmsds:
                mn = min(ex.rmsds)
                ok = "OK" if mn < 2.0 else "miss"
                print(
                    f"  {pid}: n={len(ex.rmsds)} min_rmsd={mn:.3f} "
                    f"S_top10={ok} strategy={ex.strategy_used}"
                )
            else:
                print(f"  {pid}: NA  gap={ex.gap}")

    if args.json_out:
        args.json_out = args.json_out.expanduser()
        # bootstrap load_cases_json wants either list or {cases: {pdb: {rmsds}}}
        # Keep full payload; bootstrap already accepts dict with "cases"
        bootstrap_view = {
            "cases": {
                k: {"rmsds": v["rmsds"]}
                for k, v in payload["cases"].items()
                if v["rmsds"]
            },
            "meta": {
                "arm_dir": payload["arm_dir"],
                "strategy": payload["strategy"],
                "score": payload["score"],
                "n_case_dirs": payload["n_case_dirs"],
                "n_cases_with_rmsds": payload["n_cases_with_rmsds"],
                "missing": payload["missing"],
            },
        }
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(bootstrap_view, indent=2) + "\n")
        # also write full diagnostic next to it
        diag = args.json_out.with_name(args.json_out.stem + "_full.json")
        diag.write_text(json.dumps(payload, indent=2) + "\n")
        print(f"wrote {args.json_out}")
        print(f"wrote {diag}")

    # exit 0 even if empty — caller/bootstrap reports NA
    if payload["n_cases_with_rmsds"] == 0:
        print("NOTE: no cases with RMSDs yet — bootstrap will report NA", file=sys.stderr)
        return 0
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
