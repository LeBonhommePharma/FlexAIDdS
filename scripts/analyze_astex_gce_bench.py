#!/usr/bin/env python3
"""Analyze an Astex Diverse benchmark_datasets output dir.

Reports:
  - RMSD-only top-1 / any-pose success (≤ 2.0 Å)
  - Optional PoseBusters if FLEXAIDDS_POSEBUSTERS_BIN or .venv-posebusters/bin/bust
  - GrandPartitionFunction-style competitive ranking from per-complex scores
    (diagnostic multi-ligand Ξ using best pose score as proxy for ln Z)

Does NOT claim true ΔG or experimental Ki. CF/score is a proxy.

Usage:
  python3 scripts/analyze_astex_gce_bench.py <output_dir>
"""
from __future__ import annotations

import csv
import math
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

kB = 0.001987206
T = 300.0
RMSD_CUT = 2.0


def fnum(x) -> Optional[float]:
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def find_result_csvs(root: Path) -> List[Path]:
    return sorted(root.rglob("result.csv"))


def parse_result_csv(path: Path) -> List[dict]:
    with path.open() as f:
        return list(csv.DictReader(f))


def extract_rmsd_score(rows: List[dict]) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    """Return (top1_rmsd, best_rmsd, top1_score) using lowest score as top1."""
    if not rows:
        return None, None, None

    def score(r: dict) -> float:
        for k in (
            "total_score",
            "score",
            "cf",
            "best_score",
            "energy",
            "free_energy",
            "predicted_dG",
        ):
            v = fnum(r.get(k))
            if v is not None:
                return v
        return 1e99

    def rmsd(r: dict) -> Optional[float]:
        for k, v in r.items():
            if k and "rmsd" in k.lower():
                fv = fnum(v)
                if fv is not None:
                    return fv
        return None

    ranked = sorted(rows, key=score)
    top = ranked[0]
    tr = rmsd(top)
    brs = [rmsd(r) for r in rows]
    brs = [b for b in brs if b is not None]
    br = min(brs) if brs else None
    return tr, br, score(top)


def posebusters_bin() -> Optional[Path]:
    env = os.environ.get("FLEXAIDDS_POSEBUSTERS_BIN")
    if env and Path(env).exists():
        return Path(env)
    root = Path(__file__).resolve().parents[1]
    cand = root / ".venv-posebusters" / "bin" / "bust"
    if cand.exists():
        return cand
    return None


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    root = Path(sys.argv[1]).expanduser().resolve()
    if not root.exists():
        print(f"ERROR: {root} not found")
        return 1

    csvs = find_result_csvs(root)
    print(f"=== Astex GCE / docking-power analysis ===")
    print(f"root: {root}")
    print(f"result.csv files: {len(csvs)}")
    print(f"RMSD cut: {RMSD_CUT} Å")
    print(f"PoseBusters: {posebusters_bin() or 'NOT FOUND (RMSD-only)'}")
    print()

    top1_ok = 0
    any_ok = 0
    rows_out = []
    for rc in csvs:
        code = rc.parent.name
        rows = parse_result_csv(rc)
        tr, br, sc = extract_rmsd_score(rows)
        t1 = tr is not None and tr <= RMSD_CUT
        ao = br is not None and br <= RMSD_CUT
        top1_ok += int(t1)
        any_ok += int(ao)
        rows_out.append((code, tr, br, sc, len(rows), t1, ao))
        print(
            f"  {code:6s} top1_rmsd={tr!s:>8} best_rmsd={br!s:>8} "
            f"top1_score={sc!s:>12} n={len(rows):4d} "
            f"top1_ok={t1} any_ok={ao}"
        )

    n = len(csvs)
    if n:
        print()
        print(f"TOP1 RMSD≤{RMSD_CUT}: {top1_ok}/{n} = {100*top1_ok/n:.1f}%")
        print(f"ANY  RMSD≤{RMSD_CUT}: {any_ok}/{n} = {100*any_ok/n:.1f}%")
        print(
            "NOTE: full success requires PoseBusters pass as well "
            "(skill contract). RMSD-only is interim."
        )

    # Competitive GCE diagnostic: treat each complex top1 score as F ≈ −kT ln Z
    # ⇒ ln Z ≈ −F / kT. Register as multi-ligand Ξ at equal concentration.
    print()
    print("=== Grand-canonical competitive ranking (diagnostic) ===")
    print("Proxy: log_Z_i = −score_i / (kT) from top-1 CF/score (NOT calibrated ΔG)")
    kT = kB * T
    ligands = []
    for code, tr, br, sc, nposes, t1, ao in rows_out:
        if sc is None or sc > 1e90:
            continue
        log_Z = -sc / kT
        ligands.append((code, log_Z, sc, tr))
    if not ligands:
        print("  (no scores available)")
        return 0

    # log_Xi = logsumexp(0, log_Z_i...) at c=1 M for all
    terms = [0.0] + [lz for _, lz, _, _ in ligands]
    m = max(terms)
    lxi = m + math.log(sum(math.exp(t - m) for t in terms))
    print(f"  T={T} K  ln Ξ = {lxi:.4f}  Ω = {-kT*lxi:.3f} kcal/mol")
    ranked = sorted(ligands, key=lambda x: -x[1])  # high Z first
    print("  rank  pdb     p_bound     log_Z      score      top1_rmsd")
    for i, (code, lz, sc, tr) in enumerate(ranked, 1):
        p = math.exp(lz - lxi)
        print(f"  {i:4d}  {code:6s}  {p:9.4f}  {lz:9.3f}  {sc:9.3f}  {tr}")

    # mixing entropy / collapse
    pe = math.exp(-lxi)
    probs = [pe] + [math.exp(lz - lxi) for _, lz, _, _ in ligands]
    S = -sum(p * math.log(p) for p in probs if p > 1e-300)
    p_bound = 1.0 - pe
    if p_bound > 1e-15 and len(ligands) > 1:
        Slig = 0.0
        for _, lz, _, _ in ligands:
            pt = math.exp(lz - lxi) / p_bound
            if pt > 1e-300:
                Slig -= pt * math.log(pt)
        collapse = 1.0 - Slig / math.log(len(ligands))
        collapse = max(0.0, min(1.0, collapse))
    else:
        collapse = 0.0
    print(f"  S_mix (nats)={S:.4f}  ligand_entropy_collapse={collapse:.3f}")
    print("  (collapse→1 means one complex dominates the toy multi-ligand Ξ)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
