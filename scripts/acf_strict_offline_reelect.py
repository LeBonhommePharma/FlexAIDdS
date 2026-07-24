#!/usr/bin/env python3
"""W1.1 offline ACF_STRICT election re-rank on frozen campaign heads (no re-dock).

For each target: parse ranked pose REMARK CF values; compute:
  - legacy ACF G = Emin - T ln Z over all head CFs (multiplicity-sensitive if dups)
  - strict G via unique CF collapse (mirrors SoftBetaFreeEnergy free_energy_strict)

Then rank heads by legacy soft_beta_G REMARK (if present) vs min CF (strict proxy).

This is the in-session pilot while a live full-85 holds the box: no second dock.

Uses the same free_energy math as LIB/SoftBetaFreeEnergy.h (tested against
SoftBetaStrict.ExactDuplicateInvariance vectors).
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Optional

# Import E10 helpers
sys.path.insert(0, str(Path(__file__).resolve().parent))
from e10_election_vs_scoring import (  # type: ignore
    analyze_target,
    discover_targets,
    parse_pdb_remarks,
    list_ranked_poses,
    load_result_row,
    _f,
)


def free_energy(energies: list[float], T: float) -> float:
    """Classic soft_beta::free_energy G = Emin - T ln Z (multiplicity-sensitive)."""
    finite = [e for e in energies if math.isfinite(e)]
    if not finite:
        return float("inf")
    T = max(T, 1e-12)
    Emin = min(finite)
    Z = sum(math.exp(-(e - Emin) / T) for e in finite)
    if not (Z > 0) or not math.isfinite(Z):
        return Emin
    return Emin - T * math.log(Z)


def free_energy_strict(energies: list[float], T: float) -> float:
    """UniqueGeometry: collapse exact equal CF then free_energy."""
    unique: list[float] = []
    for e in energies:
        if not math.isfinite(e):
            continue
        if e not in unique:  # exact equality as in C++
            unique.append(e)
    return free_energy(unique, T)


def reelect_target(target_dir: Path, T: float = 300.0) -> dict[str, Any]:
    row = load_result_row(target_dir)
    pdb = str(row.get("pdb_id") or target_dir.name).strip()
    ranked = list_ranked_poses(target_dir, pdb)
    cfs: list[float] = []
    head_meta: list[dict[str, Any]] = []
    for rank, path in ranked:
        rem = parse_pdb_remarks(path)
        cf = rem.get("cf", float("nan"))
        cfs.append(cf)
        head_meta.append({"rank": rank, "cf": cf, "path": str(path)})

    finite_pairs = [(h, c) for h, c in zip(head_meta, cfs) if math.isfinite(c)]
    if not finite_pairs:
        return {"pdb": pdb, "status": "no_cf"}

    # Elect by min CF among heads
    min_cf_head, min_cf_val = min(finite_pairs, key=lambda x: x[1])
    g_legacy = free_energy([c for _, c in finite_pairs], T)
    g_strict = free_energy_strict([c for _, c in finite_pairs], T)

    rmsd = _f(row.get("rmsd_hungarian"))
    if not math.isfinite(rmsd):
        rmsd = _f(row.get("rmsd_to_crystal"))
    bcr = _f(row.get("best_cluster_rmsd"))
    se = str(row.get("seed_echo", "0")).strip()
    elect_cf = _f(row.get("elected_cf"))

    # Duplicate inflation demo: clone min CF 10×
    single = [min_cf_val]
    clones = [min_cf_val] * 10
    g_legacy_clones = free_energy(clones, T)
    g_strict_clones = free_energy_strict(clones, T)
    g_single = free_energy(single, T)

    election_gap = (
        se in ("0", "0.0", "")
        and math.isfinite(bcr)
        and bcr <= 2.5
        and math.isfinite(rmsd)
        and rmsd > 2.0
    )

    return {
        "pdb": pdb,
        "status": "ok",
        "seed_echo": se,
        "rmsd_elected": rmsd,
        "bcr": bcr,
        "elected_cf": elect_cf,
        "min_cf_among_heads": min_cf_val,
        "min_cf_head_rank": min_cf_head["rank"],
        "g_legacy_heads": g_legacy,
        "g_strict_heads": g_strict,
        "strict_less_deep_than_legacy": g_strict > g_legacy + 1e-9,
        "clone_legacy_deepens": g_legacy_clones < g_single - 1.0,
        "clone_strict_invariant": abs(g_strict_clones - g_single) < 1e-6,
        "election_gap": election_gap,
        "n_heads": len(finite_pairs),
        "note": (
            "election_gap: BCR near-native but rank-0 failed — ACF_STRICT dock pilot recommended"
            if election_gap
            else ""
        ),
    }


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--campaign-dir", type=Path, required=True)
    ap.add_argument("--out-json", type=Path, required=True)
    ap.add_argument("--out-md", type=Path, default=None)
    ap.add_argument("--T", type=float, default=300.0)
    ap.add_argument(
        "--pdbs",
        nargs="*",
        default=None,
        help="Optional subset (default: all targets with result.csv)",
    )
    args = ap.parse_args(argv)

    camp = args.campaign_dir.expanduser()
    targets = discover_targets(camp)
    if args.pdbs:
        want = set(args.pdbs)
        targets = [t for t in targets if t.name in want]

    results = [reelect_target(t, T=args.T) for t in targets]
    gaps = [r for r in results if r.get("election_gap")]
    goods = [r for r in results if r.get("pdb") in ("1HNN", "1HP0", "1HQ2")]
    goods_ok = all(
        (not math.isfinite(r.get("rmsd_elected", float("nan"))) or r["rmsd_elected"] <= 2.0)
        for r in goods
        if r.get("status") == "ok" and r.get("rmsd_elected", 99) <= 2.0 or True
    )
    # goods non-regression for offline: just list current genuine status
    genuine_goods = [
        r["pdb"]
        for r in results
        if r.get("pdb") in ("1HNN", "1HP0", "1HQ2")
        and math.isfinite(r.get("rmsd_elected", float("nan")))
        and r["rmsd_elected"] <= 2.0
    ]

    payload = {
        "campaign": str(camp),
        "T": args.T,
        "n_targets": len(results),
        "n_election_gap": len(gaps),
        "genuine_goods_present": genuine_goods,
        "clone_invariance_ok": all(
            r.get("clone_strict_invariant", True)
            for r in results
            if r.get("status") == "ok"
        ),
        "results": results,
        "pilot_note": (
            "Offline re-elect only (live full-85 holds box). "
            "Dock pilot ACF_STRICT=0/1 after baseline finishes."
        ),
    }
    args.out_json.expanduser().parent.mkdir(parents=True, exist_ok=True)
    args.out_json.expanduser().write_text(json.dumps(payload, indent=2) + "\n")

    if args.out_md:
        lines = [
            "# W1.1 offline ACF_STRICT re-elect pilot",
            "",
            f"Campaign: `{camp}`",
            f"Election-gap targets: **{len(gaps)}** / {len(results)}",
            f"Genuine goods in tree: {genuine_goods}",
            f"Strict clone invariance: **{payload['clone_invariance_ok']}**",
            "",
            "| pdb | rmsd | bcr | gap | min_cf_rank | g_legacy | g_strict |",
            "|-----|------|-----|-----|-------------|----------|----------|",
        ]
        for r in results:
            if r.get("status") != "ok":
                continue
            lines.append(
                f"| {r['pdb']} | {r.get('rmsd_elected')} | {r.get('bcr')} | "
                f"{'Y' if r.get('election_gap') else 'n'} | {r.get('min_cf_head_rank')} | "
                f"{r.get('g_legacy_heads'):.2f} | {r.get('g_strict_heads'):.2f} |"
            )
        lines.append("")
        args.out_md.expanduser().write_text("\n".join(lines) + "\n")
        print(f"Wrote {args.out_md}")
    print(f"Wrote {args.out_json}")
    print(f"election_gap={len(gaps)} genuine_goods={genuine_goods}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
