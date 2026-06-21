#!/usr/bin/env python3
"""Post-run analysis for v106 Tier-2 OCD campaign.

Cross-crystal RMSD coordinate-frame fix
---------------------------------------
The benchmark runner writes ``rmsd_to_crystal`` / ``rmsd_hungarian`` by
comparing the docked pose (expressed in the ACCEPTOR receptor frame) against the
donor ligand's crystal coordinates in the DONOR crystal frame.  For the
cyclic-shift oracle cross-dock set the receptor and the ligand come from two
*different* (non-cognate) Astex crystals, so those two frames differ by the
rigid-body offset between the crystals -- up to ~50 A.  The raw RMSD is then
dominated by that coordinate-frame mismatch, not by docking quality: e.g.
1GPK<-1HNN reports 41-92 A even when the pose is correctly seated in the 1GPK
pocket (pose centroid ~(4,63,66) == 1GPK pocket; donor ligand ~(28,44,16)).

This analysis recomputes RMSD in a common frame:

1. Load the donor binding site (same frame as the donor reference ligand) and
   the acceptor binding site (== oracle_site_pdb, same frame as the docked pose).
2. Superimpose donor site -> acceptor site.  The cyclic-shift pairs are
   non-cognate and share no residue numbering, so this is a sequence-independent
   CA superposition (multi-start Kabsch/ICP with nearest-neighbour
   correspondence), not a 1:1 residue match.
3. Apply the resulting rotation+translation to the donor reference ligand to
   express it in the acceptor frame.
4. Measure symmetry-aware (Hungarian, per-element) RMSD between the docked pose
   and the transformed reference.

The binding-site fit residual is reported alongside each RMSD: because the two
pockets are different proteins, it is the physical floor below which a
cross-dock RMSD cannot meaningfully fall.

Copyright 2026 Le Bonhomme Pharma. Apache-2.0.
"""

from __future__ import annotations

import csv
import glob
import itertools
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

CAMPAIGN = Path(
    os.environ.get(
        "FLEXAIDDS_CAMPAIGN",
        "/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results/"
        "v106_20260621_ocd_noseed_consensus5r_ecb85e7",
    )
)
ORACLE_DIR = Path(
    os.environ.get(
        "FLEXAIDDS_ORACLE_SITE_DIR",
        "/Users/lp.more/Projects/FlexAIDdS/benchmarks/astex_diverse/astex_diverse",
    )
)
OUT = CAMPAIGN / "analysis"
THRESHOLD = 2.0


# --------------------------------------------------------------------------- #
# Geometry helpers
# --------------------------------------------------------------------------- #
def ca_coords(pdb: Path) -> np.ndarray:
    """CA atom coordinates from a binding-site PDB."""
    pts = []
    for line in pdb.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA":
            pts.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return np.asarray(pts, dtype=float)


def read_sdf(path: Path) -> tuple[np.ndarray, list[str]]:
    """Coordinates and element symbols from a V2000 SDF (donor reference ligand)."""
    lines = path.read_text().splitlines()
    n = int(lines[3][:3])
    xyz, els = [], []
    for i in range(n):
        rec = lines[4 + i]
        xyz.append((float(rec[0:10]), float(rec[10:20]), float(rec[20:30])))
        els.append(rec[31:34].strip().upper())
    return np.asarray(xyz, dtype=float), els


def pose_ligand(pdb: Path, elements: list[str]) -> np.ndarray | None:
    """Docked-ligand coordinates from a full-complex pose PDB.

    FlexAID appends the docked ligand as the trailing HETATM block (sentinel
    serial 90001+), after the receptor and any cofactors (e.g. HEM).  The custom
    serial overflows the fixed-width columns, so coordinates are parsed by
    whitespace (last 3 floats before occupancy/bfactor/element).  The block is
    located by matching the SDF element sequence, which also guards against
    picking up a cofactor of the same length.
    """
    het = [ln for ln in pdb.read_text().splitlines() if ln.startswith("HETATM")]
    coords, els = [], []
    for ln in het:
        toks = ln.split()
        try:
            x, y, z = (float(t) for t in toks[-6:-3])
        except ValueError:
            continue
        coords.append((x, y, z))
        els.append(toks[-1].upper())
    n = len(elements)
    target = [e.upper() for e in elements]
    # Prefer an exact contiguous element-sequence match (handles cofactors).
    for i in range(len(els) - n + 1):
        if els[i : i + n] == target:
            return np.asarray(coords[i : i + n], dtype=float)
    # Fallback: trailing n HETATM records if counts line up.
    if len(coords) >= n:
        return np.asarray(coords[-n:], dtype=float)
    return None


def kabsch(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotation R and translation t mapping P onto Q (equal length, ordered)."""
    cP, cQ = P.mean(0), Q.mean(0)
    H = (P - cP).T @ (Q - cQ)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, cQ - R @ cP


def _pca_axes(X: np.ndarray) -> np.ndarray:
    X0 = X - X.mean(0)
    _, V = np.linalg.eigh(X0.T @ X0)
    return V[:, ::-1]


def superpose_site(src: np.ndarray, dst: np.ndarray, iters: int = 80) -> tuple[np.ndarray, np.ndarray, float]:
    """Sequence-independent rigid superposition of binding-site CA clouds.

    Multi-start ICP: each start orients src principal axes onto dst principal
    axes (8 sign combinations + identity), then refines by nearest-neighbour
    Kabsch.  Returns the rotation, translation and per-atom fit residual of the
    best start.
    """
    Vs, Vd = _pca_axes(src), _pca_axes(dst)
    inits = [np.eye(3)]
    for signs in itertools.product([1, -1], repeat=3):
        R0 = Vd @ np.diag(signs) @ Vs.T
        if np.linalg.det(R0) > 0:
            inits.append(R0)
    best: tuple[float, np.ndarray, np.ndarray] | None = None
    for R0 in inits:
        R, t = R0, dst.mean(0) - R0 @ src.mean(0)
        cur = (R @ src.T).T + t
        for _ in range(iters):
            idx = np.argmin(((cur[:, None, :] - dst[None, :, :]) ** 2).sum(2), axis=1)
            R, t = kabsch(src, dst[idx])
            cur = (R @ src.T).T + t
        idx = np.argmin(((cur[:, None, :] - dst[None, :, :]) ** 2).sum(2), axis=1)
        res = float(np.sqrt(((cur - dst[idx]) ** 2).sum(1).mean()))
        if best is None or res < best[0]:
            best = (res, R, t)
    assert best is not None
    return best[1], best[2], best[0]


def rmsd_hungarian(pose: np.ndarray, ref: np.ndarray, els: list[str]) -> float:
    """Symmetry-aware RMSD via optimal within-element atom assignment."""
    n = len(pose)
    cost = np.full((n, n), 1e9)
    el = [e.upper() for e in els]
    for i in range(n):
        for j in range(n):
            if el[i] == el[j]:
                cost[i, j] = ((pose[i] - ref[j]) ** 2).sum()
    r, c = linear_sum_assignment(cost)
    return float(np.sqrt(((pose[r] - ref[c]) ** 2).sum(1).mean()))


def rmsd_direct(pose: np.ndarray, ref: np.ndarray) -> float:
    return float(np.sqrt(((pose - ref) ** 2).sum(1).mean()))


# --------------------------------------------------------------------------- #
# Per-target frame-corrected RMSD
# --------------------------------------------------------------------------- #
def corrected_rmsd(result_dir: Path) -> dict:
    """Recompute frame-corrected RMSD for one OCD target directory.

    The directory is named by the acceptor (receptor) PDB id; dock_config.json
    names the donor reference ligand.  Both binding sites live under ORACLE_DIR.
    """
    recep = result_dir.name
    cfg = json.loads((result_dir / "dock_config.json").read_text())
    donor_sdf = Path(cfg["reference_ligand"]["file"])
    donor = donor_sdf.parent.name

    accep_site_pdb = ORACLE_DIR / recep / f"{recep}_binding_site.pdb"
    donor_site_pdb = ORACLE_DIR / donor / f"{donor}_binding_site.pdb"
    out: dict = {"target": recep, "donor": donor, "pair": f"{recep}<-{donor}"}

    for p in (donor_sdf, accep_site_pdb, donor_site_pdb):
        if not p.exists():
            out["error"] = f"missing {p}"
            return out

    donor_ca, accep_ca = ca_coords(donor_site_pdb), ca_coords(accep_site_pdb)
    if len(donor_ca) < 3 or len(accep_ca) < 3:
        out["error"] = "insufficient CA atoms for superposition"
        return out

    R, t, fit = superpose_site(donor_ca, accep_ca)
    lig_xyz, els = read_sdf(donor_sdf)
    ref = (R @ lig_xyz.T).T + t  # donor reference ligand in acceptor frame
    out["pocket_fit_residual"] = round(fit, 3)
    out["n_heavy"] = len(els)

    pose_files = sorted(
        f
        for f in glob.glob(str(result_dir / f"{recep}_*.pdb"))
        if not f.endswith(("_INI.pdb", "_prepped.pdb"))
    )
    per_pose = []
    for f in pose_files:
        pose = pose_ligand(Path(f), els)
        if pose is None or len(pose) != len(ref):
            continue
        per_pose.append(
            {
                "file": Path(f).name,
                "rank": _rank_of(Path(f).name, recep),
                "rmsd_hungarian": round(rmsd_hungarian(pose, ref, els), 3),
                "rmsd_direct": round(rmsd_direct(pose, ref), 3),
            }
        )
    if not per_pose:
        out["error"] = "no parseable docked-ligand poses"
        return out

    rank0 = min(per_pose, key=lambda d: (d["rank"] if d["rank"] is not None else 1e9))
    best = min(per_pose, key=lambda d: d["rmsd_hungarian"])
    out["corrected_rmsd_rank0"] = rank0["rmsd_hungarian"]
    out["corrected_best_pose_rmsd"] = best["rmsd_hungarian"]
    out["best_pose_file"] = best["file"]
    out["n_poses_scored"] = len(per_pose)
    return out


def _rank_of(fname: str, recep: str) -> int | None:
    stem = fname[len(recep) + 1 : -4]  # "<recep>_<n>.pdb" -> "<n>"
    return int(stem) if stem.isdigit() else None


# --------------------------------------------------------------------------- #
def main() -> int:
    ad = CAMPAIGN / "astex_diverse"
    rows = []
    for rc in sorted(ad.glob("*/result.csv")):
        result_dir = rc.parent
        raw = list(csv.DictReader(rc.open()))[0]
        raw_hung = float(raw["rmsd_hungarian"])
        tds = float(raw.get("predicted_TdS", raw.get("TdS_shannon", 0)) or 0)
        corr = corrected_rmsd(result_dir)
        sel = corr.get("corrected_rmsd_rank0")
        best = corr.get("corrected_best_pose_rmsd")
        rows.append(
            {
                **corr,
                "raw_rmsd_hungarian": round(raw_hung, 3),
                "predicted_TdS": tds,
                "success": sel is not None and sel < THRESHOLD,
                "search_found": best is not None and best < THRESHOLD,
                "select_fail": (
                    best is not None and best < THRESHOLD and (sel is None or sel >= THRESHOLD)
                ),
            }
        )

    scored = [r for r in rows if "corrected_best_pose_rmsd" in r]
    n = len(scored)
    succ = sum(1 for r in scored if r["success"])
    search = sum(1 for r in scored if r["search_found"])
    sel_fail = sum(1 for r in scored if r["select_fail"])
    floors = [r["pocket_fit_residual"] for r in scored if "pocket_fit_residual" in r]

    report = {
        "generated_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "tier": "TIER-2",
        "campaign": str(CAMPAIGN),
        "rmsd_metric": "frame_corrected_oracle_crossdock_hungarian",
        "frame_fix": (
            "Donor reference ligand superimposed into the acceptor receptor frame "
            "via sequence-independent CA superposition of the donor binding site "
            "onto the acceptor (oracle) binding site before RMSD."
        ),
        "n_complete": len(rows),
        "n_scored": n,
        "success": succ,
        "success_rate": round(succ / n, 4) if n else 0.0,
        "search_found_lt_2A": search,
        "selection_failures": sel_fail,
        "search_failures": n - search,
        "pocket_fit_residual_floor": {
            "min": round(min(floors), 3) if floors else None,
            "median": round(float(np.median(floors)), 3) if floors else None,
            "max": round(max(floors), 3) if floors else None,
            "note": (
                "Non-cognate pockets: cross-dock RMSD cannot meaningfully fall "
                "below the binding-site superposition residual."
            ),
        },
        "comparators": {
            "v88_published_seeded": "79/85",
            "v102_autonomous": "8/85",
        },
        "per_target": rows,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "v106_ocd_summary.json").write_text(json.dumps(report, indent=2) + "\n")

    md = [
        "# v106 Tier-2 OCD Analysis (frame-corrected RMSD)",
        "",
        f"Generated: {report['generated_at']}",
        "",
        "RMSD is computed after superimposing each donor reference ligand into the",
        "acceptor receptor frame (sequence-independent CA superposition of the donor",
        "binding site onto the oracle binding site). Raw `rmsd_hungarian` from the",
        "runner compares mismatched crystal frames (up to ~50 A apart) and is shown",
        "for reference only.",
        "",
        "## Headline",
        "",
        "| Metric | Value |",
        "|--------|------:|",
        f"| Scored | {n}/{len(rows)} |",
        f"| **Corrected BCR (best pose < 2A)** | **{search}/{n} ({100*search/n if n else 0:.1f}%)** |",
        f"| Selected pose < 2A | {succ}/{n} |",
        f"| Selection failures | {sel_fail} |",
        f"| Pocket-fit floor (min/median/max A) | "
        f"{report['pocket_fit_residual_floor']['min']}/"
        f"{report['pocket_fit_residual_floor']['median']}/"
        f"{report['pocket_fit_residual_floor']['max']} |",
        "",
        "## Per-target",
        "",
        "| Pair (recv<-donor) | raw hung (A) | corrected sel (A) | corrected best (A) | pocket fit (A) |",
        "|--------------------|-------------:|------------------:|-------------------:|---------------:|",
    ]
    for r in rows:
        md.append(
            f"| {r.get('pair', r.get('target', '?'))} "
            f"| {r.get('raw_rmsd_hungarian', 'NA')} "
            f"| {r.get('corrected_rmsd_rank0', r.get('error', 'NA'))} "
            f"| {r.get('corrected_best_pose_rmsd', 'NA')} "
            f"| {r.get('pocket_fit_residual', 'NA')} |"
        )
    md += [
        "",
        "Note: with non-cognate (cyclic-shift) pockets the pocket-fit residual is the",
        "physical floor for cross-dock RMSD; a 'correctly pocketed' pose is bounded",
        "below by it, so sub-2A hits are not expected from this oracle metric.",
        "",
    ]
    (OUT / "v106_ocd_report.md").write_text("\n".join(md))
    print((OUT / "v106_ocd_report.md").read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
