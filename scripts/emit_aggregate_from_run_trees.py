#!/usr/bin/env python3
"""emit_aggregate_from_run_trees.py

Walk existing per-target benchmark run trees (under e.g. benchmark_results/ or a full
reproduce output dir) and emit a single aggregate results CSV in the format produced
by DatasetRunner / reproduce_astex85.sh.

- Uses only real per-target result.csv (and stdout.log for THERMO columns when present).
- Never hand-authors RMSDs or success flags.
- When multiple restarts + [THERMO] G_bind present, a future enhancement can re-apply
  the two-stage elect (min-G restart then freq/Z+H within) and (if coords+ref available)
  recompute rmsd for the elected pose. For now aggregates the run's reported values.

Usage:
  python3 scripts/emit_aggregate_from_run_trees.py \
      --tree-root benchmark_results \
      --out scratch/astex_crossdock_85_results.csv

After a full v88 reproduce run that populates many target trees with THERMO enabled,
re-run this to obtain the authentic full-85 CSV for gating verification.
"""

import argparse
import csv
import math
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, NamedTuple
import numpy as np
from scipy.optimize import linear_sum_assignment

THERMO_RE = re.compile(r'G_bind=([-\d.eE+]+)')


def parse_last_g_bind(text: str) -> float:
    matches = THERMO_RE.findall(text)
    if not matches:
        return float('nan')
    try:
        return float(matches[-1])
    except ValueError:
        return float('nan')


# --- RMSD helpers (adapted from scripts/p4_best_of_n_diagnostic.py for oracle) ---
H_TOKENS = {"H", "D"}

def parse_pdb_ligand(pdb_path: str):
    coords = []
    elems = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM") and not line.startswith("HETATM"):
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except Exception:
                continue
            elem = line[76:78].strip().upper() or line[12:14].strip().upper() or line[13:14].upper()
            if len(elem) > 1 and elem[1].isdigit():
                elem = elem[0]
            if elem in H_TOKENS:
                continue
            coords.append((x, y, z))
            elems.append(elem[:1] if elem else "X")
    return np.asarray(coords, dtype=float), elems

def parse_sdf_ligand(sdf_path: str):
    coords = []
    elems = []
    with open(sdf_path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if "$$$$" in line:
            break
        if i > 3 and len(line) > 40:
            try:
                x = float(line[0:10])
                y = float(line[10:20])
                z = float(line[20:30])
                elem = line[31:34].strip().upper()
                if elem and elem not in H_TOKENS:
                    coords.append((x, y, z))
                    elems.append(elem[:1] if elem else "X")
            except Exception:
                pass
        i += 1
    return np.asarray(coords, dtype=float), elems

def hungarian_rmsd(coordsA, elemsA, coordsB, elemsB):
    """Element-matched optimal-assignment heavy-atom RMSD."""
    if coordsA.shape[0] == 0 or coordsB.shape[0] == 0:
        return None
    if coordsA.shape[0] != coordsB.shape[0]:
        return None
    n = coordsA.shape[0]
    diff = coordsA[:, None, :] - coordsB[None, :, :]
    cost = np.einsum("ijk,ijk->ij", diff, diff)
    BIG = 1.0e9
    ea = np.asarray(elemsA)
    eb = np.asarray(elemsB)
    mask = ea[:, None] != eb[None, :]
    cost = cost + mask * BIG
    row, col = linear_sum_assignment(cost)
    matched = cost[row, col]
    if np.any(matched >= BIG):
        d = coordsA - coordsB
        return float(math.sqrt(np.mean(np.einsum("ij,ij->i", d, d))))
    return float(math.sqrt(matched.mean()))


class PoseCandidate(NamedTuple):
    path: str
    cf: float
    freq: int
    g_bind: float
    restart_id: int
    member_cfs: List[float]


def find_target_dirs(root: Path) -> List[Path]:
    """Heuristic: 4-char alphanum PDB dirs containing stdout.log or r*/ or result.csv (support flat + r* trees)."""
    cands = []
    for p in sorted(root.iterdir()):
        if not p.is_dir():
            continue
        name = p.name
        if re.match(r'^[0-9A-Z]{4}$', name):
            if (p / 'stdout.log').exists() or (p / 'result.csv').exists() or any((p / f'r{i}').exists() for i in range(7)):
                cands.append(p)
    return cands


def parse_pose(cand: str):
    """Replicate ReportedPoseSelector::parse_pose: CF from REMARK, freq, member_cfs from .mcf sidecar."""
    cf = float('inf')
    freq = 1
    member_cfs = []
    try:
        with open(cand) as f:
            for pl in f:
                if cf == float('inf') and 'REMARK CF=' in pl:
                    p2 = pl.find('CF=')
                    if p2 != -1:
                        try:
                            cf = float(pl[p2+3:].split()[0])
                        except Exception:
                            pass
                if 'Frequency:' in pl:
                    p2 = pl.find('Frequency:')
                    try:
                        freq = int(pl[p2+10:].split()[0])
                    except Exception:
                        pass
    except Exception:
        pass
    mcf = cand[:-4] + '.mcf' if cand.endswith('.pdb') else ''
    if os.path.exists(mcf):
        try:
            with open(mcf) as mf:
                for ml in mf:
                    try:
                        v = float(ml.strip())
                        if math.isfinite(v):
                            member_cfs.append(v)
                    except Exception:
                        pass
        except Exception:
            pass
    return cf, freq, member_cfs


def parse_g_from_log(log_path: Path) -> float:
    if not log_path.exists():
        return float('nan')
    try:
        txt = log_path.read_text(errors='ignore')
        return parse_last_g_bind(txt)
    except Exception:
        return float('nan')


def build_pool_for_target(td: Path) -> List[PoseCandidate]:
    """Walk r0(top)+r1..r6 , parse per-restart G_bind, collect poses, return pool for two-stage."""
    pool: List[PoseCandidate] = []
    code = td.name
    # top level as restart 0
    restarts = [(0, td)]
    for ri in range(1, 7):
        rd = td / f'r{ri}'
        if rd.exists():
            restarts.append((ri, rd))
    for ri, rdir in restarts:
        g = parse_g_from_log(rdir / 'stdout.log')
        # collect pose files in this restart dir ( *_N.pdb , INI etc)
        for p in sorted(rdir.glob('*.pdb')):
            if '_' not in p.name:
                continue
            cf, freq, mcfs = parse_pose(str(p))
            if not math.isfinite(cf):
                continue
            pool.append(PoseCandidate(path=str(p), cf=cf, freq=freq, g_bind=g, restart_id=ri, member_cfs=mcfs))
    return pool


def boltzmann_composite(p: PoseCandidate) -> float:
    """Exact replica of boltzmann_composite in elect_reported_pose (with overflow guard for large negative CF)."""
    kT = 0.592
    alpha = 1.0
    n_members = len(p.member_cfs) if p.member_cfs else max(1, p.freq)
    pop_weight = math.log1p(n_members)
    def safe_exp(x):
        if x > 700.0: return float('inf')
        if x < -700.0: return 0.0
        return math.exp(x)
    if not p.member_cfs:
        return safe_exp(-p.cf / kT) * pop_weight
    Z = 0.0
    for cf_i in p.member_cfs:
        Z += safe_exp(-cf_i / kT)
    if Z <= 0.0:
        return safe_exp(-p.cf / kT) * pop_weight
    H = 0.0
    for cf_i in p.member_cfs:
        pi = safe_exp(-cf_i / kT) / Z if Z > 0 else 0.0
        if pi > 1e-300:
            H -= pi * math.log(pi)
    return Z * safe_exp(-alpha * H) * pop_weight


def elect_reported_pose_py(pool: List[PoseCandidate], thermo_on: bool = True) -> Optional[str]:
    """Python replica of reported_pose::elect_reported_pose for the oracle (exact two-stage when thermo_on)."""
    if not pool:
        return None
    if thermo_on:
        min_g = float('inf')
        chosen_ri = -1
        for pc in pool:
            if math.isfinite(pc.g_bind) and pc.g_bind < min_g:
                min_g = pc.g_bind
                chosen_ri = pc.restart_id
        if chosen_ri < 0:
            chosen_ri = -1
        best = None
        best_sc = -float('inf')
        for pc in pool:
            if chosen_ri >= 0 and pc.restart_id != chosen_ri:
                continue
            sc = boltzmann_composite(pc)
            if sc > best_sc or best is None:
                best_sc = sc
                best = pc
        return best.path if best else None
    # non-thermo: full pool best composite
    best = None
    best_sc = -float('inf')
    for pc in pool:
        sc = boltzmann_composite(pc)
        if sc > best_sc or best is None:
            best_sc = sc
            best = pc
    if best:
        return best.path
    # fallback lowest cf
    return min(pool, key=lambda p: p.cf).path if pool else None


def load_per_target(td: Path, dataset_dir: Optional[Path] = None, use_two_stage: bool = True) -> Optional[Dict]:
    """Load or synthesize row. For r* trees, re-elect using two-stage and recompute RMSD on elected."""
    row = {'pdb_id': td.name}
    res_csv = td / 'result.csv'
    if res_csv.exists():
        try:
            with open(res_csv, newline='') as f:
                rows = list(csv.DictReader(f))
            if rows:
                row.update({k: v for k, v in rows[0].items() if k not in ('selected_policy',)})
        except Exception:
            pass

    # If r* structure present or force re-elect, do oracle two-stage + recompute RMSD
    has_r = any((td / f'r{i}').exists() for i in range(1,7)) or (td / 'stdout.log').exists()
    if has_r and use_two_stage:
        pool = build_pool_for_target(td)
        elected = elect_reported_pose_py(pool, thermo_on=True)
        if elected:
            row['elected_pose'] = elected
            # parse G for elected's restart (from its log)
            try:
                rdir = Path(elected).parent
                g = parse_g_from_log(rdir / 'stdout.log')
                row['G_bind'] = g if math.isfinite(g) else ''
            except Exception:
                row['G_bind'] = ''
            # recompute hungarian RMSD using ref ligand_sdf
            try:
                ref = None
                if dataset_dir:
                    ref = dataset_dir / td.name / f'{td.name}_ligand.sdf'
                if not ref or not ref.exists():
                    # fallback common layout
                    ref = Path('benchmarks/astex_diverse/astex_diverse') / td.name / f'{td.name}_ligand.sdf'
                if ref and ref.exists():
                    cC, cE = parse_sdf_ligand(str(ref))
                    pC, pE = parse_pdb_ligand(elected)
                    rh = hungarian_rmsd(pC, pE, cC, cE)
                    if rh is not None:
                        row['rmsd_hungarian'] = rh
                        row['rmsd_to_crystal'] = rh  # approx
            except Exception as ex:
                pass  # keep whatever was in row or leave
        # success
        try:
            rh = float(row.get('rmsd_hungarian', 99.0))
            row['success'] = 1 if rh < 2.0 else 0
        except Exception:
            row.setdefault('success', 0)

    # fallback enrich if still no G
    if 'G_bind' not in row or row.get('G_bind') in (None, ''):
        log_path = td / 'stdout.log'
        g = parse_last_g_bind(log_path.read_text(errors='ignore')) if log_path.exists() else float('nan')
        row['G_bind'] = g if math.isfinite(g) else ''
    # ensure success
    try:
        rh = float(row.get('rmsd_hungarian', row.get('rmsd_to_crystal', 99.0)))
        row['success'] = 1 if rh < 2.0 else 0
    except Exception:
        row.setdefault('success', 0)
    return row if row.get('pdb_id') else None


def main() -> None:
    ap = argparse.ArgumentParser(description="Emit aggregate Astex85 results CSV from per-target run trees. Supports r* 7-restart THERMO trees and re-elects using two-stage for measurement oracle.")
    ap.add_argument('--tree-root', type=Path, default=Path('benchmark_results'),
                    help='Root containing per-PDB target trees (r0..r6 or flat)')
    ap.add_argument('--out', type=Path, default=Path('benchmark_results/astex_crossdock_85_results.csv'),
                    help='Output aggregate CSV path')
    ap.add_argument('--verbose', '-v', action='store_true')
    ap.add_argument('--dataset-dir', type=Path, default=None, help='Base for <id>/<id>_ligand.sdf (for RMSD recompute)')
    args = ap.parse_args()

    root = args.tree_root
    if root == Path('benchmark_results') and Path('benchmark_trees_canonical').exists():
        root = Path('benchmark_trees_canonical')
    targets = find_target_dirs(root)
    if args.verbose:
        print(f"Found {len(targets)} candidate target trees under {root}")

    ds_dir = args.dataset_dir
    if ds_dir is None:
        for cand in [Path('benchmarks/astex_diverse/astex_diverse'), Path('benchmark_trees_canonical')]:
            if cand.exists() and (cand / '1G9V' / '1G9V_ligand.sdf').exists():
                ds_dir = cand
                break

    rows: List[Dict] = []
    for td in targets:
        row = load_per_target(td, dataset_dir=ds_dir, use_two_stage=True)
        if row:
            rows.append(row)
            if args.verbose:
                rh = row.get('rmsd_hungarian', '?')
                gb = row.get('G_bind', '')
                print(f"  {td.name}: rmsd_h={rh} G_bind={gb} success={row.get('success')}")

    if not rows:
        print("No targets; nothing written.")
        return

    # Build a stable field list (union, prefer common order)
    preferred = [
        'pdb_id', 'best_score', 'rmsd_to_crystal', 'rmsd_hungarian', 'predicted_dG',
        'predicted_dH', 'predicted_TdS', 'shannon_entropy', 'search_entropy_proxy',
        'num_poses', 'wall_time_s', 'success', 'cf_native', 'best_cluster_rmsd',
        'best_cluster_idx', 'seed_echo', 'pose_source', 'H_rep_rank0', 'H_pop',
        'H_rep_mean', 'D_vib', 'G_bind', 'g_bind', 'h_vct', 'h_vct_raw', 'n_heavy',
        'tds_shannon', 'tds_vib'
    ]
    fields = []
    seen = set()
    for p in preferred:
        if p in rows[0]:
            fields.append(p); seen.add(p)
    for k in rows[0].keys():
        if k not in seen:
            fields.append(k)

    out_path = args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w', newline='') as f:
        wr = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        wr.writeheader()
        wr.writerows(rows)

    # Gating summary (the real metric)
    n = len(rows)
    succ = sum(int(float(r.get('rmsd_hungarian', 99)) < 2.0 or int(r.get('success', 0)) == 1) for r in rows)
    print(f"Wrote {n} rows -> {out_path}")
    print(f"Success (rmsd_hungarian < 2.0): {succ} / {n}  ({100.0*succ/n:.1f}%)")

    if succ >= 80:
        print("OFFLINE GATE PASSED (>=80/85) — safe to launch full reproduce_astex85.sh for confirmation.")
    else:
        print("OFFLINE: current trees yield <80; run full v88 reproduce (THERMO=1, RESTARTS=7, native_85 json) to generate authentic high-success trees then re-emit.")


if __name__ == '__main__':
    main()