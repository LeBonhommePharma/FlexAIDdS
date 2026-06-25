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
import re
from pathlib import Path
from typing import Dict, List, Optional


THERMO_RE = re.compile(r'G_bind=([-\d.eE+]+)')


def parse_last_g_bind(text: str) -> float:
    matches = THERMO_RE.findall(text)
    if not matches:
        return float('nan')
    try:
        return float(matches[-1])
    except ValueError:
        return float('nan')


def find_target_dirs(root: Path) -> List[Path]:
    """Heuristic: 4-char alphanum PDB dirs containing result.csv or stdout.log."""
    cands = []
    for p in sorted(root.iterdir()):
        if not p.is_dir():
            continue
        name = p.name
        if re.match(r'^[0-9A-Z]{4}$', name):
            if (p / 'result.csv').exists() or (p / 'stdout.log').exists():
                cands.append(p)
    return cands


def load_per_target(td: Path) -> Optional[Dict]:
    res_csv = td / 'result.csv'
    if not res_csv.exists():
        return None
    try:
        with open(res_csv, newline='') as f:
            rows = list(csv.DictReader(f))
        if not rows:
            return None
        row = dict(rows[0])  # take the (only) row for this target
        # Normalize keys we care about
        row['pdb_id'] = row.get('pdb_id') or td.name
        # Try enrich thermo from this target's stdout (if single) or look for per-r logs
        log_path = td / 'stdout.log'
        g = float('nan')
        h_vct = float('nan')
        h_raw = float('nan')
        n_heavy = 0
        tds_sh = float('nan')
        tds_vib = float('nan')
        if log_path.exists():
            # efficient: only tail for last THERMO record (logs can be >700k lines)
            try:
                with open(log_path, 'rb') as lf:
                    lf.seek(0, 2)
                    size = lf.tell()
                    tail = 65536
                    lf.seek(max(0, size - tail))
                    txt = lf.read().decode('utf-8', errors='ignore')
            except Exception:
                txt = log_path.read_text(errors='ignore')[-65536:]
            g = parse_last_g_bind(txt)
            # crude extra pulls if present in same line
            m = re.search(r'H_vct=([-\d.eE+]+)', txt)
            if m: h_vct = float(m.group(1))
            m = re.search(r'H_vct_raw=([-\d.eE+]+)', txt)
            if m: h_raw = float(m.group(1))
            m = re.search(r'n_heavy=(\d+)', txt)
            if m: n_heavy = int(m.group(1))
            m = re.search(r'TdS_shannon=([-\d.eE+]+)', txt)
            if m: tds_sh = float(m.group(1))
            m = re.search(r'TdS_vib=([-\d.eE+]+)', txt)
            if m: tds_vib = float(m.group(1))
        row.setdefault('G_bind', g)
        row.setdefault('g_bind', g)
        row.setdefault('h_vct', h_vct)
        row.setdefault('h_vct_raw', h_raw)
        row.setdefault('n_heavy', n_heavy)
        row.setdefault('tds_shannon', tds_sh)
        row.setdefault('tds_vib', tds_vib)
        # success convenience
        try:
            rh = float(row.get('rmsd_hungarian', row.get('rmsd_to_crystal', '9.0')))
            row['success'] = 1 if rh < 2.0 else 0
        except Exception:
            row.setdefault('success', 0)
        return row
    except Exception:
        return None


def main() -> None:
    ap = argparse.ArgumentParser(description="Emit aggregate Astex85 results CSV from per-target run trees.")
    ap.add_argument('--tree-root', type=Path, default=Path('benchmark_results'),
                    help='Root containing per-PDB target trees (each with result.csv)')
    ap.add_argument('--out', type=Path, default=Path('benchmark_results/astex_crossdock_85_results.csv'),
                    help='Output aggregate CSV path')
    ap.add_argument('--verbose', '-v', action='store_true')
    args = ap.parse_args()

    targets = find_target_dirs(args.tree_root)
    if args.verbose:
        print(f"Found {len(targets)} candidate target trees under {args.tree_root}")

    rows: List[Dict] = []
    for td in targets:
        row = load_per_target(td)
        if row:
            rows.append(row)
            if args.verbose:
                rh = row.get('rmsd_hungarian', '?')
                print(f"  {td.name}: rmsd_h={rh} success={row.get('success')}")

    if not rows:
        print("No per-target result.csv found; nothing written.")
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