#!/usr/bin/env python3
"""
election_cf_diagnostic.py — CF-per-cluster-representative election artifact
================================================================================
FlexAIDdS elects the docked pose it reports (rank-0 head) as the single
lowest-CF cluster representative. The dominant remaining Astex-Diverse failure
mode is the *election gap*: a populous native basin whose head is not the
single lowest-CF pose loses the election to a deep-CF singleton from a small
cluster.

The raw ingredients of that election — each cluster's single-pose CF, its
population (`Frequency`) and its average CF — are written by the engine into
the cluster PDB REMARK block:

    REMARK CF=-189.85613
    REMARK Cluster 0: Rank (top):0 Average CF:-49.30565 Frequency:4

but they are *not* surfaced in the DatasetRunner `result.csv`, so the election
gap cannot be diagnosed from the tracked benchmark output. This tool extracts
them into a CSV and prototypes the proposed frequency x cluster-size soft-beta
reweighted election **offline**, so the fix can be validated before touching
the C++ election path.

The reweighted election scores each cluster by a Boltzmann soft-CF weight over
the cluster's representative energy (its average CF — a cluster-level, not
singleton, quantity) multiplied by its population:

    w_i = frequency_i * exp(-beta * (avgCF_i - avgCF_min))

Stdlib only (no numpy/pandas/scipy) so it runs in CI and on any Python 3.8+.

Usage:
    python3 scripts/election_cf_diagnostic.py <results_dir> [--out election_cf.csv]
                                              [--beta 0.02] [--n-clusters 10]
"""

import os
import re
import csv
import sys
import glob
import math
import argparse

# REMARK Cluster 0: Rank (top):0 Average CF:-49.30565 Frequency:4
_CLUSTER_RE = re.compile(
    r'REMARK Cluster\s+(?P<cluster>\d+):\s*Rank\s*\(top\):\s*(?P<rank>\d+)\s*'
    r'Average CF:\s*(?P<avgcf>[-\d.eE+]+)\s*Frequency:\s*(?P<freq>\d+)'
)


def parse_pose_remarks(pdb_path):
    """Extract (single CF, cluster avg CF, frequency, cluster top-rank) from a
    cluster PDB. Any field the engine did not write comes back as None."""
    cf = avg_cf = freq = top_rank = None
    try:
        with open(pdb_path) as fh:
            for line in fh:
                if cf is None and line.startswith('REMARK CF='):
                    try:
                        cf = float(line.split('=', 1)[1].strip())
                    except ValueError:
                        pass
                elif line.startswith('REMARK Cluster'):
                    m = _CLUSTER_RE.match(line.strip())
                    if m:
                        top_rank = int(m.group('rank'))
                        avg_cf = float(m.group('avgcf'))
                        freq = int(m.group('freq'))
                    break
    except OSError:
        pass
    return cf, avg_cf, freq, top_rank


def collect_target(target_dir, n_clusters=10):
    """Return ordered per-rank rows for one target, or [] if no cluster PDBs."""
    pdb_id = os.path.basename(target_dir.rstrip('/'))
    pattern = os.path.join(target_dir, f'{pdb_id}_*.pdb')

    def rank_of(p):
        try:
            return int(os.path.basename(p).replace(f'{pdb_id}_', '').replace('.pdb', ''))
        except ValueError:
            return -1

    pdbs = sorted(
        (p for p in glob.glob(pattern) if '_INI' not in os.path.basename(p)
         and rank_of(p) >= 0),
        key=rank_of,
    )[:n_clusters]

    rows = []
    for p in pdbs:
        cf, avg_cf, freq, top_rank = parse_pose_remarks(p)
        rows.append({
            'pdb_id': pdb_id,
            'rank': rank_of(p),
            'cf': cf,
            'cluster_avg_cf': avg_cf,
            'cluster_frequency': freq,
            'cluster_top_rank': top_rank,
        })
    return rows


def reweighted_rank(rows, beta):
    """Rank the frequency x cluster-size soft-beta election would elect.

    Uses each cluster's average CF (cluster-level representative energy),
    falling back to the head's single-pose CF when the cluster REMARK is
    absent. Returns None if no row carries usable weights."""
    def energy(r):
        return r['cluster_avg_cf'] if r['cluster_avg_cf'] is not None else r['cf']

    usable = [r for r in rows if energy(r) is not None and r['cluster_frequency']]
    if not usable:
        return None
    e_min = min(energy(r) for r in usable)
    best, best_w = None, -1.0
    for r in usable:
        w = float(r['cluster_frequency']) * math.exp(-beta * (energy(r) - e_min))
        if w > best_w:
            best_w, best = w, r['rank']
    return best


def summarize_target(rows, beta):
    """One summary row per target with the election-gap diagnostics."""
    head = next((r for r in rows if r['rank'] == 0), rows[0])
    cf_valid = [(r['rank'], r['cf']) for r in rows if r['cf'] is not None]
    cf_argmin = min(cf_valid, key=lambda x: x[1])[0] if cf_valid else None
    rw = reweighted_rank(rows, beta)
    return {
        'pdb_id': head['pdb_id'],
        'n_clusters': len(rows),
        'CF_rep_rank0': head['cf'],
        'cluster_freq_rank0': head['cluster_frequency'],
        'cluster_avgcf_rank0': head['cluster_avg_cf'],
        'CF_rep_argmin_rank': cf_argmin,
        'reweighted_rank0': rw,
        'cf_election_changed': int(rw is not None and rw != 0),
        'CF_rep_list': [r['cf'] for r in rows],
        'cluster_freq_list': [r['cluster_frequency'] for r in rows],
        'cluster_avgcf_list': [r['cluster_avg_cf'] for r in rows],
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('results_dir', help='DatasetRunner output directory')
    ap.add_argument('--out', '-o', default='election_cf.csv',
                    help='Output summary CSV (default: election_cf.csv)')
    ap.add_argument('--per-rep', default=None,
                    help='Optional per-cluster-rep CSV path')
    ap.add_argument('--beta', type=float, default=0.02,
                    help='Soft-beta temperature for the reweighted election '
                         '(default: 0.02)')
    ap.add_argument('--n-clusters', type=int, default=10,
                    help='Max cluster reps to read per target (default: 10)')
    args = ap.parse_args(argv)

    root = os.path.expanduser(args.results_dir)
    if not os.path.isdir(root):
        print(f"ERROR: {root} is not a directory", file=sys.stderr)
        return 1

    target_dirs = sorted(
        os.path.join(root, d) for d in os.listdir(root)
        if os.path.isdir(os.path.join(root, d))
    )

    summaries, per_rep = [], []
    for td in target_dirs:
        rows = collect_target(td, n_clusters=args.n_clusters)
        if not rows:
            continue
        per_rep.extend(rows)
        summaries.append(summarize_target(rows, args.beta))

    if not summaries:
        print(f"No targets with cluster PDBs found in {root}")
        return 0

    fields = ['pdb_id', 'n_clusters', 'CF_rep_rank0', 'cluster_freq_rank0',
              'cluster_avgcf_rank0', 'CF_rep_argmin_rank', 'reweighted_rank0',
              'cf_election_changed', 'CF_rep_list', 'cluster_freq_list',
              'cluster_avgcf_list']
    with open(args.out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for s in summaries:
            w.writerow({k: (v if not isinstance(v, list) else str(v))
                        for k, v in s.items()})
    print(f"Election-CF summary written -> {args.out}  ({len(summaries)} targets)")

    if args.per_rep:
        rep_fields = ['pdb_id', 'rank', 'cf', 'cluster_avg_cf',
                      'cluster_frequency', 'cluster_top_rank']
        with open(args.per_rep, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=rep_fields)
            w.writeheader()
            w.writerows(per_rep)
        print(f"Per-rep detail written -> {args.per_rep}")

    flips = [s for s in summaries if s['cf_election_changed']]
    print(f"\nfreq x soft-beta (beta={args.beta}) would flip the elected head "
          f"on {len(flips)}/{len(summaries)} targets:")
    for s in flips:
        cf0 = s['CF_rep_rank0']
        cf0s = f'{cf0:8.2f}' if cf0 is not None else '   n/a  '
        print(f"  {s['pdb_id']:6s}  elected rank0 CF={cf0s} "
              f"freq={s['cluster_freq_rank0']}  ->  would elect rank "
              f"{s['reweighted_rank0']}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
