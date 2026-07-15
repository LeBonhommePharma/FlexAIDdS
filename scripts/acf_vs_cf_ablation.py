#!/usr/bin/env python3
"""ACF-vs-CF rank-0 ablation (classic FlexAID entropy ranking).

Mirrors the C++ emission gate in LIB/cluster.cpp:

  classic (default, T>0, !force_cf_rank_emission):
      rank-0 = lowest ACF  (soft-β cluster free energy)
  force_cf / T==0 (P3b rollback):
      rank-0 = lowest representative CF

Offline-only: no docking, no re-score. Reads an existing run's ``*.cad``
(and optionally per-rank ``*_N.pdb`` REMARK CF) and reports who would win
under each policy.

Canonical exhibit is Astex 1HNN (pre-classic-entropy live run elected the
CF champion as rank-0 while the densest ACF basin sat at rank 3).

Usage:
  # Built-in 1HNN synthetic exhibit (no files needed)
  python3 scripts/acf_vs_cf_ablation.py --synthetic-1hnn

  # Live target directory with 1HNN.cad
  python3 scripts/acf_vs_cf_ablation.py results/.../1HNN

  # Explicit .cad
  python3 scripts/acf_vs_cf_ablation.py --cad path/to/1HNN.cad --json

Exit code 0 always on successful parse; prints a clear flip/no-flip verdict.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

# Cluster N: TOP=... TCF=... ACF=... freq=...
_CAD_RE = re.compile(
    r"Cluster\s+(?P<idx>\d+)\s*:\s*"
    r"TOP=(?P<top>-?\d+)\s+"
    r"TCF=(?P<tcf>[-\d.eE+]+)\s+"
    r"ACF=(?P<acf>[-\d.eE+]+)\s+"
    r"freq=(?P<freq>\d+)"
)

# REMARK CF=-189.85613
_CF_RE = re.compile(r"^REMARK\s+CF\s*=\s*(?P<cf>[-\d.eE+]+)")

# REMARK Cluster 0: Rank (top):0 Average CF:-49.30565 Frequency:4
# NOTE: engine writes Clus_ACF into the "Average CF" field (historical label).
_CLUSTER_REMARK_RE = re.compile(
    r"REMARK\s+Cluster\s+(?P<rank>\d+)\s*:\s*"
    r"Rank\s*\(top\):\s*(?P<top>-?\d+)\s+"
    r"Average CF:\s*(?P<acf>[-\d.eE+]+)\s+"
    r"Frequency:\s*(?P<freq>\d+)",
    re.IGNORECASE,
)

# Live pre-fix 1HNN.cad numbers (results/.../20260708.../1HNN).
# CF values come from rank PDBs REMARK CF= (representative evalue).
SYNTHETIC_1HNN: List[dict] = [
    {"cluster": 0, "top": 0, "acf": -49.305648, "cf": -189.85613, "freq": 4},
    {"cluster": 1, "top": 2, "acf": -83.392147, "cf": -81.80491, "freq": 8},
    {"cluster": 2, "top": 6, "acf": -48.878832, "cf": -72.93359, "freq": 5},
    {"cluster": 3, "top": 7, "acf": -263.427453, "cf": -72.05566, "freq": 29},
    {"cluster": 4, "top": 9, "acf": -221.231029, "cf": -71.55236, "freq": 24},
    {"cluster": 5, "top": 14, "acf": -51.210015, "cf": -65.03617, "freq": 6},
    {"cluster": 6, "top": 18, "acf": -38.943764, "cf": -54.69740, "freq": 4},
    {"cluster": 7, "top": 22, "acf": -105.236160, "cf": -51.35090, "freq": 14},
    {"cluster": 8, "top": 26, "acf": -90.972385, "cf": -46.08170, "freq": 10},
    {"cluster": 9, "top": 28, "acf": -10.699287, "cf": -43.31451, "freq": 1},
]


@dataclass
class ClusterRow:
    cluster: int
    top: int
    acf: float
    cf: float
    freq: int


def elect_rank0(
    rows: Sequence[ClusterRow],
    *,
    force_cf_rank_emission: bool,
    temperature: int = 300,
) -> int:
    """Return emission index (position in *rows*) for rank-0 under policy.

    Mirrors tests/test_classic_entropy_ranking.cpp::elect_rank0 and cluster.cpp.
    """
    n = len(rows)
    if n == 0:
        return -1
    order = list(range(n))
    if temperature > 0:
        order.sort(key=lambda i: rows[i].acf)
    classic = (temperature > 0) and (not force_cf_rank_emission)
    if not classic:
        order.sort(key=lambda i: rows[i].cf)
    return order[0]


def parse_cad(text: str) -> List[ClusterRow]:
    rows: List[ClusterRow] = []
    for line in text.splitlines():
        m = _CAD_RE.search(line)
        if not m:
            continue
        # TCF is soft-β top-member free energy contribution, not rep CF.
        # Use TCF as CF fallback only when no PDB REMARK is available.
        tcf = float(m.group("tcf"))
        rows.append(
            ClusterRow(
                cluster=int(m.group("idx")),
                top=int(m.group("top")),
                acf=float(m.group("acf")),
                cf=tcf,  # provisional; may be overwritten from PDB REMARK CF
                freq=int(m.group("freq")),
            )
        )
    return rows


def parse_rank_pdb_cf(pdb_path: Path) -> Tuple[Optional[float], Optional[float], Optional[int]]:
    """Return (rep_CF, ACF_from_remark, frequency) from a cluster PDB."""
    cf = acf = freq = None
    try:
        with pdb_path.open() as fh:
            for line in fh:
                if cf is None:
                    m = _CF_RE.match(line.strip())
                    if m:
                        cf = float(m.group("cf"))
                        continue
                m = _CLUSTER_REMARK_RE.search(line)
                if m:
                    acf = float(m.group("acf"))
                    freq = int(m.group("freq"))
                    break
    except OSError:
        pass
    return cf, acf, freq


def load_target_dir(target_dir: Path, pdb_id: Optional[str] = None) -> List[ClusterRow]:
    target_dir = target_dir.resolve()
    code = pdb_id or target_dir.name
    cad = target_dir / f"{code}.cad"
    if not cad.is_file():
        cads = sorted(target_dir.glob("*.cad"))
        if not cads:
            raise FileNotFoundError(f"no .cad under {target_dir}")
        cad = cads[0]
        code = cad.stem

    rows = parse_cad(cad.read_text(errors="replace"))
    if not rows:
        raise ValueError(f"no Cluster lines in {cad}")

    # Overlay representative CF from rank PDBs when present (true evalue).
    for r in rows:
        pdb = target_dir / f"{code}_{r.cluster}.pdb"
        if not pdb.is_file():
            continue
        cf, acf_rem, freq_rem = parse_rank_pdb_cf(pdb)
        if cf is not None:
            r.cf = cf
        if acf_rem is not None and abs(acf_rem - r.acf) > 1e-3:
            # Prefer .cad ACF (authoritative); remark is the same field historically.
            pass
        if freq_rem is not None and freq_rem != r.freq:
            pass
    return rows


def rows_from_synthetic(records: Sequence[dict]) -> List[ClusterRow]:
    return [
        ClusterRow(
            cluster=int(d["cluster"]),
            top=int(d.get("top", d["cluster"])),
            acf=float(d["acf"]),
            cf=float(d["cf"]),
            freq=int(d["freq"]),
        )
        for d in records
    ]


def ablation_report(
    rows: Sequence[ClusterRow],
    *,
    label: str,
    temperature: int = 300,
) -> dict:
    classic_i = elect_rank0(rows, force_cf_rank_emission=False, temperature=temperature)
    force_i = elect_rank0(rows, force_cf_rank_emission=True, temperature=temperature)
    t0_i = elect_rank0(rows, force_cf_rank_emission=False, temperature=0)

    def pack(i: int) -> Optional[dict]:
        if i < 0 or i >= len(rows):
            return None
        r = rows[i]
        return {
            "emission_index": i,
            "cluster": r.cluster,
            "top": r.top,
            "acf": r.acf,
            "cf": r.cf,
            "freq": r.freq,
        }

    classic = pack(classic_i)
    force = pack(force_i)
    flipped = classic is not None and force is not None and classic["cluster"] != force["cluster"]

    return {
        "label": label,
        "n_clusters": len(rows),
        "temperature": temperature,
        "classic_entropy_rank0": classic,
        "force_cf_rank0": force,
        "temperature_zero_rank0": pack(t0_i),
        "election_flips": bool(flipped),
        "clusters": [asdict(r) for r in rows],
        "contract": {
            "classic": "rank-0 = lowest ACF (soft-β H−T·S over cluster; CF re-sort OFF)",
            "force_cf": "rank-0 = lowest representative CF (P3b / cd9004d; CF re-sort ON)",
            "rollback": "thermodynamics.force_cf_rank_emission=true or "
            "FLEXAIDDS_FORCE_CF_RANK_EMISSION=1 or classic_entropy_ranking=false",
        },
    }


def format_text(report: dict) -> str:
    lines = [
        f"=== ACF-vs-CF ablation: {report['label']} ===",
        f"clusters={report['n_clusters']}  T={report['temperature']}",
        "",
        f"{'idx':>4}  {'cluster':>7}  {'ACF':>12}  {'CF':>12}  {'freq':>5}",
    ]
    for r in report["clusters"]:
        lines.append(
            f"{r['cluster']:4d}  {r['cluster']:7d}  {r['acf']:12.4f}  "
            f"{r['cf']:12.4f}  {r['freq']:5d}"
        )
    c = report["classic_entropy_rank0"]
    f = report["force_cf_rank0"]
    lines += [
        "",
        f"classic_entropy rank-0: cluster={c['cluster']}  ACF={c['acf']:.4f}  "
        f"CF={c['cf']:.4f}  freq={c['freq']}",
        f"force_cf        rank-0: cluster={f['cluster']}  ACF={f['acf']:.4f}  "
        f"CF={f['cf']:.4f}  freq={f['freq']}",
        "",
        "VERDICT: "
        + (
            f"ELECTION FLIP — classic elects dense ACF basin (cluster {c['cluster']}, "
            f"freq {c['freq']}); CF re-sort elects CF champion (cluster {f['cluster']}, "
            f"CF={f['cf']:.2f})."
            if report["election_flips"]
            else "no flip — ACF and CF agree on rank-0 for this ensemble."
        ),
        "",
        "Rollback: force_cf_rank_emission=true | FLEXAIDDS_FORCE_CF_RANK_EMISSION=1",
        "Docs: docs/classic_entropy_ranking.md",
    ]
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("target", nargs="?", help="Target result directory containing <PDB>.cad")
    p.add_argument("--cad", type=Path, help="Explicit .cad path")
    p.add_argument("--pdb-id", help="PDB code (default: directory name / .cad stem)")
    p.add_argument("--synthetic-1hnn", action="store_true", help="Use built-in pre-fix 1HNN numbers")
    p.add_argument("--temperature", type=int, default=300)
    p.add_argument("--json", action="store_true", help="Emit JSON report")
    p.add_argument("--csv", type=Path, help="Write cluster table CSV")
    p.add_argument("--out", type=Path, help="Write full JSON report to path")
    args = p.parse_args(argv)

    if args.synthetic_1hnn:
        rows = rows_from_synthetic(SYNTHETIC_1HNN)
        label = "synthetic-1HNN (pre-classic live exhibit)"
    elif args.cad:
        rows = parse_cad(args.cad.read_text(errors="replace"))
        label = str(args.cad)
    elif args.target:
        tdir = Path(args.target)
        rows = load_target_dir(tdir, args.pdb_id)
        label = str(tdir)
    else:
        # Default: synthetic exhibit so `python3 scripts/acf_vs_cf_ablation.py` is useful.
        rows = rows_from_synthetic(SYNTHETIC_1HNN)
        label = "synthetic-1HNN (default; pass a target dir for live data)"

    if not rows:
        print("error: no clusters loaded", file=sys.stderr)
        return 2

    report = ablation_report(rows, label=label, temperature=args.temperature)

    if args.csv:
        with args.csv.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["cluster", "top", "acf", "cf", "freq"])
            w.writeheader()
            for r in report["clusters"]:
                w.writerow(r)

    if args.out:
        args.out.write_text(json.dumps(report, indent=2) + "\n")

    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print(format_text(report))
    return 0


if __name__ == "__main__":
    sys.exit(main())
