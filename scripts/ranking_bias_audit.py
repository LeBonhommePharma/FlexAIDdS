#!/usr/bin/env python3
"""Reproducible CF ranking-bias audit for real FlexAIDdS multi-pose docks.

For each complex directory (pose PDBs + optional INI + optional result.csv):

1. Rank poses by REMARK CF (lower = better).
2. Report RMSD under **explicit references** (never mix silently):
   - ``ini``: Kabsch RMSD of ligand heavy atoms to ``*_INI.pdb``
   - ``crystal``: Kabsch RMSD to ``--crystal-ligand`` SDF (or
     ``benchmarks/astex_diverse/astex_diverse/<PDB>/<PDB>_ligand.sdf`` if present)
   - ``result_csv``: ``rmsd_to_crystal`` / ``rmsd_hungarian`` for the **elected**
     pose only (not per-cluster), from result.csv when present
3. SMFREE / total-score: scan REMARKs for TOTAL_SCORE, SMFREE, ENTROPY,
   G_bind, free energy. If absent → explicit N/A (no invented fitness).

Stdlib only (no numpy). Exit 0 for successful audits.

Usage:
  python3 scripts/ranking_bias_audit.py results/astex_jcim2015_fair_20260708_0002
  python3 scripts/ranking_bias_audit.py results/smoke_classic_vs_current/1GPK_current \\
      --json /tmp/out.json
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from collections import defaultdict, deque
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


# ── Geometry (pure Python Kabsch via quaternion) ─────────────────────────────

def _kabsch_rmsd(P: List[Tuple[float, float, float]],
                 Q: List[Tuple[float, float, float]]) -> float:
    n = len(P)
    if n < 3:
        raise ValueError("need >=3 atoms")
    pc = [sum(p[i] for p in P) / n for i in range(3)]
    qc = [sum(q[i] for q in Q) / n for i in range(3)]
    Pc = [[p[i] - pc[i] for i in range(3)] for p in P]
    Qc = [[q[i] - qc[i] for i in range(3)] for q in Q]
    S = [[0.0] * 3 for _ in range(3)]
    for a, b in zip(Pc, Qc):
        for i in range(3):
            for j in range(3):
                S[i][j] += a[i] * b[j]
    S00, S01, S02 = S[0]
    S10, S11, S12 = S[1]
    S20, S21, S22 = S[2]
    K = [
        [S00 + S11 + S22, S12 - S21, S20 - S02, S01 - S10],
        [S12 - S21, S00 - S11 - S22, S01 + S10, S02 + S20],
        [S20 - S02, S01 + S10, -S00 + S11 - S22, S12 + S21],
        [S01 - S10, S02 + S20, S12 + S21, -S00 - S11 + S22],
    ]
    v = [1.0, 0.0, 0.0, 0.0]
    for _ in range(80):
        nv = [sum(K[i][j] * v[j] for j in range(4)) for i in range(4)]
        norm = math.sqrt(sum(x * x for x in nv)) or 1.0
        v = [x / norm for x in nv]
    q0, q1, q2, q3 = v
    R = [
        [q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3, 2 * (q1 * q2 - q0 * q3), 2 * (q1 * q3 + q0 * q2)],
        [2 * (q1 * q2 + q0 * q3), q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3, 2 * (q2 * q3 - q0 * q1)],
        [2 * (q1 * q3 - q0 * q2), 2 * (q2 * q3 + q0 * q1), q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3],
    ]
    s = 0.0
    for a, b in zip(Pc, Qc):
        ra = [sum(R[i][j] * a[j] for j in range(3)) for i in range(3)]
        s += sum((ra[i] - b[i]) ** 2 for i in range(3))
    return math.sqrt(s / n)


# ── PDB / SDF parsers ────────────────────────────────────────────────────────

_CF_RE = re.compile(
    r"REMARK\s+CF\s*[=:]\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)",
    re.I,
)
_TERM_RE = {
    "CF.com": re.compile(r"REMARK\s+CF\.com\s*=\s*([+-]?\d+\.?\d*)", re.I),
    "CF.wal": re.compile(r"REMARK\s+CF\.wal\s*=\s*([+-]?\d+\.?\d*)", re.I),
    "CF.sas": re.compile(r"REMARK\s+CF\.sas\s*=\s*([+-]?\d+\.?\d*)", re.I),
}
# Per-pose SMFREE / total fitness candidates (explicit scan)
_SMFREE_KEYS = re.compile(
    r"REMARK\s+(TOTAL_SCORE|SMFREE|ENTROPY|G_bind|G_BIND|FREE_ENERGY|F_bound|TdS|TOTAL)\s*[=:]\s*"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)",
    re.I,
)

_PROT = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE", "HID", "HIE", "HIP", "HOH", "WAT", "SO4", "PO4", "NA", "CL", "MG",
}


@dataclass
class LigandAtom:
    name: str
    element: str
    xyz: Tuple[float, float, float]


def _element_from_pdb(line: str, name: str) -> str:
    if len(line) >= 78 and line[76:78].strip():
        return line[76:78].strip().upper()
    # strip digits from atom name
    letters = "".join(c for c in name if c.isalpha())
    return (letters[:1] or "C").upper()


def parse_ligand_atoms_pdb(path: Path) -> List[LigandAtom]:
    atoms: List[LigandAtom] = []
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        res = line[17:20].strip()
        if res in _PROT:
            continue
        name = line[12:16].strip()
        if name.startswith("H") or (len(name) > 0 and name[0] == "H"):
            continue
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        el = _element_from_pdb(line, name)
        if el == "H":
            continue
        atoms.append(LigandAtom(name=name, element=el, xyz=xyz))
    return atoms


def parse_sdf_heavy_atoms(path: Path) -> List[LigandAtom]:
    """Parse first molecule V2000 SDF heavy atoms (order preserved)."""
    lines = path.read_text(errors="replace").splitlines()
    atoms: List[LigandAtom] = []
    i = 0
    while i < len(lines):
        if "V2000" in lines[i] or "V3000" in lines[i]:
            # counts line is this line for V2000
            parts = lines[i].split()
            try:
                nat = int(parts[0])
            except (ValueError, IndexError):
                i += 1
                continue
            for j in range(1, nat + 1):
                if i + j >= len(lines):
                    break
                toks = lines[i + j].split()
                if len(toks) < 4:
                    continue
                try:
                    x, y, z = float(toks[0]), float(toks[1]), float(toks[2])
                except ValueError:
                    continue
                el = toks[3].upper()
                if el == "H":
                    continue
                atoms.append(LigandAtom(name=f"{el}{j}", element=el, xyz=(x, y, z)))
            break
        i += 1
    return atoms


def rmsd_match_by_name(a: Sequence[LigandAtom], b: Sequence[LigandAtom]) -> Optional[float]:
    da = {x.name: x.xyz for x in a}
    db = {x.name: x.xyz for x in b}
    keys = sorted(set(da) & set(db))
    if len(keys) < 3:
        return None
    return _kabsch_rmsd([da[k] for k in keys], [db[k] for k in keys])


def rmsd_match_by_element_order(a: Sequence[LigandAtom], b: Sequence[LigandAtom]) -> Optional[float]:
    """Match by sorted element then original order (for SDF vs PDB)."""
    def key_el(xs: Sequence[LigandAtom]):
        # group by element preserving order within element
        return sorted(xs, key=lambda x: (x.element, x.name))

    aa, bb = key_el(a), key_el(b)
    n = min(len(aa), len(bb))
    if n < 3:
        return None
    # only pair if elements match at positions
    P, Q = [], []
    ia = ib = 0
    # greedy by element multiset order
    qa: Dict[str, deque] = defaultdict(deque)
    qb: Dict[str, deque] = defaultdict(deque)
    for x in aa:
        qa[x.element].append(x.xyz)
    for x in bb:
        qb[x.element].append(x.xyz)
    for el in sorted(set(qa) & set(qb)):
        while qa[el] and qb[el]:
            P.append(qa[el].popleft())
            Q.append(qb[el].popleft())
    if len(P) < 3:
        return None
    return _kabsch_rmsd(P, Q)


def parse_pose_remarks(path: Path) -> Dict[str, Any]:
    text = path.read_text(errors="replace")
    out: Dict[str, Any] = {"CF": None, "CF.com": None, "CF.wal": None, "smfree_fields": {}}
    m = _CF_RE.search(text)
    if m:
        out["CF"] = float(m.group(1))
    for k, rx in _TERM_RE.items():
        m = rx.search(text)
        if m:
            out[k] = float(m.group(1))
    for m in _SMFREE_KEYS.finditer(text):
        out["smfree_fields"][m.group(1).upper()] = float(m.group(2))
    return out


def spearman(xs: Sequence[float], ys: Sequence[float]) -> Optional[float]:
    n = len(xs)
    if n < 3:
        return None

    def ranks(vals: Sequence[float]) -> List[float]:
        order = sorted(range(n), key=lambda i: vals[i])
        r = [0.0] * n
        for rank, i in enumerate(order):
            r[i] = float(rank + 1)
        return r

    rx, ry = ranks(xs), ranks(ys)
    mx = sum(rx) / n
    my = sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    denx = math.sqrt(sum((rx[i] - mx) ** 2 for i in range(n)))
    deny = math.sqrt(sum((ry[i] - my) ** 2 for i in range(n)))
    if denx == 0.0 or deny == 0.0:
        return None
    return num / (denx * deny)


@dataclass
class PoseRow:
    file: str
    cf: float
    cf_com: Optional[float]
    cf_wal: Optional[float]
    rmsd_ini: Optional[float]
    rmsd_crystal: Optional[float]
    smfree_fields: Dict[str, float] = field(default_factory=dict)


def discover_complex_dirs(root: Path) -> List[Path]:
    root = root.resolve()
    if not root.exists():
        return []
    # If root itself has poses
    if list(root.glob("*_[0-9]*.pdb")) or list(root.glob("*_0.pdb")):
        return [root]
    return sorted(
        p for p in root.iterdir()
        if p.is_dir() and (
            list(p.glob("*_[0-9]*.pdb")) or (p / "result.csv").is_file()
        )
    )


def infer_pdb_id(d: Path) -> str:
    for p in sorted(d.glob("*_0.pdb")):
        return p.name[: -len("_0.pdb")]
    for p in sorted(d.glob("*_INI.pdb")):
        return p.name[: -len("_INI.pdb")]
    return d.name


def resolve_crystal_sdf(pdb_id: str, repo: Path, explicit: Optional[Path]) -> Optional[Path]:
    if explicit and explicit.is_file():
        return explicit
    cand = repo / "benchmarks" / "astex_diverse" / "astex_diverse" / pdb_id / f"{pdb_id}_ligand.sdf"
    if cand.is_file():
        return cand
    # case variants
    for p in (repo / "benchmarks" / "astex_diverse").rglob(f"{pdb_id}_ligand.sdf"):
        return p
    return None


def audit_complex(
    d: Path,
    repo: Path,
    crystal_sdf: Optional[Path] = None,
) -> Dict[str, Any]:
    pdb_id = infer_pdb_id(d)
    ini_path = d / f"{pdb_id}_INI.pdb"
    if not ini_path.is_file():
        inis = list(d.glob("*_INI.pdb"))
        ini_path = inis[0] if inis else None

    pose_paths = sorted(d.glob(f"{pdb_id}_[0-9]*.pdb"))
    if not pose_paths:
        pose_paths = sorted(
            p for p in d.glob("*_[0-9]*.pdb") if "INI" not in p.name
        )

    ini_atoms = parse_ligand_atoms_pdb(ini_path) if ini_path and ini_path.is_file() else []
    cryst_path = resolve_crystal_sdf(pdb_id, repo, crystal_sdf)
    cryst_atoms = parse_sdf_heavy_atoms(cryst_path) if cryst_path else []

    rows: List[PoseRow] = []
    for pp in pose_paths:
        rem = parse_pose_remarks(pp)
        if rem["CF"] is None:
            continue
        pose_atoms = parse_ligand_atoms_pdb(pp)
        r_ini = rmsd_match_by_name(ini_atoms, pose_atoms) if ini_atoms else None
        if r_ini is None and ini_atoms:
            r_ini = rmsd_match_by_element_order(ini_atoms, pose_atoms)
        r_xtal = None
        if cryst_atoms:
            r_xtal = rmsd_match_by_element_order(cryst_atoms, pose_atoms)
        rows.append(
            PoseRow(
                file=pp.name,
                cf=float(rem["CF"]),
                cf_com=rem.get("CF.com"),
                cf_wal=rem.get("CF.wal"),
                rmsd_ini=r_ini,
                rmsd_crystal=r_xtal,
                smfree_fields=dict(rem["smfree_fields"]),
            )
        )

    if not rows:
        return {
            "pdb_id": pdb_id,
            "path": str(d),
            "error": "no_cf_poses",
            "rmsd_references": {
                "ini": str(ini_path) if ini_path else None,
                "crystal_sdf": str(cryst_path) if cryst_path else None,
                "result_csv": None,
            },
            "smfree_per_pose": "N/A — no TOTAL_SCORE/SMFREE/ENTROPY/G_bind REMARKs on pose PDBs",
        }

    by_cf = sorted(rows, key=lambda r: r.cf)
    top1 = by_cf[0]

    def rank_best(rmsd_attr: str) -> Optional[int]:
        with_r = [r for r in rows if getattr(r, rmsd_attr) is not None]
        if not with_r:
            return None
        best = min(with_r, key=lambda r: getattr(r, rmsd_attr))
        for i, r in enumerate(by_cf):
            if r.file == best.file:
                return i + 1
        return None

    def rho(rmsd_attr: str) -> Optional[float]:
        pairs = [(r.cf, getattr(r, rmsd_attr)) for r in rows if getattr(r, rmsd_attr) is not None]
        if len(pairs) < 3:
            return None
        return spearman([a for a, _ in pairs], [b for _, b in pairs])

    # result.csv elected-only crystal RMSD
    result_csv_rmsd = None
    result_csv_path = d / "result.csv"
    cf_native_csv = None
    best_score_csv = None
    if result_csv_path.is_file():
        r0 = list(csv.DictReader(result_csv_path.open()))[0]
        for k in ("rmsd_to_crystal", "rmsd_hungarian"):
            if r0.get(k) not in (None, "", "NA"):
                try:
                    result_csv_rmsd = float(r0[k])
                    break
                except ValueError:
                    pass
        try:
            if r0.get("cf_native") not in (None, "", "NA"):
                cf_native_csv = float(r0["cf_native"])
            if r0.get("best_score") not in (None, "", "NA"):
                best_score_csv = float(r0["best_score"])
        except ValueError:
            pass

    ini_cf = None
    if ini_path and ini_path.is_file():
        ini_cf = parse_pose_remarks(ini_path).get("CF")

    best_ini = min((r for r in rows if r.rmsd_ini is not None), key=lambda r: r.rmsd_ini, default=None)
    best_xtal = min((r for r in rows if r.rmsd_crystal is not None), key=lambda r: r.rmsd_crystal, default=None)

    oracle_ini = any(r.rmsd_ini is not None and r.rmsd_ini <= 2.0 for r in rows)
    oracle_xtal = any(r.rmsd_crystal is not None and r.rmsd_crystal <= 2.0 for r in rows)
    top1_ini_hit = top1.rmsd_ini is not None and top1.rmsd_ini <= 2.0
    top1_xtal_hit = top1.rmsd_crystal is not None and top1.rmsd_crystal <= 2.0

    smfree_keys = sorted({k for r in rows for k in r.smfree_fields})
    smfree_note = (
        "present on some poses: " + ",".join(smfree_keys)
        if smfree_keys
        else "N/A — no TOTAL_SCORE/SMFREE/ENTROPY/G_bind REMARKs on pose PDBs (CF-only ranking)"
    )

    return {
        "pdb_id": pdb_id,
        "path": str(d),
        "n_poses": len(rows),
        "rmsd_references": {
            "ini": {
                "path": str(ini_path) if ini_path else None,
                "method": "Kabsch ligand heavy atoms; match atom name then element multiset",
                "note": "INI is the engine initial/seed structure — may NOT equal crystal",
            },
            "crystal": {
                "path": str(cryst_path) if cryst_path else None,
                "method": "Kabsch heavy atoms vs SDF; element-order multiset match",
            },
            "result_csv_elected_only": {
                "path": str(result_csv_path) if result_csv_path.is_file() else None,
                "rmsd_to_crystal": result_csv_rmsd,
                "note": "Applies only to elected/top pose, not full cluster ranking",
            },
        },
        "smfree_per_pose": smfree_note,
        "ini_CF": ini_cf,
        "top1_file": top1.file,
        "top1_CF": top1.cf,
        "top1_rmsd_ini": top1.rmsd_ini,
        "top1_rmsd_crystal": top1.rmsd_crystal,
        "top1_CF.com": top1.cf_com,
        "top1_CF.wal": top1.cf_wal,
        "gap_top1_minus_ini": (top1.cf - ini_cf) if ini_cf is not None else None,
        "gap_top1_minus_cf_native_csv": (
            (best_score_csv - cf_native_csv)
            if best_score_csv is not None and cf_native_csv is not None
            else None
        ),
        "best_rmsd_ini": best_ini.rmsd_ini if best_ini else None,
        "best_rmsd_ini_file": best_ini.file if best_ini else None,
        "rank_of_best_rmsd_ini_under_CF": rank_best("rmsd_ini"),
        "best_rmsd_crystal": best_xtal.rmsd_crystal if best_xtal else None,
        "best_rmsd_crystal_file": best_xtal.file if best_xtal else None,
        "rank_of_best_rmsd_crystal_under_CF": rank_best("rmsd_crystal"),
        "oracle_ini_le_2": oracle_ini,
        "oracle_crystal_le_2": oracle_xtal,
        "nn_miss_ini": oracle_ini and not top1_ini_hit,
        "nn_miss_crystal": oracle_xtal and not top1_xtal_hit,
        "spearman_CF_rmsd_ini": rho("rmsd_ini"),
        "spearman_CF_rmsd_crystal": rho("rmsd_crystal"),
        "result_csv_elected_rmsd_crystal": result_csv_rmsd,
        "poses": [asdict(r) for r in by_cf],
    }


def format_table(results: List[Dict[str, Any]]) -> str:
    lines = [
        "=== CF ranking bias audit (dual RMSD references) ===",
        "RMSD_ini  = Kabsch vs *_INI.pdb (seed; may ≠ crystal)",
        "RMSD_xtal = Kabsch vs crystal ligand SDF when available",
        "result_csv rmsd_to_crystal = elected pose only",
        "SMFREE/total: see smfree_per_pose field (N/A if REMARKs absent)",
        "",
        f"{'pdb':8} {'n':>3} {'top1_CF':>10} {'rmsd_ini':>9} {'rmsd_xtal':>9} "
        f"{'csv_xtal':>9} {'rk_ini':>6} {'rk_xt':>5} {'ρ_ini':>7} {'ρ_xt':>7} {'gap_INI':>9}",
    ]
    for r in results:
        if r.get("error"):
            lines.append(f"{r.get('pdb_id','?'):8} ERROR {r['error']}")
            continue
        def f(v, nd=3):
            if v is None:
                return "n/a"
            if isinstance(v, float):
                return f"{v:.{nd}f}"
            return str(v)
        lines.append(
            f"{r['pdb_id']:8} {r['n_poses']:3d} {r['top1_CF']:10.2f} "
            f"{f(r.get('top1_rmsd_ini')):>9} {f(r.get('top1_rmsd_crystal')):>9} "
            f"{f(r.get('result_csv_elected_rmsd_crystal')):>9} "
            f"{f(r.get('rank_of_best_rmsd_ini_under_CF'),0):>6} "
            f"{f(r.get('rank_of_best_rmsd_crystal_under_CF'),0):>5} "
            f"{f(r.get('spearman_CF_rmsd_ini')):>7} {f(r.get('spearman_CF_rmsd_crystal')):>7} "
            f"{f(r.get('gap_top1_minus_ini'),2):>9}"
        )
        lines.append(f"         smfree: {r.get('smfree_per_pose')}")
    return "\n".join(lines) + "\n"


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("results_dir", type=Path, help="Complex dir or parent of complexes")
    ap.add_argument("--json", type=Path, default=None, help="Write JSON summary")
    ap.add_argument("--crystal-ligand", type=Path, default=None, help="Override crystal SDF")
    args = ap.parse_args(argv)

    root = repo_root()
    target = args.results_dir
    if not target.is_absolute():
        # try cwd then repo
        if not target.exists():
            alt = root / target
            if alt.exists():
                target = alt
    dirs = discover_complex_dirs(target)
    if not dirs:
        print(f"No complex dirs under {target}", file=sys.stderr)
        return 0

    results = [audit_complex(d, root, args.crystal_ligand) for d in dirs]
    print(format_table(results))
    if args.json:
        args.json.write_text(json.dumps(results, indent=2) + "\n")
        print(f"Wrote {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
