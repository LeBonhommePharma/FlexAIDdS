#!/usr/bin/env python3
"""Native CF audit: elected pose CF vs crystal/native CF and the gap.

For each complex under a results directory, report:

  pdb_id | elected_pose | CF_pose | CF_native | gap | rmsd | pathology

where:
  gap        = CF_pose - CF_native   (negative ⇒ decoy scores better; lower CF is better)
  pathology  = gap < -pathology_gap  (default 5.0) — decoy much better than native

Data sources (first available wins):
  CF_pose:   REMARK CF= on elected pose PDB, else result.csv best_score
  CF_native: REMARK CF= on *_INI.pdb, else result.csv cf_native
  rmsd:      result.csv rmsd_hungarian / rmsd_to_crystal, else REMARK RMSD on pose PDB
  elected:   result.csv elected_pose path, else {pdb}_0.pdb, else best-CF numbered pose

Usage:
  python3 scripts/audit_native_cf.py results/astex_jcim2015_fair_20260708_0002
  python3 scripts/audit_native_cf.py results/.../1GPK --json audit.json

Exit code is always 0 for successful audits (missing data is reported, not fatal).
Paths are resolved relative to the repo root or the caller's CWD — no hard-coded
machine paths.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


# ── Patterns ──────────────────────────────────────────────────────────────────

# FlexAID writes: "REMARK CF=-73.72798" or "REMARK CF= 1.02993" (space after =)
_CF_RE = re.compile(
    r"REMARK\s+CF\s*[=:]\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)",
    re.IGNORECASE,
)
# Alternate: "REMARK CF -73.7" (space delimiter, no =)
_CF_SPACE_RE = re.compile(
    r"REMARK\s+CF\s+([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*$",
    re.IGNORECASE,
)
_RMSD_RE = re.compile(
    r"REMARK\s+([+-]?(?:\d+(?:\.\d*)?|\.\d+))\s+RMSD\s+to\s+ref",
    re.IGNORECASE,
)
_POSE_NUM_RE = re.compile(r"^(.+)_(\d+)$")

DEFAULT_PATHOLOGY_GAP = 5.0  # CF units; lower CF is better


# ── Pure parsers (unit-tested) ────────────────────────────────────────────────

def parse_cf_from_remark_lines(lines: Iterable[str]) -> Optional[float]:
    """Extract the first REMARK CF value from PDB remark lines.

    Accepts common FlexAID forms:
      REMARK CF=-12.5
      REMARK CF= -12.5
      REMARK CF: -12.5
      REMARK CF -12.5
    """
    for line in lines:
        if "cf" not in line.lower():
            continue
        m = _CF_RE.search(line)
        if m:
            try:
                return float(m.group(1))
            except ValueError:
                continue
        m = _CF_SPACE_RE.search(line)
        if m:
            try:
                return float(m.group(1))
            except ValueError:
                continue
    return None


def parse_cf_from_pdb_text(text: str) -> Optional[float]:
    """Parse CF from full PDB file text."""
    return parse_cf_from_remark_lines(text.splitlines())


def parse_rmsd_from_remark_lines(lines: Iterable[str]) -> Optional[float]:
    """Extract first 'REMARK <x> RMSD to ref. structure' value if present."""
    for line in lines:
        if "RMSD" not in line.upper():
            continue
        m = _RMSD_RE.search(line)
        if m:
            try:
                return float(m.group(1))
            except ValueError:
                continue
    return None


def is_pathology(
    cf_pose: Optional[float],
    cf_native: Optional[float],
    gap_threshold: float = DEFAULT_PATHOLOGY_GAP,
) -> bool:
    """True when decoy CF is better than native by more than *gap_threshold*.

    Lower CF is better, so pathology when CF_pose < CF_native - gap_threshold
    (equivalently gap = CF_pose - CF_native < -gap_threshold).
    """
    if cf_pose is None or cf_native is None:
        return False
    return (cf_pose - cf_native) < -float(gap_threshold)


# ── Filesystem helpers ────────────────────────────────────────────────────────

def repo_root() -> Path:
    """Resolve repository root from this script's location."""
    return Path(__file__).resolve().parents[1]


def resolve_path(path: str | Path) -> Path:
    """Resolve a user path; prefer absolute, else CWD, else repo-relative."""
    p = Path(path).expanduser()
    if p.is_absolute():
        return p.resolve()
    cwd_candidate = (Path.cwd() / p).resolve()
    if cwd_candidate.exists():
        return cwd_candidate
    repo_candidate = (repo_root() / p).resolve()
    if repo_candidate.exists():
        return repo_candidate
    return cwd_candidate


def read_cf_from_pdb(path: Path) -> Optional[float]:
    try:
        text = path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return None
    return parse_cf_from_pdb_text(text)


def read_rmsd_from_pdb(path: Path) -> Optional[float]:
    try:
        text = path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return None
    return parse_rmsd_from_remark_lines(text.splitlines())


def _float_or_none(raw: Any) -> Optional[float]:
    if raw is None:
        return None
    s = str(raw).strip()
    if s in ("", "N/A", "NA", "nan", "None", "."):
        return None
    try:
        v = float(s)
    except (TypeError, ValueError):
        return None
    if v != v:  # NaN
        return None
    return v


def load_result_csv(path: Path) -> Optional[Dict[str, str]]:
    """Load the first data row of a per-complex result.csv as a dict."""
    if not path.is_file():
        return None
    try:
        with path.open(newline="", encoding="utf-8", errors="ignore") as fh:
            rows = list(csv.DictReader(fh))
    except OSError:
        return None
    if not rows:
        return None
    return {k: (v if v is not None else "") for k, v in rows[0].items()}


def discover_complex_dirs(root: Path) -> List[Path]:
    """Return directories that look like per-complex docking outputs.

    Accepts either:
      - a single complex dir (contains *_0.pdb / result.csv / dock_config.json)
      - a parent run dir with child complex dirs
    """
    root = root.resolve()
    if not root.is_dir():
        return []

    def looks_like_complex(d: Path) -> bool:
        if not d.is_dir():
            return False
        if (d / "result.csv").is_file() or (d / "dock_config.json").is_file():
            return True
        # Numbered pose PDBs or INI
        for p in d.glob("*_*.pdb"):
            if p.name.upper().endswith("_INI.PDB"):
                return True
            stem = p.stem
            if _POSE_NUM_RE.match(stem):
                return True
        return False

    if looks_like_complex(root):
        return [root]

    children = sorted(
        (c for c in root.iterdir() if c.is_dir() and not c.name.startswith(".")),
        key=lambda p: p.name.upper(),
    )
    complexes = [c for c in children if looks_like_complex(c)]
    return complexes


def guess_pdb_id(complex_dir: Path, csv_row: Optional[Dict[str, str]] = None) -> str:
    if csv_row:
        for key in ("pdb_id", "target", "complex", "id", "code"):
            v = (csv_row.get(key) or "").strip()
            if v:
                return v.upper()
    name = complex_dir.name
    # Strip common suffixes like 1GPK_cleft0
    m = re.match(r"^([0-9][A-Za-z0-9]{3})", name)
    if m:
        return m.group(1).upper()
    return name.upper()


def find_ini_pdb(complex_dir: Path, pdb_id: str) -> Optional[Path]:
    candidates = [
        complex_dir / f"{pdb_id}_INI.pdb",
        complex_dir / f"{pdb_id}_ini.pdb",
        complex_dir / f"{pdb_id}_native.pdb",
        complex_dir / "INI.pdb",
        complex_dir / "native.pdb",
    ]
    for p in candidates:
        if p.is_file():
            return p
    # Case-insensitive glob fallback
    for p in complex_dir.glob("*_INI.pdb"):
        return p
    for p in complex_dir.glob("*_ini.pdb"):
        return p
    return None


def list_pose_pdbs(complex_dir: Path, pdb_id: str) -> List[Tuple[int, Path]]:
    """Return (rank_index, path) for numbered pose PDBs, sorted by index."""
    poses: List[Tuple[int, Path]] = []
    for p in complex_dir.glob("*.pdb"):
        name = p.name
        upper = name.upper()
        if upper.endswith("_INI.PDB") or upper.endswith("_PRUNED.PDB"):
            continue
        if "NATIVE" in upper and not _POSE_NUM_RE.match(p.stem):
            continue
        m = _POSE_NUM_RE.match(p.stem)
        if not m:
            continue
        prefix, idx_s = m.group(1), m.group(2)
        # Prefer poses matching pdb_id prefix when available
        if pdb_id and not prefix.upper().startswith(pdb_id.upper()[:4]):
            # still accept if directory only has one family
            pass
        try:
            idx = int(idx_s)
        except ValueError:
            continue
        poses.append((idx, p))
    poses.sort(key=lambda t: t[0])
    # Prefer same-prefix poses if any match pdb_id
    same = [t for t in poses if t[1].stem.upper().startswith(pdb_id.upper())]
    return same if same else poses


def elect_pose(
    complex_dir: Path,
    pdb_id: str,
    csv_row: Optional[Dict[str, str]],
) -> Tuple[Optional[Path], str]:
    """Pick the elected (reported best) pose path and a short label.

    Returns (path_or_None, label) where label is a basename or csv path string.
    """
    if csv_row:
        elected = (csv_row.get("elected_pose") or "").strip()
        if elected:
            ep = Path(elected)
            # Relative to CWD / repo / complex dir
            candidates = [
                ep if ep.is_absolute() else None,
                (Path.cwd() / ep).resolve() if not ep.is_absolute() else None,
                (repo_root() / ep).resolve() if not ep.is_absolute() else None,
                complex_dir / ep.name,
                complex_dir / ep,
            ]
            for c in candidates:
                if c is not None and c.is_file():
                    return c, c.name
            # Keep label even if file missing
            return None, Path(elected).name

    poses = list_pose_pdbs(complex_dir, pdb_id)
    if not poses:
        return None, ""

    # Prefer rank 0
    for idx, path in poses:
        if idx == 0:
            return path, path.name

    # Else lowest CF among numbered poses
    best_path: Optional[Path] = None
    best_cf: Optional[float] = None
    for _, path in poses:
        cf = read_cf_from_pdb(path)
        if cf is None:
            continue
        if best_cf is None or cf < best_cf:
            best_cf = cf
            best_path = path
    if best_path is not None:
        return best_path, best_path.name
    return poses[0][1], poses[0][1].name


# ── Audit core ────────────────────────────────────────────────────────────────

@dataclass
class AuditRow:
    pdb_id: str
    elected_pose: str
    cf_pose: Optional[float]
    cf_native: Optional[float]
    gap: Optional[float]
    rmsd: Optional[float]
    pathology: bool
    source_pose: str = ""
    source_native: str = ""
    complex_dir: str = ""
    notes: str = ""

    def as_jsonable(self) -> Dict[str, Any]:
        d = asdict(self)
        # Keep None as null; format floats for stability
        for k in ("cf_pose", "cf_native", "gap", "rmsd"):
            v = d[k]
            if isinstance(v, float):
                d[k] = round(v, 6)
        return d


def audit_complex(
    complex_dir: Path,
    pathology_gap: float = DEFAULT_PATHOLOGY_GAP,
) -> AuditRow:
    complex_dir = complex_dir.resolve()
    csv_row = load_result_csv(complex_dir / "result.csv")
    pdb_id = guess_pdb_id(complex_dir, csv_row)

    pose_path, pose_label = elect_pose(complex_dir, pdb_id, csv_row)
    ini_path = find_ini_pdb(complex_dir, pdb_id)

    cf_pose: Optional[float] = None
    source_pose = ""
    if pose_path is not None:
        cf_pose = read_cf_from_pdb(pose_path)
        if cf_pose is not None:
            source_pose = f"REMARK:{pose_path.name}"
    if cf_pose is None and csv_row is not None:
        cf_pose = _float_or_none(csv_row.get("best_score"))
        if cf_pose is not None:
            source_pose = "result.csv:best_score"

    cf_native: Optional[float] = None
    source_native = ""
    if ini_path is not None:
        cf_native = read_cf_from_pdb(ini_path)
        if cf_native is not None:
            source_native = f"REMARK:{ini_path.name}"
    if cf_native is None and csv_row is not None:
        cf_native = _float_or_none(csv_row.get("cf_native"))
        if cf_native is not None:
            source_native = "result.csv:cf_native"

    rmsd: Optional[float] = None
    if csv_row is not None:
        for key in ("rmsd_hungarian", "rmsd_to_crystal", "best_cluster_rmsd", "rmsd"):
            rmsd = _float_or_none(csv_row.get(key))
            if rmsd is not None:
                break
    if rmsd is None and pose_path is not None:
        rmsd = read_rmsd_from_pdb(pose_path)

    gap: Optional[float] = None
    if cf_pose is not None and cf_native is not None:
        gap = cf_pose - cf_native

    notes_parts: List[str] = []
    if pose_path is None and not pose_label:
        notes_parts.append("no_elected_pose")
    if cf_pose is None:
        notes_parts.append("missing_cf_pose")
    if cf_native is None:
        notes_parts.append("missing_cf_native")
    if csv_row is None:
        notes_parts.append("no_result_csv")

    return AuditRow(
        pdb_id=pdb_id,
        elected_pose=pose_label or (pose_path.name if pose_path else ""),
        cf_pose=cf_pose,
        cf_native=cf_native,
        gap=gap,
        rmsd=rmsd,
        pathology=is_pathology(cf_pose, cf_native, pathology_gap),
        source_pose=source_pose,
        source_native=source_native,
        complex_dir=str(complex_dir),
        notes=";".join(notes_parts),
    )


def audit_results_dir(
    results_dir: Path,
    pathology_gap: float = DEFAULT_PATHOLOGY_GAP,
) -> List[AuditRow]:
    complexes = discover_complex_dirs(results_dir)
    return [audit_complex(c, pathology_gap=pathology_gap) for c in complexes]


# ── Presentation ──────────────────────────────────────────────────────────────

def _fmt(v: Optional[float], width: int = 10, prec: int = 3) -> str:
    if v is None:
        return f"{'n/a':>{width}}"
    return f"{v:>{width}.{prec}f}"


def format_table(rows: Sequence[AuditRow]) -> str:
    headers = (
        "pdb_id",
        "elected_pose",
        "CF_pose",
        "CF_native",
        "gap",
        "rmsd",
        "pathology",
    )
    # Column widths
    col_pdb = max(6, max((len(r.pdb_id) for r in rows), default=6))
    col_pose = max(12, max((len(r.elected_pose or "-") for r in rows), default=12))
    col_pose = min(col_pose, 28)

    lines: List[str] = []
    hdr = (
        f"{'pdb_id':<{col_pdb}}  "
        f"{'elected_pose':<{col_pose}}  "
        f"{'CF_pose':>10}  "
        f"{'CF_native':>10}  "
        f"{'gap':>10}  "
        f"{'rmsd':>8}  "
        f"{'pathology':>9}"
    )
    lines.append(hdr)
    lines.append("-" * len(hdr))
    for r in rows:
        pose = (r.elected_pose or "-")[:col_pose]
        path_flag = "YES" if r.pathology else ""
        lines.append(
            f"{r.pdb_id:<{col_pdb}}  "
            f"{pose:<{col_pose}}  "
            f"{_fmt(r.cf_pose, 10, 3)}  "
            f"{_fmt(r.cf_native, 10, 3)}  "
            f"{_fmt(r.gap, 10, 3)}  "
            f"{_fmt(r.rmsd, 8, 3)}  "
            f"{path_flag:>9}"
        )
    n_path = sum(1 for r in rows if r.pathology)
    n_have = sum(1 for r in rows if r.gap is not None)
    lines.append("")
    lines.append(
        f"complexes={len(rows)}  with_gap={n_have}  "
        f"pathology={n_path}  (gap < -pathology_threshold; lower CF is better)"
    )
    # Silence unused headers binding (kept for documentation)
    _ = headers
    return "\n".join(lines)


def build_summary(rows: Sequence[AuditRow], results_dir: Path, pathology_gap: float) -> Dict[str, Any]:
    return {
        "results_dir": str(results_dir),
        "pathology_gap_threshold": pathology_gap,
        "n_complexes": len(rows),
        "n_with_gap": sum(1 for r in rows if r.gap is not None),
        "n_pathology": sum(1 for r in rows if r.pathology),
        "rows": [r.as_jsonable() for r in rows],
    }


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Audit CF_pose vs CF_native for elected poses under a results directory. "
            "Exit 0 always for successful audit runs."
        )
    )
    p.add_argument(
        "results_dir",
        help="Per-complex dir or parent run dir (e.g. results/astex_.../1GPK or parent)",
    )
    p.add_argument(
        "--json",
        dest="json_out",
        default=None,
        metavar="PATH",
        help="Optional path to write JSON summary",
    )
    p.add_argument(
        "--pathology-gap",
        type=float,
        default=DEFAULT_PATHOLOGY_GAP,
        help=(
            f"Flag pathology when CF_pose < CF_native - GAP "
            f"(default: {DEFAULT_PATHOLOGY_GAP})"
        ),
    )
    p.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress table; still write JSON if requested",
    )
    return p.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    results_dir = resolve_path(args.results_dir)
    if not results_dir.is_dir():
        print(f"WARNING: results_dir not found or not a directory: {results_dir}", file=sys.stderr)
        if args.json_out:
            out = resolve_path(args.json_out) if not Path(args.json_out).is_absolute() else Path(args.json_out)
            # Still write empty summary for tooling
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text(
                json.dumps(
                    {
                        "results_dir": str(results_dir),
                        "error": "not_a_directory",
                        "n_complexes": 0,
                        "rows": [],
                    },
                    indent=2,
                )
                + "\n",
                encoding="utf-8",
            )
        return 0  # audit always exits 0

    rows = audit_results_dir(results_dir, pathology_gap=args.pathology_gap)
    if not rows:
        print(f"WARNING: no complex directories found under {results_dir}", file=sys.stderr)

    if not args.quiet:
        print(f"# Native CF audit: {results_dir}")
        print(f"# pathology when gap = CF_pose - CF_native < -{args.pathology_gap:g}")
        print(format_table(rows))

    if args.json_out:
        out_path = Path(args.json_out).expanduser()
        if not out_path.is_absolute():
            # Prefer CWD for explicit relative outs
            out_path = (Path.cwd() / out_path).resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        summary = build_summary(rows, results_dir, args.pathology_gap)
        out_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        if not args.quiet:
            print(f"\nJSON summary written to {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
