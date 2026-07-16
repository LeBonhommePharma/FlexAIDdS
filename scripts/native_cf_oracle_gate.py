#!/usr/bin/env python3
"""Native CF oracle gate — fail closed when crystal-like CF loses to decoys.

Scientific contract
-------------------
The GA ranks poses with the **CF/contact-function scoring proxy** (not true ΔG).
If CF_native (crystal LIG_ref / INI) is far worse than the best unseeded GA CF,
the scoring function rejects the native under the current prep/matrix. Ranking
experiments (Softβ, entropy election, etc.) are **forbidden** until this gate
passes on canary targets (pilot8: 1P62, 1T40, …).

Gate (lower CF is better):
  PASS  when CF_native <= best_ga_cf + tolerance
  FAIL  when CF_native >  best_ga_cf + tolerance   (native not competitive)
  FAIL  when CF_native is a sentinel (e.g. 10000) or missing when required

Data sources (first available wins):
  CF_native:
    1. REMARK CF= on *_INI.pdb under work or results
    2. result.csv cf_native (if finite and not sentinel)
  best_ga_cf:
    1. min REMARK CF= over numbered pose PDBs (exclude INI)
    2. result.csv best_score / elected_cf / cf_best_cluster

Usage:
  python3 scripts/native_cf_oracle_gate.py --work $Q/work/B0/1P62
  python3 scripts/native_cf_oracle_gate.py --work $WORK --results $OUT/1P62
  python3 scripts/native_cf_oracle_gate.py --ini 1P62_INI.pdb --poses-dir $OUT/1P62
  python3 scripts/native_cf_oracle_gate.py --work $W --tolerance 5.0 --json gate.json

Exit codes:
  0  native competitive
  1  native not competitive (pathology)
  2  usage / missing required inputs
  3  CF_native missing or sentinel (cannot evaluate fairness)
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

# Reuse CF parsers from audit_native_cf when available.
_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

try:
    from audit_native_cf import (  # type: ignore
        parse_cf_from_pdb_text,
        read_cf_from_pdb,
    )
except ImportError:  # pragma: no cover
    _CF_RE = re.compile(
        r"REMARK\s+CF\s*[=:]\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)",
        re.IGNORECASE,
    )

    def parse_cf_from_pdb_text(text: str) -> Optional[float]:
        for line in text.splitlines():
            m = _CF_RE.search(line)
            if m:
                try:
                    return float(m.group(1))
                except ValueError:
                    continue
        return None

    def read_cf_from_pdb(path: Path) -> Optional[float]:
        try:
            return parse_cf_from_pdb_text(
                path.read_text(encoding="utf-8", errors="ignore")
            )
        except OSError:
            return None


# CF values at or above this are treated as "unscored / sentinel" (FlexAID INI often 10000).
DEFAULT_SENTINEL_CF = 9999.0
DEFAULT_TOLERANCE = 0.0  # strict: native must be ≤ best GA (+0)

EXIT_PASS = 0
EXIT_FAIL_PATHOLOGY = 1
EXIT_USAGE = 2
EXIT_MISSING_NATIVE = 3

_POSE_NUM_RE = re.compile(r"^(.+)_(\d+)$")


@dataclass
class OracleGateResult:
    ok: bool
    exit_code: int
    cf_native: Optional[float] = None
    best_ga_cf: Optional[float] = None
    gap: Optional[float] = None  # best_ga_cf - cf_native (negative ⇒ decoy better)
    tolerance: float = DEFAULT_TOLERANCE
    source_native: str = ""
    source_ga: str = ""
    work_dir: str = ""
    results_dir: str = ""
    n_poses_scored: int = 0
    messages: List[str] = field(default_factory=list)
    ranking_forbidden: bool = False

    def as_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        for k in ("cf_native", "best_ga_cf", "gap", "tolerance"):
            v = d.get(k)
            if isinstance(v, float):
                d[k] = round(v, 6)
        return d


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
    if v != v:
        return None
    return v


def is_sentinel_cf(cf: Optional[float], sentinel: float = DEFAULT_SENTINEL_CF) -> bool:
    if cf is None:
        return True
    return cf >= float(sentinel)


def load_result_csv(path: Path) -> Optional[Dict[str, str]]:
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


def find_ini(search_roots: Sequence[Path], pdb_hint: str = "") -> Optional[Path]:
    candidates: List[Path] = []
    for root in search_roots:
        if not root.is_dir():
            continue
        if pdb_hint:
            candidates.extend(
                [
                    root / f"{pdb_hint}_INI.pdb",
                    root / f"{pdb_hint}_ini.pdb",
                ]
            )
        candidates.extend(
            [
                root / "INI.pdb",
                *sorted(root.glob("*_INI.pdb")),
                *sorted(root.glob("*_ini.pdb")),
            ]
        )
        # one-level restart / r* subdirs
        for sub in sorted(root.iterdir()):
            if not sub.is_dir():
                continue
            if sub.name.startswith("restart_") or re.fullmatch(r"r\d+", sub.name):
                candidates.extend(sorted(sub.glob("*_INI.pdb")))
                candidates.append(sub / "INI.pdb")
    for p in candidates:
        if p.is_file():
            return p
    return None


def list_pose_pdbs(root: Path) -> List[Path]:
    """Numbered pose PDBs under root (not INI, not receptor-only)."""
    if not root.is_dir():
        return []
    poses: List[Path] = []
    for p in sorted(root.glob("*.pdb")):
        upper = p.name.upper()
        if upper.endswith("_INI.PDB") or "NATIVE" in upper and not _POSE_NUM_RE.match(p.stem):
            continue
        if upper.endswith("_PRUNED.PDB"):
            continue
        if _POSE_NUM_RE.match(p.stem):
            poses.append(p)
    # Also search r*/restart_* one level
    for sub in sorted(root.iterdir()):
        if not sub.is_dir():
            continue
        if not (sub.name.startswith("restart_") or re.fullmatch(r"r\d+", sub.name)):
            continue
        for p in sorted(sub.glob("*.pdb")):
            if p.name.upper().endswith("_INI.PDB"):
                continue
            if _POSE_NUM_RE.match(p.stem):
                poses.append(p)
    return poses


def best_cf_from_poses(poses: Sequence[Path]) -> Tuple[Optional[float], str, int]:
    best: Optional[float] = None
    best_src = ""
    n = 0
    for p in poses:
        cf = read_cf_from_pdb(p)
        if cf is None or is_sentinel_cf(cf):
            continue
        n += 1
        if best is None or cf < best:
            best = cf
            best_src = f"REMARK:{p.name}"
    return best, best_src, n


def resolve_cf_native(
    ini: Optional[Path],
    csv_row: Optional[Dict[str, str]],
    sentinel: float,
) -> Tuple[Optional[float], str]:
    ini_sentinel_seen = False
    if ini is not None and ini.is_file():
        cf = read_cf_from_pdb(ini)
        if cf is not None and not is_sentinel_cf(cf, sentinel):
            return cf, f"REMARK:{ini.name}"
        if cf is not None and is_sentinel_cf(cf, sentinel):
            ini_sentinel_seen = True
            # Fall through to CSV for a real native CF if present
    if csv_row:
        for key in ("cf_native", "native_cf", "CF_native"):
            cf = _float_or_none(csv_row.get(key))
            if cf is None or is_sentinel_cf(cf, sentinel):
                continue
            # Placeholder 0.0 with unscored INI is not a real native CF
            if ini_sentinel_seen and abs(cf) < 1e-9:
                continue
            return cf, f"result.csv:{key}"
    if ini_sentinel_seen:
        return None, "INI_sentinel_no_valid_cf_native"
    return None, ""


def resolve_best_ga_cf(
    pose_dirs: Sequence[Path],
    csv_row: Optional[Dict[str, str]],
) -> Tuple[Optional[float], str, int]:
    all_poses: List[Path] = []
    for d in pose_dirs:
        all_poses.extend(list_pose_pdbs(d))
    # Dedup by resolve
    seen = set()
    uniq: List[Path] = []
    for p in all_poses:
        rp = p.resolve()
        if rp in seen:
            continue
        seen.add(rp)
        uniq.append(p)

    best, src, n = best_cf_from_poses(uniq)
    if best is not None:
        return best, src, n

    if csv_row:
        for key in (
            "best_score",
            "elected_cf",
            "cf_best_cluster",
            "best_cf",
            "CF_best",
        ):
            cf = _float_or_none(csv_row.get(key))
            if cf is not None and not is_sentinel_cf(cf):
                return cf, f"result.csv:{key}", 0
    return None, "", 0


def evaluate_gate(
    *,
    cf_native: Optional[float],
    best_ga_cf: Optional[float],
    tolerance: float,
    source_native: str = "",
    source_ga: str = "",
    n_poses: int = 0,
    work_dir: str = "",
    results_dir: str = "",
    sentinel: float = DEFAULT_SENTINEL_CF,
    require_poses: bool = True,
) -> OracleGateResult:
    msgs: List[str] = []

    if cf_native is None:
        msgs.append("CF_native missing — cannot evaluate native competitiveness")
        return OracleGateResult(
            ok=False,
            exit_code=EXIT_MISSING_NATIVE,
            cf_native=None,
            best_ga_cf=best_ga_cf,
            tolerance=tolerance,
            source_native=source_native,
            source_ga=source_ga,
            work_dir=work_dir,
            results_dir=results_dir,
            n_poses_scored=n_poses,
            messages=msgs,
            ranking_forbidden=True,
        )

    if is_sentinel_cf(cf_native, sentinel):
        msgs.append(
            f"CF_native={cf_native} is sentinel (>= {sentinel}); "
            "INI was not scored in pocket. Re-score crystal or fix prep."
        )
        return OracleGateResult(
            ok=False,
            exit_code=EXIT_MISSING_NATIVE,
            cf_native=cf_native,
            best_ga_cf=best_ga_cf,
            tolerance=tolerance,
            source_native=source_native,
            source_ga=source_ga,
            work_dir=work_dir,
            results_dir=results_dir,
            n_poses_scored=n_poses,
            messages=msgs,
            ranking_forbidden=True,
        )

    if best_ga_cf is None:
        if require_poses:
            msgs.append(
                "best_ga_cf missing — no scored poses / result.csv; "
                "gate cannot confirm competitiveness"
            )
            return OracleGateResult(
                ok=False,
                exit_code=EXIT_USAGE,
                cf_native=cf_native,
                best_ga_cf=None,
                tolerance=tolerance,
                source_native=source_native,
                source_ga=source_ga,
                work_dir=work_dir,
                results_dir=results_dir,
                n_poses_scored=n_poses,
                messages=msgs,
                ranking_forbidden=True,
            )
        msgs.append("No GA poses yet; only CF_native available (preflight incomplete)")
        return OracleGateResult(
            ok=False,
            exit_code=EXIT_USAGE,
            cf_native=cf_native,
            best_ga_cf=None,
            tolerance=tolerance,
            source_native=source_native,
            source_ga=source_ga,
            work_dir=work_dir,
            results_dir=results_dir,
            n_poses_scored=n_poses,
            messages=msgs,
            ranking_forbidden=True,
        )

    # gap = best_ga - native; negative means decoy better (pathology)
    gap = best_ga_cf - cf_native
    # PASS: native <= best_ga + tol  ⇔  cf_native - best_ga_cf <= tol
    # ⇔  -gap <= tol  ⇔ gap >= -tol when... wait:
    #   CF_native <= best_ga_cf + tolerance
    competitive = cf_native <= (best_ga_cf + tolerance)

    if competitive:
        msgs.append(
            f"PASS: CF_native={cf_native:.4f} <= best_ga_cf={best_ga_cf:.4f} "
            f"+ tol={tolerance:.4f} (gap_ga_minus_native={gap:.4f})"
        )
        msgs.append(
            "Native is competitive under CF proxy — ranking experiments may proceed "
            "only with this documented gate pass."
        )
        return OracleGateResult(
            ok=True,
            exit_code=EXIT_PASS,
            cf_native=cf_native,
            best_ga_cf=best_ga_cf,
            gap=gap,
            tolerance=tolerance,
            source_native=source_native,
            source_ga=source_ga,
            work_dir=work_dir,
            results_dir=results_dir,
            n_poses_scored=n_poses,
            messages=msgs,
            ranking_forbidden=False,
        )

    msgs.append(
        f"FAIL: CF_native={cf_native:.4f} > best_ga_cf={best_ga_cf:.4f} "
        f"+ tol={tolerance:.4f} (gap_ga_minus_native={gap:.4f}; decoy better by {-gap:.4f})"
    )
    msgs.append(
        "RANKING EXPERIMENTS FORBIDDEN until gate passes on canary: "
        "CF rejects native under this prep/matrix. Do NOT claim Softβ fixes sampling."
    )
    return OracleGateResult(
        ok=False,
        exit_code=EXIT_FAIL_PATHOLOGY,
        cf_native=cf_native,
        best_ga_cf=best_ga_cf,
        gap=gap,
        tolerance=tolerance,
        source_native=source_native,
        source_ga=source_ga,
        work_dir=work_dir,
        results_dir=results_dir,
        n_poses_scored=n_poses,
        messages=msgs,
        ranking_forbidden=True,
    )


def run_gate(
    *,
    work: Optional[Path] = None,
    results: Optional[Path] = None,
    ini: Optional[Path] = None,
    poses_dir: Optional[Path] = None,
    pdb_id: str = "",
    tolerance: float = DEFAULT_TOLERANCE,
    sentinel: float = DEFAULT_SENTINEL_CF,
    require_poses: bool = True,
) -> OracleGateResult:
    search_roots: List[Path] = []
    pose_dirs: List[Path] = []
    csv_row: Optional[Dict[str, str]] = None

    if work is not None:
        work = work.expanduser().resolve()
        search_roots.append(work)
        pose_dirs.append(work)
        if (work / "result.csv").is_file():
            csv_row = load_result_csv(work / "result.csv")
        # Guess pdb from work name
        if not pdb_id:
            pdb_id = work.name.upper()

    if results is not None:
        results = results.expanduser().resolve()
        search_roots.append(results)
        pose_dirs.append(results)
        if (results / "result.csv").is_file():
            csv_row = load_result_csv(results / "result.csv") or csv_row
        if not pdb_id:
            pdb_id = results.name.upper()

    if poses_dir is not None:
        poses_dir = poses_dir.expanduser().resolve()
        pose_dirs.append(poses_dir)
        search_roots.append(poses_dir)

    if ini is not None:
        ini = ini.expanduser().resolve()
    else:
        ini = find_ini(search_roots, pdb_hint=pdb_id)

    cf_native, src_n = resolve_cf_native(ini, csv_row, sentinel)
    best_ga, src_g, n_poses = resolve_best_ga_cf(pose_dirs, csv_row)

    return evaluate_gate(
        cf_native=cf_native,
        best_ga_cf=best_ga,
        tolerance=tolerance,
        source_native=src_n,
        source_ga=src_g,
        n_poses=n_poses,
        work_dir=str(work) if work else "",
        results_dir=str(results) if results else "",
        sentinel=sentinel,
        require_poses=require_poses,
    )


def format_result(res: OracleGateResult) -> str:
    status = "PASS" if res.ok else "FAIL"
    lines = [
        f"{status} exit={res.exit_code} ranking_forbidden={res.ranking_forbidden}",
        f"  CF_native={res.cf_native}  ({res.source_native or 'n/a'})",
        f"  best_ga_cf={res.best_ga_cf}  ({res.source_ga or 'n/a'})  n_poses={res.n_poses_scored}",
        f"  gap(best_ga - native)={res.gap}  tolerance={res.tolerance}",
    ]
    for m in res.messages:
        lines.append(f"  {m}")
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--work", type=Path, help="Prep/work dir (CONFIG, LIG_ref, TARGET)")
    ap.add_argument("--results", type=Path, help="Results dir with poses / result.csv")
    ap.add_argument("--ini", type=Path, help="Explicit INI PDB with REMARK CF=")
    ap.add_argument("--poses-dir", type=Path, help="Directory of numbered pose PDBs")
    ap.add_argument("--pdb", default="", help="PDB id hint for file discovery")
    ap.add_argument(
        "--tolerance",
        type=float,
        default=DEFAULT_TOLERANCE,
        help="Allow CF_native <= best_ga + tolerance (default 0)",
    )
    ap.add_argument(
        "--sentinel",
        type=float,
        default=DEFAULT_SENTINEL_CF,
        help=f"CF >= this treated as unscored (default {DEFAULT_SENTINEL_CF})",
    )
    ap.add_argument(
        "--allow-missing-poses",
        action="store_true",
        help="Do not require GA poses (still fails if best_ga missing)",
    )
    ap.add_argument("--json", dest="json_out", type=Path, default=None)
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args(argv)

    if not any([args.work, args.results, args.ini, args.poses_dir]):
        print(
            "ERROR: provide at least one of --work / --results / --ini / --poses-dir",
            file=sys.stderr,
        )
        return EXIT_USAGE

    res = run_gate(
        work=args.work,
        results=args.results,
        ini=args.ini,
        poses_dir=args.poses_dir,
        pdb_id=args.pdb.upper(),
        tolerance=args.tolerance,
        sentinel=args.sentinel,
        require_poses=not args.allow_missing_poses,
    )

    if not args.quiet:
        print(format_result(res))

    if args.json_out:
        out = args.json_out.expanduser()
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(res.as_dict(), indent=2) + "\n", encoding="utf-8")
        if not args.quiet:
            print(f"JSON written to {out}")

    return int(res.exit_code)


if __name__ == "__main__":
    sys.exit(main())
