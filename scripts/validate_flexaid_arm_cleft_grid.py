#!/usr/bin/env python3
"""Validate identical cleft / grid / matrix across FlexAID arms A, B0, B.

For each PDB in pilot8 (default) or full85 / custom list, checks work trees:

  {work_root}/{A,B0,B}/{PDB}/{PDB}_spheres.pdb   → MD5 must match across arms
  {work_root}/{A,B0,B}/{PDB}/CONFIG.inp          → SPACER, LOCCLF basename,
                                                    IMATRX md5 must match
                                                  (TEMPER / CLUSTA may differ)

Optional runtime log check for:
  "will build a grid with spacing 0.375"

Does **not** touch running docks — read-only.

Usage:
  python3 scripts/validate_flexaid_arm_cleft_grid.py
  python3 scripts/validate_flexaid_arm_cleft_grid.py --panel pilot8
  python3 scripts/validate_flexaid_arm_cleft_grid.py --panel full85
  python3 scripts/validate_flexaid_arm_cleft_grid.py --check-logs
  python3 scripts/validate_flexaid_arm_cleft_grid.py --json /tmp/cleft_grid_report.json

Exit codes:
  0 — all checked targets MATCH (or no targets and --allow-empty)
  1 — one or more FAIL
  2 — usage / path errors
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

# Keep in sync with scripts/generate_flexaid_inp.py
PILOT8 = ["1G9V", "1GPK", "1MEH", "1P62", "1Q4G", "1R9O", "1T40", "2BYS"]
ARMS_DEFAULT = ("A", "B0", "B")
MATRIX_MD5_PIN = "72d7c7396702331d96ff12d18f831796"
EXPECTED_SPACER = "0.375"
GRID_LOG_RE = re.compile(
    r"will build a grid with spacing\s+([0-9]*\.?[0-9]+)", re.IGNORECASE
)

# Keys that MUST be identical across arms (cleft/grid/scoring matrix).
MUST_MATCH_KEYS = ("spheres_md5", "spacer", "locclf_basename", "imatrix_md5")
# Keys that MAY differ by design (entropy arm protocol).
MAY_DIFFER_KEYS = ("temper", "clusta")


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def default_work_root() -> Path:
    env = os.environ.get("FLEXAID_WORK_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    local = os.environ.get("FLEXAIDDS_LOCAL_ROOT", str(Path.home() / "flexaidds_results"))
    return (Path(local).expanduser() / "three_engine_entropy_q1" / "work").resolve()


def default_log_candidates() -> List[Path]:
    """Likely arm runtime logs (local-first layout)."""
    local = Path(
        os.environ.get("FLEXAIDDS_LOCAL_ROOT", str(Path.home() / "flexaidds_results"))
    ).expanduser()
    logdir = os.environ.get("FLEXAIDDS_LOCAL_LOGDIR")
    cands: List[Path] = []
    if logdir:
        cands.append(Path(logdir).expanduser())
    cands.append(local / "logs" / "three_engine")
    q = os.environ.get("FLEXAIDDS_QUEUE_ROOT")
    if q:
        cands.append(Path(q).expanduser() / "logs")
    out: List[Path] = []
    seen = set()
    for d in cands:
        try:
            r = d.resolve()
        except OSError:
            continue
        if r in seen:
            continue
        seen.add(r)
        out.append(r)
    return out


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_config(path: Path) -> Dict[str, str]:
    """Extract SPACER, LOCCLF, IMATRX, CLUSTA, TEMPER from CONFIG.inp."""
    out: Dict[str, str] = {}
    text = path.read_text(errors="replace")
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        # RNGOPT LOCCLF /path/to/spheres.pdb
        if line.upper().startswith("RNGOPT"):
            parts = line.split()
            # ["RNGOPT", "LOCCLF", path] or ["RNGOPT", "LOCCLF", path, ...]
            if len(parts) >= 3 and parts[1].upper() == "LOCCLF":
                out["locclf_path"] = parts[2]
                out["locclf_basename"] = Path(parts[2]).name
            continue
        toks = line.split(None, 1)
        key = toks[0].upper()
        val = toks[1].strip() if len(toks) > 1 else ""
        if key == "SPACER":
            out["spacer"] = val.split()[0] if val else ""
        elif key == "IMATRX":
            out["imatrix_path"] = val.split()[0] if val else ""
        elif key == "CLUSTA":
            out["clusta"] = val.split()[0] if val else ""
        elif key == "TEMPER":
            out["temper"] = val.split()[0] if val else ""
    return out


def resolve_config(work_pdb: Path) -> Optional[Path]:
    """Prefer top-level CONFIG.inp; fall back to restart_0."""
    top = work_pdb / "CONFIG.inp"
    if top.is_file():
        return top
    r0 = work_pdb / "restart_0" / "CONFIG.inp"
    if r0.is_file():
        return r0
    return None


def load_pdb_list(
    panel: str,
    pdb: Optional[str],
    list_file: Optional[str],
    work_root: Path,
    arms: Sequence[str],
) -> List[str]:
    if pdb:
        return [pdb.upper()]
    if list_file:
        return [
            ln.strip().upper()
            for ln in Path(list_file).read_text().splitlines()
            if ln.strip() and not ln.startswith("#")
        ]
    if panel == "pilot8":
        return list(PILOT8)
    if panel == "full85":
        # Prefer union of PDB dirs present under any arm work tree; else inputs/astex_diverse.
        found: set[str] = set()
        for arm in arms:
            arm_dir = work_root / arm
            if not arm_dir.is_dir():
                continue
            for p in arm_dir.iterdir():
                if p.is_dir() and re.fullmatch(r"[0-9A-Za-z]{4}", p.name):
                    found.add(p.name.upper())
        if found:
            return sorted(found)
        # Fallback: queue inputs if FLEXAIDDS_QUEUE_ROOT set
        q = os.environ.get("FLEXAIDDS_QUEUE_ROOT")
        if q:
            d = Path(q).expanduser() / "inputs" / "astex_diverse"
            if d.is_dir():
                return sorted(
                    p.name.upper()
                    for p in d.iterdir()
                    if p.is_dir() and re.fullmatch(r"[0-9A-Za-z]{4}", p.name)
                )
        # Repo Astex Diverse tree
        repo_astex = (
            repo_root() / "benchmarks" / "astex_diverse" / "astex_diverse"
        )
        if repo_astex.is_dir():
            return sorted(
                p.name.upper()
                for p in repo_astex.iterdir()
                if p.is_dir() and re.fullmatch(r"[0-9A-Za-z]{4}", p.name)
            )
        raise SystemExit(
            "full85: no PDB dirs under work arms and no Astex tree found"
        )
    raise SystemExit(f"unknown panel: {panel}")


@dataclass
class ArmSnapshot:
    arm: str
    present: bool
    spheres_path: Optional[str] = None
    spheres_md5: Optional[str] = None
    config_path: Optional[str] = None
    spacer: Optional[str] = None
    locclf_basename: Optional[str] = None
    locclf_path: Optional[str] = None
    imatrix_path: Optional[str] = None
    imatrix_md5: Optional[str] = None
    temper: Optional[str] = None
    clusta: Optional[str] = None
    errors: List[str] = field(default_factory=list)


@dataclass
class TargetReport:
    pdb: str
    status: str  # MATCH | FAIL | SKIP
    arms: Dict[str, dict] = field(default_factory=dict)
    mismatches: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)


def snapshot_arm(work_root: Path, arm: str, pdb: str) -> ArmSnapshot:
    wdir = work_root / arm / pdb
    snap = ArmSnapshot(arm=arm, present=wdir.is_dir())
    if not snap.present:
        snap.errors.append(f"missing work dir {wdir}")
        return snap

    spheres = wdir / f"{pdb}_spheres.pdb"
    if spheres.is_file():
        snap.spheres_path = str(spheres)
        try:
            snap.spheres_md5 = md5_file(spheres)
        except OSError as e:
            snap.errors.append(f"spheres md5 failed: {e}")
    else:
        snap.errors.append(f"missing spheres {spheres}")

    cfg = resolve_config(wdir)
    if cfg is None:
        snap.errors.append(f"missing CONFIG.inp under {wdir}")
        return snap
    snap.config_path = str(cfg)
    try:
        fields = parse_config(cfg)
    except OSError as e:
        snap.errors.append(f"CONFIG read failed: {e}")
        return snap

    snap.spacer = fields.get("spacer")
    snap.locclf_basename = fields.get("locclf_basename")
    snap.locclf_path = fields.get("locclf_path")
    snap.imatrix_path = fields.get("imatrix_path")
    snap.temper = fields.get("temper")
    snap.clusta = fields.get("clusta")

    if not snap.spacer:
        snap.errors.append("CONFIG missing SPACER")
    if not snap.locclf_basename:
        snap.errors.append("CONFIG missing RNGOPT LOCCLF")
    if not snap.imatrix_path:
        snap.errors.append("CONFIG missing IMATRX")
    else:
        ip = Path(snap.imatrix_path)
        if ip.is_file():
            try:
                snap.imatrix_md5 = md5_file(ip)
            except OSError as e:
                snap.errors.append(f"IMATRX md5 failed: {e}")
        else:
            snap.errors.append(f"IMATRX file missing: {ip}")

    return snap


def compare_target(
    pdb: str, snaps: Dict[str, ArmSnapshot], require_all_arms: bool
) -> TargetReport:
    present = {a: s for a, s in snaps.items() if s.present}
    rep = TargetReport(pdb=pdb, status="MATCH")
    rep.arms = {a: asdict(s) for a, s in snaps.items()}

    if len(present) < 2:
        rep.status = "FAIL" if require_all_arms else "SKIP"
        rep.mismatches.append(
            f"need ≥2 arms with work trees; present={sorted(present)}"
        )
        return rep

    if require_all_arms and len(present) < len(snaps):
        missing = sorted(set(snaps) - set(present))
        rep.status = "FAIL"
        rep.mismatches.append(f"missing arms: {missing}")

    for a, s in present.items():
        if s.errors:
            rep.status = "FAIL"
            for e in s.errors:
                rep.mismatches.append(f"{a}: {e}")

    # Cross-arm must-match
    for key in MUST_MATCH_KEYS:
        vals = {a: getattr(present[a], key) for a in present}
        unique = {v for v in vals.values() if v is not None}
        if not unique:
            rep.status = "FAIL"
            rep.mismatches.append(f"{key}: no values across arms")
            continue
        if len(unique) > 1:
            rep.status = "FAIL"
            rep.mismatches.append(f"{key} differs: {vals}")
        elif None in vals.values():
            rep.status = "FAIL"
            rep.mismatches.append(f"{key} incomplete: {vals}")

    # Expected spacer value (when present)
    spacers = {getattr(s, "spacer") for s in present.values() if s.spacer}
    if spacers and spacers != {EXPECTED_SPACER}:
        # Not automatically FAIL if arms agree but non-canonical; flag as note
        # unless they disagree (already FAIL above).
        if len(spacers) == 1:
            rep.notes.append(
                f"SPACER={next(iter(spacers))} (protocol expects {EXPECTED_SPACER})"
            )

    # Matrix pin advisory
    imds = {s.imatrix_md5 for s in present.values() if s.imatrix_md5}
    if imds and imds != {MATRIX_MD5_PIN}:
        rep.notes.append(
            f"IMATRX md5={sorted(imds)} (protocol pin {MATRIX_MD5_PIN})"
        )

    # LOCCLF basename should match {PDB}_spheres.pdb
    expected_loc = f"{pdb}_spheres.pdb"
    for a, s in present.items():
        if s.locclf_basename and s.locclf_basename != expected_loc:
            rep.notes.append(
                f"{a}: LOCCLF basename {s.locclf_basename!r} != {expected_loc!r}"
            )

    # Record may-differ values for the report (informational)
    for key in MAY_DIFFER_KEYS:
        vals = {a: getattr(present[a], key) for a in present}
        if len({v for v in vals.values() if v is not None}) > 1:
            rep.notes.append(f"{key} differs by design: {vals}")

    return rep


def scan_logs_for_grid(
    log_paths: Iterable[Path],
    max_bytes_per_file: int = 80_000_000,
) -> Dict[str, object]:
    """Scan log files for 'will build a grid with spacing X' lines."""
    hits: List[dict] = []
    files_scanned = 0
    files_missing = 0
    spacings: Dict[str, int] = {}

    expanded: List[Path] = []
    for p in log_paths:
        if p.is_dir():
            for name in (
                "run_3dsig_red_pair_serial.log",
                "run_A_pilot8.log",
                "run_B0_pilot8.log",
                "run_B_pilot8.log",
            ):
                cand = p / name
                if cand.is_file():
                    expanded.append(cand)
            # Also any run_*_pilot8.log / *3dsig*.log
            for cand in sorted(p.glob("run_*pilot8*.log")) + sorted(
                p.glob("*3dsig*.log")
            ):
                if cand not in expanded:
                    expanded.append(cand)
        elif p.is_file():
            expanded.append(p)
        else:
            files_missing += 1

    # de-dupe
    seen = set()
    uniq: List[Path] = []
    for p in expanded:
        r = p.resolve()
        if r in seen:
            continue
        seen.add(r)
        uniq.append(p)

    for path in uniq:
        files_scanned += 1
        try:
            size = path.stat().st_size
            # Read from start in chunks; grid message is early in each docking run
            # but serial log is multi-target — scan whole file up to cap.
            remaining = min(size, max_bytes_per_file)
            with path.open("rb") as f:
                buf = b""
                line_no = 0
                while remaining > 0:
                    chunk = f.read(min(1 << 20, remaining))
                    if not chunk:
                        break
                    remaining -= len(chunk)
                    buf += chunk
                    while b"\n" in buf:
                        raw, buf = buf.split(b"\n", 1)
                        line_no += 1
                        try:
                            line = raw.decode("utf-8", errors="replace")
                        except Exception:
                            continue
                        m = GRID_LOG_RE.search(line)
                        if m:
                            sp = m.group(1)
                            spacings[sp] = spacings.get(sp, 0) + 1
                            if len(hits) < 50:
                                hits.append(
                                    {
                                        "file": str(path),
                                        "line": line_no,
                                        "spacing": sp,
                                        "text": line.strip()[:200],
                                    }
                                )
        except OSError as e:
            hits.append({"file": str(path), "error": str(e)})

    ok = bool(spacings) and set(spacings) == {EXPECTED_SPACER}
    status = "MATCH" if ok else ("FAIL" if spacings else "SKIP")
    return {
        "status": status,
        "expected_spacing": EXPECTED_SPACER,
        "spacings_seen": spacings,
        "files_scanned": files_scanned,
        "files_missing_or_empty_dirs": files_missing,
        "sample_hits": hits,
        "note": (
            "MATCH means every 'will build a grid with spacing X' hit used "
            f"{EXPECTED_SPACER}; SKIP if no hits found"
        ),
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Validate same cleft/grid/matrix across FlexAID arms A/B0/B"
    )
    ap.add_argument(
        "--work-root",
        type=Path,
        default=None,
        help="Work root containing A/, B0/, B/ (default: FLEXAID_WORK_ROOT or local)",
    )
    ap.add_argument(
        "--arms",
        default="A,B0,B",
        help="Comma-separated arms (default: A,B0,B)",
    )
    ap.add_argument(
        "--panel",
        choices=("pilot8", "full85"),
        default="pilot8",
        help="PDB panel (default: pilot8)",
    )
    ap.add_argument("--pdb", help="Single PDB id")
    ap.add_argument("--list-file", help="File with one PDB id per line")
    ap.add_argument(
        "--require-all-arms",
        action="store_true",
        help="FAIL if any listed arm is missing for a PDB",
    )
    ap.add_argument(
        "--check-logs",
        action="store_true",
        help="Also scan runtime logs for 'will build a grid with spacing 0.375'",
    )
    ap.add_argument(
        "--log",
        action="append",
        default=[],
        help="Log file or directory (repeatable). Default: local three_engine logs",
    )
    ap.add_argument("--json", type=Path, help="Write full JSON report")
    ap.add_argument(
        "--allow-empty",
        action="store_true",
        help="Exit 0 when no PDBs to check",
    )
    args = ap.parse_args(argv)

    work_root = (args.work_root or default_work_root()).expanduser().resolve()
    arms = tuple(a.strip() for a in args.arms.split(",") if a.strip())
    if not arms:
        print("ERROR: no arms", file=sys.stderr)
        return 2

    if not work_root.is_dir():
        print(f"ERROR: work root not found: {work_root}", file=sys.stderr)
        return 2

    try:
        pdbs = load_pdb_list(
            args.panel, args.pdb, args.list_file, work_root, arms
        )
    except SystemExit as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2

    if not pdbs:
        print("ERROR: empty PDB list", file=sys.stderr)
        return 0 if args.allow_empty else 2

    print(f"=== validate_flexaid_arm_cleft_grid ===")
    print(f"work_root={work_root}")
    print(f"arms={list(arms)} panel={args.panel} n_pdb={len(pdbs)}")
    print(f"must_match={list(MUST_MATCH_KEYS)}")
    print(f"may_differ={list(MAY_DIFFER_KEYS)}")
    print()

    reports: List[TargetReport] = []
    n_match = n_fail = n_skip = 0

    for pdb in pdbs:
        snaps = {a: snapshot_arm(work_root, a, pdb) for a in arms}
        rep = compare_target(pdb, snaps, require_all_arms=args.require_all_arms)
        reports.append(rep)
        if rep.status == "MATCH":
            n_match += 1
            tag = "MATCH"
        elif rep.status == "SKIP":
            n_skip += 1
            tag = "SKIP"
        else:
            n_fail += 1
            tag = "FAIL"

        # Compact one-liner
        present_arms = [a for a, s in snaps.items() if s.present]
        md5s = {
            a: (snaps[a].spheres_md5 or "?")[:8]
            for a in present_arms
        }
        spacers = {a: snaps[a].spacer for a in present_arms}
        tempers = {a: snaps[a].temper for a in present_arms}
        clustas = {a: snaps[a].clusta for a in present_arms}
        print(
            f"{tag:5} {pdb}  arms={present_arms}  "
            f"spheres={md5s}  SPACER={spacers}  "
            f"TEMPER={tempers}  CLUSTA={clustas}"
        )
        for m in rep.mismatches:
            print(f"      ! {m}")
        for n in rep.notes:
            print(f"      · {n}")

    log_report: Optional[dict] = None
    if args.check_logs:
        log_paths: List[Path] = (
            [Path(p).expanduser() for p in args.log]
            if args.log
            else default_log_candidates()
        )
        print()
        print("--- optional log grid-spacing check ---")
        log_report = scan_logs_for_grid(log_paths)
        print(
            f"logs: status={log_report['status']} "
            f"spacings_seen={log_report['spacings_seen']} "
            f"files_scanned={log_report['files_scanned']}"
        )
        for h in (log_report.get("sample_hits") or [])[:5]:
            if "error" in h:
                print(f"  log error {h.get('file')}: {h['error']}")
            else:
                print(
                    f"  {h.get('file')}:{h.get('line')}: spacing={h.get('spacing')}"
                )
        if log_report["status"] == "FAIL":
            n_fail += 1

    print()
    print(
        f"SUMMARY  MATCH={n_match}  FAIL={n_fail}  SKIP={n_skip}  "
        f"total={len(pdbs)}"
    )
    overall = "PASS" if n_fail == 0 else "FAIL"
    print(f"OVERALL  {overall}")

    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "work_root": str(work_root),
        "arms": list(arms),
        "panel": args.panel,
        "pdbs": pdbs,
        "must_match": list(MUST_MATCH_KEYS),
        "may_differ": list(MAY_DIFFER_KEYS),
        "expected_spacer": EXPECTED_SPACER,
        "matrix_md5_pin": MATRIX_MD5_PIN,
        "summary": {
            "match": n_match,
            "fail": n_fail,
            "skip": n_skip,
            "total": len(pdbs),
            "overall": overall,
        },
        "targets": [asdict(r) for r in reports],
        "log_grid_check": log_report,
    }

    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(payload, indent=2) + "\n")
        print(f"wrote {args.json}")

    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
