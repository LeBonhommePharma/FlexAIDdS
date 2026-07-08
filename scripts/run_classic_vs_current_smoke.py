#!/usr/bin/env python3
"""Run FlexAID-classic vs current FlexAIDdS smoke docks on one Astex target.

Compares steric/scoring stacks with equal GA budget. Does not invent RMSD —
parses REMARK CF from elected pose PDBs and optional INI, plus result.csv if
the engine wrote one.

Usage (from repo root):
  python3 scripts/run_classic_vs_current_smoke.py \\
    --binary /path/to/FlexAIDdS \\
    --pdb-id 1GPK \\
    --out-dir results/smoke_classic_vs_current

Environment:
  FLEXAIDDS_BIN  optional default for --binary
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def parse_cf_remarks(pdb_path: Path) -> Dict[str, Optional[float]]:
    text = pdb_path.read_text(errors="replace")
    out: Dict[str, Optional[float]] = {
        "CF": None,
        "CF.com": None,
        "CF.wal": None,
        "CF.sas": None,
        "CF.hbond": None,
    }
    m = re.search(r"REMARK CF=([-\d.]+)", text)
    if m:
        out["CF"] = float(m.group(1))
    for key in ("CF.com", "CF.wal", "CF.sas", "CF.hbond", "CF.app"):
        m = re.search(rf"REMARK {re.escape(key)}=([-\d.]+)", text)
        if m:
            out[key] = float(m.group(1))
    return out


def find_elected_pose(out_prefix: Path) -> Optional[Path]:
    """Prefer *_0.pdb (rank 0 cluster) next to output prefix."""
    parent = out_prefix.parent
    stem = out_prefix.name
    candidates = [
        parent / f"{stem}_0.pdb",
        parent / f"{stem}.pdb",
    ]
    # Also FlexAIDdS sometimes writes PDBID_0.pdb based on complex id
    candidates.extend(sorted(parent.glob("*_0.pdb")))
    candidates.extend(sorted(parent.glob("*.pdb")))
    seen = set()
    for c in candidates:
        if c.is_file() and c.name.upper() != "INI" and "INI" not in c.name:
            if c.resolve() in seen:
                continue
            seen.add(c.resolve())
            # Skip pure INI
            if c.name.endswith("_INI.pdb"):
                continue
            return c
    return None


def find_ini(out_dir: Path) -> Optional[Path]:
    for p in out_dir.glob("*_INI.pdb"):
        return p
    for p in out_dir.glob("*INI*.pdb"):
        return p
    return None


def read_result_csv(out_dir: Path) -> Optional[Dict[str, str]]:
    p = out_dir / "result.csv"
    if not p.is_file():
        # search one level
        hits = list(out_dir.glob("**/result.csv"))
        if not hits:
            return None
        p = hits[0]
    with p.open() as f:
        rows = list(csv.DictReader(f))
    return rows[0] if rows else None


def run_arm(
    *,
    binary: Path,
    receptor: Path,
    ligand: Path,
    config: Path,
    out_prefix: Path,
    timeout_s: int,
) -> Dict[str, Any]:
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    log_path = out_prefix.parent / f"{out_prefix.name}.log"
    cmd = [
        str(binary),
        str(receptor),
        str(ligand),
        "-c",
        str(config),
        "-o",
        str(out_prefix),
    ]
    t0 = time.monotonic()
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(repo_root()),
            capture_output=True,
            text=True,
            timeout=timeout_s,
        )
        elapsed = time.monotonic() - t0
        log_path.write_text(
            f"CMD: {' '.join(cmd)}\nEXIT: {proc.returncode}\nELAPSED_S: {elapsed:.2f}\n\n"
            f"=== STDOUT ===\n{proc.stdout}\n\n=== STDERR ===\n{proc.stderr}\n"
        )
        return {
            "ok": proc.returncode == 0,
            "returncode": proc.returncode,
            "elapsed_s": elapsed,
            "log": str(log_path),
            "cmd": cmd,
        }
    except subprocess.TimeoutExpired as exc:
        elapsed = time.monotonic() - t0
        log_path.write_text(
            f"CMD: {' '.join(cmd)}\nTIMEOUT after {elapsed:.1f}s\n"
            f"stdout partial:\n{exc.stdout}\nstderr partial:\n{exc.stderr}\n"
        )
        return {
            "ok": False,
            "returncode": -1,
            "elapsed_s": elapsed,
            "log": str(log_path),
            "cmd": cmd,
            "error": "timeout",
        }


def summarize_arm(name: str, out_dir: Path, out_prefix: Path, run_meta: Dict[str, Any]) -> Dict[str, Any]:
    pose = find_elected_pose(out_prefix)
    ini = find_ini(out_dir)
    # also look for INI next to prefix stem
    if ini is None:
        cand = out_prefix.parent / f"{out_prefix.name}_INI.pdb"
        if cand.is_file():
            ini = cand
    row: Dict[str, Any] = {
        "arm": name,
        "run_ok": run_meta.get("ok"),
        "elapsed_s": run_meta.get("elapsed_s"),
        "returncode": run_meta.get("returncode"),
        "log": run_meta.get("log"),
        "elected_pose": str(pose) if pose else None,
        "ini_path": str(ini) if ini else None,
    }
    if pose:
        row.update({f"pose_{k}": v for k, v in parse_cf_remarks(pose).items()})
    if ini:
        row.update({f"ini_{k}": v for k, v in parse_cf_remarks(ini).items()})
    if pose and ini and row.get("pose_CF") is not None and row.get("ini_CF") is not None:
        row["gap_pose_minus_ini"] = row["pose_CF"] - row["ini_CF"]
        # pathology: pose much more favorable (more negative) than native/INI
        row["decoy_beats_native"] = row["pose_CF"] < (row["ini_CF"] - 5.0)
    rc = read_result_csv(out_dir)
    if rc:
        row["result_csv"] = {k: rc.get(k) for k in (
            "pdb_id", "best_score", "rmsd_to_crystal", "cf_native",
            "predicted_dG", "num_poses", "success", "G_bind", "H_vct",
        ) if k in rc}
        try:
            if rc.get("best_score") and rc.get("cf_native"):
                row["gap_best_minus_cf_native"] = float(rc["best_score"]) - float(rc["cf_native"])
        except (TypeError, ValueError):
            pass
    return row


def side_by_side(classic: Dict[str, Any], current: Dict[str, Any]) -> str:
    lines = [
        "# FlexAID-classic vs current FlexAIDdS — single-target smoke",
        "",
        "| Quantity | Classic | Current |",
        "|----------|---------|---------|",
    ]

    def cell(d: Dict[str, Any], *keys: str) -> str:
        for k in keys:
            if d.get(k) is not None:
                v = d[k]
                if isinstance(v, float):
                    return f"{v:.4f}"
                return str(v)
        return "—"

    rows = [
        ("run_ok", "run_ok", "run_ok"),
        ("elapsed_s", "elapsed_s", "elapsed_s"),
        ("pose CF", "pose_CF", "pose_CF"),
        ("pose CF.com", "pose_CF.com", "pose_CF.com"),
        ("pose CF.wal", "pose_CF.wal", "pose_CF.wal"),
        ("INI CF (native seed)", "ini_CF", "ini_CF"),
        ("gap pose−INI", "gap_pose_minus_ini", "gap_pose_minus_ini"),
        ("decoy_beats_native", "decoy_beats_native", "decoy_beats_native"),
    ]
    for label, ck, cu in rows:
        lines.append(
            f"| {label} | {cell(classic, ck)} | {cell(current, cu)} |"
        )

    # result.csv RMSD if present
    for arm_name, d in (("classic", classic), ("current", current)):
        rc = d.get("result_csv") or {}
        if rc:
            lines.append("")
            lines.append(f"## {arm_name} result.csv")
            for k, v in rc.items():
                lines.append(f"- **{k}**: {v}")

    lines.append("")
    lines.append("## Config deltas (classic ← FlexAID-like steric stack)")
    lines.append("- soft_wall_cutoff: 0 (classic hard) vs 0.4 (current soft)")
    lines.append("- permeability: 1.0 vs 0.9")
    lines.append("- intramolecular: false vs true")
    lines.append("- normalize_area: true on both (false explodes CF scale; do not use)")
    lines.append("- hbond scoring off vs on (rank still off)")
    lines.append("- fitness_model: CF vs SMFREE; thermo_engine off vs on")
    lines.append("- Equal GA budget in both JSON configs (see configs/smoke_classic_vs_current/)")
    lines.append("")
    lines.append("## Logs")
    lines.append(f"- classic: `{classic.get('log')}`")
    lines.append(f"- current: `{current.get('log')}`")
    return "\n".join(lines) + "\n"


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--binary",
        default=os.environ.get("FLEXAIDDS_BIN", ""),
        help="Path to FlexAIDdS executable",
    )
    ap.add_argument("--pdb-id", default="1GPK")
    ap.add_argument(
        "--out-dir",
        default="results/smoke_classic_vs_current",
        help="Output directory under repo root",
    )
    ap.add_argument("--timeout", type=int, default=1800, help="Per-arm timeout seconds")
    ap.add_argument(
        "--skip-run",
        action="store_true",
        help="Only summarize existing outputs (no dock)",
    )
    args = ap.parse_args(argv)
    root = repo_root()
    binary = Path(args.binary) if args.binary else None
    if not args.skip_run:
        if not binary or not binary.is_file():
            # try common locations without resolve_build pin
            for cand in (
                root / "build" / "FlexAIDdS",
                root / "build_lto" / "FlexAIDdS",
                Path("/Users/lp.more/Projects/FlexAIDdS/build_lto/FlexAIDdS"),
                Path("/Users/lp.more/Projects/FlexAIDdS/build/FlexAIDdS"),
            ):
                if cand.is_file() and os.access(cand, os.X_OK):
                    binary = cand
                    break
        if not binary or not binary.is_file():
            print("ERROR: FlexAIDdS binary not found. Pass --binary.", file=sys.stderr)
            return 2

    pdb_id = args.pdb_id
    struct = root / "benchmarks" / "astex_diverse" / "astex_diverse" / pdb_id
    receptor = struct / f"{pdb_id}.pdb"
    ligand = struct / f"{pdb_id}_ligand.sdf"
    if not receptor.is_file() or not ligand.is_file():
        print(f"ERROR: missing structures under {struct}", file=sys.stderr)
        return 2

    cfg_dir = root / "configs" / "smoke_classic_vs_current"
    classic_cfg = cfg_dir / f"{pdb_id}_classic.json"
    current_cfg = cfg_dir / f"{pdb_id}_current.json"
    if not classic_cfg.is_file() or not current_cfg.is_file():
        print(f"ERROR: configs missing in {cfg_dir}", file=sys.stderr)
        return 2

    out_root = root / args.out_dir
    classic_dir = out_root / f"{pdb_id}_classic"
    current_dir = out_root / f"{pdb_id}_current"
    classic_prefix = classic_dir / pdb_id
    current_prefix = current_dir / pdb_id

    metas: Dict[str, Dict[str, Any]] = {}
    if not args.skip_run:
        assert binary is not None
        print(f"Binary: {binary}")
        print(f"Target: {pdb_id}")
        # Sequential arms to avoid disk/CPU thrash on laptop; still fair budget
        for name, cfg, pref in (
            ("classic", classic_cfg, classic_prefix),
            ("current", current_cfg, current_prefix),
        ):
            print(f"=== Running {name} arm ===")
            metas[name] = run_arm(
                binary=binary,
                receptor=receptor,
                ligand=ligand,
                config=cfg,
                out_prefix=pref,
                timeout_s=args.timeout,
            )
            print(f"  ok={metas[name]['ok']} elapsed={metas[name]['elapsed_s']:.1f}s")
    else:
        metas = {
            "classic": {"ok": None, "elapsed_s": None, "returncode": None, "log": None},
            "current": {"ok": None, "elapsed_s": None, "returncode": None, "log": None},
        }

    classic_sum = summarize_arm("classic", classic_dir, classic_prefix, metas["classic"])
    current_sum = summarize_arm("current", current_dir, current_prefix, metas["current"])

    report = side_by_side(classic_sum, current_sum)
    report_path = out_root / f"{pdb_id}_classic_vs_current.md"
    json_path = out_root / f"{pdb_id}_classic_vs_current.json"
    report_path.write_text(report)
    json_path.write_text(json.dumps({"classic": classic_sum, "current": current_sum}, indent=2) + "\n")
    print(report)
    print(f"Wrote {report_path}")
    print(f"Wrote {json_path}")
    return 0 if (args.skip_run or (metas["classic"].get("ok") and metas["current"].get("ok"))) else 1


if __name__ == "__main__":
    sys.exit(main())
