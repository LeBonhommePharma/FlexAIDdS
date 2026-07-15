#!/usr/bin/env python3
"""
compare_astex_2015_fair_vs_full_ds.py

Clean side-by-side comparison of:
  - 2015-fair: original FlexAID-like (no entropy/thermo, minimal modern features, oracle sites)
  - full-dS: modern FlexAIDdS with all enhancements (THERMO, multiple restarts, elitism, etc.)

Defaults output directories to iCloud Drive when possible (via FLEXAIDDS_ICLOUD
or standard ~/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks/).
Uses timestamped paths under iCloud/results/ (or local results/) .

After successful --run-both (or analysis paths), automatically invokes
scripts/safe_archive_to_icoud.py on the result directories (with --keep-local)
to perform safer-than-safe copy (rsync + SHA verify + atomic markers + no-evict).
This respects iCloud File Provider best practices; no direct risky writes for
the archive step (staging+verification delegated to safe archiver).

Usage:
  # Dry-run (validate commands, no long docking) — uses iCloud defaults if available
  python3 scripts/compare_astex_2015_fair_vs_full_ds.py --dry-run

  # Run both configurations (will take time); auto-archives on success
  python3 scripts/compare_astex_2015_fair_vs_full_ds.py --run-both --workers 4

  # Force iCloud for defaults (or auto-detected when FLEXAIDDS_ICLOUD is set)
  python3 scripts/compare_astex_2015_fair_vs_full_ds.py --icloud --run-both --dry-run

  # Analyze existing result dirs (may be under iCloud)
  python3 scripts/compare_astex_2015_fair_vs_full_ds.py \
      --fair-dir "$FLEXAIDDS_ICLOUD/results/astex_2015_fair_20260708_123456" \
      --full-dir "$FLEXAIDDS_ICLOUD/results/astex_full_ds_20260708_123456"

  # Or provide local and let analysis trigger archive of copies
  python3 scripts/compare_astex_2015_fair_vs_full_ds.py \
      --fair-dir results/astex_jcim2015_fair_XXX \
      --full-dir results/astex_full_ds_YYY

After runs, it invokes summarize_astex_single_ga.py (or equivalent parsing)
and prints a nice comparison table. Full git SHA, binary SHA, env, and
provenance captured via capture_provenance() into comparison_provenance.json .

Follows AGENTS.md: reproducibility, clean env, zero assumptions, verify runs.
Source of truth for iCloud logic: scripts/safe_archive_to_icoud.py:get_icoud_base()
"""

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Any, Optional
import statistics

REPO_ROOT = Path(__file__).resolve().parents[1]
SUMMARIZER = REPO_ROOT / "summarize_astex_single_ga.py"
DATASET_RUNNER = REPO_ROOT / ".grok/skills/flexaidds/scripts/dataset_runner.py"
MANIFEST = "crossdock_json:benchmarks/datasets/benchmark_astex_native_85.json"

# 2015-fair config: disable everything modern
FAIR_ENV = {
    "FLEXAIDDS_THERMO": "0",
    "FLEXAIDDS_SEED_ELITISM": "0",
    "FLEXAIDDS_RESTARTS": "1",
    "FLEXAIDDS_PARALLEL_RESTARTS": "0",
    "FLEXAIDDS_CONSENSUS_SCORER": "0",
    "FLEXAIDDS_RECEPTOR_ROTAMER_PREP": "0",
    "FLEXAIDDS_NATIVE_SEED_FRAC": "0",
    "FLEXAIDDS_SOFTCORE_WAL": "0",
    "FLEXAIDDS_SOFTCORE_FLOOR": "0",
    "FLEXAIDDS_T_EFF": "1.0",
    "FLEXAIDDS_TENCOM_SCALE": "0",
}

# Full modern dS (from REPRODUCIBILITY.md / reproduce_astex85.sh)
FULL_ENV = {
    "FLEXAIDDS_THERMO": "1",
    "FLEXAIDDS_T_EFF": "0.596",
    "FLEXAIDDS_TENCOM_SCALE": "1.0",
    "FLEXAIDDS_RESTARTS": "7",
    "FLEXAIDDS_PARALLEL_RESTARTS": "1",
    "FLEXAIDDS_EVAL_SCALE_DIHEDRAL": "1",
    "FLEXAIDDS_CONSENSUS_SCORER": "1",
    "FLEXAIDDS_SEED_ELITISM": "1",
    "FLEXAIDDS_N_ELITE": "1",
    "FLEXAIDDS_BUDGET_SCALE": "1",
    "FLEXAIDDS_SOFTCORE_WAL": "1",
    "FLEXAIDDS_SOFTCORE_FLOOR": "0.5",
    "FLEXAIDDS_T_HOT": "500",
    "FLEXAIDDS_NATIVE_SEED_FRAC": "0.90",
    "FLEXAIDDS_RECEPTOR_ROTAMER_PREP": "1",
}


# --- iCloud support (matches safer-than-safe contract in safe_archive_to_icoud.py) ---
# Do NOT bake machine-specific /Users/... paths. Resolve at runtime via env or Path.home().
ICLOUD_STANDARD = "Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks"

def get_icoud_base() -> Path:
    """Return iCloud base dir for FlexAIDdS benchmarks.

    Honors FLEXAIDDS_ICLOUD (if set, appends /FlexAIDdS_benchmarks if needed).
    Falls back to standard macOS iCloud container + /FlexAIDdS_benchmarks .
    """
    env = os.environ.get("FLEXAIDDS_ICLOUD")
    if env:
        p = Path(env).expanduser()
        if p.name != "FlexAIDdS_benchmarks":
            p = p / "FlexAIDdS_benchmarks"
        return p
    return Path.home() / ICLOUD_STANDARD


def should_use_icloud(force_flag: bool) -> bool:
    """Auto-detect or respect --icloud.

    Prefers explicit env, then flag, then on darwin presence of iCloud container.
    """
    if force_flag:
        return True
    if os.environ.get("FLEXAIDDS_ICLOUD"):
        return True
    if sys.platform == "darwin":
        std = Path.home() / ICLOUD_STANDARD
        # Heuristic: container parent or the dir itself signals iCloud Drive is active
        container_parent = std.parent.parent  # .../com~apple~CloudDocs
        if std.exists() or container_parent.exists():
            return True
    return False


def write_json_atomic(path: Path, obj: dict) -> None:
    """Atomic write (safer for iCloud File Provider)."""
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True), encoding="utf-8")
    tmp.replace(path)


def resolve_default_dir(kind: str, use_icloud: bool, ts: str) -> Path:
    """Return timestamped default for fair/full dirs under iCloud/results/ or local results/."""
    if use_icloud:
        base = get_icoud_base() / "results"
    else:
        base = REPO_ROOT / "results"
    if kind == "fair":
        return base / f"astex_2015_fair_{ts}"
    else:
        return base / f"astex_full_ds_{ts}"


# Reproducibility helpers (preserved and now integrated; were previously dead code after __main__)
def capture_provenance(results_dir: Path, config_name: str, env: Dict[str, str], cmd: list, binary: str):
    """Capture full reproducibility info for objectivity. Writes comparison_provenance.json .

    Includes git SHA, clean status, platform, FLEXAIDDS_* env, binary sha256, full cmd.
    Uses atomic write for iCloud safety.
    """
    prov = {
        "config": config_name,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "git_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True).strip(),
        "git_status_clean": subprocess.call(["git", "diff", "--quiet", "HEAD"], cwd=REPO_ROOT) == 0,
        "platform": os.uname()._asdict() if hasattr(os, "uname") else str(os.name),
        "binary": binary,
        "binary_sha256": None,
        "env_vars": {k: v for k, v in env.items() if k.startswith("FLEXAIDDS_")},
        "full_command": " ".join(map(str, cmd)),
        "cwd": str(REPO_ROOT),
    }
    try:
        if binary and Path(binary).exists():
            import hashlib
            prov["binary_sha256"] = hashlib.sha256(Path(binary).read_bytes()).hexdigest()
    except Exception:
        pass
    out_path = results_dir / "comparison_provenance.json"
    write_json_atomic(out_path, prov)
    return prov


def print_comparison_table(fair: Dict, full: Dict):
    print("\n# Objective Side-by-Side: 2015-fair vs full-dS Astex Diverse")
    print("## Metrics (using summarize_astex_single_ga.py logic: RMSD_hungarian < 2.0 success)")
    print("| Metric              | 2015-fair          | full-dS            | Delta          |")
    print("|---------------------|--------------------|--------------------|----------------|")
    for k in ["success_rate", "success_count", "mean_rmsd", "median_rmsd"]:
        fv = fair.get(k, "N/A")
        uv = full.get(k, "N/A")
        if isinstance(fv, float) and isinstance(uv, float):
            d = uv - fv
            ds = f"{d:+.4f}"
        else:
            ds = "N/A"
        print(f"| {k:<19} | {str(fv):<18} | {str(uv):<18} | {ds:<14} |")
    print("\n## Reproducibility notes")
    print("- Both runs should use identical input manifest (oracle binding sites).")
    print("- full-dS enables entropy/thermo/restarts/elitism as in REPRODUCIBILITY.md (94.1% claimed).")
    print("- 2015-fair disables them for closer match to original JCIM 2015 FlexAID (no entropy).")
    print("- Run with same binary, same data checkout, capture git SHA + binary SHA.")
    print("- For maximum objectivity, also run the exact reproduce_astex85.sh for full-dS reference.")


def run_cmd(cmd: list, env: Optional[Dict[str, str]] = None, cwd: Optional[Path] = None, capture: bool = False) -> subprocess.CompletedProcess:
    full_env = os.environ.copy()
    if env:
        full_env.update(env)
    print(f"[RUN] {' '.join(map(str, cmd))}")
    return subprocess.run(cmd, env=full_env, cwd=cwd or REPO_ROOT, capture_output=capture, text=True)

def launch_config(name: str, env: Dict[str, str], results_dir: Path, workers: int = 4, binary: Optional[str] = None, dry_run: bool = False):
    """Launch dataset runner for one config. Returns (proc, cmd) so caller can capture provenance."""
    results_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, str(DATASET_RUNNER),
        "--dataset", "astex_diverse",
        "--tier", "2",
        "--results-dir", str(results_dir),
        "--workers", str(workers),
    ]
    if binary:
        cmd += ["--binary", binary]
    if dry_run:
        cmd.append("--dry-run")
    else:
        cmd.append("--resume")

    proc = run_cmd(cmd, env=env)
    return proc, cmd

def analyze_dir(res_dir: Path) -> Dict[str, Any]:
    """Use the project's summarizer or parse results."""
    if not res_dir.exists():
        return {"error": "dir not found"}

    # Try the standard summarizer first
    if SUMMARIZER.exists():
        try:
            out = subprocess.check_output([sys.executable, str(SUMMARIZER), str(res_dir)], text=True, cwd=REPO_ROOT)
            # Parse key lines
            metrics = {}
            for line in out.splitlines():
                if "successes=" in line:
                    # e.g. successes=38/85
                    parts = line.split()
                    for p in parts:
                        if "successes=" in p:
                            val = p.split("=")[1]
                            s, tot = val.split("/")
                            metrics["success_count"] = int(s)
                            metrics["total"] = int(tot)
                            metrics["success_rate"] = int(s) / int(tot)
                if "mean_rmsd=" in line:
                    metrics["mean_rmsd"] = float(line.split("=")[1])
                if "median_rmsd=" in line:
                    metrics["median_rmsd"] = float(line.split("=")[1])
            if metrics:
                return metrics
        except Exception as e:
            print(f"[WARN] summarizer failed: {e}")

    # Fallback: look for common CSV names
    candidates = list(res_dir.glob("*astex*results*.csv")) + list(res_dir.glob("**/astex_diverse_results.csv"))
    for csv in candidates:
        try:
            import csv as pycsv
            rows = list(pycsv.DictReader(open(csv)))
            valid = [r for r in rows if 0.0 <= float(r.get("rmsd_hungarian", -1)) < 900]
            successes = [r for r in valid if float(r.get("rmsd_hungarian", 99)) < 2.0]
            rmsds = [float(r["rmsd_hungarian"]) for r in valid]
            return {
                "success_count": len(successes),
                "total": len(rows),
                "success_rate": len(successes) / len(rows) if rows else 0,
                "mean_rmsd": statistics.mean(rmsds) if rmsds else 0,
                "median_rmsd": statistics.median(rmsds) if rmsds else 0,
                "source": str(csv),
            }
        except Exception:
            pass

    # Look for per-target JSONs or other
    jsons = list(res_dir.glob("**/*_holo.json")) + list(res_dir.glob("**/*result*.json"))
    if jsons:
        # simplistic
        return {"note": f"found {len(jsons)} json results, manual inspection recommended"}

    return {"error": "no parseable results found"}

def main():
    parser = argparse.ArgumentParser(description="Side-by-side 2015-fair vs full-dS Astex Diverse comparison")
    parser.add_argument("--run-both", action="store_true", help="Launch both configurations")
    parser.add_argument("--dry-run", action="store_true", help="Dry run only (no actual docking)")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--binary", help="Path to benchmark_datasets binary")
    parser.add_argument("--fair-dir", help="Existing/target fair (2015-like) results dir. Defaults to timestamped path under iCloud/results/ or local results/ when omitted.")
    parser.add_argument("--full-dir", help="Existing/target full-dS results dir. Defaults to timestamped path under iCloud/results/ or local results/ when omitted.")
    parser.add_argument("--no-launch", action="store_true", help="Only analyze, do not launch")
    parser.add_argument("--icloud", action="store_true", help="Force default --fair-dir/--full-dir under iCloud (auto-detected if FLEXAIDDS_ICLOUD set or iCloud container visible)")
    args = parser.parse_args()

    # Resolve iCloud preference and timestamped defaults (only when not explicitly provided)
    use_icloud = should_use_icloud(args.icloud)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")

    if args.fair_dir:
        fair_dir = Path(args.fair_dir)
    else:
        fair_dir = resolve_default_dir("fair", use_icloud, ts)

    if args.full_dir:
        full_dir = Path(args.full_dir)
    else:
        full_dir = resolve_default_dir("full", use_icloud, ts)

    if use_icloud:
        print(f"[iCloud] Using iCloud base for defaults: {get_icoud_base()}")
    print(f"[defaults] fair_dir={fair_dir}")
    print(f"[defaults] full_dir={full_dir}")

    launched = False
    if args.run_both and not args.no_launch:
        print("=== Launching 2015-fair ===")
        proc_f, cmd_f = launch_config("2015-fair", FAIR_ENV, fair_dir, args.workers, args.binary, args.dry_run)

        print("\n=== Launching full-dS ===")
        proc_u, cmd_u = launch_config("full-dS", FULL_ENV, full_dir, args.workers, args.binary, args.dry_run)
        launched = True

        # Always capture provenance (preserves git SHA, binary SHA, env, command, etc.)
        try:
            capture_provenance(fair_dir, "2015-fair", FAIR_ENV, cmd_f, args.binary or "")
            capture_provenance(full_dir, "full-dS", FULL_ENV, cmd_u, args.binary or "")
        except Exception as e:
            print(f"[WARN] provenance capture failed: {e}")

        if not args.dry_run:
            # Auto-archive to iCloud using safer-than-safe archiver (with --keep-local)
            # This delegates rsync+SHA-verify+atomic+.verified marker; never direct risky write here.
            print("\n=== Auto-archiving to iCloud (safer-than-safe, --keep-local) ===")
            archiver = REPO_ROOT / "scripts" / "safe_archive_to_icoud.py"
            icloud_arch_dest = get_icoud_base() / "archived"
            for label, d in [("2015-fair", fair_dir), ("full-dS", full_dir)]:
                if d.exists() and d.is_dir():
                    try:
                        arch_cmd = [
                            sys.executable, str(archiver),
                            "--source", str(d),
                            "--dest", str(icloud_arch_dest),
                            "--keep-local",
                        ]
                        run_cmd(arch_cmd)
                    except Exception as e:
                        print(f"[WARN] safe_archive step for {label} ({d}) had issue (data preserved locally): {e}")
            if proc_f.returncode == 0 and proc_u.returncode == 0:
                print("Launches + provenance + iCloud archive steps completed.")
            # Fall through to analysis below so user gets immediate table (results are ready post-subprocess)
            # (Previously returned early; now we analyze + report the exact dirs used.)

    # Analyze (always reached for --dry-run, --no-launch, or after --run-both)
    print("\n" + "="*70)
    print("SIDE-BY-SIDE COMPARISON: 2015-fair vs full-dS (Astex Diverse 85)")
    print("="*70)

    fair_metrics = analyze_dir(fair_dir)
    full_metrics = analyze_dir(full_dir)

    print(f"\n2015-fair dir : {fair_dir}")
    print(f"full-dS dir   : {full_dir}")

    # Also call the preserved print_comparison_table for the detailed reproducibility output
    try:
        print_comparison_table(fair_metrics, full_metrics)
    except Exception:
        pass

    print("\n| Metric                  | 2015-fair          | full-dS            | Delta (full - fair) |")
    print("|-------------------------|--------------------|--------------------|---------------------|")

    def fmt(val):
        if isinstance(val, float):
            if val < 1:
                return f"{val:.4f}"
            return f"{val:.2f}"
        return str(val)

    for key in ["success_rate", "success_count", "mean_rmsd", "median_rmsd"]:
        f = fair_metrics.get(key, "N/A")
        fu = full_metrics.get(key, "N/A")
        if isinstance(f, (int, float)) and isinstance(fu, (int, float)):
            delta = fu - f
            delta_str = f"{delta:+.4f}" if key in ("success_rate", "mean_rmsd", "median_rmsd") else f"{int(delta):+d}"
        else:
            delta_str = "N/A"
        print(f"| {key:23} | {fmt(f):18} | {fmt(fu):18} | {delta_str:19} |")

    print("\nReference (from repo REPRODUCIBILITY.md for full-dS enhanced):")
    print("  Published (with all dS): 80/85 = 94.1%, mean RMSD 0.81, median 0.33")
    print("  (This used 7 restarts, THERMO=1, seed elitism, etc.)")

    print("\nOriginal 2015 paper context (from repo docs):")
    print("  FlexAID (no entropy) baseline ~55-62% on Astex Diverse.")
    print("  Use --package on a real run for full provenance.")

    if "error" in fair_metrics or "error" in full_metrics:
        print("\n[INFO] Some metrics could not be auto-parsed. Run the summarizer manually:")
        print(f"  python3 {SUMMARIZER} <dir>")

    # Auto-archive also on pure analysis paths (per requirements: after --run-both or analysis)
    # Only when not dry-run and the dirs contain plausible result artifacts.
    if not args.dry_run and (not args.run_both or args.no_launch or not launched):
        # Trigger only if it looks like we have data (avoid archiving empty defaults)
        if (fair_dir.exists() or full_dir.exists()):
            try:
                print("\n=== Analysis path: auto-archive of provided/analyzed dirs to iCloud (safer, --keep-local) ===")
                archiver = REPO_ROOT / "scripts" / "safe_archive_to_icoud.py"
                icloud_arch_dest = get_icoud_base() / "archived"
                for label, d in [("2015-fair", fair_dir), ("full-dS", full_dir)]:
                    if d.exists() and d.is_dir():
                        arch_cmd = [sys.executable, str(archiver), "--source", str(d), "--dest", str(icloud_arch_dest), "--keep-local"]
                        run_cmd(arch_cmd)
            except Exception as e:
                print(f"[WARN] analysis-triggered archive had issue: {e}")

if __name__ == "__main__":
    main()
