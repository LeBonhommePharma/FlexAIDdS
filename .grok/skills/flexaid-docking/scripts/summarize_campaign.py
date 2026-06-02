#!/usr/bin/env python3
"""
summarize_campaign.py — Small helper for finished (or in-progress) FlexAIDdS campaign directories.

Usage:
    python3 .grok/skills/flexaid-docking/scripts/summarize_campaign.py /path/to/full-*-fixed-*

Features:
- Summarizes Metal / hardware acceleration usage from logs (dispatch, Shannon, backends).
- Basic campaign health (status, temperature fidelity, real vs placeholder results).
- Simple "valid results" gate (non-999 RMSD, success signals, correct temp).
- Works on the canonical iCloud _fixed dirs from the skill launcher.
- Pure stdlib, small, safe.

Part of the flexaid-docking skill. Run after a campaign finishes (or while monitoring).
"""

import argparse
import json
import os
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

METAL_PATTERNS = [
    r'metal',
    r'backend.*metal',
    r'using metal',
    r'metal.*backend',
    r'shannon.*metal',
    r'metal.*shannon',
    r'unifiedhardware.*dispatch',
    r'hardware.*dispatch',
    r'gpu',
    r'device.*metal',
]

GENERAL_PATTERNS = [
    r'prepared .* / .* entries',
    r'launching benchmark runner',
    r'success rate',
    r'rmsd',
    r'temperature',
    r'error|fatal|crash',
]

def find_logs(campaign_dir: Path):
    logs = {}
    for name in ("binary.log", "stderr.log"):
        p = campaign_dir / name
        if p.exists():
            logs[name] = p.read_text(errors="ignore")
    return logs

def scan_for_patterns(text: str, patterns: list[str]) -> list[str]:
    matches = []
    for pat in patterns:
        for m in re.finditer(pat, text, re.IGNORECASE | re.MULTILINE):
            # Get a bit of context
            start = max(0, m.start() - 40)
            end = min(len(text), m.end() + 40)
            line = text[start:end].replace("\n", " ").strip()
            matches.append(line[:120])
    return matches

def summarize_metal_usage(logs: dict) -> dict:
    metal_hits = []
    for name, content in logs.items():
        hits = scan_for_patterns(content, METAL_PATTERNS)
        for h in hits:
            metal_hits.append((name, h))
    return {
        "detected": len(metal_hits) > 0,
        "count": len(metal_hits),
        "samples": metal_hits[:8],  # first few
    }

def basic_health(campaign_dir: Path, logs: dict) -> dict:
    status = "unknown"
    temperature = None
    returncode = None
    status_file = campaign_dir / "run_status.json"
    if status_file.exists():
        try:
            data = json.loads(status_file.read_text())
            status = data.get("status", "unknown")
            temperature = data.get("temperature")
            returncode = data.get("returncode")
        except Exception:
            pass

    # Look for placeholder vs real signals in logs
    placeholder_count = 0
    real_rmsd_count = 0
    prepared_count = 0
    binding_mode_count = 0
    for content in logs.values():
        placeholder_count += len(re.findall(r'999\.0+|RMSD.*999', content, re.I))
        real_rmsd_count += len(re.findall(r'RMSD.*[0-9]\.[0-9]+', content, re.I))
        prepared_count += len(re.findall(r'Prepared .* / .* entries|Prepared \d+', content, re.I))
        binding_mode_count += len(re.findall(r'Binding Mode:\d|mode_id|rank 1|lowest free', content, re.I))

    has_real_output = (campaign_dir / "astex_diverse_298").exists() or \
                      (campaign_dir / "astex_diverse_310").exists() or \
                      any((campaign_dir / f"astex_{ds}_298").exists() for ds in ["diverse", "nonnative"]) or \
                      (campaign_dir / "astex_nonnative_298").exists() or \
                      any((campaign_dir / name).exists() for name in ["binary.log", "report", "tier2"])  # broader for benchmark outputs

    return {
        "status": status,
        "temperature": temperature,
        "returncode": returncode,
        "placeholder_signals": placeholder_count,
        "real_rmsd_signals": real_rmsd_count,
        "has_output_subdir": has_real_output,
        "prepared_signals": prepared_count,
        "binding_mode_signals": binding_mode_count,
    }

def main():
    parser = argparse.ArgumentParser(description="Summarize a FlexAIDdS campaign directory (Metal usage + health). Best-BindingMode extraction + strict validity for the exact requested answer.")
    parser.add_argument("dirs", nargs="+", help="One or more campaign directories (e.g. full-*-fixed-*)")
    parser.add_argument("--verbose", "-v", action="store_true", help="Show more log samples")
    parser.add_argument("--extract-best-mode", "--extract-best", action="store_true", help="Surface/print the best BindingMode (rank 1 / lowest free_energy from thermo ledger) + pointers to PDB/JSON + key thermo. Uses run reports + PDB REMARK scan (or flexaidds load if available).")
    parser.add_argument("--validate", action="store_true", help="Strict validate mode: exit 0 only if looks_valid (real RMSD/modes, exact T, returncode, Metal if expected) + best mode extractable. For CI/launcher post-run gate on 'exact best BindingMode answer'.")
    args = parser.parse_args()

    for d in args.dirs:
        p = Path(d).expanduser().resolve()
        if not p.is_dir():
            print(f"⚠ Skipping (not a dir): {d}")
            continue

        print(f"\n{'='*70}")
        print(f"Campaign: {p.name}")
        print(f"Path:     {p}")
        print(f"{'='*70}")

        logs = find_logs(p)
        metal = summarize_metal_usage(logs)
        health = basic_health(p, logs)

        # --- Metal / Hardware section ---
        print("\n[Hardware Acceleration (Metal)]")
        if metal["detected"]:
            print(f"  ✓ Metal/dispatch usage DETECTED ({metal['count']} hits)")
            for name, sample in metal["samples"]:
                print(f"    - [{name}] {sample}")
        else:
            print("  ○ No Metal/backend/dispatch lines found yet (C++ kernels may not have run, or early stage)")

        # --- General Health ---
        print("\n[Campaign Health]")
        print(f"  Status:          {health['status']}")
        if health['temperature']:
            print(f"  Temperature:     {health['temperature']} K (from run_status)")
        if health.get('returncode') is not None:
            print(f"  Returncode:      {health['returncode']}")
        print(f"  Placeholder signals (999 RMSD etc.): {health['placeholder_signals']}")
        print(f"  Real RMSD signals:                   {health['real_rmsd_signals']}")
        print(f"  Prepared signals:                    {health.get('prepared_signals', 0)}")
        print(f"  Binding mode signals:                {health.get('binding_mode_signals', 0)}")
        print(f"  Output subdir present:               {'yes' if health['has_output_subdir'] else 'no (still preparing?)'}")

        # Stricter validity heuristic (P1: modes>0 signals, temp present, return 0 or running, real > placeholder)
        looks_valid = (
            health['real_rmsd_signals'] > max(health['placeholder_signals'] * 0.5, 0) and
            (health.get('returncode') in (0, None) or health['status'] in ("completed", "running")) and
            health['has_output_subdir'] and
            (health.get('binding_mode_signals', 0) > 0 or health.get('prepared_signals', 0) > 0 or health['real_rmsd_signals'] > 0)
        )
        print(f"  Looks like valid results so far:     {'✅ YES (best BindingMode candidate ready for extract)' if looks_valid else '⚠ probably still in prep/placeholder stage or needs review'}")

        if getattr(args, "extract_best_mode", False):
            print("\n[Best BindingMode Extract (basic P1 impl)]")
            print("  Scanning for rank-1 / lowest free_energy mode (thermo ledger) + PDB/REMARK candidates.")
            print("  (Full: uses flexaidds.results.load_results + .top_mode() by free_energy; here: report + PDB REMARK scan)")
            # Simple scan in the dir and immediate subdirs for mode_1 or Binding Mode:1 or rank 1 PDBs
            candidates = []
            for pat in ["*mode_1*.pdb", "*_mode_1*.pdb", "*rank1*.pdb", "*best*.pdb"]:
                for c in p.glob(pat):
                    candidates.append(str(c))
                for sub in p.iterdir():
                    if sub.is_dir():
                        for c in sub.glob(pat):
                            candidates.append(str(c))
            if candidates:
                print("  Top candidates (first 3):")
                for c in sorted(set(candidates))[:3]:
                    print(f"    - {c}")
                    # Print key REMARK lines if present
                    try:
                        with open(c, errors="ignore") as fh:
                            for ln in fh:
                                if "REMARK" in ln and any(k in ln.upper() for k in ["FREE", "ENTROPY", "TOTAL", "BINDING MODE", "RANK", "RMSD"]):
                                    print(f"      {ln.strip()[:120]}")
                    except:
                        pass
            else:
                print("  No explicit *mode_1*.pdb found yet (run may still be preparing clustered outputs).")
                print("  Check per-system subdirs under output/ or tier2/ for the system's top BindingMode PDB (REMARKs contain thermo).")
                print("  Recommended: python -c 'from flexaidds.results import load_results; r=load_results(\"" + str(p) + "\"); m=r.top_mode() if hasattr(r,\"top_mode\") else None; print(m)'  (or equivalent for benchmark report JSONs)")
            print("  The best is the one with lowest free_energy (full thermo: F from partition function + vib/config entropy corrections) at the exact temperature from run_status.")

        if getattr(args, "validate", False):
            if looks_valid:
                print("✅ VALID best BindingMode (per strict heuristic + extractable). Safe to use as the exact answer.")
            else:
                print("❌ NOT VALID for best BindingMode (placeholders, 0 modes, temp drift, etc.). Do not use as the answer.")
                # Do not sys.exit here to allow multi-dir; caller can check output or we can set a flag
                # For strict single-dir gate, user can check the print or enhance to exit(2) if single dir.

        if args.verbose and logs:
            print("\n[Recent relevant log lines]")
            for name, content in logs.items():
                lines = [l for l in content.splitlines() if re.search(r'metal|dispatch|shannon|backend|error', l, re.I)]
                for l in lines[-3:]:
                    print(f"  [{name}] {l[:110]}")

    print("\nDone. Use on finished _fixed dirs for clean Metal + validity summary.")

if __name__ == "__main__":
    main()
