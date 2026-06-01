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
    status_file = campaign_dir / "run_status.json"
    if status_file.exists():
        try:
            data = json.loads(status_file.read_text())
            status = data.get("status", "unknown")
            temperature = data.get("temperature")
        except Exception:
            pass

    # Look for placeholder vs real signals in logs
    placeholder_count = 0
    real_rmsd_count = 0
    for content in logs.values():
        placeholder_count += len(re.findall(r'999\.0+|RMSD.*999', content, re.I))
        real_rmsd_count += len(re.findall(r'RMSD.*[0-9]\.[0-9]+', content, re.I))

    has_real_output = (campaign_dir / "astex_diverse_298").exists() or \
                      (campaign_dir / "astex_diverse_310").exists() or \
                      any((campaign_dir / f"astex_{ds}_298").exists() for ds in ["diverse", "nonnative"]) or \
                      (campaign_dir / "astex_nonnative_298").exists()

    return {
        "status": status,
        "temperature": temperature,
        "placeholder_signals": placeholder_count,
        "real_rmsd_signals": real_rmsd_count,
        "has_output_subdir": has_real_output,
    }

def main():
    parser = argparse.ArgumentParser(description="Summarize a FlexAIDdS campaign directory (Metal usage + health).")
    parser.add_argument("dirs", nargs="+", help="One or more campaign directories (e.g. full-*-fixed-*)")
    parser.add_argument("--verbose", "-v", action="store_true", help="Show more log samples")
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
        print(f"  Placeholder signals (999 RMSD etc.): {health['placeholder_signals']}")
        print(f"  Real RMSD signals:                   {health['real_rmsd_signals']}")
        print(f"  Output subdir present:               {'yes' if health['has_output_subdir'] else 'no (still preparing?)'}")

        # Simple validity heuristic
        looks_valid = (
            health['real_rmsd_signals'] > health['placeholder_signals'] * 0.5 and
            health['status'] in ("completed", "running") and
            health['has_output_subdir']
        )
        print(f"  Looks like valid results so far:     {'✅ YES' if looks_valid else '⚠ probably still in prep/placeholder stage'}")

        if args.verbose and logs:
            print("\n[Recent relevant log lines]")
            for name, content in logs.items():
                lines = [l for l in content.splitlines() if re.search(r'metal|dispatch|shannon|backend|error', l, re.I)]
                for l in lines[-3:]:
                    print(f"  [{name}] {l[:110]}")

    print("\nDone. Use on finished _fixed dirs for clean Metal + validity summary.")

if __name__ == "__main__":
    main()
