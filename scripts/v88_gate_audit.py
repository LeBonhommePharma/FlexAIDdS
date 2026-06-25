#!/usr/bin/env python3
"""Pure parser for v88 gate audit. No hacks, no synthetic."""
import re
import json
import csv
import sys
import os
from pathlib import Path

def parse_v88_dict(repro_path: Path):
    txt = repro_path.read_text()
    m = re.search(r"published = \{(.+?)\n\}", txt, re.DOTALL)
    if not m:
        raise ValueError("no published dict")
    d = {}
    for line in m.group(1).split(","):
        mm = re.search(r'"(\w+)":([0-9.]+)', line)
        if mm:
            d[mm.group(1)] = float(mm.group(2))
    n = len(d)
    loose = sum(1 for v in d.values() if v < 2.0)
    strict = sum(1 for v in d.values() if 0 < v < 2.0)
    zeros = sum(1 for v in d.values() if v == 0.0)
    return {"n": n, "loose_lt2": loose, "strict_0_lt2": strict, "seed_echo": zeros}

def parse_current_csv(csv_path: Path):
    if not csv_path.exists():
        return {"n": 0, "succ": 0, "rate": 0.0, "note": "no csv"}
    rows = list(csv.DictReader(csv_path.open()))
    n = len(rows)
    succ = sum(1 for r in rows if 0 < float(r.get("rmsd_hungarian") or 99) < 2.0)
    rate = round(100 * succ / n, 1) if n else 0
    return {"n": n, "succ": succ, "rate": rate}

def main():
    repro = Path("REPRODUCIBILITY.md")
    scr = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/var/folders/8b/tgtvwb_j6zd_g03vl1w4ykfw0000gn/T/grok-goal-86a3d1efec00/implementer")
    v88 = parse_v88_dict(repro)
    curr_csv = scr / "astex_crossdock_85_results.csv"
    curr = parse_current_csv(curr_csv)
    audit = {
        "v88_dict": v88,
        "historical_claim": {"rate": 91.4, "n": 78, "note": "v88 with NATIVE_SEED_FRAC=0.90 (seed-echo)"},
        "current_strict": curr,
        "gate": "0 < rmsd_hungarian < 2.0",
        "strict_note": "78/85 not reproduced under strict gate; seed-echo zeros (33) explain historical claim vs dict strict 46"
    }
    out = scr / "gate_audit.json"
    out.write_text(json.dumps(audit, indent=2))
    print(json.dumps(audit, indent=2))

if __name__ == "__main__":
    main()
