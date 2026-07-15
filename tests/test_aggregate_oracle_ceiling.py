#!/usr/bin/env python3
"""Drive scripts/aggregate_oracle_ceiling.py on real v43 campaign CSV."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "aggregate_oracle_ceiling.py"
V43 = Path.home() / "flexaidds_results" / "v43_20260613_softcore_natural"


def test_v43_bcr_ceiling_exceeds_90():
    if not V43.is_dir() or not (V43 / "astex_diverse_results.csv").is_file():
        print("SKIP: v43 campaign not present on this machine")
        return
    out = subprocess.check_output(
        [sys.executable, str(SCRIPT), str(V43)], text=True
    )
    report = json.loads(out)
    assert report["N"] == 85, report
    assert report["ceiling_n"] == 78, report
    assert report["ceiling_rate"] > 0.90, report
    assert report["exceeds_90_ceiling"] is True
    # top-1 historical lower than BCR (election gap)
    assert report["top1_success_rmsd_n"] == 69, report
    print(
        f"PASS v43 ceiling={report['ceiling_n']}/{report['N']} "
        f"rate={report['ceiling_rate']:.4f} top1={report['top1_success_rmsd_n']}"
    )


if __name__ == "__main__":
    test_v43_bcr_ceiling_exceeds_90()
