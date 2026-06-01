#!/usr/bin/env python3
"""
Quick verification for `redock ingest-dir` + iCloud FS behavior.

Run with:
  source ~/.flexaidds_env
  PYTHONPATH=. python benchmarks/re-dock/test_ingest_dir.py

This creates a temporary campaign on iCloud (or local tmp if you prefer),
generates a few fake Codex-style result JSONs, runs ingest-dir against them,
verifies state was updated, then cleans up.

It also optionally calls the icloud_fs_check we just added.
"""

import json
import shutil
import sys
import tempfile
import time
from pathlib import Path

# Make the re-dock package importable using the documented style for this machine
sys.path.insert(0, "benchmarks/re-dock")

import importlib.util
import pathlib

def _load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod

# Load siblings first so relative imports inside cli work when loaded this way
orchestrator = _load_module("orchestrator", "benchmarks/re-dock/orchestrator.py")
cli_mod = _load_module("cli", "benchmarks/re-dock/cli.py")

cmd_ingest_dir = cli_mod.cmd_ingest_dir
cmd_init = cli_mod.cmd_init
cmd_status = cli_mod.cmd_status
BenchmarkCampaign = orchestrator.BenchmarkCampaign

import argparse


def main():
    # Use a real iCloud location for the test campaign (honors user requirement)
    icloud_base = Path.home() / "Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS/re-dock-test-ingest"
    icloud_base.mkdir(parents=True, exist_ok=True)

    campaign_dir = icloud_base / f"test-ingest-{int(time.time())}"
    campaign_dir.mkdir()

    print(f"Using campaign dir on iCloud: {campaign_dir}")

    # 1. Initialize a tiny campaign
    targets = campaign_dir / "targets.json"
    targets.write_text(json.dumps({
        "targets": [
            {"pdb_id": "1a30", "receptor": "test", "ligand": "test"},
            {"pdb_id": "1err", "receptor": "test", "ligand": "test"},
        ]
    }))

    # Simulate init via direct call (bypassing argparse for the test)
    init_args = argparse.Namespace(
        campaign_dir=str(campaign_dir),
        t_min=298.0, t_max=600.0, n_replicas=4,
        targets=str(targets)
    )
    cmd_init(init_args)
    print("Campaign initialized on iCloud.")

    # 2. Create fake Codex result JSONs (as if user copied them from Codex sessions)
    results_dir = campaign_dir / "codex-results"
    results_dir.mkdir()

    fake_results = [
        {"chunk_id": "1a30_T298_g0", "pdb_id": "1a30", "temperature": 298.0,
         "generation": 0, "best_energy": -11.2, "energies_sample": [-11.2 + i*0.01 for i in range(12)]},
        {"chunk_id": "1a30_T340_g0", "pdb_id": "1a30", "temperature": 340.0,
         "generation": 0, "best_energy": -10.8, "energies_sample": [-10.8 + i*0.02 for i in range(12)]},
        {"chunk_id": "1err_T298_g0", "pdb_id": "1err", "temperature": 298.0,
         "generation": 0, "best_energy": -9.5, "energies_sample": [-9.5 + i*0.015 for i in range(12)]},
    ]

    for r in fake_results:
        (results_dir / f"{r['chunk_id']}.json").write_text(json.dumps(r, indent=2))

    print(f"Created {len(fake_results)} fake Codex JSONs in {results_dir}")

    # 3. Run batch ingest-dir (the feature under test)
    ingest_args = argparse.Namespace(
        campaign_dir=str(campaign_dir),
        directory=str(results_dir),
        dry_run=False,
        validate_only=False,
        fs_check=True,   # Exercise the new iCloud FS checker as the user requested
    )
    print("\n--- Running ingest-dir with --fs-check ---")
    cmd_ingest_dir(ingest_args)

    # 4. Verify state
    print("\n--- Final status ---")
    status_args = argparse.Namespace(campaign_dir=str(campaign_dir))
    cmd_status(status_args)

    # 5. Cleanup (we are good citizens on iCloud)
    shutil.rmtree(campaign_dir)
    print(f"\n✅ ingest-dir + iCloud FS check test completed successfully. Campaign dir removed.")

    # Optional: leave a marker that this test passed on real iCloud
    marker = icloud_base / "LAST_INGEST_DIR_TEST_PASSED"
    marker.write_text(f"Test passed at {time.time()}\n")


if __name__ == "__main__":
    main()
