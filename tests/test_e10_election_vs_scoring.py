#!/usr/bin/env python3
"""Tests for scripts/e10_election_vs_scoring.py (offline E10 diagnostic)."""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "e10_election_vs_scoring.py"


def _load():
    spec = importlib.util.spec_from_file_location("e10_election_vs_scoring", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _write_target(tmp: Path, pdb: str, *, rmsd: float, bcr: float, cf: float, soft_g: float, freq: int):
    d = tmp / pdb
    d.mkdir(parents=True)
    fields = [
        "pdb_id",
        "rmsd_hungarian",
        "best_cluster_rmsd",
        "elected_cf",
        "seed_echo",
        "pose_source",
    ]
    with (d / "result.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerow(
            {
                "pdb_id": pdb,
                "rmsd_hungarian": f"{rmsd:.4f}",
                "best_cluster_rmsd": f"{bcr:.4f}",
                "elected_cf": f"{cf:.6f}",
                "seed_echo": "0",
                "pose_source": "ga_cluster",
            }
        )
    pdb_text = (
        f"REMARK CF={cf:.5f}\n"
        f"REMARK soft_beta_G = {soft_g:.6f}\n"
        f"REMARK frequency = {freq}\n"
        f"REMARK CF.com=-100.0\n"
        "ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n"
    )
    (d / f"{pdb}_0.pdb").write_text(pdb_text)
    # second head with better CF
    (d / f"{pdb}_1.pdb").write_text(
        f"REMARK CF={cf - 5.0:.5f}\n"
        "ATOM      1  C   LIG A   1       1.000   0.000   0.000  1.00  0.00           C\n"
    )
    return d


def test_election_gap_detected(tmp_path: Path):
    mod = _load()
    d = _write_target(
        tmp_path, "1G9V", rmsd=4.5, bcr=2.06, cf=-28.0, soft_g=-2800.0, freq=13226
    )
    r = mod.analyze_target(d)
    assert r.election_gap is True
    assert r.size_bias_suspect is True
    assert r.n_heads >= 2
    assert r.n_heads_cf_better_than_elected >= 1


def test_genuine_no_gap(tmp_path: Path):
    mod = _load()
    d = _write_target(
        tmp_path, "1HNN", rmsd=1.58, bcr=1.42, cf=-50.0, soft_g=-50.0, freq=5
    )
    r = mod.analyze_target(d)
    assert r.election_gap is False
    assert r.size_bias_suspect is False


def test_cli_writes_json(tmp_path: Path):
    mod = _load()
    d = _write_target(
        tmp_path, "1G9V", rmsd=10.0, bcr=2.3, cf=-31.0, soft_g=-2000.0, freq=100
    )
    out = tmp_path / "e10.json"
    rc = mod.main(["--target-dir", str(d), "--out-json", str(out)])
    assert rc == 0
    data = json.loads(out.read_text())
    assert data["n_targets"] == 1
    assert data["n_election_gap"] == 1
