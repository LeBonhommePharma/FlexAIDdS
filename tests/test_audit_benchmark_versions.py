#!/usr/bin/env python3
"""Unit tests for scripts/audit_benchmark_versions.py."""

from __future__ import annotations

import csv
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "audit_benchmark_versions.py"


def _load():
    spec = importlib.util.spec_from_file_location("audit_benchmark_versions", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_row_is_genuine_rejects_seed_echo():
    mod = _load()
    row = {
        "seed_echo": "1",
        "pose_source": "ga_cluster",
        "rmsd_hungarian": "0.5",
        "rmsd_to_crystal": "0.5",
    }
    assert not mod.row_is_genuine(row)


def test_row_is_genuine_rejects_seed_pose_source():
    mod = _load()
    row = {
        "seed_echo": "0",
        "pose_source": "ini_elitism",
        "rmsd_hungarian": "0.4",
        "rmsd_to_crystal": "0.4",
    }
    assert not mod.row_is_genuine(row)


def test_row_is_genuine_accepts_ga_under_2a():
    mod = _load()
    row = {
        "seed_echo": "0",
        "pose_source": "ga_cluster",
        "rmsd_hungarian": "1.5",
        "rmsd_to_crystal": "1.8",
    }
    assert mod.row_is_genuine(row)


def test_row_s1_uses_crystal_rmsd():
    mod = _load()
    # ordered >2 but hungarian ≤2 → S1 false
    row = {
        "rmsd_to_crystal": "5.0",
        "rmsd_hungarian": "0.5",
    }
    assert not mod.row_s1(row)
    row2 = {"rmsd_to_crystal": "1.2", "rmsd_hungarian": "5.0"}
    assert mod.row_s1(row2)


def test_discover_skips_clouddocs_and_noise(tmp_path: Path):
    mod = _load()
    root = tmp_path / "results"
    root.mkdir()
    (root / "canary_foo").mkdir()
    (root / "C0_full85_bar").mkdir()
    (root / "v_comcap_softbeta_x").mkdir()
    (root / "workorders").mkdir()
    (root / "data_9dc9").mkdir()
    camps = root / "campaigns"
    camps.mkdir()
    (camps / "C0_full85_nested").mkdir()
    found = {p.name for p in mod.discover_campaigns([root])}
    assert "canary_foo" in found
    assert "C0_full85_bar" in found
    assert "v_comcap_softbeta_x" in found
    assert "C0_full85_nested" in found
    assert "workorders" not in found
    assert "data_9dc9" not in found


def test_audit_canary_from_summary_json(tmp_path: Path):
    mod = _load()
    camp = tmp_path / "canary_test_20260724"
    camp.mkdir()
    summary = {
        "n_genuine": 1,
        "n_total": 3,
        "rows": [
            {
                "pdb": "1AAA",
                "rmsd_hungarian": 1.1,
                "rmsd_to_crystal": 1.2,
                "seed_echo": 0,
                "pose_source": "ga_cluster",
                "genuine": True,
            },
            {
                "pdb": "1BBB",
                "rmsd_hungarian": 4.0,
                "rmsd_to_crystal": 5.0,
                "seed_echo": 0,
                "pose_source": "ga_cluster",
                "genuine": False,
            },
            {
                "pdb": "1CCC",
                "rmsd_hungarian": 0.3,
                "rmsd_to_crystal": 0.4,
                "seed_echo": 1,
                "pose_source": "ini_elitism",
                "genuine": False,
            },
        ],
    }
    (camp / "summary.json").write_text(json.dumps(summary))
    (camp / "provenance.txt").write_text(
        "root /tmp/x\n"
        "commit abcdef1234567890abcdef1234567890abcdef12\n"
        "matrix_md5 72d7c7396702331d96ff12d18f831796\n"
    )
    a = mod.audit_campaign_dir(camp)
    assert a.kind == "canary"
    assert a.n_targets == 3
    assert a.n_genuine == 1
    assert a.genuine_rate == pytest.approx(1 / 3)
    assert a.git_commit.startswith("abcdef12")
    assert a.matrix_md5.startswith("72d7c739")


def test_audit_campaign_from_result_csv(tmp_path: Path):
    mod = _load()
    camp = tmp_path / "C0_full85_fixture"
    camp.mkdir()
    tdir = camp / "1G9V"
    tdir.mkdir()
    fields = [
        "pdb_id",
        "rmsd_to_crystal",
        "rmsd_hungarian",
        "seed_echo",
        "pose_source",
        "best_cluster_rmsd",
    ]
    with (tdir / "result.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerow(
            {
                "pdb_id": "1G9V",
                "rmsd_to_crystal": "5.5",
                "rmsd_hungarian": "5.2",
                "seed_echo": "0",
                "pose_source": "ga_cluster",
                "best_cluster_rmsd": "2.2",
            }
        )
    (camp / "RUN_RECEIPT.json").write_text(
        json.dumps(
            {
                "matrix_md5": "9dc93717dfed0698006d88dd6a9627bc",
                "git_commit": "deadbeefcafebabe0123456789abcdef01234567",
                "seed_elitism": "0",
            }
        )
    )
    a = mod.audit_campaign_dir(camp)
    assert a.kind == "campaign"
    assert a.n_targets == 1
    assert a.n_genuine == 0
    assert a.matrix_md5.startswith("9dc93717")
    assert a.git_commit.startswith("deadbeef")


def test_audit_final_md_fallback(tmp_path: Path):
    mod = _load()
    camp = tmp_path / "canary_only_final"
    camp.mkdir()
    (camp / "FINAL.md").write_text("# Canary\n\nGenuine: 1 / 3\n")
    a = mod.audit_campaign_dir(camp)
    assert a.n_targets == 3
    assert a.n_genuine == 1
    assert "FINAL.md" in a.notes


def test_cli_writes_json_and_md(tmp_path: Path):
    root = tmp_path / "results"
    root.mkdir()
    camp = root / "canary_cli"
    camp.mkdir()
    (camp / "FINAL.md").write_text("Genuine: 0 / 3\n")
    out = tmp_path / "out"
    proc = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--roots",
            str(root),
            "--out-dir",
            str(out),
            "--repo",
            str(ROOT),
        ],
        capture_output=True,
        text=True,
        cwd=str(ROOT),
    )
    assert proc.returncode == 0, proc.stderr
    jpath = out / "audit_benchmark_versions.json"
    mpath = out / "audit_benchmark_versions.md"
    assert jpath.is_file()
    assert mpath.is_file()
    data = json.loads(jpath.read_text())
    assert data["n_campaigns"] == 1
    assert data["campaigns"][0]["n_genuine"] == 0
    assert "version_deltas" in data
    assert "genuine" in mpath.read_text().lower()


def test_build_version_deltas_empty_without_commits():
    mod = _load()
    # No commits → empty deltas
    deltas = mod.build_version_deltas([], ROOT)
    assert deltas == []
