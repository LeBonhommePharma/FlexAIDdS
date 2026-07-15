"""Security and completion tests for the Bonhomme Fleet dashboard."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def _load(name: str, relative: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


server = _load("fleet_status_server_test", "benchmarks/m3pro/dashboard/fleet_status_server.py")
monitor = _load("fleet_monitor_test", "benchmarks/m3pro/dashboard/fleet_monitor.py")


def test_pose_file_does_not_commit_target(tmp_path: Path) -> None:
    target = tmp_path / "1GPK"
    target.mkdir()
    (target / "1GPK_0.pdb").write_text("ATOM\n", encoding="utf-8")

    assert server._classify_target(str(target))["state"] == "in_progress"
    assert not server._classify_target(str(target))["has_result"]
    assert monitor.scan_target_dir(str(target)).state == "in_progress"
    assert not monitor.scan_target_dir(str(target)).has_result


def test_only_nonempty_result_csv_commits_target(tmp_path: Path) -> None:
    target = tmp_path / "1GPK"
    target.mkdir()
    (target / "result.csv").touch()
    assert server._classify_target(str(target))["state"] != "done"
    assert monitor.scan_target_dir(str(target)).state != "done"

    (target / "result.csv").write_text("pdb_id,success_pb\n1GPK,1\n", encoding="utf-8")
    assert server._classify_target(str(target))["state"] == "done"
    assert monitor.scan_target_dir(str(target)).state == "done"


def test_progress_components_and_symlinks_cannot_escape_roots(tmp_path: Path) -> None:
    results = tmp_path / "results" / "tier2"
    results.mkdir(parents=True)
    with pytest.raises(ValueError, match="invalid dataset"):
        server.deep_scan_run([str(results)], "..", "run1")

    outside = tmp_path / "outside"
    outside.mkdir()
    link = results / "astex"
    try:
        link.symlink_to(outside, target_is_directory=True)
    except OSError:
        pytest.skip("symlinks are unavailable")
    with pytest.raises(ValueError, match="traversal"):
        server.deep_scan_run([str(results)], "astex", "run1")


def test_server_is_loopback_and_cors_off_by_default() -> None:
    args = server.build_argparser().parse_args([])
    assert args.host == "127.0.0.1"
    assert args.cors_origin is None


def test_merged_status_matches_pwa_root_contract(tmp_path: Path) -> None:
    payload = {
        "runner": "worker-one",
        "devices": [{"deviceID": "worker-one", "model": "Mac"}],
        "activeChunks": [{"id": "chunk-0000"}],
        "metrics": {
            "totalChunks": 2,
            "completedChunks": 1,
            "failedChunks": 0,
            "orphanedChunks": 0,
            "activeDevices": 1,
            "totalTFLOPS": 0.0,
        },
    }
    (tmp_path / "fleet_status_worker.json").write_text(json.dumps(payload), encoding="utf-8")
    merged = server.merge_fleet_status(str(tmp_path))

    assert merged["devices"] == payload["devices"]
    assert merged["activeChunks"] == payload["activeChunks"]
    assert merged["metrics"]["completedChunks"] == 1
    assert merged["metrics"]["totalChunks"] == 2


def test_campaign_status_files_are_discovered_and_parsed(tmp_path: Path) -> None:
    campaign = tmp_path / "results" / "campaigns" / "astex-fleet"
    campaign.mkdir(parents=True)
    payload = {
        "campaign_id": "astex-fleet",
        "dataset": "astex-diverse",
        "totalChunks": 85,
        "completedChunks": 4,
        "states": {"completed": 4, "running": 1, "failed": 0},
        "metrics": {"failedChunks": 0, "estimatedRemainingSeconds": None},
    }
    status = campaign / "status.json"
    status.write_text(json.dumps(payload), encoding="utf-8")

    assert str(status) in server.discover_fleet_status_files(str(tmp_path))
    assert str(status) in monitor.discover_fleet_status_files(str(tmp_path))
    parsed = monitor.parse_runner(str(status), [])
    assert parsed.name == "astex-fleet"
    assert parsed.active_dataset == "astex-diverse"
    assert parsed.datasets["astex-diverse"].status == "running"
