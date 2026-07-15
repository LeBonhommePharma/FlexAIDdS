"""Bonhomme Fleet control-plane tests."""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "flexaidds_fleet_module", ROOT / "python" / "flexaidds" / "fleet.py"
)
assert SPEC and SPEC.loader
fleet = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fleet)


def _binary(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    path.chmod(0o755)
    return path


def _plan_args(tmp_path: Path, **overrides: object) -> argparse.Namespace:
    runner = _binary(tmp_path / "benchmark_datasets", "runner-v1\n")
    engine = _binary(tmp_path / "FlexAIDdS", "engine-v1\n")
    posebusters = _binary(tmp_path / "bust", "posebusters-v1\n")
    codes = tmp_path / "codes.txt"
    codes.write_text("1GPK\n1HNN\n1P62\n1Q4G\n", encoding="utf-8")
    values = {
        "campaign": str(tmp_path / "campaign"),
        "campaign_id": "fleet-test",
        "runner": str(runner),
        "engine": str(engine),
        "posebusters_bin": str(posebusters),
        "benchmark": "astex",
        "dataset": "astex-diverse",
        "mode": "defined-cleft-redock",
        "codes_file": str(codes),
        "chunks": 2,
        "compute_root": str(tmp_path / "compute"),
        "cache": None,
        "threads": 1,
        "omp_threads": 4,
        "ga_population": 1000,
        "ga_generations": 6000,
        "temperature": 298.0,
        "grid_spacing": 0.375,
        "job_timeout_seconds": 10800,
        "clustering": "CF",
        "gpu": None,
        "lease_seconds": 300,
        "max_attempts": 3,
        "min_free_gb": 0.0,
        "env": ["FLEXAIDDS_RESTARTS=5", "FLEXAIDDS_SEED_ELITISM=0"],
    }
    values.update(overrides)
    return argparse.Namespace(**values)


def test_plan_is_deterministic_and_pins_binary_hashes(tmp_path: Path) -> None:
    args = _plan_args(tmp_path)
    result = fleet.plan_campaign(args)
    campaign = Path(args.campaign)
    manifest, manifest_hash = fleet.load_manifest(campaign)

    assert result["manifest_sha256"] == manifest_hash
    assert manifest["chunks"][0]["codes"] == ["1GPK", "1P62"]
    assert manifest["chunks"][1]["codes"] == ["1HNN", "1Q4G"]
    assert manifest["runner"]["sha256"] == fleet.sha256_file(Path(args.runner))
    assert manifest["engine"]["sha256"] == fleet.sha256_file(Path(args.engine))
    assert manifest["validators"]["posebusters"]["sha256"] == fleet.sha256_file(Path(args.posebusters_bin))
    assert manifest["environment"]["FLEXAIDDS_HVIB"] == "1"
    assert manifest["environment"]["FLEXAIDDS_POSEBUST_BACKEND"] == "bust_cli"
    assert manifest["environment"]["FLEXAIDDS_POSEBUSTERS_BIN"] == str(
        Path(args.posebusters_bin).resolve()
    )
    assert manifest["environment"]["FLEXAIDDS_SEED_ELITISM"] == "0"


def test_astex_defined_cleft_manifest_is_portable_and_complete() -> None:
    manifest_path = ROOT / "benchmarks/datasets/benchmark_astex_native_85.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    pairs = manifest["pairs"]

    assert manifest["n_pairs"] == len(pairs) == 85
    assert len({pair["receptor_id"] for pair in pairs}) == 85
    for pair in pairs:
        assert pair["receptor_id"] == pair["ligand_id"]
        for field in ("receptor_pdb", "ligand_sdf", "oracle_site_pdb"):
            value = Path(pair[field])
            assert not value.is_absolute(), (pair["receptor_id"], field, value)
            assert (manifest_path.parent / value).resolve().is_file(), (
                pair["receptor_id"], field, value
            )


def test_claims_are_exclusive_and_stale_leases_are_fenced(tmp_path: Path) -> None:
    args = _plan_args(tmp_path, chunks=1, lease_seconds=30)
    fleet.plan_campaign(args)
    campaign = Path(args.campaign)
    manifest, _ = fleet.load_manifest(campaign)

    chunk1, claim1 = fleet.claim_next(campaign, manifest, "worker-one")
    assert chunk1["id"] == "chunk-0000"
    assert fleet.claim_next(campaign, manifest, "worker-two") is None

    heartbeat = fleet._heartbeat_path(campaign, claim1)
    data = fleet.read_json(heartbeat)
    data["timestamp"] = "2000-01-01T00:00:00Z"
    fleet.atomic_write_json(heartbeat, data)
    chunk2, claim2 = fleet.claim_next(campaign, manifest, "worker-two")

    assert chunk2["id"] == chunk1["id"]
    assert claim2["lease_epoch"] == claim1["lease_epoch"] + 1
    assert claim2["lease_token"] != claim1["lease_token"]
    assert not fleet.lease_owned(campaign, claim1)
    assert fleet.lease_owned(campaign, claim2)


def test_runner_command_contains_attempt_and_provenance_pins(tmp_path: Path) -> None:
    args = _plan_args(tmp_path, chunks=1)
    fleet.plan_campaign(args)
    campaign = Path(args.campaign)
    manifest, manifest_hash = fleet.load_manifest(campaign)
    chunk, claim = fleet.claim_next(campaign, manifest, "worker-one")
    command = fleet.build_runner_command(
        manifest, manifest_hash, chunk, claim, tmp_path / "local"
    )

    assert command[command.index("--manifest-sha256") + 1] == manifest_hash
    assert command[command.index("--runner-sha256") + 1] == manifest["runner"]["sha256"]
    assert command[command.index("--engine-sha256") + 1] == manifest["engine"]["sha256"]
    assert command[command.index("--attempt-id") + 1] == claim["attempt_id"]
    assert command[command.index("--mode") + 1] == "defined-cleft-redock"


def test_chunk_summary_requires_exact_assigned_target_set(tmp_path: Path) -> None:
    args = _plan_args(tmp_path, chunks=1)
    fleet.plan_campaign(args)
    campaign = Path(args.campaign)
    manifest, manifest_hash = fleet.load_manifest(campaign)
    chunk, claim = fleet.claim_next(campaign, manifest, "worker-one")
    summary = {
        "schema_version": 1,
        "type": "flexaidds-fleet-chunk-result",
        "campaign_id": manifest["campaign_id"],
        "chunk_id": chunk["id"],
        "attempt_id": claim["attempt_id"],
        "worker_id": claim["worker_id"],
        "dataset": manifest["benchmark"],
        "provenance": {
            "manifest_sha256": manifest_hash,
            "runner_sha256": manifest["runner"]["sha256"],
            "engine_sha256": manifest["engine"]["sha256"],
        },
        "summary": {"target_count": len(chunk["codes"])},
        "targets": [{"pdb_id": code} for code in chunk["codes"]],
    }

    assert fleet.summary_matches_chunk(summary, manifest, manifest_hash, chunk, claim)
    summary["targets"][-1]["pdb_id"] = summary["targets"][0]["pdb_id"]
    assert not fleet.summary_matches_chunk(summary, manifest, manifest_hash, chunk, claim)


def test_verified_archive_is_immutable(tmp_path: Path) -> None:
    source = tmp_path / "source"
    source.mkdir()
    (source / "result.csv").write_text("ok\n", encoding="utf-8")
    destination = tmp_path / "archive" / "attempt-1"

    expected = fleet.tree_sha256(source)
    assert fleet.copy_tree_verified(source, destination) == expected
    assert fleet.tree_sha256(destination) == expected
    with pytest.raises(fleet.FleetError, match="already exists"):
        fleet.copy_tree_verified(source, destination)


def test_aggregate_reports_only_strict_scientific_counters(tmp_path: Path) -> None:
    args = _plan_args(tmp_path, chunks=1)
    fleet.plan_campaign(args)
    campaign = Path(args.campaign)
    manifest, manifest_hash = fleet.load_manifest(campaign)
    chunk = manifest["chunks"][0]
    accepted = {
        "schema_version": 1,
        "status": "completed",
        "campaign_id": manifest["campaign_id"],
        "chunk_id": chunk["id"],
        "attempt_id": "attempt-fixed",
        "worker_id": "worker-one",
        "lease_epoch": 1,
        "manifest_sha256": manifest_hash,
        "accepted_at": fleet.utc_now(),
        "archive": "artifacts/chunk-0000/attempt-fixed",
        "archive_tree_sha256": "a" * 64,
        "summary": {
            "targets": [
                {
                    "pdb_id": "1GPK",
                    "execution_completed": True,
                    "success_rmsd": True,
                    "posebusters_ran": True,
                    "posebusters_pass": False,
                    "success_pb": False,
                    "tencom_status": "ok",
                    "eigen_status": "ok",
                    "validators_complete": True,
                    "protocol_claim_eligible": True,
                    "claim_ready": False,
                },
                {
                    "pdb_id": "1HNN",
                    "execution_completed": True,
                    "success_rmsd": True,
                    "posebusters_ran": True,
                    "posebusters_pass": True,
                    "success_pb": True,
                    "tencom_status": "ok",
                    "eigen_status": "ok",
                    "validators_complete": True,
                    "protocol_claim_eligible": True,
                    "claim_ready": True,
                },
            ]
        },
    }
    fleet.exclusive_write_json(
        fleet._chunk_dir(campaign, chunk["id"]) / "result.json", accepted
    )

    summary = fleet.aggregate_campaign(campaign, manifest, manifest_hash)
    assert summary["execution_completed"] == 2
    assert summary["success_rmsd"] == 2
    assert summary["success_pb"] == 1
    assert summary["claim_ready"] == 1
    assert summary["success_pb_rate"] == 0.5


def test_secret_or_unknown_environment_keys_are_rejected() -> None:
    with pytest.raises(fleet.FleetError, match="not allowed"):
        fleet.parse_environment(["OPENAI_API_KEY=secret"])
    with pytest.raises(fleet.FleetError, match="not allowed"):
        fleet.parse_environment(["UNREVIEWED_KNOB=1"])


def test_campaign_paths_reject_traversal(tmp_path: Path) -> None:
    with pytest.raises(fleet.FleetError, match="outside campaign"):
        fleet.confined_path(tmp_path / "campaign", "chunks", "..", "..", "escape")


def test_oracle_campaign_is_rejected_as_production_fleet(tmp_path: Path) -> None:
    args = _plan_args(tmp_path, mode="oracle-ceiling")
    with pytest.raises(fleet.FleetError, match="diagnostic-only"):
        fleet.plan_campaign(args)


def test_end_to_end_worker_archives_then_commits_and_aggregates(tmp_path: Path) -> None:
    args = _plan_args(tmp_path, chunks=1)
    runner = Path(args.runner)
    runner.write_text(
        f"""#!{sys.executable}
import json, os, sys

def value(flag):
    return sys.argv[sys.argv.index(flag) + 1]

codes = value('--only-codes').split(',')
output_dir = value('--output')
for code in codes:
    target_dir = os.path.join(output_dir, code)
    os.makedirs(target_dir, exist_ok=True)
    with open(os.path.join(target_dir, 'result.csv'), 'x', encoding='utf-8') as stream:
        stream.write('pdb_id,success_pb\\n' + code + ',1\\n')
payload = {{
    'schema_version': 1,
    'type': 'flexaidds-fleet-chunk-result',
    'campaign_id': value('--campaign-id'),
    'chunk_id': value('--chunk-id'),
    'attempt_id': value('--attempt-id'),
    'worker_id': value('--worker-id'),
    'dataset': value('--benchmark'),
    'provenance': {{
        'manifest_sha256': value('--manifest-sha256'),
        'runner_sha256': value('--runner-sha256'),
        'engine_sha256': value('--engine-sha256'),
    }},
    'summary': {{'target_count': len(codes)}},
    'targets': [{{
        'pdb_id': target_code,
        'execution_completed': True,
        'success_rmsd': True,
        'posebusters_ran': True,
        'posebusters_pass': True,
        'success_pb': True,
        'tencom_status': 'ok',
        'eigen_status': 'ok',
        'validators_complete': True,
        'protocol_claim_eligible': True,
        'claim_ready': True,
    }} for target_code in codes],
}}
with open(value('--output-json'), 'x', encoding='utf-8') as stream:
    json.dump(payload, stream)
""",
        encoding="utf-8",
    )
    runner.chmod(0o755)
    fleet.plan_campaign(args)
    campaign = Path(args.campaign)

    worker_args = argparse.Namespace(
        campaign=str(campaign), worker_id="worker-one", once=True,
        max_chunks=1, keep_local=False,
    )
    assert fleet.run_worker(worker_args) == 0

    manifest, manifest_hash = fleet.load_manifest(campaign)
    accepted = fleet.read_json(campaign / "chunks" / "chunk-0000" / "result.json")
    archive = campaign / accepted["archive"]
    assert archive.is_dir()
    assert fleet.tree_sha256(archive) == accepted["archive_tree_sha256"]
    assert not Path(manifest["compute_root"]).joinpath(
        manifest["campaign_id"], "chunk-0000", accepted["attempt_id"]
    ).exists()

    summary = fleet.aggregate_campaign(campaign, manifest, manifest_hash)
    assert summary["complete"] is True
    assert summary["success_pb"] == len(manifest["chunks"][0]["codes"])
    assert summary["claim_ready"] == len(manifest["chunks"][0]["codes"])
