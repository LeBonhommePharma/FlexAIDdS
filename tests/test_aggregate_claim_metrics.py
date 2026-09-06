"""Receipt consistency, frozen observation units, and CLI rejection regressions.

Synthetic CSV controls exercise admission logic. They do not stand in for real
validator execution; every report must state that limitation.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/aggregate_claim_metrics.py"
DEFAULT_PIN = "9dc93717dfed0698006d88dd6a9627bc"
POSE = hashlib.sha256(b"synthetic elected pose").hexdigest()
PB_INPUT = hashlib.sha256(b"synthetic PB input SDF").hexdigest()


def _load():
    spec = importlib.util.spec_from_file_location("aggregate_claim_metrics_test", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def agg():
    return _load()


def _row(pdb_id="1G9V", **changes):
    row = dict(pdb_id=pdb_id, seed="12345", endpoint="argmin(ACF)",
               seed_echo="0", native_pose_seeded="0", protocol_claim_eligible="1",
               matrix_md5=DEFAULT_PIN, rmsd_to_crystal="1.0", rmsd_hungarian="0.4",
               success_rmsd="1", success_pb="1", pb_pass="1", claim_ready="1",
               score_pose_consistent="1", score_pose_delta="0", pb_backend="bust_cli",
               tencom_status="ok", eigen_status="ok", eigen_n_modes="1", elected_H_vib="-1.5",
               pb_ran="1", pb_n_pass="27", pb_n_fail="0", pb_n_checks="27",
               pose_sha256=POSE, rmsd_pose_sha256=POSE, posebusters_pose_sha256=POSE,
               tencom_pose_sha256=POSE, posebusters_input_sha256=PB_INPUT,
               best_cluster_rmsd="0.5")
    row.update({k: str(v) for k, v in changes.items()})
    return row


def _write(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return path


def _cli(path, *args, script=SCRIPT):
    result = subprocess.run([sys.executable, "-B", str(script), "--csv", str(path), "--quiet", *args],
                            text=True, capture_output=True, timeout=15)
    report = json.loads(result.stdout) if result.stdout.strip() else None
    return result, report


def _aggregate(agg, rows, **kwargs):
    return agg.aggregate_rows(rows, DEFAULT_PIN, "test", **kwargs)


def test_valid_receipt_control_cli(tmp_path):
    result, report = _cli(_write(tmp_path / "result.csv", [_row()]))
    assert result.returncode == 0, result.stderr
    assert report["metrics"]["STRICT"]["n"] == 1
    assert report["headline"]["rate"] == 1 / 85
    assert report["evidence_level"] == "validated_receipt_fields"
    assert report["artifacts_verified"] is False
    assert report["aggregation"]["mode"] == "single_observation"


@pytest.mark.parametrize("changes", [
    {"pb_pass": "0", "success_pb": "0"},
    {"pb_pass": "0", "success_pb": "1"},
    {"rmsd_to_crystal": "-1"}, {"rmsd_to_crystal": "9"},
    {"rmsd_to_crystal": "nan"}, {"rmsd_to_crystal": "inf"},
    {"score_pose_consistent": "0"}, {"score_pose_delta": "0.0001001"},
    {"score_pose_delta": "nan"}, {"score_pose_delta": "inf"},
    {"score_pose_delta": ""}, {"tencom_status": ""}, {"eigen_status": "fail"},
    {"pb_backend": "internal"}, {"pose_sha256": ""}, {"rmsd_pose_sha256": ""},
    {"tencom_pose_sha256": "c" * 64}, {"posebusters_pose_sha256": "missing"},
    {"posebusters_input_sha256": ""}, {"matrix_md5": ""}, {"matrix_md5": "0" * 32},
    {"protocol_claim_eligible": ""}, {"protocol_claim_eligible": "0"},
    {"seed_echo": ""}, {"seed_echo": "1"}, {"native_pose_seeded": "1"},
    {"native_pose_seed_fraction": "0.1"}, {"claim_ready": "0"},
    {"success_rmsd": "0"}, {"success_pb": "0"},
    {"docking_completed": "0"}, {"docking_completed": ""},
    {"docking_exit_code": "-1"}, {"docking_exit_code": ""},
    {"num_poses": "0"}, {"num_poses": "-1"}, {"num_poses": "1.5"},
    {"pb_ran": "0"}, {"pb_n_fail": "1"}, {"pb_n_checks": "0"},
    {"pb_n_pass": "0"}, {"pb_n_pass": "26", "pb_n_checks": "27", "pb_n_fail": "0"},
    {"eigen_n_modes": "0"}, {"eigen_n_modes": "nan"},
    {"elected_H_vib": "nan"}, {"native_pose_seed_fraction": ""},
    {"pb_n_checks": "1", "pb_n_pass": "1"},
    {"pb_failed_keys": "bond_lengths"}, {"rmsd_fail_reason": "ref_empty"},
    {"pb_ran": ""}, {"pb_n_checks": ""}, {"elected_H_vib": ""},
])
def test_stale_strict_flag_cannot_override_bad_evidence(tmp_path, changes):
    result, report = _cli(_write(tmp_path / "result.csv", [_row(**changes)]))
    assert result.returncode == 1, result.stderr
    assert report["metrics"]["STRICT"]["n"] == 0
    assert report["dropped_rows"][0]["reasons"]


def test_four_field_self_attestation_is_not_strict(tmp_path):
    minimal = dict(pdb_id="1G9V", seed_echo="0", native_pose_seeded="0", claim_ready="1")
    result, report = _cli(_write(tmp_path / "result.csv", [minimal]))
    assert result.returncode == 1
    assert report["metrics"]["STRICT"]["n"] == 0


def test_matching_malformed_hashes_not_receipts(agg):
    row = _row(**{key: "abc" for key in (
        "pose_sha256", "rmsd_pose_sha256", "posebusters_pose_sha256", "tencom_pose_sha256", "posebusters_input_sha256")})
    assert not agg.is_strict_success(row)


def test_score_tolerance_and_rmsd_boundary(agg):
    assert agg.is_strict_success(_row(score_pose_delta="-0.0001", rmsd_to_crystal="2.0"))
    assert not agg.is_strict_success(_row(score_pose_delta="-0.0001000001"))


def test_diagnostics_survive_strict_failure(agg):
    report = _aggregate(agg, [_row(claim_ready="0", tencom_status="fail")])
    assert report["N_claim"] == 0
    assert report["N_protocol_eligible"] == 1
    assert {key: metric["n"] for key, metric in report["metrics"].items()} == dict(S1=1, S2=1, STRICT=0, S3=1)


def test_nonfinite_rmsd_does_not_fall_back_to_stale_success_flags(agg):
    assert not agg.is_s1(_row(rmsd_to_crystal="nan"))
    assert not agg.is_s3(_row(best_cluster_rmsd="nan", success_s3="1"))


def test_hungarian_never_substitutes_for_serial(agg):
    row = _row(rmsd_to_crystal="5.0", rmsd_hungarian="0.1")
    assert not agg.is_s1(row)
    assert not agg.is_strict_success(row)


def test_s2_recomputes_pb_and_requires_pose_identity(agg):
    assert not agg.is_s2(_row(pb_pass="0", success_pb="1"), True)
    assert not agg.is_s2(_row(posebusters_pose_sha256=""), True)
    # Engine success_pb includes serial RMSD; it must not veto a graph-symmetry diagnostic.
    assert agg.is_s2(_row(success_pb="0"), True)


def test_seeded_rows_excluded_from_diagnostics(agg):
    report = _aggregate(agg, [_row(seed_echo="1")])
    assert all(metric["n"] == 0 for metric in report["metrics"].values())


def test_off_manifest_targets_never_inflate_any_metric(agg):
    report = _aggregate(agg, [_row(), _row("9XXX")])
    assert all(metric["n"] == 1 for metric in report["metrics"].values())
    assert report["off_manifest_targets"] == ["9XXX"]


@pytest.mark.parametrize("count", [2, 86])
def test_duplicate_observations_rejected_cli(tmp_path, count):
    result, report = _cli(_write(tmp_path / "result.csv", [_row()] * count))
    assert result.returncode == 2
    assert report is None
    assert "duplicate observation" in result.stderr


def test_repeated_targets_require_frozen_seed_list(agg):
    rows = [_row(seed=seed) for seed in ("1", "2", "3")]
    with pytest.raises(ValueError, match="expected-seeds"):
        _aggregate(agg, rows)


def test_majority_of_expected_seeds_not_union(agg):
    rows = [_row(seed="1"), _row(seed="2", claim_ready="0"), _row(seed="3", claim_ready="0")]
    report = _aggregate(agg, rows, expected_seeds=["1", "2", "3"])
    assert report["metrics"]["STRICT"]["n"] == 0
    assert report["metrics"]["S1"]["n"] == 1
    report = _aggregate(agg, [_row(seed="1"), _row(seed="2")], expected_seeds=["1", "2", "3"])
    assert report["metrics"]["STRICT"]["n"] == 1
    assert dict(pdb_id="1G9V", seed="3") in report["aggregation"]["missing_expected_observations"]


def test_missing_seeds_and_even_seed_ties_fail(agg):
    report = _aggregate(agg, [_row(seed="1")], expected_seeds=["1", "2", "3"])
    assert report["metrics"]["STRICT"]["n"] == 0
    report = _aggregate(agg, [_row(seed="1"), _row(seed="2")], expected_seeds=["1", "2", "3", "4"])
    assert report["metrics"]["STRICT"]["n"] == 0


@pytest.mark.parametrize("seeds", [["1", "1"], [], [""]])
def test_invalid_expected_seed_lists_rejected(agg, seeds):
    with pytest.raises(ValueError, match="expected seeds"):
        _aggregate(agg, [_row()], expected_seeds=seeds)


def test_unexpected_seed_rejected(agg):
    with pytest.raises(ValueError, match="unexpected seed"):
        _aggregate(agg, [_row(seed="4")], expected_seeds=["1", "2", "3"])


def test_mixed_arms_rejected_and_explicit_selection_counted(agg):
    rows = [_row(arm="A"), _row(arm="B")]
    with pytest.raises(ValueError, match="mixed arm"):
        _aggregate(agg, rows)
    report = _aggregate(agg, rows, arm="A")
    assert report["aggregation"]["N_filtered_by_arm"] == 1
    assert report["aggregation"]["arm"] == "A"
    assert report["metrics"]["STRICT"]["n"] == 1


def test_mixed_endpoints_rejected(agg):
    with pytest.raises(ValueError, match="mixed endpoint"):
        _aggregate(agg, [_row(), _row("1GM8", endpoint="argmin(CF)")])


@pytest.mark.parametrize("text", [
    "pdb_id,claim_ready,claim_ready\n1G9V,0,1\n",
    "pdb_id,claim_ready\n1G9V,1,extra\n",
    "pdb_id,claim_ready\n1G9V\n", "pdb_id,\n1G9V,1\n",
    "pdb_id,target\n1G9V,1GM8\n", "claim_ready\n1\n",
])
def test_malformed_csv_rejected_cli(tmp_path, text):
    path = tmp_path / "result.csv"
    path.write_text(text)
    result, report = _cli(path)
    assert result.returncode == 2
    assert report is None


def test_stale_tree_cannot_override_summary(agg, tmp_path):
    _write(tmp_path / "old" / "result.csv", [_row()])
    _write(tmp_path / "summary.csv", [_row(claim_ready="0")])
    with pytest.raises(ValueError, match="ambiguous campaign sources"):
        agg.load_campaign_rows(tmp_path)


def test_multiple_flat_summaries_require_explicit_selection(agg, tmp_path):
    _write(tmp_path / "results.csv", [_row()])
    _write(tmp_path / "summary.csv", [_row()])
    with pytest.raises(ValueError, match="ambiguous"):
        agg.load_campaign_rows(tmp_path)


def test_per_target_file_cannot_discard_extra_rows(agg, tmp_path):
    _write(tmp_path / "one" / "result.csv", [_row(), _row("1GM8")])
    with pytest.raises(ValueError, match="exactly one"):
        agg.load_campaign_rows(tmp_path)


def test_manifest_integrity_and_fixed_denominator(agg):
    codes, digest = agg.load_target_manifest()
    assert len(codes) == 85
    assert digest == hashlib.sha256(",".join(codes).encode()).hexdigest()
    report = _aggregate(agg, [_row()])
    assert report["N_denominator"] == 85
    assert report["N_missing_from_manifest"] == 84


@pytest.mark.parametrize("mutation", ["digest", "duplicate", "count", "empty"])
def test_corrupt_manifest_fails(agg, tmp_path, mutation):
    codes, digest = agg.load_target_manifest()
    data = dict(schema="flexaidds.astex.target_manifest/v1", N=85, targets=codes, sha256_of_sorted_codes=digest)
    if mutation == "digest": data["sha256_of_sorted_codes"] = "0" * 64
    if mutation == "duplicate": data["targets"][-1] = data["targets"][0]
    if mutation == "count": data["N"] = 84
    if mutation == "empty": data["targets"] = []
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        agg.load_target_manifest(path)


def test_missing_manifest_never_silently_shrinks_denominator(tmp_path):
    script = tmp_path / "standalone" / "aggregate_claim_metrics.py"
    script.parent.mkdir()
    script.write_bytes(SCRIPT.read_bytes())
    source = _write(tmp_path / "result.csv", [_row()])
    result, report = _cli(source, script=script)
    assert result.returncode == 2
    assert report is None
    assert "manifest missing" in result.stderr
    result, report = _cli(source, "--legacy-observed-denominator", "--diagnostic-only", "--headline", "s1", script=script)
    assert result.returncode == 0
    assert report["headline"]["diagnostic_only"]
    assert report["metrics"]["STRICT"]["rate"] is None
    assert report["metrics"]["STRICT"]["n"] == 0


def test_legacy_denominator_cannot_emit_strict_headline(tmp_path):
    result, report = _cli(_write(tmp_path / "result.csv", [_row()]), "--legacy-observed-denominator", "--diagnostic-only")
    assert result.returncode == 2
    assert report is None


def test_matrix_receipt_conflicts_are_errors(agg, tmp_path):
    (tmp_path / "RUN_RECEIPT.json").write_text(json.dumps(dict(matrix_md5="0" * 32)))
    with pytest.raises(ValueError, match="conflicting matrix pins"):
        agg.load_matrix_pin(tmp_path, DEFAULT_PIN)


@pytest.mark.parametrize("sha", ["", "abc", "c" * 64])
def test_sidecar_requires_matching_valid_pose_hash(agg, sha):
    row = _row()
    sidecar = {"1G9V": dict(pdb_id="1G9V", rmsd_symmcorr="0.5", pose_sha256=sha)}
    report = agg.join_symmcorr([row], sidecar)
    assert report["joined"] == 0
    assert report["refused_sha_mismatch"] == ["1G9V"]
    assert "rmsd_symmcorr" not in row


def test_sidecar_same_target_multiple_poses_joins_by_hash(agg, tmp_path):
    second = hashlib.sha256(b"second elected pose").hexdigest()
    side = [dict(pdb_id="1G9V", rmsd_symmcorr="0.5", pose_sha256=POSE, status="ok"),
            dict(pdb_id="1G9V", rmsd_symmcorr="1.5", pose_sha256=second, status="ok")]
    rows = [_row(), _row(pose_sha256=second, seed="2")]
    report = agg.join_symmcorr(rows, agg.load_symmcorr_sidecar(_write(tmp_path / "side.csv", side)))
    assert report["joined"] == 2
    assert [r["rmsd_symmcorr"] for r in rows] == ["0.5", "1.5"]


def test_sidecar_duplicate_identity_rejected(agg, tmp_path):
    row = dict(pdb_id="1G9V", rmsd_symmcorr="0.5", pose_sha256=POSE, status="ok")
    with pytest.raises(ValueError, match="duplicate sidecar"):
        agg.load_symmcorr_sidecar(_write(tmp_path / "side.csv", [row, row]))


def test_s3_primary_requires_diagnostic_flag(tmp_path):
    source = _write(tmp_path / "result.csv", [_row()])
    result, _ = _cli(source, "--headline", "s3")
    assert result.returncode == 2
    result, report = _cli(source, "--headline", "s3", "--diagnostic-only")
    assert result.returncode == 0
    assert report["headline"]["diagnostic_only"]


def test_cloud_docs_requires_staging(agg):
    with pytest.raises(ValueError, match="stage CloudDocs"):
        agg.load_rows_from_csv(Path("/tmp/Mobile Documents/com~apple~CloudDocs/result.csv"))


def test_seed_echo_explicit_zero_float_accepted(agg):
    assert agg.is_strict_success(_row(seed_echo="0.0", native_pose_seeded="0.0"))


def test_explicit_completed_runtime_receipt_accepted(agg):
    assert agg.is_strict_success(_row(docking_completed="1", docking_exit_code="0", num_poses="2",
                                      pb_ran="1", pb_n_checks="27", pb_n_fail="0", pb_n_pass="27",
                                      eigen_n_modes="3", elected_H_vib="-1.5"))


def test_cli_exit_uses_target_majority_not_individual_seed_success(tmp_path):
    result, report = _cli(_write(tmp_path / "result.csv", [_row(seed="1")]), "--expected-seeds", "1,2,3")
    assert report["N_claim"] == 1
    assert report["metrics"]["STRICT"]["n"] == 0
    assert result.returncode == 1



def test_decimal_seed_aliases_cannot_manufacture_a_majority(agg):
    with pytest.raises(ValueError, match="duplicate observation"):
        _aggregate(agg, [_row(seed="1"), _row(seed="01")], expected_seeds=["1", "01", "2"])
    with pytest.raises(ValueError, match="unique list"):
        _aggregate(agg, [_row(seed="1")], expected_seeds=["1", "01", "2"])


@pytest.mark.parametrize("seed", ["abc", "-1", "1e3", str(2**64)])
def test_invalid_seed_domain_rejected(agg, seed):
    with pytest.raises(ValueError, match="unsigned 64-bit"):
        _aggregate(agg, [_row(seed=seed)])


def test_valid_decimal_seed_alias_normalizes_identity(agg):
    report = _aggregate(agg, [_row(seed="0001")], expected_seeds=["1"])
    assert report["metrics"]["STRICT"]["n"] == 1
    assert report["aggregation"]["expected_seeds"] == ["1"]



def test_explicit_local_results_root_preserved(agg, monkeypatch, tmp_path):
    monkeypatch.setenv("FLEXAIDDS_RESULTS", str(tmp_path / "chosen"))
    monkeypatch.setenv("FLEXAIDDS_LOCAL_ROOT", str(tmp_path / "other"))
    assert agg.resolve_c0_full85_dir() == tmp_path / "chosen" / agg.C0_FULL85_REL


@pytest.mark.parametrize("field", ["pb_ran", "pb_n_checks", "pb_n_pass", "pb_n_fail", "eigen_n_modes", "elected_H_vib"])
def test_current_strict_requires_complete_pb_eigen_receipts(agg, field):
    row = _row()
    del row[field]
    report = _aggregate(agg, [row])
    assert report["metrics"]["STRICT"]["n"] == 0
    assert report["metrics"]["S1"]["n"] == 1



@pytest.mark.parametrize("changes", [dict(pb_ran="0"), dict(pb_n_checks="1"),
                                    dict(pb_n_fail="1"), dict(pb_failed_keys="bond_lengths"),
                                    dict(pb_backend="internal")])
def test_s2_cannot_ignore_its_own_pb_contradictions(agg, changes):
    report = _aggregate(agg, [_row(**changes)])
    assert report["metrics"]["S1"]["n"] == 1
    assert report["metrics"]["S2"]["n"] == 0
