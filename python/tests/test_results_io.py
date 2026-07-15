import csv
import json
from pathlib import Path

import pytest

from flexaidds.models import DockingResult
from flexaidds.results import load_results


def _write_pdb(path: Path, remarks: list[str]) -> None:
    lines = [f"REMARK {line}\n" for line in remarks]
    lines.extend(
        [
            "ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n",
            "END\n",
        ]
    )
    path.write_text("".join(lines), encoding="utf-8")


def test_load_results_groups_poses_by_mode(tmp_path: Path) -> None:
    _write_pdb(
        tmp_path / "binding_mode_1_pose_1.pdb",
        [
            "binding_mode = 1",
            "pose_rank = 1",
            "CF = -42.5",
            "free_energy = -41.0",
            "enthalpy = -40.0",
            "entropy = 0.0033",
            "temperature = 300.0",
        ],
    )
    _write_pdb(
        tmp_path / "binding_mode_1_pose_2.pdb",
        [
            "binding_mode = 1",
            "pose_rank = 2",
            "CF = -39.0",
            "temperature = 300.0",
        ],
    )
    _write_pdb(
        tmp_path / "binding_mode_2_pose_1.pdb",
        [
            "binding_mode = 2",
            "pose_rank = 1",
            "CF = -35.0",
            "free_energy = -34.2",
            "temperature = 300.0",
        ],
    )

    result = load_results(tmp_path)

    assert result.n_modes == 2
    assert result.temperature == 300.0
    assert result.binding_modes[0].mode_id == 1
    assert result.binding_modes[0].n_poses == 2
    assert result.binding_modes[0].best_cf == -42.5
    assert result.binding_modes[0].free_energy == -41.0
    assert result.binding_modes[1].mode_id == 2
    assert result.binding_modes[1].best_cf == -35.0


def test_load_results_cluster_remark_mode_identity(tmp_path: Path) -> None:
    """Real production Cluster N REMARK lines must not collapse to mode_id=1."""
    for mid, cf in ((0, -86.1), (1, -69.5), (2, -2.9)):
        _write_pdb(
            tmp_path / f"complex_{mid}.pdb",
            [
                f"CF = {cf}",
                f"Cluster {mid}: Rank (top):{mid} Average CF:{cf - 10} Frequency:10",
            ],
        )
    result = load_results(tmp_path)
    mode_ids = sorted(m.mode_id for m in result.binding_modes)
    assert mode_ids == [0, 1, 2]
    # CF-only files should recompute F (labelled)
    for mode in result.binding_modes:
        assert mode.free_energy is not None
        assert mode.metadata.get("ledger_source") == "ensemble_estimate_from_cf"


def test_load_results_engine_ledger_not_recomputed(tmp_path: Path) -> None:
    _write_pdb(
        tmp_path / "mode_1_pose_1.pdb",
        [
            "binding_mode = 1",
            "pose_rank = 1",
            "CF = -42.5",
            "free_energy = -41.0",
            "enthalpy = -40.0",
            "entropy = 0.0033",
            "temperature = 300.0",
        ],
    )
    result = load_results(tmp_path)
    mode = result.binding_modes[0]
    assert mode.free_energy == pytest.approx(-41.0)
    assert mode.metadata.get("ledger_source") == "engine_remark"


def test_load_results_uses_filename_heuristics_when_remarks_missing(tmp_path: Path) -> None:
    _write_pdb(
        tmp_path / "mode_7_pose_3.pdb",
        [
            "CF = -11.5",
            "temperature = 298.15",
        ],
    )
    _write_pdb(
        tmp_path / "cluster-12_conformer-4.pdb",
        [
            "CF = -9.25",
            "temperature = 298.15",
        ],
    )

    result = load_results(tmp_path)
    first_mode = result.binding_modes[0]
    first_pose = first_mode.poses[0]
    second_mode = result.binding_modes[1]
    second_pose = second_mode.poses[0]

    assert first_mode.mode_id == 7
    assert first_pose.pose_rank == 3
    assert first_mode.best_cf == -11.5
    assert second_mode.mode_id == 12
    assert second_pose.pose_rank == 4
    assert second_mode.best_cf == -9.25


def test_from_json_round_trip_string(tmp_path: Path) -> None:
    """to_json() → from_json() round-trips mode-level scalars via string."""
    _write_pdb(
        tmp_path / "binding_mode_1_pose_1.pdb",
        [
            "binding_mode = 1",
            "pose_rank = 1",
            "CF = -42.5",
            "free_energy = -41.0",
            "enthalpy = -40.0",
            "entropy = 0.0033",
            "heat_capacity = 0.012",
            "std_energy = 1.5",
            "temperature = 300.0",
        ],
    )
    _write_pdb(
        tmp_path / "binding_mode_2_pose_1.pdb",
        [
            "binding_mode = 2",
            "pose_rank = 1",
            "CF = -35.0",
            "free_energy = -34.2",
            "temperature = 300.0",
        ],
    )

    original = load_results(tmp_path)
    json_text = original.to_json()
    restored = DockingResult.from_json(json_text)

    assert restored.n_modes == original.n_modes
    assert restored.temperature == original.temperature
    assert str(restored.source_dir) == str(original.source_dir)

    for orig_mode, rest_mode in zip(original.binding_modes, restored.binding_modes):
        assert rest_mode.mode_id == orig_mode.mode_id
        assert rest_mode.rank == orig_mode.rank
        assert rest_mode.free_energy == orig_mode.free_energy
        assert rest_mode.enthalpy == orig_mode.enthalpy
        assert rest_mode.entropy == orig_mode.entropy
        assert rest_mode.heat_capacity == orig_mode.heat_capacity
        assert rest_mode.std_energy == orig_mode.std_energy
        assert rest_mode.best_cf == orig_mode.best_cf
        assert rest_mode.temperature == orig_mode.temperature


def test_from_json_round_trip_file(tmp_path: Path) -> None:
    """to_json(path) → from_json(path) round-trips via a file on disk."""
    _write_pdb(
        tmp_path / "mode_1_pose_1.pdb",
        [
            "binding_mode = 1",
            "pose_rank = 1",
            "CF = -10.0",
            "free_energy = -9.5",
            "temperature = 310.0",
        ],
    )

    original = load_results(tmp_path)
    json_path = tmp_path / "results.json"
    original.to_json(json_path)

    restored = DockingResult.from_json(json_path)

    assert restored.n_modes == 1
    assert restored.temperature == 310.0
    assert restored.binding_modes[0].free_energy == -9.5
    assert restored.binding_modes[0].best_cf == -10.0


def test_from_json_empty_modes() -> None:
    """from_json handles results with no binding modes."""
    payload = json.dumps({
        "source_dir": "/tmp/empty",
        "temperature": 300.0,
        "n_modes": 0,
        "metadata": {},
        "binding_modes": [],
    })
    restored = DockingResult.from_json(payload)
    assert restored.n_modes == 0
    assert restored.temperature == 300.0


# ---------------------------------------------------------------------------
# Backward compatibility: old REMARK formats
# ---------------------------------------------------------------------------

def test_old_uppercase_remark_keys(tmp_path: Path) -> None:
    """Upper-case REMARK keys (older engine versions) are parsed correctly."""
    _write_pdb(
        tmp_path / "binding_mode_1_pose_1.pdb",
        [
            "BINDING_MODE = 1",
            "POSE_RANK = 1",
            "CF = -30.0",
            "FREE_ENERGY = -29.5",
            "TEMPERATURE = 300.0",
        ],
    )
    result = load_results(tmp_path)
    assert result.n_modes == 1
    mode = result.binding_modes[0]
    assert mode.best_cf == -30.0


def test_missing_temperature_remark_falls_back(tmp_path: Path) -> None:
    """When the REMARK Temperature line is absent, temperature is None or falls back."""
    _write_pdb(
        tmp_path / "binding_mode_1_pose_1.pdb",
        [
            "binding_mode = 1",
            "pose_rank = 1",
            "CF = -22.0",
            # NO temperature line
        ],
    )
    result = load_results(tmp_path)
    # Should not raise; temperature may be None when not provided
    assert result.n_modes == 1
    mode = result.binding_modes[0]
    assert mode.best_cf == -22.0


def test_mixed_remark_versions_same_directory(tmp_path: Path) -> None:
    """Directory mixing old-style and new-style REMARK formats loads all modes."""
    _write_pdb(
        tmp_path / "binding_mode_1_pose_1.pdb",
        [
            "binding_mode = 1",
            "pose_rank = 1",
            "CF = -20.0",
            "temperature = 300.0",
        ],
    )
    _write_pdb(
        tmp_path / "binding_mode_2_pose_1.pdb",
        [
            "BINDING_MODE = 2",  # old-style upper-case
            "POSE_RANK = 1",
            "CF = -18.0",
            "TEMPERATURE = 300.0",
        ],
    )
    result = load_results(tmp_path)
    assert result.n_modes == 2
    cf_values = {m.best_cf for m in result.binding_modes}
    assert -20.0 in cf_values
    assert -18.0 in cf_values


def test_no_pdb_files_returns_empty_result(tmp_path: Path) -> None:
    """A directory with no PDB files returns an empty DockingResult gracefully."""
    result = load_results(tmp_path)
    assert result.n_modes == 0


def test_partial_remarks_only_cf(tmp_path: Path) -> None:
    """Files with only CF REMARK still load; F is recomputed from CF (labelled)."""
    _write_pdb(
        tmp_path / "binding_mode_1_pose_1.pdb",
        [
            "binding_mode = 1",
            "pose_rank = 1",
            "CF = -11.0",
            "temperature = 298.15",
        ],
    )
    result = load_results(tmp_path)
    assert result.n_modes == 1
    mode = result.binding_modes[0]
    assert mode.best_cf == -11.0
    # Smoothing: recompute ensemble estimate from CF when ledger REMARKs missing
    assert mode.free_energy is not None
    assert mode.free_energy == pytest.approx(-11.0)  # single sample → F ≈ E
    assert mode.metadata.get("ledger_source") == "ensemble_estimate_from_cf"


def test_multiple_poses_per_mode_aggregated(tmp_path: Path) -> None:
    """Multiple poses from one mode return the best CF as best_cf."""
    for rank, cf in [(1, -25.0), (2, -18.0), (3, -12.0)]:
        _write_pdb(
            tmp_path / f"binding_mode_1_pose_{rank}.pdb",
            [
                "binding_mode = 1",
                f"pose_rank = {rank}",
                f"CF = {cf}",
                "temperature = 300.0",
            ],
        )
    result = load_results(tmp_path)
    assert result.n_modes == 1
    assert result.binding_modes[0].best_cf == -25.0
    assert result.binding_modes[0].n_poses == 3


def test_to_csv_round_trip(tmp_path: Path) -> None:
    """to_csv() produces valid CSV that round-trips mode IDs and CF values."""
    for mode_id, cf in [(1, -30.0), (2, -25.5)]:
        _write_pdb(
            tmp_path / f"binding_mode_{mode_id}_pose_1.pdb",
            [
                f"binding_mode = {mode_id}",
                "pose_rank = 1",
                f"CF = {cf}",
                "temperature = 300.0",
            ],
        )
    result = load_results(tmp_path)
    csv_text = result.to_csv()

    assert csv_text.strip()
    lines = [l for l in csv_text.splitlines() if l.strip()]
    # Header + 2 data rows
    assert len(lines) == 3

    reader = csv.DictReader(csv_text.splitlines())
    rows = list(reader)
    assert len(rows) == 2
    cfs = {float(r["best_cf"]) for r in rows}
    assert -30.0 in cfs
    assert -25.5 in cfs


def test_subdirectory_pdb_files_loaded(tmp_path: Path) -> None:
    """PDB files in sub-directories are discovered recursively."""
    subdir = tmp_path / "run1"
    subdir.mkdir()
    _write_pdb(
        subdir / "binding_mode_1_pose_1.pdb",
        [
            "binding_mode = 1",
            "pose_rank = 1",
            "CF = -17.0",
            "temperature = 300.0",
        ],
    )
    result = load_results(tmp_path)
    assert result.n_modes == 1
    assert result.binding_modes[0].best_cf == -17.0
