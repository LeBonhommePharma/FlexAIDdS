"""Tests for flexaidds.results – helper internals and load_results error paths.

Priority 3 coverage.  No C++ extension needed.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from flexaidds.models import BindingModeResult, PoseResult
from flexaidds.results import (
    _collect_pose_files,
    _mode_frequency,
    _mode_metadata,
    _mode_metric,
    _mode_temperature,
    _build_mode,
    load_results,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_pose(
    *,
    path: str = "x.pdb",
    mode_id: int = 1,
    pose_rank: int = 1,
    cf: float | None = None,
    temperature: float | None = None,
    free_energy: float | None = None,
    enthalpy: float | None = None,
    entropy: float | None = None,
    heat_capacity: float | None = None,
    std_energy: float | None = None,
    remarks: dict | None = None,
) -> PoseResult:
    return PoseResult(
        path=Path(path),
        mode_id=mode_id,
        pose_rank=pose_rank,
        cf=cf,
        temperature=temperature,
        free_energy=free_energy,
        enthalpy=enthalpy,
        entropy=entropy,
        heat_capacity=heat_capacity,
        std_energy=std_energy,
        remarks=remarks or {},
    )


def _write_pdb(path: Path, remarks: list[str]) -> None:
    lines = [f"REMARK {r}\n" for r in remarks]
    lines += [
        "ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n",
        "END\n",
    ]
    path.write_text("".join(lines), encoding="utf-8")


# ===========================================================================
# _collect_pose_files
# ===========================================================================

class TestCollectPoseFiles:
    def test_finds_pdb_files(self, tmp_path):
        (tmp_path / "a.pdb").write_text("END\n")
        (tmp_path / "b.pdb").write_text("END\n")
        files = _collect_pose_files(tmp_path)
        assert len(files) == 2

    def test_finds_ent_extension(self, tmp_path):
        (tmp_path / "pose.ent").write_text("END\n")
        files = _collect_pose_files(tmp_path)
        assert len(files) == 1

    def test_ignores_non_pdb_files(self, tmp_path):
        (tmp_path / "data.txt").write_text("hello\n")
        (tmp_path / "config.inp").write_text("PDBNAM receptor.pdb\n")
        files = _collect_pose_files(tmp_path)
        assert files == []

    def test_recurses_into_subdirectories(self, tmp_path):
        sub = tmp_path / "mode_1"
        sub.mkdir()
        (sub / "pose_1.pdb").write_text("END\n")
        (sub / "pose_2.pdb").write_text("END\n")
        (tmp_path / "top.pdb").write_text("END\n")
        files = _collect_pose_files(tmp_path)
        assert len(files) == 3

    def test_extension_case_insensitive(self, tmp_path):
        (tmp_path / "POSE.PDB").write_text("END\n")
        files = _collect_pose_files(tmp_path)
        assert len(files) == 1

    def test_returns_sorted_paths(self, tmp_path):
        for name in ["c.pdb", "a.pdb", "b.pdb"]:
            (tmp_path / name).write_text("END\n")
        files = _collect_pose_files(tmp_path)
        assert files == sorted(files)


# ===========================================================================
# _mode_temperature
# ===========================================================================

class TestModeTemperature:
    def test_returns_first_non_none(self):
        poses = [
            _make_pose(temperature=None),
            _make_pose(temperature=300.0),
            _make_pose(temperature=310.0),
        ]
        assert _mode_temperature(poses) == 300.0

    def test_returns_none_when_all_none(self):
        poses = [_make_pose(), _make_pose()]
        assert _mode_temperature(poses) is None

    def test_empty_list_returns_none(self):
        assert _mode_temperature([]) is None


# ===========================================================================
# _mode_metric
# ===========================================================================

class TestModeMetric:
    def test_returns_first_non_none(self):
        poses = [
            _make_pose(free_energy=None),
            _make_pose(free_energy=-9.5),
            _make_pose(free_energy=-8.0),
        ]
        assert _mode_metric(poses, "free_energy") == -9.5

    def test_returns_none_when_all_none(self):
        poses = [_make_pose(), _make_pose()]
        assert _mode_metric(poses, "free_energy") is None

    def test_all_fields(self):
        pose = _make_pose(
            enthalpy=-8.0, entropy=0.001, heat_capacity=0.05, std_energy=0.3
        )
        assert _mode_metric([pose], "enthalpy") == -8.0
        assert _mode_metric([pose], "entropy") == 0.001
        assert _mode_metric([pose], "heat_capacity") == 0.05
        assert _mode_metric([pose], "std_energy") == 0.3


# ===========================================================================
# _mode_frequency
# ===========================================================================

class TestModeFrequency:
    def test_reads_frequency_from_remarks(self):
        pose = _make_pose(remarks={"frequency": 15})
        assert _mode_frequency([pose]) == 15

    def test_reads_nposes_alias(self):
        pose = _make_pose(remarks={"nposes": 8})
        assert _mode_frequency([pose]) == 8

    def test_reads_population_alias(self):
        pose = _make_pose(remarks={"population": 4})
        assert _mode_frequency([pose]) == 4

    def test_reads_cluster_size_alias(self):
        pose = _make_pose(remarks={"cluster_size": 12})
        assert _mode_frequency([pose]) == 12

    def test_reads_size_alias(self):
        pose = _make_pose(remarks={"size": 6})
        assert _mode_frequency([pose]) == 6

    def test_falls_back_to_pose_count(self):
        poses = [_make_pose(), _make_pose(), _make_pose()]
        assert _mode_frequency(poses) == 3

    def test_empty_list_returns_none(self):
        assert _mode_frequency([]) is None

    def test_non_int_remark_is_ignored(self):
        # String values should be skipped; fallback to count
        poses = [_make_pose(remarks={"frequency": "many"}), _make_pose()]
        assert _mode_frequency(poses) == 2


# ===========================================================================
# _mode_metadata
# ===========================================================================

class TestModeMetadata:
    def test_returns_shared_constant_keys(self):
        poses = [
            _make_pose(remarks={"run_id": 42, "cf": -10.0}),
            _make_pose(remarks={"run_id": 42, "cf": -9.5}),
        ]
        meta = _mode_metadata(poses)
        assert meta["run_id"] == 42
        assert "cf" not in meta  # differs across poses

    def test_non_shared_keys_excluded(self):
        poses = [
            _make_pose(remarks={"a": 1, "b": 2}),
            _make_pose(remarks={"b": 2}),  # 'a' missing
        ]
        meta = _mode_metadata(poses)
        assert "a" not in meta
        assert meta["b"] == 2

    def test_empty_poses_returns_empty_dict(self):
        assert _mode_metadata([]) == {}

    def test_single_pose_all_keys_shared(self):
        pose = _make_pose(remarks={"x": 1, "y": 2})
        meta = _mode_metadata([pose])
        assert meta == {"x": 1, "y": 2}


# ===========================================================================
# _build_mode
# ===========================================================================

class TestBuildMode:
    def test_poses_sorted_by_rank_then_name(self):
        poses = [
            _make_pose(path="b.pdb", pose_rank=2),
            _make_pose(path="a.pdb", pose_rank=1),
            _make_pose(path="c.pdb", pose_rank=2),
        ]
        mode = _build_mode(1, poses)
        assert mode.poses[0].path.name == "a.pdb"
        assert mode.poses[1].path.name == "b.pdb"
        assert mode.poses[2].path.name == "c.pdb"

    def test_best_cf_is_most_negative(self):
        poses = [
            _make_pose(cf=-8.0),
            _make_pose(cf=-12.0),
            _make_pose(cf=-10.0),
        ]
        mode = _build_mode(1, poses)
        assert mode.best_cf == pytest.approx(-12.0)

    def test_best_cf_none_when_no_scored_poses(self):
        poses = [_make_pose(), _make_pose()]
        mode = _build_mode(1, poses)
        assert mode.best_cf is None

    def test_mode_id_set(self):
        mode = _build_mode(7, [_make_pose(cf=-5.0)])
        assert mode.mode_id == 7

    def test_free_energy_taken_from_first_available(self):
        poses = [_make_pose(free_energy=None), _make_pose(free_energy=-9.0)]
        mode = _build_mode(1, poses)
        assert mode.free_energy == pytest.approx(-9.0)

    def test_temperature_taken_from_first_available(self):
        poses = [_make_pose(temperature=None), _make_pose(temperature=300.0)]
        mode = _build_mode(1, poses)
        assert mode.temperature == pytest.approx(300.0)

    def test_frequency_defaults_to_pose_count(self):
        poses = [_make_pose(), _make_pose(), _make_pose()]
        mode = _build_mode(1, poses)
        assert mode.frequency == 3

    def test_n_poses(self):
        poses = [_make_pose(), _make_pose()]
        mode = _build_mode(1, poses)
        assert mode.n_poses == 2


# ===========================================================================
# load_results – error paths
# ===========================================================================

class TestLoadResultsErrors:
    def test_nonexistent_path_raises_file_not_found(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="does not exist"):
            load_results(tmp_path / "nonexistent_dir")

    def test_file_path_raises_not_a_directory(self, tmp_path):
        f = tmp_path / "file.txt"
        f.write_text("hello\n")
        with pytest.raises(NotADirectoryError, match="directory"):
            load_results(f)

    def test_empty_directory_returns_empty_result(self, tmp_path):
        result = load_results(tmp_path)
        assert result.n_modes == 0
        assert result.metadata["n_pose_files"] == 0

    def test_directory_with_no_pdb_files_returns_empty_result(self, tmp_path):
        (tmp_path / "data.txt").write_text("hello\n")
        result = load_results(tmp_path)
        assert result.n_modes == 0
        assert result.metadata["n_pose_files"] == 0


# ===========================================================================
# load_results – grouping and metadata
# ===========================================================================

class TestLoadResultsGrouping:
    def test_multiple_modes_with_multiple_poses(self, tmp_path):
        for mode in (1, 2, 3):
            for pose in (1, 2):
                _write_pdb(
                    tmp_path / f"mode_{mode}_pose_{pose}.pdb",
                    [
                        f"binding_mode = {mode}",
                        f"pose_rank = {pose}",
                        f"CF = {-10.0 - mode}",
                        "temperature = 300.0",
                    ],
                )
        result = load_results(tmp_path)
        assert result.n_modes == 3
        for bm in result.binding_modes:
            assert bm.n_poses == 2

    def test_modes_ordered_by_id(self, tmp_path):
        for mode in (3, 1, 2):
            _write_pdb(
                tmp_path / f"mode_{mode}_pose_1.pdb",
                [f"binding_mode = {mode}", "CF = -10.0", "temperature = 300.0"],
            )
        result = load_results(tmp_path)
        ids = [m.mode_id for m in result.binding_modes]
        assert ids == sorted(ids)

    def test_metadata_contains_n_pose_files(self, tmp_path):
        _write_pdb(tmp_path / "mode_1_pose_1.pdb", ["CF = -5.0"])
        _write_pdb(tmp_path / "mode_1_pose_2.pdb", ["CF = -4.0"])
        result = load_results(tmp_path)
        assert result.metadata["n_pose_files"] == 2

    def test_temperature_propagated_to_docking_result(self, tmp_path):
        _write_pdb(
            tmp_path / "mode_1_pose_1.pdb",
            ["CF = -5.0", "temperature = 310.0"],
        )
        result = load_results(tmp_path)
        assert result.temperature == pytest.approx(310.0)

    def test_temperature_none_when_no_remark(self, tmp_path):
        _write_pdb(tmp_path / "mode_1_pose_1.pdb", ["CF = -5.0"])
        result = load_results(tmp_path)
        assert result.temperature is None

    def test_str_path_accepted(self, tmp_path):
        _write_pdb(tmp_path / "mode_1_pose_1.pdb", ["CF = -5.0"])
        result = load_results(str(tmp_path))
        assert result.n_modes == 1
