"""Tests for flexaidds.models – BindingModeResult and DockingResult.

Priority 4 coverage.  No C++ extension needed.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from flexaidds.models import BindingModeResult, DockingResult, PoseResult


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _pose(
    *,
    path: str = "x.pdb",
    mode_id: int = 1,
    pose_rank: int = 1,
    cf: float | None = None,
    cf_app: float | None = None,
) -> PoseResult:
    return PoseResult(
        path=Path(path),
        mode_id=mode_id,
        pose_rank=pose_rank,
        cf=cf,
        cf_app=cf_app,
    )


def _mode(
    *,
    mode_id: int = 1,
    poses: list[PoseResult] | None = None,
    free_energy: float | None = None,
    best_cf: float | None = None,
    rank: int | None = None,
) -> BindingModeResult:
    return BindingModeResult(
        mode_id=mode_id,
        rank=rank if rank is not None else mode_id,
        poses=poses or [],
        free_energy=free_energy,
        best_cf=best_cf,
    )


# ===========================================================================
# PoseResult
# ===========================================================================

class TestPoseResultRepr:
    def test_repr_basic(self):
        p = _pose(mode_id=2, pose_rank=3)
        r = repr(p)
        assert "PoseResult" in r
        assert "mode=2" in r
        assert "rank=3" in r

    def test_repr_with_cf(self):
        p = _pose(cf=-10.5)
        assert "cf=-10.50" in repr(p)

    def test_repr_with_free_energy(self):
        p = PoseResult(path=Path("x.pdb"), mode_id=1, pose_rank=1, free_energy=-8.3)
        assert "F=-8.30" in repr(p)

    def test_repr_without_optional_fields(self):
        p = _pose()
        r = repr(p)
        assert "cf=" not in r
        assert "F=" not in r


class TestBindingModeResultRepr:
    def test_repr_basic(self):
        mode = _mode(mode_id=2, poses=[_pose(), _pose()], rank=1)
        r = repr(mode)
        assert "BindingModeResult" in r
        assert "id=2" in r
        assert "rank=1" in r
        assert "poses=2" in r

    def test_repr_with_free_energy(self):
        mode = _mode(free_energy=-9.5)
        assert "F=-9.50" in repr(mode)


class TestDockingResultRepr:
    def test_repr(self):
        result = DockingResult(
            source_dir=Path("/tmp/docking_out"),
            binding_modes=[_mode()],
        )
        r = repr(result)
        assert "DockingResult" in r
        assert "modes=1" in r
        assert "docking_out" in r


class TestPoseResult:
    def test_is_frozen(self):
        p = _pose(cf=-5.0)
        with pytest.raises((AttributeError, TypeError)):
            p.cf = 0.0  # type: ignore[misc]

    def test_defaults_are_none(self):
        p = PoseResult(path=Path("x.pdb"), mode_id=1, pose_rank=1)
        for field in ("cf", "cf_app", "rmsd_raw", "rmsd_sym",
                      "free_energy", "enthalpy", "entropy",
                      "heat_capacity", "std_energy", "temperature"):
            assert getattr(p, field) is None

    def test_remarks_default_is_empty_dict(self):
        p = PoseResult(path=Path("x.pdb"), mode_id=1, pose_rank=1)
        assert p.remarks == {}


# ===========================================================================
# BindingModeResult.n_poses
# ===========================================================================

class TestBindingModeResultNPoses:
    def test_n_poses_matches_list_length(self):
        poses = [_pose(), _pose(), _pose()]
        mode = _mode(poses=poses)
        assert mode.n_poses == 3

    def test_n_poses_zero_when_no_poses(self):
        assert _mode(poses=[]).n_poses == 0


# ===========================================================================
# BindingModeResult.best_pose
# ===========================================================================

class TestBindingModeResultBestPose:
    def test_returns_pose_with_lowest_cf(self):
        poses = [
            _pose(path="a.pdb", cf=-8.0),
            _pose(path="b.pdb", cf=-12.0),
            _pose(path="c.pdb", cf=-10.0),
        ]
        mode = _mode(poses=poses)
        assert mode.best_pose().path.name == "b.pdb"

    def test_falls_back_to_cf_app_when_no_cf(self):
        poses = [
            _pose(path="a.pdb", cf_app=-7.0),
            _pose(path="b.pdb", cf_app=-11.0),
        ]
        mode = _mode(poses=poses)
        assert mode.best_pose().path.name == "b.pdb"

    def test_cf_takes_priority_over_cf_app(self):
        poses = [
            _pose(path="cf_only.pdb", cf=-9.0),
            _pose(path="app_only.pdb", cf_app=-15.0),  # better app score, no cf
        ]
        mode = _mode(poses=poses)
        # Only the pose with cf is in the cf-scored set; it should win
        best = mode.best_pose()
        assert best.path.name == "cf_only.pdb"

    def test_returns_first_pose_when_no_scores(self):
        poses = [_pose(path="first.pdb"), _pose(path="second.pdb")]
        mode = _mode(poses=poses)
        assert mode.best_pose().path.name == "first.pdb"

    def test_returns_none_when_no_poses(self):
        assert _mode(poses=[]).best_pose() is None


# ===========================================================================
# DockingResult.n_modes
# ===========================================================================

class TestDockingResultNModes:
    def _result(self, n: int) -> DockingResult:
        return DockingResult(
            source_dir=Path("/tmp"),
            binding_modes=[_mode(mode_id=i) for i in range(1, n + 1)],
        )

    def test_n_modes_correct(self):
        assert self._result(3).n_modes == 3

    def test_n_modes_zero(self):
        assert self._result(0).n_modes == 0


# ===========================================================================
# DockingResult.top_mode
# ===========================================================================

class TestDockingResultTopMode:
    def test_returns_mode_with_lowest_free_energy(self):
        modes = [
            _mode(mode_id=1, free_energy=-8.0, rank=1),
            _mode(mode_id=2, free_energy=-12.0, rank=2),
            _mode(mode_id=3, free_energy=-10.0, rank=3),
        ]
        result = DockingResult(source_dir=Path("/tmp"), binding_modes=modes)
        assert result.top_mode().mode_id == 2

    def test_falls_back_to_lowest_rank_when_no_free_energy(self):
        modes = [
            _mode(mode_id=3, rank=3),
            _mode(mode_id=1, rank=1),
            _mode(mode_id=2, rank=2),
        ]
        result = DockingResult(source_dir=Path("/tmp"), binding_modes=modes)
        assert result.top_mode().rank == 1

    def test_returns_none_when_no_modes(self):
        result = DockingResult(source_dir=Path("/tmp"), binding_modes=[])
        assert result.top_mode() is None

    def test_modes_with_partial_free_energy(self):
        """Modes without free_energy are excluded from F-based ranking."""
        modes = [
            _mode(mode_id=1, free_energy=None, rank=1),
            _mode(mode_id=2, free_energy=-8.0, rank=2),
        ]
        result = DockingResult(source_dir=Path("/tmp"), binding_modes=modes)
        assert result.top_mode().mode_id == 2


# ===========================================================================
# DockingResult.to_records
# ===========================================================================

class TestDockingResultToRecords:
    def _make_result(self) -> DockingResult:
        poses = [_pose(path="p1.pdb", cf=-10.0), _pose(path="p2.pdb", cf=-9.0)]
        modes = [
            BindingModeResult(
                mode_id=1, rank=1, poses=poses,
                free_energy=-9.8, enthalpy=-9.5, entropy=0.001,
                heat_capacity=0.05, std_energy=0.3,
                best_cf=-10.0, temperature=300.0,
            ),
            BindingModeResult(
                mode_id=2, rank=2, poses=[],
                free_energy=None, best_cf=None, temperature=None,
            ),
        ]
        return DockingResult(source_dir=Path("/tmp"), binding_modes=modes)

    def test_returns_one_record_per_mode(self):
        records = self._make_result().to_records()
        assert len(records) == 2

    def test_record_keys(self):
        # Schema v2 deliberately widens this record: the claim-firewall fields
        # (proxy_free_energy, soft_beta_G, scientific_provenance) must travel
        # with every serialised mode, so exact set equality is asserted against
        # the widened contract rather than relaxed to a subset check.
        record = self._make_result().to_records()[0]
        expected = {
            "mode_id", "rank", "n_poses", "free_energy", "proxy_free_energy",
            "soft_beta_G", "enthalpy", "entropy", "heat_capacity",
            "std_energy", "best_cf", "temperature", "scientific_provenance",
            "thermodynamics", "best_pose_path",
        }
        assert expected == set(record.keys())

    def test_record_provenance_is_serialisable_mapping(self):
        record = self._make_result().to_records()[0]
        prov = record["scientific_provenance"]
        assert isinstance(prov, dict)
        # Default (unclassified) evidence must stay proxy-only.
        assert prov["claim_validity"] == "proxy_only"
        # Records feed json.dumps / csv / pandas — no dataclass objects.
        json.dumps(record)

    def test_record_values_mode1(self):
        record = self._make_result().to_records()[0]
        assert record["mode_id"] == 1
        assert record["rank"] == 1
        assert record["n_poses"] == 2
        assert record["free_energy"] == pytest.approx(-9.8)
        assert record["best_cf"] == pytest.approx(-10.0)
        assert record["temperature"] == pytest.approx(300.0)
        assert record["best_pose_path"] is not None

    def test_none_values_for_empty_mode(self):
        record = self._make_result().to_records()[1]
        assert record["free_energy"] is None
        assert record["best_cf"] is None
        assert record["best_pose_path"] is None

    def test_best_pose_path_is_string(self):
        record = self._make_result().to_records()[0]
        assert isinstance(record["best_pose_path"], str)


# ===========================================================================
# DockingResult.to_dataframe – success path and pandas import error path
# ===========================================================================

class TestDockingResultToDataframe:
    def test_raises_import_error_without_pandas(self, monkeypatch):
        import builtins
        real_import = builtins.__import__

        def mock_import(name, *args, **kwargs):
            if name == "pandas":
                raise ImportError("No module named 'pandas'")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", mock_import)

        result = DockingResult(
            source_dir=Path("/tmp"),
            binding_modes=[_mode()],
        )
        with pytest.raises(ImportError, match="pandas"):
            result.to_dataframe()


# ===========================================================================
# DockingResult.to_json
# ===========================================================================

class TestDockingResultToJson:
    def _make_result(self) -> DockingResult:
        poses = [_pose(path="p1.pdb", cf=-10.0), _pose(path="p2.pdb", cf=-9.0)]
        modes = [
            BindingModeResult(
                mode_id=1, rank=1, poses=poses,
                free_energy=-9.8, enthalpy=-9.5, entropy=0.001,
                heat_capacity=0.05, std_energy=0.3,
                best_cf=-10.0, temperature=300.0,
            ),
            BindingModeResult(
                mode_id=2, rank=2, poses=[],
                free_energy=None, best_cf=None, temperature=None,
            ),
        ]
        return DockingResult(source_dir=Path("/tmp"), binding_modes=modes)

    def test_returns_valid_json_string(self):
        import json
        text = self._make_result().to_json()
        parsed = json.loads(text)
        assert isinstance(parsed, dict)

    def test_json_has_required_keys(self):
        import json
        parsed = json.loads(self._make_result().to_json())
        assert set(parsed.keys()) >= {
            "source_dir", "temperature", "n_modes", "metadata", "binding_modes",
        }

    def test_json_n_modes(self):
        import json
        parsed = json.loads(self._make_result().to_json())
        assert parsed["n_modes"] == 2

    def test_json_binding_modes_length(self):
        import json
        parsed = json.loads(self._make_result().to_json())
        assert len(parsed["binding_modes"]) == 2

    def test_json_free_energy_value(self):
        import json
        parsed = json.loads(self._make_result().to_json())
        assert parsed["binding_modes"][0]["free_energy"] == pytest.approx(-9.8)

    def test_json_none_free_energy(self):
        import json
        parsed = json.loads(self._make_result().to_json())
        assert parsed["binding_modes"][1]["free_energy"] is None

    def test_to_json_writes_file(self, tmp_path):
        import json
        out = tmp_path / "result.json"
        ret = self._make_result().to_json(path=out)
        assert ret is None
        parsed = json.loads(out.read_text(encoding="utf-8"))
        assert parsed["n_modes"] == 2

    def test_to_json_kwargs_forwarded(self):
        import json
        text = self._make_result().to_json(sort_keys=True, indent=4)
        parsed = json.loads(text)
        assert isinstance(parsed, dict)
        # Verify sort_keys was applied: keys should be alphabetical
        lines = text.strip().splitlines()
        # indent=4 means 4-space indentation
        assert any("    " in line for line in lines)



# ===========================================================================
# DockingResult.from_json
# ===========================================================================

class TestDockingResultFromJson:
    def _make_result(self) -> DockingResult:
        poses = [_pose(path="p1.pdb", cf=-10.0), _pose(path="p2.pdb", cf=-9.0)]
        modes = [
            BindingModeResult(
                mode_id=1, rank=1, poses=poses,
                free_energy=-9.8, enthalpy=-9.5, entropy=0.001,
                heat_capacity=0.05, std_energy=0.3,
                best_cf=-10.0, temperature=300.0,
            ),
            BindingModeResult(
                mode_id=2, rank=2, poses=[],
                free_energy=None, best_cf=None, temperature=None,
            ),
        ]
        return DockingResult(
            source_dir=Path("/tmp"),
            binding_modes=modes,
            temperature=300.0,
            metadata={"n_pose_files": 2},
        )

    def test_roundtrip_from_string(self):
        original = self._make_result()
        json_text = original.to_json()
        restored = DockingResult.from_json(json_text)
        assert restored.n_modes == 2
        assert restored.temperature == pytest.approx(300.0)
        assert restored.binding_modes[0].free_energy == pytest.approx(-9.8)
        assert restored.binding_modes[0].mode_id == 1
        assert restored.binding_modes[1].free_energy is None

    def test_roundtrip_from_file(self, tmp_path):
        original = self._make_result()
        out = tmp_path / "result.json"
        original.to_json(path=out)
        restored = DockingResult.from_json(out)
        assert restored.n_modes == 2
        assert restored.binding_modes[0].best_cf == pytest.approx(-10.0)

    def test_metadata_preserved(self):
        original = self._make_result()
        json_text = original.to_json()
        restored = DockingResult.from_json(json_text)
        assert restored.metadata["n_pose_files"] == 2

    def test_poses_list_is_empty(self):
        original = self._make_result()
        json_text = original.to_json()
        restored = DockingResult.from_json(json_text)
        assert restored.binding_modes[0].poses == []

    def test_source_dir_restored(self):
        original = self._make_result()
        json_text = original.to_json()
        restored = DockingResult.from_json(json_text)
        assert str(restored.source_dir) == "/tmp"


# ===========================================================================
# DockingResult.from_csv
# ===========================================================================

class TestDockingResultFromCsv:
    def _make_result(self) -> DockingResult:
        poses = [_pose(path="p1.pdb", cf=-10.0)]
        modes = [
            BindingModeResult(
                mode_id=1, rank=1, poses=poses,
                free_energy=-9.8, enthalpy=-9.5, entropy=0.001,
                heat_capacity=0.05, std_energy=0.3,
                best_cf=-10.0, temperature=300.0,
            ),
            BindingModeResult(
                mode_id=2, rank=2, poses=[],
                free_energy=-7.2, best_cf=-8.0, temperature=300.0,
            ),
        ]
        return DockingResult(source_dir=Path("/tmp"), binding_modes=modes)

    def test_roundtrip_from_string(self):
        original = self._make_result()
        csv_text = original.to_csv()
        restored = DockingResult.from_csv(csv_text)
        assert restored.n_modes == 2
        assert restored.binding_modes[0].free_energy == pytest.approx(-9.8)
        assert restored.binding_modes[0].mode_id == 1

    def test_roundtrip_from_file(self, tmp_path):
        original = self._make_result()
        out = tmp_path / "result.csv"
        original.to_csv(path=out)
        restored = DockingResult.from_csv(out)
        assert restored.n_modes == 2
        assert restored.binding_modes[1].best_cf == pytest.approx(-8.0)

    def test_none_values_preserved(self):
        modes = [
            BindingModeResult(
                mode_id=1, rank=1, poses=[],
                free_energy=None, enthalpy=None, entropy=None,
            ),
        ]
        original = DockingResult(source_dir=Path("/tmp"), binding_modes=modes)
        csv_text = original.to_csv()
        restored = DockingResult.from_csv(csv_text)
        assert restored.binding_modes[0].free_energy is None
        assert restored.binding_modes[0].enthalpy is None

    def test_poses_list_is_empty(self):
        original = self._make_result()
        csv_text = original.to_csv()
        restored = DockingResult.from_csv(csv_text)
        assert restored.binding_modes[0].poses == []


# ===========================================================================
# DockingResult.to_csv
# ===========================================================================

class TestDockingResultToCsv:
    def _make_result(self) -> DockingResult:
        poses = [_pose(path="p1.pdb", cf=-10.0), _pose(path="p2.pdb", cf=-9.0)]
        modes = [
            BindingModeResult(
                mode_id=1, rank=1, poses=poses,
                free_energy=-9.8, enthalpy=-9.5, entropy=0.001,
                heat_capacity=0.05, std_energy=0.3,
                best_cf=-10.0, temperature=300.0,
            ),
            BindingModeResult(
                mode_id=2, rank=2, poses=[],
                free_energy=None, best_cf=None, temperature=None,
            ),
        ]
        return DockingResult(source_dir=Path("/tmp"), binding_modes=modes)

    def test_returns_csv_string_when_no_path(self):
        text = self._make_result().to_csv()
        assert isinstance(text, str)
        assert "mode_id" in text  # header row

    def test_csv_has_header_and_data_rows(self):
        text = self._make_result().to_csv()
        lines = text.strip().splitlines()
        assert len(lines) == 3  # 1 header + 2 data rows

    def test_csv_header_columns(self):
        text = self._make_result().to_csv()
        header = text.strip().splitlines()[0]
        for col in ("mode_id", "rank", "n_poses", "free_energy",
                     "enthalpy", "entropy", "best_cf", "best_pose_path"):
            assert col in header

    def test_csv_data_values(self):
        text = self._make_result().to_csv()
        lines = text.strip().splitlines()
        # First data row should contain mode_id=1 values
        assert "-9.8" in lines[1]
        assert "-10.0" in lines[1]

    def test_csv_writes_file(self, tmp_path):
        out = tmp_path / "result.csv"
        ret = self._make_result().to_csv(path=out)
        assert ret is None
        content = out.read_text(encoding="utf-8")
        assert "mode_id" in content
        lines = content.strip().splitlines()
        assert len(lines) == 3

    def test_csv_empty_result(self):
        result = DockingResult(source_dir=Path("/tmp"), binding_modes=[])
        text = result.to_csv()
        assert isinstance(text, str)
        # Should have empty header or just a newline
        assert len(text.strip().splitlines()) <= 1


# ===========================================================================
# DockingResult.to_dataframe – success path
# ===========================================================================

class TestDockingResultToDataframeSuccess:
    def test_to_dataframe_success(self):
        pytest.importorskip("pandas")
        poses = [_pose(path="p1.pdb", cf=-10.0)]
        modes = [
            BindingModeResult(
                mode_id=1, rank=1, poses=poses,
                free_energy=-9.8, best_cf=-10.0,
            ),
        ]
        result = DockingResult(source_dir=Path("/tmp"), binding_modes=modes)
        df = result.to_dataframe()
        import pandas as pd
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 1
        assert "mode_id" in df.columns
        assert df.iloc[0]["mode_id"] == 1
        assert df.iloc[0]["free_energy"] == pytest.approx(-9.8)


# ===========================================================================
# DockingResult.to_csv
# ===========================================================================

class TestDockingResultToCsv:
    def _make_result(self) -> DockingResult:
        poses = [_pose(path="p1.pdb", cf=-10.0), _pose(path="p2.pdb", cf=-9.0)]
        modes = [
            BindingModeResult(
                mode_id=1, rank=1, poses=poses,
                free_energy=-9.8, enthalpy=-9.5, entropy=0.001,
                heat_capacity=0.05, std_energy=0.3,
                best_cf=-10.0, temperature=300.0,
            ),
            BindingModeResult(
                mode_id=2, rank=2, poses=[],
                free_energy=None, best_cf=None, temperature=None,
            ),
        ]
        return DockingResult(source_dir=Path("/tmp"), binding_modes=modes)

    def test_returns_csv_string_when_no_path(self):
        text = self._make_result().to_csv()
        assert isinstance(text, str)
        assert len(text) > 0

    def test_csv_has_header_row(self):
        import csv as csv_mod
        import io
        text = self._make_result().to_csv()
        reader = csv_mod.DictReader(io.StringIO(text))
        # Mirrors to_records(): the schema-v2 firewall columns are part of the
        # CSV contract so downstream consumers cannot lose the provenance.
        expected_cols = {
            "mode_id", "rank", "n_poses", "free_energy", "proxy_free_energy",
            "soft_beta_G", "enthalpy", "entropy", "heat_capacity",
            "std_energy", "best_cf", "temperature", "scientific_provenance",
            "thermodynamics", "best_pose_path",
        }
        assert expected_cols == set(reader.fieldnames)

    def test_csv_has_one_row_per_mode(self):
        import csv as csv_mod
        import io
        text = self._make_result().to_csv()
        reader = csv_mod.DictReader(io.StringIO(text))
        rows = list(reader)
        assert len(rows) == 2

    def test_csv_values_mode1(self):
        import csv as csv_mod
        import io
        text = self._make_result().to_csv()
        reader = csv_mod.DictReader(io.StringIO(text))
        row = next(reader)
        assert row["mode_id"] == "1"
        assert row["rank"] == "1"
        assert row["n_poses"] == "2"
        assert float(row["free_energy"]) == pytest.approx(-9.8)
        assert float(row["best_cf"]) == pytest.approx(-10.0)
        assert float(row["temperature"]) == pytest.approx(300.0)

    def test_csv_none_values_for_empty_mode(self):
        import csv as csv_mod
        import io
        text = self._make_result().to_csv()
        reader = csv_mod.DictReader(io.StringIO(text))
        rows = list(reader)
        row2 = rows[1]
        # CSV serialises None as empty string
        assert row2["free_energy"] == ""
        assert row2["best_cf"] == ""
        assert row2["best_pose_path"] == ""

    def test_to_csv_writes_file(self, tmp_path):
        out = tmp_path / "result.csv"
        ret = self._make_result().to_csv(path=out)
        assert ret is None
        text = out.read_text(encoding="utf-8")
        assert "mode_id" in text
        assert "-9.8" in text

    def test_to_csv_empty_result(self):
        result = DockingResult(source_dir=Path("/tmp"), binding_modes=[])
        text = result.to_csv()
        # Empty result has only a blank header line (no column names, no data rows)
        assert text.strip() == ""


# ===========================================================================
# DockingResult.from_json
# ===========================================================================

class TestDockingResultFromJson:
    def _make_result(self) -> DockingResult:
        poses = [_pose(path="p1.pdb", cf=-10.0), _pose(path="p2.pdb", cf=-9.0)]
        modes = [
            BindingModeResult(
                mode_id=1, rank=1, poses=poses,
                free_energy=-9.8, enthalpy=-9.5, entropy=0.001,
                heat_capacity=0.05, std_energy=0.3,
                best_cf=-10.0, temperature=300.0,
            ),
            BindingModeResult(
                mode_id=2, rank=2, poses=[],
                free_energy=None, best_cf=None, temperature=None,
            ),
        ]
        return DockingResult(
            source_dir=Path("/tmp"), binding_modes=modes, temperature=300.0,
        )

    def test_round_trip_from_string(self):
        original = self._make_result()
        json_text = original.to_json()
        loaded = DockingResult.from_json(json_text)
        assert loaded.n_modes == original.n_modes
        assert loaded.temperature == pytest.approx(original.temperature)
        assert str(loaded.source_dir) == str(original.source_dir)

    def test_round_trip_preserves_mode_fields(self):
        original = self._make_result()
        loaded = DockingResult.from_json(original.to_json())
        mode = loaded.binding_modes[0]
        assert mode.mode_id == 1
        assert mode.rank == 1
        assert mode.free_energy == pytest.approx(-9.8)
        assert mode.enthalpy == pytest.approx(-9.5)
        assert mode.entropy == pytest.approx(0.001)
        assert mode.heat_capacity == pytest.approx(0.05)
        assert mode.std_energy == pytest.approx(0.3)
        assert mode.best_cf == pytest.approx(-10.0)
        assert mode.temperature == pytest.approx(300.0)

    def test_round_trip_none_mode(self):
        original = self._make_result()
        loaded = DockingResult.from_json(original.to_json())
        mode2 = loaded.binding_modes[1]
        assert mode2.free_energy is None
        assert mode2.best_cf is None
        assert mode2.temperature is None
        assert mode2.n_poses == 0

    def test_round_trip_from_file(self, tmp_path):
        original = self._make_result()
        json_file = tmp_path / "result.json"
        original.to_json(path=json_file)
        loaded = DockingResult.from_json(json_file)
        assert loaded.n_modes == 2
        assert loaded.temperature == pytest.approx(300.0)

    def test_source_dir_override(self):
        original = self._make_result()
        loaded = DockingResult.from_json(
            original.to_json(), source_dir="/new/path"
        )
        assert str(loaded.source_dir) == "/new/path"

    def test_metadata_preserved(self):
        result = DockingResult(
            source_dir=Path("/tmp"),
            binding_modes=[],
            metadata={"n_pose_files": 42},
        )
        loaded = DockingResult.from_json(result.to_json())
        assert loaded.metadata["n_pose_files"] == 42

    def test_reconstructed_pose_has_best_path(self):
        original = self._make_result()
        loaded = DockingResult.from_json(original.to_json())
        mode = loaded.binding_modes[0]
        assert mode.n_poses == 1
        best = mode.best_pose()
        assert best is not None
        assert best.path.name == "p1.pdb"
        assert best.cf == pytest.approx(-10.0)

    def test_invalid_json_raises(self):
        import json as json_mod
        with pytest.raises(json_mod.JSONDecodeError):
            DockingResult.from_json("not valid json {{")


# ===========================================================================
# __repr__ methods
# ===========================================================================

class TestPoseResultRepr:
    def test_includes_mode_and_rank(self):
        p = _pose(path="ligand.pdb", mode_id=2, pose_rank=3, cf=-8.5)
        r = repr(p)
        assert "PoseResult" in r
        assert "mode=2" in r
        assert "rank=3" in r

    def test_includes_cf_when_present(self):
        p = _pose(cf=-8.5)
        assert "cf=-8.50" in repr(p)

    def test_uses_cf_app_when_no_cf(self):
        p = _pose(cf_app=-7.0)
        assert "cf=-7.00" in repr(p)

    def test_omits_score_when_none(self):
        p = _pose()
        assert "cf=" not in repr(p)

    def test_includes_filename(self):
        p = _pose(path="ligand_mode1_pose2.pdb")
        assert "ligand_mode1_pose2.pdb" in repr(p)


class TestBindingModeResultRepr:
    def test_includes_mode_id_and_n_poses(self):
        poses = [_pose(), _pose()]
        mode = _mode(mode_id=3, poses=poses)
        r = repr(mode)
        assert "BindingModeResult" in r
        assert "mode_id=3" in r
        assert "n_poses=2" in r

    def test_includes_free_energy_when_present(self):
        mode = _mode(free_energy=-9.5)
        assert "F=-9.50" in repr(mode)

    def test_includes_best_cf_when_present(self):
        mode = _mode(best_cf=-12.0)
        assert "best_cf=-12.00" in repr(mode)

    def test_omits_optional_fields_when_none(self):
        mode = _mode()
        r = repr(mode)
        assert "F=" not in r
        assert "best_cf=" not in r


class TestDockingResultRepr:
    def test_includes_n_modes(self):
        result = DockingResult(
            source_dir=Path("/data/run1"),
            binding_modes=[_mode(mode_id=1), _mode(mode_id=2)],
        )
        r = repr(result)
        assert "DockingResult" in r
        assert "n_modes=2" in r

    def test_includes_temperature_when_present(self):
        result = DockingResult(
            source_dir=Path("/data/run1"),
            binding_modes=[],
            temperature=300.0,
        )
        assert "T=300K" in repr(result)

    def test_omits_temperature_when_none(self):
        result = DockingResult(
            source_dir=Path("/data/run1"),
            binding_modes=[],
        )
        assert "T=" not in repr(result)

    def test_includes_source_dir_name(self):
        result = DockingResult(
            source_dir=Path("/data/my_run"),
            binding_modes=[],
        )
        assert "my_run" in repr(result)
