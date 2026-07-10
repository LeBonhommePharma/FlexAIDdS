"""Tests for Phase 3 PyMOL plugin features.

Tests cover:
- Entropy heatmap spatial entropy computation (pure Python & NumPy)
- CGO and pseudoatom renderers
- Mode animation coordinate interpolation and Kabsch alignment
- ITC comparison data parsing (single & multi-ligand CSV)
- Interactive docking configuration generation and async mode
- Docking cancellation callback

All tests mock PyMOL to avoid requiring it as a dependency.
"""

import csv
import math
import os
import sys
import tempfile
import types
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

try:
    import numpy as np
    _HAS_NUMPY = True
except ImportError:
    _HAS_NUMPY = False


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def tmp_dir():
    """Create a temporary directory for test files."""
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        yield Path(d)


@pytest.fixture(autouse=True)
def mock_pymol():
    """Mock PyMOL for all tests in this module."""
    pymol_mod = types.ModuleType("pymol")
    pymol_mod.cmd = MagicMock()
    pymol_mod.stored = MagicMock()
    with patch.dict(sys.modules, {"pymol": pymol_mod}):
        yield pymol_mod


@pytest.fixture
def sample_pdb(tmp_dir):
    """Create a minimal PDB file with ATOM records."""
    pdb_path = tmp_dir / "pose1.pdb"
    pdb_path.write_text(
        "REMARK Binding Mode:1 Best CF in Binding Mode:-5.0 Binding Mode Frequency:3\n"
        "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       3.000   4.000   5.000  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   1       4.000   5.000   6.000  1.00  0.00           O\n"
        "ATOM      5  H   ALA A   1       1.500   2.500   3.500  1.00  0.00           H\n"
        "END\n"
    )
    return pdb_path


@pytest.fixture
def sample_pdb2(tmp_dir):
    """Create a second PDB file with different coordinates."""
    pdb_path = tmp_dir / "pose2.pdb"
    pdb_path.write_text(
        "REMARK Binding Mode:1 Best CF in Binding Mode:-4.5 Binding Mode Frequency:3\n"
        "ATOM      1  N   ALA A   1       5.000   6.000   7.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       6.000   7.000   8.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       7.000   8.000   9.000  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   1       8.000   9.000  10.000  1.00  0.00           O\n"
        "ATOM      5  H   ALA A   1       5.500   6.500   7.500  1.00  0.00           H\n"
        "END\n"
    )
    return pdb_path


@pytest.fixture
def sample_itc_csv(tmp_dir):
    """Create a sample ITC experimental data CSV."""
    csv_path = tmp_dir / "itc_data.csv"
    csv_path.write_text("dG,dH,TdS\n-8.5,-12.3,3.8\n")
    return csv_path


# ---------------------------------------------------------------------------
# Entropy Heatmap Tests
# ---------------------------------------------------------------------------

class TestEntropyHeatmapHelpers:
    """Test entropy heatmap computation."""

    def test_read_pose_coords_filters_hydrogens(self, sample_pdb):
        """Hydrogen atoms should be excluded from coordinate extraction."""
        from pymol_plugin.entropy_heatmap import _read_pose_coords
        coords = _read_pose_coords(str(sample_pdb), ligand_only=False)
        # Should have 4 heavy atoms (N, CA, C, O), not H
        assert len(coords) == 4

    def test_read_pose_coords_returns_xyz(self, sample_pdb):
        """Coordinates should be (x, y, z) tuples."""
        from pymol_plugin.entropy_heatmap import _read_pose_coords
        coords = _read_pose_coords(str(sample_pdb), ligand_only=False)
        assert coords[0] == pytest.approx((1.0, 2.0, 3.0))
        assert coords[1] == pytest.approx((2.0, 3.0, 4.0))

    def test_read_pose_coords_missing_file(self):
        """Non-existent file should return empty list."""
        from pymol_plugin.entropy_heatmap import _read_pose_coords
        coords = _read_pose_coords("/nonexistent/file.pdb")
        assert coords == []

    def test_compute_grid_bounds(self, sample_pdb):
        """Grid bounds should encompass all coordinates with padding."""
        from pymol_plugin.entropy_heatmap import (
            _read_pose_coords, _compute_grid_bounds
        )
        coords = _read_pose_coords(str(sample_pdb), ligand_only=False)
        lo, hi = _compute_grid_bounds(coords, padding=2.0)
        assert lo[0] < 1.0  # x_min - padding
        assert hi[0] > 4.0  # x_max + padding

    def test_entropy_color_low(self):
        """Low entropy should be burgundy red."""
        from pymol_plugin.entropy_heatmap import _entropy_color
        r, g, b = _entropy_color(0.0)
        assert r == pytest.approx(0.502)
        assert g == pytest.approx(0.0)
        assert b == pytest.approx(0.125)

    def test_entropy_color_high(self):
        """High entropy should be purple blue."""
        from pymol_plugin.entropy_heatmap import _entropy_color
        r, g, b = _entropy_color(1.0)
        assert r == pytest.approx(0.294)
        assert g == pytest.approx(0.0)
        assert b == pytest.approx(0.510)

    def test_entropy_color_mid(self):
        """Mid entropy should be white."""
        from pymol_plugin.entropy_heatmap import _entropy_color
        r, g, b = _entropy_color(0.5)
        assert r == pytest.approx(1.0)
        assert g == pytest.approx(1.0)
        assert b == pytest.approx(1.0)

    def test_ligand_only_fallback(self, sample_pdb):
        """ligand_only=True falls back to all atoms when no HETATM found."""
        from pymol_plugin.entropy_heatmap import _read_pose_coords
        # sample_pdb has only ATOM records, so ligand_only should fall back
        coords = _read_pose_coords(str(sample_pdb), ligand_only=True)
        assert len(coords) == 4  # Falls back to all heavy atoms


# ---------------------------------------------------------------------------
# Mode Animation Tests
# ---------------------------------------------------------------------------

class TestModeAnimationHelpers:
    """Test mode animation coordinate interpolation."""

    def test_read_atom_coords(self, sample_pdb):
        """Should read all ATOM coordinates including hydrogens."""
        from pymol_plugin.mode_animation import _read_atom_coords
        coords = _read_atom_coords(str(sample_pdb))
        assert len(coords) == 5  # All atoms including H

    def test_interpolate_coords_t0(self, sample_pdb, sample_pdb2):
        """t=0 should return coords1."""
        from pymol_plugin.mode_animation import (
            _read_atom_coords, _interpolate_coords
        )
        c1 = _read_atom_coords(str(sample_pdb))
        c2 = _read_atom_coords(str(sample_pdb2))
        result = _interpolate_coords(c1, c2, 0.0)
        for i in range(len(result)):
            assert result[i] == pytest.approx(c1[i])

    def test_interpolate_coords_t1(self, sample_pdb, sample_pdb2):
        """t=1 should return coords2."""
        from pymol_plugin.mode_animation import (
            _read_atom_coords, _interpolate_coords
        )
        c1 = _read_atom_coords(str(sample_pdb))
        c2 = _read_atom_coords(str(sample_pdb2))
        result = _interpolate_coords(c1, c2, 1.0)
        for i in range(len(result)):
            assert result[i] == pytest.approx(c2[i])

    def test_interpolate_coords_midpoint(self, sample_pdb, sample_pdb2):
        """t=0.5 should return the midpoint."""
        from pymol_plugin.mode_animation import (
            _read_atom_coords, _interpolate_coords
        )
        c1 = _read_atom_coords(str(sample_pdb))
        c2 = _read_atom_coords(str(sample_pdb2))
        result = _interpolate_coords(c1, c2, 0.5)
        # First atom: (1,2,3) and (5,6,7) -> midpoint (3,4,5)
        assert result[0] == pytest.approx([3.0, 4.0, 5.0])

    def test_interpolate_unequal_lengths(self):
        """Unequal coordinate lists should use the shorter length."""
        from pymol_plugin.mode_animation import _interpolate_coords
        c1 = [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]]
        c2 = [[10.0, 10.0, 10.0], [11.0, 11.0, 11.0]]
        result = _interpolate_coords(c1, c2, 0.5)
        assert len(result) == 2


# ---------------------------------------------------------------------------
# ITC Comparison Tests
# ---------------------------------------------------------------------------

class TestITCComparison:
    """Test ITC data parsing and comparison helpers."""

    def test_parse_itc_csv(self, sample_itc_csv):
        """Should parse ITC experimental data from CSV."""
        from pymol_plugin.itc_comparison import _parse_itc_csv
        data = _parse_itc_csv(str(sample_itc_csv))
        assert data["dG"] == pytest.approx(-8.5)
        assert data["dH"] == pytest.approx(-12.3)
        assert data["TdS"] == pytest.approx(3.8)

    def test_parse_itc_csv_missing_file(self, tmp_dir):
        """Missing CSV file should raise an error."""
        from pymol_plugin.itc_comparison import _parse_itc_csv
        with pytest.raises((FileNotFoundError, OSError)):
            _parse_itc_csv(str(tmp_dir / "nonexistent.csv"))


# ---------------------------------------------------------------------------
# Interactive Docking Config Tests
# ---------------------------------------------------------------------------

class TestInteractiveDockingHelpers:
    """Test interactive docking configuration generation."""

    def test_write_minimal_config(self, tmp_dir):
        """Should generate a valid FlexAID configuration file."""
        from pymol_plugin.interactive_docking import _write_minimal_config
        config_path = str(tmp_dir / "test.inp")
        _write_minimal_config(
            config_path,
            receptor_pdb="receptor.pdb",
            ligand_file="ligand.mol2",
            center=(10.0, 20.0, 30.0),
            radius=15.0,
            temperature=310,
            n_results=5,
        )
        text = Path(config_path).read_text()
        assert "PDBNAM receptor.pdb" in text
        assert "INPLIG ligand.mol2" in text
        assert "TEMPER 310" in text
        assert "NRGOUT 5" in text
        assert "COMPLF 10.000 20.000 30.000" in text
        assert "SPACER 15.0" in text

    def test_write_minimal_config_defaults(self, tmp_dir):
        """Default temperature and n_results should be used."""
        from pymol_plugin.interactive_docking import _write_minimal_config
        config_path = str(tmp_dir / "defaults.inp")
        _write_minimal_config(
            config_path,
            receptor_pdb="rec.pdb",
            ligand_file="lig.mol2",
            center=(0.0, 0.0, 0.0),
            radius=10.0,
        )
        text = Path(config_path).read_text()
        assert "TEMPER 300" in text
        assert "NRGOUT 10" in text


# ---------------------------------------------------------------------------
# NumPy Entropy Computation Tests
# ---------------------------------------------------------------------------

class TestEntropyNumPy:
    """Test NumPy-accelerated entropy computation."""

    @pytest.mark.skipif(not _HAS_NUMPY, reason="NumPy not installed")
    def test_numpy_pure_agreement(self, sample_pdb):
        """NumPy and pure-Python backends should produce matching results."""
        from pymol_plugin.entropy_heatmap import (
            _read_pose_coords,
            _compute_grid_bounds,
            _compute_spatial_entropy_numpy,
            _compute_spatial_entropy_pure,
        )
        coords = _read_pose_coords(str(sample_pdb), ligand_only=False)
        # Simulate 2 poses with slightly offset coords
        pose1 = coords
        pose2 = [(x + 0.5, y + 0.5, z + 0.5) for x, y, z in coords]
        all_c = pose1 + pose2
        bounds = _compute_grid_bounds(all_c)

        result_np = _compute_spatial_entropy_numpy(
            [pose1, pose2], bounds, grid_spacing=2.0, sigma=2.0,
        )
        result_py = _compute_spatial_entropy_pure(
            [pose1, pose2], bounds, grid_spacing=2.0, sigma=2.0,
        )

        assert len(result_np) == len(result_py)
        for (x1, y1, z1, s1), (x2, y2, z2, s2) in zip(result_np, result_py):
            assert x1 == pytest.approx(x2, abs=1e-6)
            assert y1 == pytest.approx(y2, abs=1e-6)
            assert s1 == pytest.approx(s2, abs=1e-4)

    @pytest.mark.skipif(not _HAS_NUMPY, reason="NumPy not installed")
    def test_numpy_has_numpy_flag(self):
        """Module should detect NumPy availability."""
        from pymol_plugin.entropy_heatmap import _HAS_NUMPY as flag
        assert flag is True


# ---------------------------------------------------------------------------
# Kabsch Alignment Tests
# ---------------------------------------------------------------------------

class TestKabschAlignment:
    """Test Kabsch RMSD-optimal alignment."""

    @pytest.mark.skipif(not _HAS_NUMPY, reason="NumPy not installed")
    def test_kabsch_identity(self):
        """Aligning identical coordinates should return them unchanged."""
        from pymol_plugin.mode_animation import _kabsch_align
        coords = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
        aligned = _kabsch_align(coords, coords)
        for i in range(len(coords)):
            assert aligned[i] == pytest.approx(coords[i], abs=1e-10)

    @pytest.mark.skipif(not _HAS_NUMPY, reason="NumPy not installed")
    def test_kabsch_translation(self):
        """Should align a translated copy back to the target."""
        from pymol_plugin.mode_animation import _kabsch_align
        target = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
        shifted = [[x + 10, y + 20, z + 30] for x, y, z in target]
        aligned = _kabsch_align(shifted, target)
        for i in range(len(target)):
            assert aligned[i] == pytest.approx(target[i], abs=1e-8)

    @pytest.mark.skipif(not _HAS_NUMPY, reason="NumPy not installed")
    def test_kabsch_reduces_rmsd(self, sample_pdb, sample_pdb2):
        """Alignment should not increase RMSD compared to unaligned."""
        from pymol_plugin.mode_animation import _read_atom_coords, _kabsch_align

        c1 = _read_atom_coords(str(sample_pdb))
        c2 = _read_atom_coords(str(sample_pdb2))
        n = min(len(c1), len(c2))

        def rmsd(a, b):
            return (sum(
                (a[i][j] - b[i][j]) ** 2
                for i in range(n) for j in range(3)
            ) / n) ** 0.5

        rmsd_before = rmsd(c1, c2)
        c2_aligned = _kabsch_align(c2, c1)
        rmsd_after = rmsd(c1, c2_aligned)
        assert rmsd_after <= rmsd_before + 1e-10

    def test_kabsch_no_numpy_passthrough(self):
        """Without NumPy, _kabsch_align should return input unchanged."""
        from pymol_plugin.mode_animation import _kabsch_align, _HAS_NUMPY
        coords = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
        if not _HAS_NUMPY:
            result = _kabsch_align(coords, coords)
            assert result is coords


# ---------------------------------------------------------------------------
# Multi-Ligand ITC CSV Tests
# ---------------------------------------------------------------------------

class TestITCMultiLigand:
    """Test multi-row ITC CSV parsing."""

    def test_parse_multi_row(self, tmp_dir):
        """Should parse multiple ligand rows."""
        csv_path = tmp_dir / "multi_itc.csv"
        csv_path.write_text(
            "ligand,dG,dH,TdS\n"
            "aspirin,-8.5,-12.3,3.8\n"
            "ibuprofen,-7.2,-10.1,2.9\n"
        )
        from pymol_plugin.itc_comparison import _parse_itc_csv_multi
        rows = _parse_itc_csv_multi(str(csv_path))
        assert len(rows) == 2
        assert rows[0]["ligand"] == "aspirin"
        assert rows[0]["dG"] == pytest.approx(-8.5)
        assert rows[1]["ligand"] == "ibuprofen"
        assert rows[1]["dH"] == pytest.approx(-10.1)

    def test_parse_multi_row_no_ligand_column(self, tmp_dir):
        """Should work without a ligand name column."""
        csv_path = tmp_dir / "no_ligand.csv"
        csv_path.write_text("dG,dH,TdS\n-8.0,-11.0,3.0\n-7.0,-9.0,2.0\n")
        from pymol_plugin.itc_comparison import _parse_itc_csv_multi
        rows = _parse_itc_csv_multi(str(csv_path))
        assert len(rows) == 2
        assert "ligand" not in rows[0]

    def test_single_row_backward_compat(self, sample_itc_csv):
        """_parse_itc_csv should still return first row dict."""
        from pymol_plugin.itc_comparison import _parse_itc_csv
        data = _parse_itc_csv(str(sample_itc_csv))
        assert data["dG"] == pytest.approx(-8.5)


# ---------------------------------------------------------------------------
# Async Docking & Cancel Tests
# ---------------------------------------------------------------------------

class TestDockingAsync:
    """Test async docking infrastructure."""

    def test_callback_cancel(self):
        """DockingProgressCallback.cancel() should set cancelled flag."""
        from pymol_plugin.interactive_docking import DockingProgressCallback
        cb = DockingProgressCallback()
        assert cb.cancelled is False
        cb.cancel()
        assert cb.cancelled is True

    def test_active_docking_initially_none(self):
        """Module-level _active_docking should start as None."""
        from pymol_plugin.interactive_docking import _active_docking
        assert _active_docking is None


# ---------------------------------------------------------------------------
# CGO Renderer Tests
# ---------------------------------------------------------------------------

class TestCGORenderer:
    """Test CGO heatmap rendering."""

    def test_render_cgo_calls_load_cgo(self, mock_pymol):
        """_render_cgo should call cmd.load_cgo with sphere data."""
        from pymol_plugin.entropy_heatmap import _render_cgo
        points = [(1.0, 2.0, 3.0, 0.5), (4.0, 5.0, 6.0, 0.8)]
        _render_cgo(points, "test_cgo", sphere_scale=0.4, transparency=0.3)
        mock_pymol.cmd.load_cgo.assert_called_once()
        args = mock_pymol.cmd.load_cgo.call_args
        cgo_list = args[0][0]
        obj_name = args[0][1]
        assert obj_name == "test_cgo"
        # Should contain ALPHA, COLOR, SPHERE primitives
        assert 25.0 in cgo_list  # ALPHA
        assert 6.0 in cgo_list   # COLOR
        assert 7.0 in cgo_list   # SPHERE


# ---------------------------------------------------------------------------
# Shared session + load-path sync
# ---------------------------------------------------------------------------

class TestPluginSession:
    """Shared SESSION state used by all plugin modules."""

    def test_object_and_group_naming(self, mock_pymol):
        from pymol_plugin.session import object_name, group_name, mode_label
        assert object_name("flexaids", 1, 2) == "flexaids_mode1_pose2"
        assert group_name("dock", 3) == "dock_mode3"
        assert mode_label(7) == "mode7"

    def test_session_clear(self, mock_pymol):
        from pymol_plugin.session import SESSION
        SESSION.prefix = "tmp"
        SESSION.objects[1] = ["a"]
        SESSION.clear()
        assert SESSION.result is None
        assert SESSION.objects == {}
        assert SESSION.prefix == "flexaids"

    def test_get_mode(self, mock_pymol, sample_pdb, tmp_dir):
        from pymol_plugin.session import SESSION
        from flexaidds.models import PoseResult, BindingModeResult, DockingResult

        pose = PoseResult(path=sample_pdb, mode_id=1, pose_rank=1, cf=-5.0)
        mode = BindingModeResult(
            mode_id=1, rank=1, poses=[pose], free_energy=-8.0, best_cf=-5.0,
        )
        SESSION.clear()
        SESSION.result = DockingResult(
            source_dir=tmp_dir, binding_modes=[mode], temperature=298.0,
        )
        assert SESSION.get_mode(1) is mode
        assert SESSION.get_mode(99) is None
        SESSION.clear()


class TestLoadPathSync:
    """flexaids_load and flexaids_load_results share object naming + state."""

    def _write_mode_pdb(self, path: Path, mode_id: int, cf: float) -> None:
        path.write_text(
            f"REMARK Binding Mode:{mode_id} Best CF in Binding Mode:{cf} "
            f"Binding Mode Frequency:1\n"
            f"REMARK Free Energy:{-abs(cf) - 2.0}\n"
            f"REMARK Enthalpy:{-abs(cf) - 1.0}\n"
            f"REMARK Entropy:0.01\n"
            f"REMARK Temperature:298.15\n"
            "HETATM    1  C1  LIG A   1       1.000   2.000   3.000  1.00  0.00           C\n"
            "HETATM    2  C2  LIG A   1       2.000   3.000   4.000  1.00  0.00           C\n"
            "END\n"
        )

    def test_load_docking_results_populates_session(self, mock_pymol, tmp_dir):
        from pymol_plugin.session import SESSION
        from pymol_plugin import results_adapter

        pdb = tmp_dir / "mode1_pose1.pdb"
        self._write_mode_pdb(pdb, 1, -5.5)
        SESSION.clear()
        results_adapter.load_docking_results(str(tmp_dir), prefix="flexaids")

        assert SESSION.result is not None
        assert SESSION.n_modes >= 1
        assert 1 in SESSION.objects or any(SESSION.objects)
        # Objects use canonical names
        names = [n for ns in SESSION.objects.values() for n in ns]
        assert any(n.startswith("flexaids_mode") for n in names)
        mock_pymol.cmd.load.assert_called()
        SESSION.clear()

    def test_load_binding_modes_syncs_results_adapter(self, mock_pymol, tmp_dir):
        from pymol_plugin.session import SESSION
        from pymol_plugin import visualization, results_adapter

        pdb = tmp_dir / "mode1_pose1.pdb"
        self._write_mode_pdb(pdb, 1, -4.0)
        SESSION.clear()
        visualization.load_binding_modes(str(tmp_dir), temperature=298.15)

        assert SESSION.result is not None
        # results_adapter can resolve the mode
        mode = results_adapter._get_mode(1)
        assert mode is not None
        assert mode.mode_id == 1
        # visualization records exist
        assert "mode1" in SESSION.mode_records
        SESSION.clear()

    def test_unload_results(self, mock_pymol, tmp_dir):
        from pymol_plugin.session import SESSION
        from pymol_plugin import results_adapter

        pdb = tmp_dir / "mode1_pose1.pdb"
        self._write_mode_pdb(pdb, 1, -3.0)
        results_adapter.load_docking_results(str(tmp_dir))
        assert SESSION.result is not None
        results_adapter.unload_results(delete_objects=1)
        assert SESSION.result is None
        assert SESSION.objects == {}

    def test_show_mode_details_handles_missing(self, mock_pymol, capsys):
        from pymol_plugin.session import SESSION
        from pymol_plugin import results_adapter

        SESSION.clear()
        results_adapter.show_mode_details(99)
        out = capsys.readouterr().out
        assert "ERROR" in out

    def test_boltzmann_weights_stable(self, mock_pymol):
        from pymol_plugin.visualization import _compute_boltzmann_weights
        w = _compute_boltzmann_weights([-10.0, -9.0, -5.0], 300.0)
        assert len(w) == 3
        assert abs(sum(w) - 1.0) < 1e-9
        assert w[0] > w[1] > w[2]

    def test_mode_name_accepts_integer_string(self, mock_pymol, tmp_dir):
        from pymol_plugin.session import SESSION
        from pymol_plugin import visualization
        from flexaidds.models import PoseResult, BindingModeResult, DockingResult

        pose = PoseResult(path=tmp_dir / "p.pdb", mode_id=2, pose_rank=1, cf=-1.0)
        mode = BindingModeResult(mode_id=2, rank=1, poses=[pose], best_cf=-1.0)
        SESSION.clear()
        SESSION.result = DockingResult(
            source_dir=tmp_dir, binding_modes=[mode], temperature=300.0,
        )
        SESSION.objects[2] = ["flexaids_mode2_pose1"]
        visualization.sync_from_session()
        rec = visualization._find_mode("2")
        assert rec is not None
        assert rec.mode_id == 2
        SESSION.clear()


class TestPluginHelp:
    def test_flexaids_help_prints(self, mock_pymol, capsys):
        from pymol_plugin import flexaids_help
        flexaids_help()
        out = capsys.readouterr().out
        assert "flexaids_load" in out
        assert "flexaids_entropy_heatmap" in out
        assert "CF/contact-function" in out
