"""Tests for flexaidds.visualization module.

Tests cover PyMOL unavailability handling, function interfaces, and
palette cycling logic. All tests mock PyMOL to avoid requiring it.
"""

import sys
import types
import warnings
from unittest.mock import MagicMock, patch, call

import pytest


# ---------------------------------------------------------------------------
# Helpers: mock PyMOL module
# ---------------------------------------------------------------------------

def _make_pymol_mock():
    """Create a mock pymol package with cmd submodule."""
    pymol_mod = types.ModuleType("pymol")
    cmd_mod = MagicMock()
    pymol_mod.cmd = cmd_mod
    return pymol_mod, cmd_mod


# ---------------------------------------------------------------------------
# _require_pymol — graceful degradation
# ---------------------------------------------------------------------------

class TestRequirePymol:
    def test_raises_when_pymol_unavailable(self):
        """_require_pymol raises ImportError when PyMOL is not installed."""
        # Force-reload visualization with _PYMOL_AVAILABLE = False
        import flexaidds.visualization as viz
        original = viz._PYMOL_AVAILABLE
        try:
            viz._PYMOL_AVAILABLE = False
            with pytest.raises(ImportError, match="PyMOL is required"):
                viz._require_pymol()
        finally:
            viz._PYMOL_AVAILABLE = original

    def test_no_error_when_pymol_available(self):
        """_require_pymol succeeds when PyMOL is available."""
        import flexaidds.visualization as viz
        original = viz._PYMOL_AVAILABLE
        try:
            viz._PYMOL_AVAILABLE = True
            # Should not raise
            viz._require_pymol()
        finally:
            viz._PYMOL_AVAILABLE = original


# ---------------------------------------------------------------------------
# load_binding_mode
# ---------------------------------------------------------------------------

class TestLoadBindingMode:
    def test_loads_each_pdb_file(self):
        """Each PDB path is loaded with correct object name and styling."""
        import flexaidds.visualization as viz

        mock_cmd = MagicMock()
        original_available = viz._PYMOL_AVAILABLE
        original_cmd = getattr(viz, "_cmd", None)
        try:
            viz._PYMOL_AVAILABLE = True
            viz._cmd = mock_cmd

            mode = MagicMock()
            pdb_paths = ["/tmp/pose_0.pdb", "/tmp/pose_1.pdb", "/tmp/pose_2.pdb"]
            viz.load_binding_mode(mode, pdb_paths, mode_name="test_mode",
                                  color="green", show="lines")

            assert mock_cmd.load.call_count == 3
            mock_cmd.load.assert_any_call("/tmp/pose_0.pdb", "test_mode_000")
            mock_cmd.load.assert_any_call("/tmp/pose_1.pdb", "test_mode_001")
            mock_cmd.load.assert_any_call("/tmp/pose_2.pdb", "test_mode_002")
            assert mock_cmd.color.call_count == 3
            assert mock_cmd.show.call_count == 3
        finally:
            viz._PYMOL_AVAILABLE = original_available
            if original_cmd is not None:
                viz._cmd = original_cmd

    def test_raises_without_pymol(self):
        import flexaidds.visualization as viz
        original = viz._PYMOL_AVAILABLE
        try:
            viz._PYMOL_AVAILABLE = False
            with pytest.raises(ImportError):
                viz.load_binding_mode(MagicMock(), ["/tmp/a.pdb"])
        finally:
            viz._PYMOL_AVAILABLE = original

    def test_empty_pdb_list(self):
        """No calls made when pdb_paths is empty."""
        import flexaidds.visualization as viz

        mock_cmd = MagicMock()
        original_available = viz._PYMOL_AVAILABLE
        original_cmd = getattr(viz, "_cmd", None)
        try:
            viz._PYMOL_AVAILABLE = True
            viz._cmd = mock_cmd

            viz.load_binding_mode(MagicMock(), [], mode_name="empty")
            mock_cmd.load.assert_not_called()
        finally:
            viz._PYMOL_AVAILABLE = original_available
            if original_cmd is not None:
                viz._cmd = original_cmd


# ---------------------------------------------------------------------------
# load_population — palette cycling and warnings
# ---------------------------------------------------------------------------

class TestLoadPopulation:
    def test_warns_when_no_pdbs_found(self, tmp_path):
        """Issues RuntimeWarning when no PDB files match the mode pattern."""
        import flexaidds.visualization as viz

        mock_cmd = MagicMock()
        original_available = viz._PYMOL_AVAILABLE
        original_cmd = getattr(viz, "_cmd", None)
        try:
            viz._PYMOL_AVAILABLE = True
            viz._cmd = mock_cmd

            # Population with one mode, no matching PDB files in empty dir
            population = [MagicMock()]  # one mode
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                viz.load_population(population, str(tmp_path))
                assert len(w) == 1
                assert issubclass(w[0].category, RuntimeWarning)
                assert "No PDB files found" in str(w[0].message)
        finally:
            viz._PYMOL_AVAILABLE = original_available
            if original_cmd is not None:
                viz._cmd = original_cmd

    def test_loads_matching_pdbs(self, tmp_path):
        """PDB files matching mode pattern are loaded."""
        import flexaidds.visualization as viz

        mock_cmd = MagicMock()
        original_available = viz._PYMOL_AVAILABLE
        original_cmd = getattr(viz, "_cmd", None)
        try:
            viz._PYMOL_AVAILABLE = True
            viz._cmd = mock_cmd

            # Create PDB files matching the pattern *_1_*.pdb
            (tmp_path / "result_1_pose1.pdb").write_text("ATOM")
            (tmp_path / "result_1_pose2.pdb").write_text("ATOM")
            # A non-matching file
            (tmp_path / "result_2_pose1.pdb").write_text("ATOM")

            population = [MagicMock()]  # one mode
            viz.load_population(population, str(tmp_path))

            # Should have loaded 2 PDBs (matching *_1_*.pdb)
            assert mock_cmd.load.call_count == 2
        finally:
            viz._PYMOL_AVAILABLE = original_available
            if original_cmd is not None:
                viz._cmd = original_cmd

    def test_palette_cycles(self, tmp_path):
        """Colors cycle when there are more modes than palette colors."""
        import flexaidds.visualization as viz

        mock_cmd = MagicMock()
        original_available = viz._PYMOL_AVAILABLE
        original_cmd = getattr(viz, "_cmd", None)
        try:
            viz._PYMOL_AVAILABLE = True
            viz._cmd = mock_cmd

            # Create files for 3 modes
            for mode_idx in range(1, 4):
                (tmp_path / f"result_{mode_idx}_pose1.pdb").write_text("ATOM")

            population = [MagicMock(), MagicMock(), MagicMock()]  # 3 modes

            # Palette with only 2 colors → third mode wraps
            viz.load_population(population, str(tmp_path),
                                palette=["red", "blue"])

            # 3 modes × 1 file each = 3 load calls
            assert mock_cmd.load.call_count == 3
            # First mode gets red, second blue, third red (wraps).
            # Object naming is `{prefix}_mode{id}_pose{rank}` since 2fb4eb76
            # ("mode identity"), which changed load_population but left these
            # assertions on the old `mode_NN_NNN` form.
            color_calls = mock_cmd.color.call_args_list
            assert color_calls[0] == call("red", "flexaids_mode1_pose1")
            assert color_calls[1] == call("blue", "flexaids_mode2_pose1")
            assert color_calls[2] == call("red", "flexaids_mode3_pose1")
        finally:
            viz._PYMOL_AVAILABLE = original_available
            if original_cmd is not None:
                viz._cmd = original_cmd


# ---------------------------------------------------------------------------
# color_by_energy
# ---------------------------------------------------------------------------

class TestColorByEnergy:
    def test_empty_energies_no_op(self):
        """Empty energies list does nothing."""
        import flexaidds.visualization as viz

        mock_cmd = MagicMock()
        original_available = viz._PYMOL_AVAILABLE
        original_cmd = getattr(viz, "_cmd", None)
        try:
            viz._PYMOL_AVAILABLE = True
            viz._cmd = mock_cmd

            viz.color_by_energy("obj", [])
            mock_cmd.set_color.assert_not_called()
        finally:
            viz._PYMOL_AVAILABLE = original_available
            if original_cmd is not None:
                viz._cmd = original_cmd

    def test_single_energy_no_crash(self):
        """Single energy should handle zero range gracefully."""
        import flexaidds.visualization as viz

        mock_cmd = MagicMock()
        original_available = viz._PYMOL_AVAILABLE
        original_cmd = getattr(viz, "_cmd", None)
        try:
            viz._PYMOL_AVAILABLE = True
            viz._cmd = mock_cmd

            viz.color_by_energy("obj", [-5.0])
            assert mock_cmd.set_color.call_count == 1
            assert mock_cmd.color.call_count == 1
        finally:
            viz._PYMOL_AVAILABLE = original_available
            if original_cmd is not None:
                viz._cmd = original_cmd

    def test_gradient_direction(self):
        """Lowest energy gets burgundy (frac=0), highest gets purple (frac=1)."""
        import flexaidds.visualization as viz

        mock_cmd = MagicMock()
        original_available = viz._PYMOL_AVAILABLE
        original_cmd = getattr(viz, "_cmd", None)
        try:
            viz._PYMOL_AVAILABLE = True
            viz._cmd = mock_cmd

            viz.color_by_energy("obj", [-10.0, 0.0, 10.0])

            # First call (state 1, energy=-10, frac=0): burgundy red [0.502, 0.0, 0.125]
            first_color = mock_cmd.set_color.call_args_list[0]
            assert first_color[0][0] == "_energy_1"
            rgb = first_color[0][1]
            assert pytest.approx(rgb[0], abs=0.01) == 0.502  # burgundy R
            assert pytest.approx(rgb[1], abs=0.01) == 0.0    # G always 0
            assert pytest.approx(rgb[2], abs=0.01) == 0.125  # burgundy B

            # Last call (state 3, energy=10, frac=1): purple blue [0.294, 0.0, 0.510]
            last_color = mock_cmd.set_color.call_args_list[2]
            rgb_last = last_color[0][1]
            assert pytest.approx(rgb_last[0], abs=0.01) == 0.294  # purple R
            assert pytest.approx(rgb_last[1], abs=0.01) == 0.0    # G always 0
            assert pytest.approx(rgb_last[2], abs=0.01) == 0.510  # purple B
        finally:
            viz._PYMOL_AVAILABLE = original_available
            if original_cmd is not None:
                viz._cmd = original_cmd

    def test_raises_without_pymol(self):
        import flexaidds.visualization as viz
        original = viz._PYMOL_AVAILABLE
        try:
            viz._PYMOL_AVAILABLE = False
            with pytest.raises(ImportError):
                viz.color_by_energy("obj", [-5.0, -3.0])
        finally:
            viz._PYMOL_AVAILABLE = original


# ---------------------------------------------------------------------------
# show_cleft_spheres
# ---------------------------------------------------------------------------

class TestShowCleftSpheres:
    def test_loads_and_styles_spheres(self):
        import flexaidds.visualization as viz

        mock_cmd = MagicMock()
        original_available = viz._PYMOL_AVAILABLE
        original_cmd = getattr(viz, "_cmd", None)
        try:
            viz._PYMOL_AVAILABLE = True
            viz._cmd = mock_cmd

            viz.show_cleft_spheres("/tmp/cleft.pdb", color="green",
                                   transparency=0.7)

            mock_cmd.load.assert_called_once_with("/tmp/cleft.pdb",
                                                   "cleft_spheres")
            mock_cmd.show.assert_called_once_with("spheres", "cleft_spheres")
            mock_cmd.color.assert_called_once_with("green", "cleft_spheres")
            mock_cmd.set.assert_called_once_with("sphere_transparency", 0.7,
                                                  "cleft_spheres")
        finally:
            viz._PYMOL_AVAILABLE = original_available
            if original_cmd is not None:
                viz._cmd = original_cmd


# ---------------------------------------------------------------------------
# setup_publication_view
# ---------------------------------------------------------------------------

class TestSetupPublicationView:
    def test_basic_settings_applied(self):
        import flexaidds.visualization as viz

        mock_cmd = MagicMock()
        mock_cmd.get_object_list.return_value = ["receptor", "ligand"]
        original_available = viz._PYMOL_AVAILABLE
        original_cmd = getattr(viz, "_cmd", None)
        try:
            viz._PYMOL_AVAILABLE = True
            viz._cmd = mock_cmd

            viz.setup_publication_view()

            mock_cmd.bg_color.assert_called_once_with("white")
            mock_cmd.orient.assert_called_once()
            mock_cmd.zoom.assert_called_once()
            # Receptor object should get cartoon + color
            mock_cmd.show.assert_any_call("cartoon", "receptor")
            mock_cmd.color.assert_any_call("grey80", "receptor")
        finally:
            viz._PYMOL_AVAILABLE = original_available
            if original_cmd is not None:
                viz._cmd = original_cmd

    def test_missing_receptor_no_crash(self):
        """If receptor object is not in scene, skip its styling."""
        import flexaidds.visualization as viz

        mock_cmd = MagicMock()
        mock_cmd.get_object_list.return_value = ["ligand"]  # no receptor
        original_available = viz._PYMOL_AVAILABLE
        original_cmd = getattr(viz, "_cmd", None)
        try:
            viz._PYMOL_AVAILABLE = True
            viz._cmd = mock_cmd

            viz.setup_publication_view()

            # Should not call show("cartoon", "receptor")
            for c in mock_cmd.show.call_args_list:
                assert c != call("cartoon", "receptor")
        finally:
            viz._PYMOL_AVAILABLE = original_available
            if original_cmd is not None:
                viz._cmd = original_cmd

    def test_ray_trace_and_png(self):
        import flexaidds.visualization as viz

        mock_cmd = MagicMock()
        mock_cmd.get_object_list.return_value = []
        original_available = viz._PYMOL_AVAILABLE
        original_cmd = getattr(viz, "_cmd", None)
        try:
            viz._PYMOL_AVAILABLE = True
            viz._cmd = mock_cmd

            viz.setup_publication_view(ray_trace=True,
                                       output_png="/tmp/render.png")

            mock_cmd.ray.assert_called_once_with(1200, 900)
            mock_cmd.png.assert_called_once_with("/tmp/render.png", dpi=300)
        finally:
            viz._PYMOL_AVAILABLE = original_available
            if original_cmd is not None:
                viz._cmd = original_cmd

    def test_no_ray_trace_png_only(self):
        import flexaidds.visualization as viz

        mock_cmd = MagicMock()
        mock_cmd.get_object_list.return_value = []
        original_available = viz._PYMOL_AVAILABLE
        original_cmd = getattr(viz, "_cmd", None)
        try:
            viz._PYMOL_AVAILABLE = True
            viz._cmd = mock_cmd

            viz.setup_publication_view(output_png="/tmp/render.png")

            mock_cmd.ray.assert_not_called()
            mock_cmd.png.assert_called_once()
        finally:
            viz._PYMOL_AVAILABLE = original_available
            if original_cmd is not None:
                viz._cmd = original_cmd
