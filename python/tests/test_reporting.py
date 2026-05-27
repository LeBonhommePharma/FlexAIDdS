"""
Tests for flexaidds.reporting (Task 9)
"""

import tempfile
from pathlib import Path

import pytest

from flexaidds.reporting import (
    generate_pymol_script,
    generate_markdown_report,
    generate_temperature_scan_plot,
    write_all_reports,
)
from flexaidds.models import DockingResult, BindingModeResult, PoseResult


def make_sample_result():
    """Create a minimal DockingResult for testing."""
    pose = PoseResult(path=Path("pose.pdb"), rank=1, cf=-10.0)
    mode = BindingModeResult(
        mode_id=1,
        rank=1,
        poses=[pose],
        free_energy=-9.5,
        enthalpy=-10.0,
        entropy=0.0017,
        heat_capacity=0.0005,
        temperature=300.0,
        thermodynamics={
            "G_total_kcal_mol": -9.5,
            "H_eff_kcal_mol": -10.0,
            "minus_T_S_config_kcal_mol": 0.5,
            "components_complete": True,
        },
    )
    return DockingResult(
        source_dir=Path("/tmp/test_results"),
        binding_modes=[mode],
    )


def test_generate_pymol_script(tmp_path):
    result = make_sample_result()
    out = tmp_path / "test.pml"
    script = generate_pymol_script(result, out)
    assert "G_total" in script or "free_energy" in script
    assert out.exists()
    assert "EXPERIMENTAL" not in script  # since components_complete=True


def test_generate_markdown_report(tmp_path):
    result = make_sample_result()
    out = tmp_path / "test.md"
    md = generate_markdown_report(result, out)
    assert "G_total" in md
    assert out.exists()


def test_write_all_reports(tmp_path):
    result = make_sample_result()
    written = write_all_reports(result, tmp_path, prefix="test_report")
    assert "pymol" in written
    assert "markdown" in written
    assert (tmp_path / "test_report.pml").exists()
    assert (tmp_path / "test_report.md").exists()


def test_generate_temperature_scan_plot_no_matplotlib(monkeypatch, tmp_path):
    """Should gracefully return None if matplotlib not installed."""
    monkeypatch.setitem(__import__("sys").modules, "matplotlib", None)
    monkeypatch.setitem(__import__("sys").modules, "matplotlib.pyplot", None)

    scan = [{"T_K": 300, "G_kcal_mol": -10, "H_kcal_mol": -9, "S_kcal_mol_K": 0.003}]
    out = generate_temperature_scan_plot(scan, tmp_path / "plot.png")
    assert out is None


def test_generate_temperature_scan_plot_with_data(tmp_path):
    """Basic smoke test - will only run if matplotlib is available."""
    pytest.importorskip("matplotlib")

    scan = [
        {"T_K": 280, "G_kcal_mol": -8.0, "H_kcal_mol": -9.0, "S_kcal_mol_K": 0.0036, "Cv_kcal_mol_K": 0.01},
        {"T_K": 300, "G_kcal_mol": -9.5, "H_kcal_mol": -10.0, "S_kcal_mol_K": 0.0017, "Cv_kcal_mol_K": 0.02},
        {"T_K": 320, "G_kcal_mol": -10.5, "H_kcal_mol": -10.5, "S_kcal_mol_K": 0.0, "Cv_kcal_mol_K": 0.015},
    ]
    out = generate_temperature_scan_plot(scan, tmp_path / "scan.png")
    assert out is not None
    assert Path(out).exists()
