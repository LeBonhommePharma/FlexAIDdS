"""Verify the flexaidds package imports gracefully without C++ bindings.

The read-only result loading (models, io, results) must work even when the
compiled _core extension is unavailable.  All public API types should be
importable and instantiable via pure-Python fallback implementations.
"""

import importlib
import sys
import types
from pathlib import Path
from unittest import mock


_MODULES_TO_CLEAR = [
    "flexaidds._core", "flexaidds", "flexaidds.thermodynamics",
    "flexaidds.models", "flexaidds.io", "flexaidds.results",
    "flexaidds.encom", "flexaidds._fallback_types",
    "flexaidds.__version__",
]


def _reimport_without_core():
    """Helper: force-reimport flexaidds with _core blocked."""
    saved = {k: sys.modules.pop(k, None) for k in _MODULES_TO_CLEAR}
    sys.modules["flexaidds._core"] = None  # triggers ImportError on import
    return saved


def _restore_modules(saved):
    """Helper: restore saved module state."""
    for k in list(saved.keys()):
        sys.modules.pop(k, None)
    for k, v in saved.items():
        if v is not None:
            sys.modules[k] = v


def test_package_imports_without_core_extension():
    """Importing flexaidds must not crash when _core is missing."""
    saved = _reimport_without_core()
    try:
        import flexaidds

        # Read-only classes should be available
        assert hasattr(flexaidds, "PoseResult")
        assert hasattr(flexaidds, "BindingModeResult")
        assert hasattr(flexaidds, "DockingResult")
        assert hasattr(flexaidds, "load_results")
        assert hasattr(flexaidds, "Thermodynamics")
        assert hasattr(flexaidds, "StatMechEngine")
    finally:
        _restore_modules(saved)


def test_encom_fallback_types_available():
    """ENCoMEngine, NormalMode, VibrationalEntropy should be usable without C++."""
    saved = _reimport_without_core()
    try:
        import flexaidds

        # ENCoM types should be actual classes, not None
        assert flexaidds.ENCoMEngine is not None, "ENCoMEngine should not be None"
        assert flexaidds.NormalMode is not None, "NormalMode should not be None"
        assert flexaidds.VibrationalEntropy is not None, "VibrationalEntropy should not be None"

        # Should be instantiable
        mode = flexaidds.NormalMode(index=1, eigenvalue=0.5, frequency=0.7)
        assert mode.index == 1
        assert mode.eigenvalue == 0.5

        vs = flexaidds.VibrationalEntropy(S_vib_kcal_mol_K=0.01, temperature=300.0)
        assert vs.temperature == 300.0
    finally:
        _restore_modules(saved)


def test_fallback_stub_types_available():
    """WHAMBin, TIPoint, Replica, State, BoltzmannLUT should be usable without C++."""
    saved = _reimport_without_core()
    try:
        import flexaidds

        assert flexaidds.WHAMBin is not None, "WHAMBin should not be None"
        assert flexaidds.TIPoint is not None, "TIPoint should not be None"
        assert flexaidds.Replica is not None, "Replica should not be None"
        assert flexaidds.State is not None, "State should not be None"
        assert flexaidds.BoltzmannLUT is not None, "BoltzmannLUT should not be None"

        # Should be instantiable
        wbin = flexaidds.WHAMBin(coord_center=1.0, free_energy=-5.0)
        assert wbin.coord_center == 1.0

        ti = flexaidds.TIPoint(lambda_val=0.5, dV_dlambda=-2.0)
        assert ti.lambda_val == 0.5

        rep = flexaidds.Replica(temperature=400.0)
        assert rep.temperature == 400.0

        state = flexaidds.State(energy=-10.0, count=3)
        assert state.count == 3

        lut = flexaidds.BoltzmannLUT(temperature=300.0)
        assert lut.temperature == 300.0
    finally:
        _restore_modules(saved)


def test_all_exports_non_none():
    """All types in __all__ should be non-None after import."""
    saved = _reimport_without_core()
    try:
        import flexaidds

        for name in flexaidds.__all__:
            val = getattr(flexaidds, name, "MISSING")
            assert val != "MISSING", f"{name} not found in flexaidds"
            assert val is not None, f"{name} is None in flexaidds (missing fallback)"
    finally:
        _restore_modules(saved)


def test_load_results_works_without_core(tmp_path: Path):
    """load_results() should work without C++ bindings."""
    # Write a minimal PDB file
    pdb = tmp_path / "mode_1_pose_1.pdb"
    pdb.write_text(
        "REMARK binding_mode = 1\n"
        "REMARK CF = -5.0\n"
        "ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "END\n"
    )

    from flexaidds.results import load_results

    result = load_results(tmp_path)
    assert result.n_modes == 1
    assert result.binding_modes[0].best_cf == -5.0


def test_has_core_true_when_extension_loads():
    """HAS_CORE_BINDINGS must be True whenever the *compiled* _core loads.

    Regression guard: pure-Python-only symbols (TemperatureScanPoint,
    DeltaCpFit) must not be imported from _core — that used to raise
    ImportError and force HAS_CORE_BINDINGS=False even with a working .so,
    breaking accelerated pip wheels.

    Note: ``python/tests/conftest.py`` injects a lightweight stub module named
    ``flexaidds._core`` for pure-Python CI. We must clear that stub first and
    only assert when a real extension binary is present.
    """
    import pytest

    # Drop package modules *and* any session-level stub so a real .so/.pyd can
    # be discovered. Keep a copy for restore.
    saved = {
        k: sys.modules.pop(k, None)
        for k in list(sys.modules)
        if k == "flexaidds" or k.startswith("flexaidds.")
    }
    try:
        try:
            import flexaidds._core as core  # noqa: F401
        except ImportError:
            pytest.skip("_core extension not built in this environment")

        # Reject the pure-Python test stub (no binary origin).
        core_file = getattr(core, "__file__", None) or ""
        if not any(core_file.endswith(ext) for ext in (".so", ".pyd", ".dll")):
            pytest.skip("_core extension not built in this environment")

        import flexaidds

        assert flexaidds.HAS_CORE_BINDINGS is True
        # Pure-Python-only types must still be available.
        assert flexaidds.TemperatureScanPoint is not None
        assert flexaidds.DeltaCpFit is not None
        # The public API remains the serializable Python façade while its
        # numerical backend is the compiled engine.
        thermo_module = importlib.import_module("flexaidds.thermodynamics")
        assert flexaidds.StatMechEngine is thermo_module.StatMechEngine
        assert flexaidds.Thermodynamics is thermo_module.Thermodynamics
        assert (
            flexaidds.ThermodynamicBreakdown
            is thermo_module.ThermodynamicBreakdown
        )
        assert (
            flexaidds.ScientificProvenance
            is thermo_module.ScientificProvenance
        )
        engine = flexaidds.StatMechEngine(300.0)
        assert isinstance(engine._engine, core.StatMechEngine)

        # Direct low-level bindings carry the same fail-closed predicates.
        provenance = core.ScientificProvenance()
        provenance.energy_domain = core.EnergyDomain.CALIBRATED_KCAL_PER_MOL
        provenance.ensemble_measure = core.EnsembleMeasure.ENUMERATED_MICROSTATES
        provenance.reference_state = core.ReferenceState.MATCHED_ASSOCIATION_CYCLE
        provenance.energy_provenance = "\u00a0\u2003"
        provenance.measure_provenance = "sha256:d4735e3a265e16eee03f59718b9b5d03019c07d8b6c51f90da3a666eec13ab35"
        provenance.reference_provenance = "sha256:4e07408562bedb8b60ce05c1decfe3ad16b72230967de01f640b7e4729b49fce"
        assert provenance.is_proxy_only()
        assert not provenance.allows_canonical_physical_claim()

        provenance.energy_provenance = "sha256:6b86b273ff34fce19d6b804eff5a3f5747ada4eaa22f1d49c01e52ddb7875b4b"
        native_engine = core.StatMechEngine(300.0, provenance)
        native_engine.add_sample(-10.0)
        native_result = native_engine.compute()
        native_breakdown = native_engine.compute_breakdown()
        assert native_result.allows_binding_physical_claim()
        assert not native_result.is_proxy_only()
        assert native_breakdown.allows_binding_physical_claim()
        assert not native_breakdown.is_proxy_only()
    finally:
        for k in list(sys.modules):
            if k == "flexaidds" or k.startswith("flexaidds."):
                sys.modules.pop(k, None)
        for k, v in saved.items():
            if v is not None:
                sys.modules[k] = v
