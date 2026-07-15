"""pytest configuration for FlexAIDdS Python test suite.

Installs a lightweight stub for the ``flexaidds._core`` C++ extension before
any test module is collected.  This lets pure-Python modules (io, models,
results, docking, …) be imported and tested without a compiled extension,
while tests that genuinely need the C++ engine are individually guarded with
``pytest.mark.needs_core``.
"""

from __future__ import annotations

import sys
import types
import pytest


def _make_core_stub() -> types.ModuleType:
    """Return a minimal module that satisfies ``flexaidds/__init__.py``."""
    stub = types.ModuleType("flexaidds._core")
    # Marker so the package never treats this placeholder as a real extension.
    stub._FLEXAIDDS_CORE_STUB = True

    # Physical constants the package exposes at module level
    stub.kB_kcal = 0.001987206  # kcal mol⁻¹ K⁻¹
    stub.kB_SI = 1.380649e-23   # J K⁻¹

    # Sentinel classes – construction raises so tests that need the real
    # engine will fail loudly rather than silently passing on wrong values.
    class _NeedsRealCore:
        def __init__(self, *args, **kwargs):
            raise RuntimeError(
                "C++ _core not built – mark this test with "
                "@pytest.mark.needs_core or skip it."
            )

    stub.StatMechEngine = _NeedsRealCore
    stub.BoltzmannLUT = _NeedsRealCore
    stub.Thermodynamics = _NeedsRealCore
    stub.ENCoMEngine = _NeedsRealCore
    stub.NormalMode = _NeedsRealCore
    stub.Replica = _NeedsRealCore
    stub.State = _NeedsRealCore
    stub.TIPoint = _NeedsRealCore
    stub.VibrationalEntropy = _NeedsRealCore
    stub.WHAMBin = _NeedsRealCore
    stub.TorsionalENM = _NeedsRealCore
    stub.TorsionalNormalMode = _NeedsRealCore
    stub.Conformer = _NeedsRealCore
    stub.compute_shannon_entropy = _NeedsRealCore
    stub.compute_torsional_vibrational_entropy = _NeedsRealCore
    stub.run_shannon_thermo_stack = _NeedsRealCore
    return stub


def pytest_configure(config):
    """Register the ``needs_core`` marker."""
    config.addinivalue_line(
        "markers",
        "needs_core: test requires compiled C++ _core extension",
    )


def pytest_sessionstart(session):
    """Inject the stub before any test module is imported."""
    if "flexaidds._core" not in sys.modules:
        stub = _make_core_stub()
        sys.modules["flexaidds._core"] = stub
        # Also make it accessible as an attribute on the parent package stub
        # in case the package itself is already partially imported.
        if "flexaidds" in sys.modules:
            sys.modules["flexaidds"]._core = stub


# ---------------------------------------------------------------------------
# Fixture: skip tests that need the real C++ core when it isn't available
# ---------------------------------------------------------------------------

def _core_available() -> bool:
    try:
        import flexaidds._core as c
        # The stub's StatMechEngine raises on construction; the real one
        # returns a usable object.
        c.StatMechEngine(300.0)
        return True
    except Exception:
        return False


skip_without_core = pytest.mark.skipif(
    not _core_available(),
    reason="C++ _core extension not built; skipping StatMechEngine tests",
)
