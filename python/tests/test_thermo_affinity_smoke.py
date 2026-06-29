"""Smoke test for ThermoAffinitySuite v2.1 (dry-run, no external data)."""
from __future__ import annotations

import os
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT / "python") not in os.sys.path:
    os.sys.path.insert(0, str(REPO_ROOT / "python"))

from benchmarks.runner import ThermoAffinitySuite  # type: ignore


def test_thermo_suite_imports_and_smoke():
    # pytest marker is optional; this runs under plain python -c too
    try:
        import pytest  # noqa
        _ = pytest.mark.slow
    except Exception:
        pass
    with tempfile.TemporaryDirectory() as td:
        suite = ThermoAffinitySuite(results_dir=td, dry_run=True)
        # Discover should at least see itc187.yaml
        cfgs = suite.discover_thermo_datasets()
        assert any(c.slug == "itc187" for c in cfgs)

        # Dry run on itc187 only (synthetic)
        res = suite.run(datasets=["itc187"], tier=1, dry_run=True, validate_pb=True, max_targets=1)
        assert "itc187" in res
        r = res["itc187"]
        assert r.n_targets >= 0
        assert hasattr(r, "pb_valid_rate")
        assert isinstance(r.metrics, dict)


def test_pb_wrapper_standalone():
    from benchmarks.runner import run_posebusters, HAS_POSEBUSTERS  # type: ignore
    HAS_PB = HAS_POSEBUSTERS  # compat alias for test
    # Should not crash even without files
    out = run_posebusters("/non/existent.pdb", rmsd=1.1)
    assert "success_pb" in out
    assert "failures" in out
    # When no pb installed, still produces a result based on rmsd alone
    if not HAS_PB:
        assert out["success_pb"] is True  # rmsd<2
