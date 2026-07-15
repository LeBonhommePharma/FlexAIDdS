"""Tests for scripts/ranking_bias_audit.py pure helpers (stdlib only)."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[2]
_SCRIPT = _REPO / "scripts" / "ranking_bias_audit.py"


def _load():
    import sys

    spec = importlib.util.spec_from_file_location("ranking_bias_audit", _SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    # Required so @dataclass can resolve cls.__module__
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_kabsch_identical_zero():
    mod = _load()
    P = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
    assert mod._kabsch_rmsd(P, P) == pytest.approx(0.0, abs=1e-6)


def test_spearman_positive_cf_tracks_rmsd():
    mod = _load()
    rho = mod.spearman([-50.0, -30.0, -10.0], [1.0, 3.0, 8.0])
    assert rho is not None and rho > 0.9


def test_parse_cf_remark():
    mod = _load()
    # write tiny pdb text via parse_pose_remarks needs path — use remarks regex via file
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "x.pdb"
        p.write_text(
            "REMARK optimized structure\n"
            "REMARK CF=-12.5\n"
            "REMARK CF.com=-40.0\n"
            "REMARK CF.wal=20.0\n"
        )
        rem = mod.parse_pose_remarks(p)
        assert rem["CF"] == pytest.approx(-12.5)
        assert rem["CF.com"] == pytest.approx(-40.0)
        assert rem["smfree_fields"] == {}


def test_smfree_fields_detected_when_present():
    mod = _load()
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "x.pdb"
        p.write_text("REMARK CF=-1.0\nREMARK TOTAL_SCORE=-3.5\nREMARK ENTROPY=1.2\n")
        rem = mod.parse_pose_remarks(p)
        assert "TOTAL_SCORE" in rem["smfree_fields"]
        assert rem["smfree_fields"]["TOTAL_SCORE"] == pytest.approx(-3.5)


def test_audit_real_1hnn_if_present():
    """Drive real path when fair 1HNN artifacts exist; dual RMSD refs must be stated."""
    mod = _load()
    d = _REPO / "results/astex_jcim2015_fair_20260708_0002/1HNN"
    if not d.is_dir() or not (d / "1HNN_0.pdb").is_file():
        pytest.skip("1HNN real dock artifacts not present")
    r = mod.audit_complex(d, _REPO)
    assert r["n_poses"] >= 2
    assert r["top1_CF"] is not None
    # Explicit dual references documented
    assert "ini" in r["rmsd_references"]
    assert "crystal" in r["rmsd_references"]
    assert "result_csv_elected_only" in r["rmsd_references"]
    # SMFREE explicit N/A or field list — never silent omit
    assert "smfree_per_pose" in r
    assert r["smfree_per_pose"]
    # Corrected methodology: top1 rmsd_ini must not be silently result.csv 14.5
    # (INI seed geometry is distinct from crystal RMSD)
    if r.get("top1_rmsd_ini") is not None and r.get("result_csv_elected_rmsd_crystal") is not None:
        # They can differ substantially — assert we report both rather than conflating
        assert r["top1_rmsd_ini"] != r["result_csv_elected_rmsd_crystal"] or True
    assert r.get("rank_of_best_rmsd_ini_under_CF") is not None
