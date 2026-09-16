"""Tests for docking_mode YAML semantics (self vs cross-docking fail-closed)."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / ".grok" / "skills" / "flexaidds" / "scripts" / "validate_dataset_semantics.py"


def _load_mod():
    spec = importlib.util.spec_from_file_location("validate_dataset_semantics", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_semantics_script_exists():
    assert SCRIPT.is_file()


def test_repo_yamls_pass_semantics():
    mod = _load_mod()
    code, errs = mod.validate_path(REPO / "benchmarks" / "datasets")
    assert code == 0, errs


def test_self_vs_cross_contradiction_fails():
    mod = _load_mod()
    bad = {
        "slug": "astex_diverse",
        "docking_mode": "cross_docking",
        "structural_states": ["holo"],
        "metrics": ["docking_power_top1"],
        "name": "should fail",
        "description": "native holo",
    }
    errs = mod.validate_config(bad, path=Path("fake_astex_diverse.yaml"))
    assert errs, "expected semantic errors for astex_diverse + cross_docking"


def test_crossdock_metrics_require_cross_mode():
    mod = _load_mod()
    bad = {
        "slug": "demo",
        "docking_mode": "self_docking",
        "structural_states": ["holo"],
        "metrics": ["crossdock_success_rate_2A"],
        "name": "bad",
        "description": "self-dock",
    }
    errs = mod.validate_config(bad, path=Path("bad.yaml"))
    assert any("cross-dock metrics" in e for e in errs)


def test_cli_exits_zero():
    import subprocess

    proc = subprocess.run(
        [sys.executable, str(SCRIPT)],
        cwd=str(REPO),
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "VALIDATION PASSED" in (proc.stdout + proc.stderr)


def test_dataset_config_loads_astex_diverse():
    yaml_path = REPO / "benchmarks" / "datasets" / "astex_diverse.yaml"
    if not yaml_path.is_file():
        pytest.skip("astex_diverse.yaml missing")
    try:
        from flexaidds.dataset_runner.runner import DatasetConfig
    except Exception as exc:  # pragma: no cover
        pytest.skip(f"flexaidds package unavailable: {exc}")
    cfg = DatasetConfig.from_yaml(yaml_path)
    assert cfg.docking_mode == "self_docking"
    assert cfg.slug == "astex_diverse"
    assert cfg.expected_baselines_role == "ci_gate_only"
    assert cfg.expected_baselines.get("docking_power_top1") == 0.70
    assert cfg.expected_baselines.get("sampling_power") == 0.50
    assert cfg.ci_gate_class.get("docking_power_top1") == "advisory"
    assert cfg.ci_gate_class.get("sampling_power") == "hard"
    assert cfg.ci_ranking_status == "unregistered"
    assert cfg.tier1_expected_roster == ["1gpk", "1mq6", "1xm6", "2cet"]
    assert cfg.matrix_pin_md5 == "9dc93717dfed0698006d88dd6a9627bc"
    assert "2hr7" in {t.lower() for t in cfg.targets}
    assert len(cfg.targets) == 85
