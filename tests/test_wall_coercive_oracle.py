#!/usr/bin/env python3
"""Tests for wall oracle production config resolution (shipped script)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _load_wall():
    path = ROOT / "scripts" / "wall_coercive_oracle.py"
    name = "wall_coercive_oracle"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_resolve_dock_config_uses_ops_gates_not_diagnostic():
    mod = _load_wall()
    cfg = mod.resolve_dock_config(ROOT, "1M2Z")
    assert cfg is not None, "ops/gates/configs/1M2Z_dock_config.json must exist"
    assert "ops/gates/configs" in str(cfg).replace("\\", "/")
    assert "1M2Z" in cfg.name
    assert "diagnostic/probe_config" not in str(cfg).replace("\\", "/")


def test_clean_panel_excludes_1g9v():
    mod = _load_wall()
    assert "1G9V" not in mod.CLEAN_PANEL
    assert set(mod.CLEAN_PANEL) == {"1J3J", "1K3U", "1L7F", "1N1M", "1M2Z"}


def test_manifest_loads_production_configs():
    mod = _load_wall()
    man = ROOT / "ops" / "gates" / "panel_manifest.tsv"
    assert man.is_file()
    rows = mod.load_manifest(ROOT, man)
    assert len(rows) >= 5
    by_pdb = {r["pdb"]: r for r in rows}
    assert "1M2Z" in by_pdb
    cfg = by_pdb["1M2Z"]["config"]
    assert cfg.is_file()
    assert "ops/gates/configs" in str(cfg).replace("\\", "/")
    # Manifest decoy must be diagnostic falsemin; native is crystal SDF
    assert by_pdb["1M2Z"]["native"].suffix.lower() == ".sdf"
    assert "falsemin" in by_pdb["1M2Z"]["decoy"].name.lower() or by_pdb["1M2Z"][
        "decoy"
    ].suffix.lower() == ".pdb"
