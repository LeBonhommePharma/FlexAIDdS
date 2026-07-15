"""Tests for .grok/skills/flexaidds/flexaidds_skill.py manifest handler."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SKILL_ROOT = REPO_ROOT / ".grok" / "skills" / "flexaidds"
sys.path.insert(0, str(SKILL_ROOT))
sys.path.insert(0, str(REPO_ROOT / "python"))

from flexaidds_skill import (  # noqa: E402
    FigureParameters,
    figure_parameters_from_mapping,
    generate_flexaids_figure,
    invoke_manifest_action,
)


def test_enthalpy_value_alias_maps_to_tds():
    params = figure_parameters_from_mapping(
        {"entropy_value": 0.5, "enthalpy_value": 1.1, "index_value": 0.8}
    )
    assert params.tds_value == 1.1
    assert params.entropy_value == 0.5


def test_unknown_manifest_keys_rejected():
    with pytest.raises(ValueError, match="Unknown figure parameter"):
        figure_parameters_from_mapping({"entropy_value": 0.5, "evil": 1})


def test_deferred_generation_returns_metadata():
    result = generate_flexaids_figure(
        params={"entropy_value": 0.93, "enthalpy_value": 1.4, "index_value": 0.92}
    )
    assert "prompt" in result and len(result["prompt"]) > 100
    assert "metadata" in result and isinstance(result["metadata"], dict)
    assert "params" in result
    assert result["path"] is None


def test_prompt_override_injection_blocked():
    with pytest.raises(ValueError, match="prompt override"):
        generate_flexaids_figure(
            params=FigureParameters(),
            prompt_overrides={"title": "ok\ninjected"},
        )


def test_invoke_manifest_action():
    manifest = json.loads((SKILL_ROOT / "flexaidds_skill_manifest.json").read_text())
    assert manifest["skill"] == "flexaidds"
    result = invoke_manifest_action(
        "generate_flexaids_figure",
        {"params": {"enthalpy_value": 1.4, "index_value": 0.92}},
    )
    assert result["metadata"]


def test_bin_validate_skill_wrapper_runs():
    wrapper = SKILL_ROOT / "bin" / "validate-skill"
    assert wrapper.is_file() and not wrapper.is_symlink()
    proc = subprocess.run([str(wrapper)], cwd=REPO_ROOT, capture_output=True, text=True, timeout=30)
    assert proc.returncode == 0, proc.stderr
    assert "VALIDATION PASSED" in proc.stdout