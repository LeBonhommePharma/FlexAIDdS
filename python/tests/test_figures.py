"""Pure-Python tests for flexaidds.figures (imagine / publication cover integration).

These tests exercise:
- Gate 6 detection from multiple JSON shapes (audit + repro)
- Best-mode summary extraction (with fallback)
- Prompt construction with real-value injection + required branding elements
- prepare_publication_figures layout + metadata + gate enforcement
- No C++ / _core required (uses the pure results loader path + mocks)

Run:
    python -m pytest python/tests/test_figures.py -q --tb=line
"""

from __future__ import annotations

import json
from pathlib import Path

try:
    import pytest
except ImportError:
    pytest = None  # only needed when running under pytest

# The module under test (must import cleanly without core bindings)
from flexaidds.figures import (
    build_imagine_cover_prompt,
    build_imagine_animation_prompt,
    check_gate6_passed,
    extract_best_mode_summary,
    prepare_publication_figures,
    BANNER,
)


def _write_json(p: Path, obj: dict) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, indent=2), encoding="utf-8")


def _make_minimal_results(tmp: Path, *, gate6_pass: bool = True, mode_id: int = 1) -> Path:
    """Create a tiny fake results dir that load_results can partially use + gate JSONs."""
    rd = tmp / "fake_run"
    rd.mkdir(parents=True, exist_ok=True)

    # Minimal pose file (load_results will discover it; best_pose logic tolerates it)
    poses_dir = rd / "poses"
    poses_dir.mkdir()
    pose = poses_dir / "mode1_pose1.pdb"
    pose.write_text(
        "REMARK 999 FlexAIDdS best scoring mode\n"
        "ATOM      1  CA  LIG L   1       1.000   2.000   3.000  1.00 10.00           C\n"
        "END\n",
        encoding="utf-8",
    )

    # Reproducibility JSON with gate + numbers (used by extract + gate checker)
    repro = {
        "provenance": {
            "timestamp": "2026-06-01T12:00:00Z",
            "flexaidds_git": {"commit": "deadbeef1234567", "branch": "master"},
            "gate_results": {
                "gate6_crosscheck": {"passed": gate6_pass, "max_dev": 0.008}
            },
        },
        "top_mode": {
            "mode_id": mode_id,
            "free_energy": -11.42,
            "enthalpy": -8.17,
            "entropy": -3.25,  # treated as approx -TΔS in summary
        },
        "input_hashes": {"receptor": "abc", "ligand": "def"},
    }
    _write_json(rd / "reproducibility.json", repro)

    # Also drop an audit-style JSON (tests multiple discovery paths)
    audit = {
        "total_sampled": {"F_config_kcal_mol": -11.4},
        "provenance": {
            "temperature_K": 298.15,
            "gate_results": {"gate6_crosscheck": {"passed": gate6_pass}},
        },
    }
    _write_json(rd / "thermo_audit.json", audit)

    # Minimal ligand hint
    (rd / "ligand.mol2").write_text("@<TRIPOS>MOLECULE\nbiotin\n", encoding="utf-8")

    return rd


def test_check_gate6_passed_true_and_false(tmp_path: Path) -> None:
    rd_pass = _make_minimal_results(tmp_path / "p", gate6_pass=True)
    assert check_gate6_passed(rd_pass) is True

    rd_fail = _make_minimal_results(tmp_path / "f", gate6_pass=False)
    assert check_gate6_passed(rd_fail) is False

    # Empty dir -> conservative False
    empty = tmp_path / "empty"
    empty.mkdir()
    assert check_gate6_passed(empty) is False


def test_extract_best_mode_summary_injects_real_values(tmp_path: Path) -> None:
    rd = _make_minimal_results(tmp_path, gate6_pass=True)
    s = extract_best_mode_summary(rd)
    assert s["gate6_passed"] is True
    assert s["delta_g"] == pytest.approx(-11.42, abs=0.01)
    assert s["delta_h"] == pytest.approx(-8.17, abs=0.01)
    assert s["git_sha"].startswith("deadbeef") or s["git_sha"] != "unknown"
    assert "best_pose_pdb" in s or s.get("best_pose_pdb") is not None  # may be None if loader strict, but summary must not crash


def test_build_prompts_contain_required_elements_and_injected_numbers(tmp_path: Path) -> None:
    rd = _make_minimal_results(tmp_path, gate6_pass=True)
    s = extract_best_mode_summary(rd)

    cover = build_imagine_cover_prompt(s)
    anim = build_imagine_animation_prompt(s)

    # Injected real thermo
    assert "ΔG = -11.42" in cover or "ΔG = -11.42" in cover.replace("−", "-")
    assert "-11.42" in cover and "-8.17" in cover

    # Exact banner + footer style
    assert BANNER in cover
    assert BANNER in anim
    assert "gate6:PASS" in cover or "gate6:PASS" in cover.upper()
    assert "git:deadbeef" in cover or "git:dead" in cover  # short sha

    # Aesthetic / NRDD requirements (from the approved spec)
    assert "Nature Reviews Drug Discovery" in cover
    assert "SwitchCraft" in cover or "high-end" in cover.lower() or "cinematic" in cover.lower()
    assert "blue-to-red" in cover.lower() or "entropy heatmap" in cover.lower()
    assert "induced-fit" in cover.lower()
    assert "#22D3EE" in cover or "teal" in cover.lower() or "cyan" in cover.lower()

    # New P3 polish: critical interactions are injected and called out
    assert "key_interactions" in cover.lower() or "critical" in cover.lower() or "H-bond" in cover or "hydrophobic" in cover.lower()
    assert "critical" in anim.lower() or "interaction" in anim.lower() or "H-bond" in anim or "label" in anim.lower()

    # PLIP nice figs + fonts + CF priority (latest P3)
    assert "PLIP" in cover or "plip" in cover.lower() or "PLIP-style" in cover
    assert "JetBrains Mono" in cover or "thebonhomme.com" in cover.lower()
    assert "CF" in cover or "favourable" in cover.lower() or "most favourable" in cover.lower()

# --- NRDD / FlexAID∆S journal cover generation (the 5-point /flexaidds integration) ---

def test_nrdd_cover_params_validation():
    from flexaidds.figures import NRDDCoverParams, generate_flexaids_nrdd_cover
    # good
    p = NRDDCoverParams(entropy_value=0.93, enthalpy_value=1.4, index_value=0.92)
    assert p.index_value == 0.92
    # bad ranges (robust validation)
    try:
        NRDDCoverParams(index_value=3.0)
        assert False, "should have raised"
    except ValueError:
        pass
    try:
        NRDDCoverParams(entropy_value=99)
        assert False, "should have raised"
    except ValueError:
        pass
    # generation
    res = generate_flexaids_nrdd_cover(entropy_value=0.93, enthalpy_value=1.4, index_value=0.92, style="dramatic_faces")
    assert "prompt" in res and "metadata" in res
    prompt = res["prompt"]
    assert "JetBrains Mono" in prompt or "thebonhomme.com" in prompt.lower()
    assert "I_E-E" in prompt or "Enthalpy-Entropy Index" in prompt or "E–E" in prompt or "ENTROPY-ENTHALPY INDEX" in prompt.upper()
    assert "0.93" in prompt and "1.4" in prompt and "0.92" in prompt
    assert "FlexAID∆S" in prompt or "FLEXAID∆S" in prompt
    assert res["metadata"]["scientific_note"].startswith("Visualisation only")
    # from results_dir (graceful)
    res2 = generate_flexaids_nrdd_cover(results_dir=".", style="molecular_gauge")  # will use defaults or fail softly
    assert "prompt" in res2


def test_prepare_writes_layout_and_respects_gate(tmp_path: Path) -> None:
    # Gate pass path
    rd_pass = _make_minimal_results(tmp_path / "pass", gate6_pass=True)
    res = prepare_publication_figures(rd_pass, visualize=True, require_gate6=True)
    assert res["proceeded"] is True
    assert res["gate6_passed"] is True
    fd = Path(res["figures_dir"])
    assert fd.exists()
    assert (fd / "prompt_cover.txt").exists()
    assert (fd / "prompt_animation.txt").exists()
    assert (fd / "figure_metadata.json").exists()

    meta = json.loads((fd / "figure_metadata.json").read_text())
    assert meta["gate6_passed"] is True
    assert "required_elements" in meta and "/flexaids-docking" in str(meta["required_elements"])

    # Gate fail + require -> skip
    rd_fail = _make_minimal_results(tmp_path / "fail", gate6_pass=False)
    res2 = prepare_publication_figures(rd_fail, visualize=True, require_gate6=True)
    assert res2.get("proceeded") is False
    assert res2.get("skipped") == "gate6"

    # Force still proceeds
    res3 = prepare_publication_figures(rd_fail, visualize=True, require_gate6=True, force=True)
    assert res3["proceeded"] is True


def test_prepare_does_not_create_figures_without_visualize(tmp_path: Path) -> None:
    rd = _make_minimal_results(tmp_path, gate6_pass=True)
    res = prepare_publication_figures(rd, visualize=False)
    assert res["proceeded"] is False
    assert not (rd / "figures").exists()
