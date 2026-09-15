"""P0/P1 reproducibility hardening source-audit (no C++ build required).

Guards:
  1) DualAssembly growth is fail-closed on default docking (natural_deltaG==0
     unless --natural / FLEXAIDDS_NATURAL=1 / advanced.enable_natural=true).
  2) predicted_dG emits has_free_energy + cf_fallback; CF is never experimental ΔG.
  3) Astex denominator N=85 including 2HR7-as-failure; no astex85_codes_84.txt;
     expected_baselines 0.70 is ci_gate_only.
  4) PoseHelix / PoseLocal stay default OFF.
  5) Zhao JPCB 2011 is not cited with Werner DOI 10.1021/jp109255g as a live
     attribution.
  6) tencm_cpu_fallback OpenMP teams are num_threads-clamped.
  7) SUPPORT_MATRIX compiler floors match CMake (GCC 14 / Clang 18 / AppleClang 16).
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _read(rel: str) -> str:
    path = ROOT / rel
    assert path.is_file(), f"missing required file: {rel}"
    return path.read_text(encoding="utf-8", errors="replace")


class TestNaturalFailClosed:
    def test_gaboom_uses_docking_natural_growth_enabled(self):
        text = _read("LIB/gaboom.cpp")
        assert "docking_natural_growth_enabled" in text
        assert "natural_deltaG = 0.0" in text
        assert "FLEXAIDDS_NATURAL" in text

    def test_enable_natural_defaults_off(self):
        defaults = _read("LIB/config_defaults.h")
        assert re.search(r'"enable_natural"\s*,\s*V\(false\)', defaults)
        fa = _read("LIB/flexaid.h")
        assert "enable_natural" in fa
        top = _read("LIB/top.cpp")
        assert "FA->enable_natural=0" in top
        assert 'arg == "--natural"' in top

    def test_pose_helix_local_default_off(self):
        hdr = _read("LIB/NATURaL/NATURaLDualAssembly.h")
        assert "enable_pose_helix_rewrite = false" in hdr
        assert "pose_local_thermo_rewrite = false" in hdr


class TestPredictedDGFlags:
    def test_contract_header_exists(self):
        text = _read("LIB/predicted_dg_contract.h")
        assert "has_free_energy" in text
        assert "cf_fallback" in text
        assert "make_predicted_dg" in text
        assert "experimental" in text.lower()

    def test_dataset_runner_emits_flags(self):
        text = _read("LIB/DatasetRunner.cpp")
        assert "make_predicted_dg" in text
        assert "has_free_energy,cf_fallback" in text
        hdr = _read("LIB/DatasetRunner.h")
        assert "has_free_energy" in hdr
        assert "cf_fallback" in hdr

    def test_thermodynamics_docs_name_flags(self):
        text = _read("docs/thermodynamics.md")
        assert "has_free_energy" in text
        assert "cf_fallback" in text
        assert "never treat" in text.lower()


class TestAstexDenominator:
    def test_no_astex85_codes_84_file(self):
        forbidden = [
            ROOT / "state" / "astex85_codes_84.txt",
            ROOT / "benchmarks" / "protocols" / "astex85_codes_84.txt",
            ROOT / "benchmarks" / "astex_repro" / "astex85_codes_84.txt",
        ]
        hits = [str(p) for p in forbidden if p.exists()]
        assert hits == [], f"astex85_codes_84.txt must not exist: {hits}"

    def test_science_exclusions_keeps_2hr7_in_n85(self):
        text = _read("benchmarks/protocols/science_exclusions.md")
        assert "N=85" in text or "N = 85" in text
        assert "2HR7" in text
        assert "failure" in text.lower()
        assert "astex85_codes_84.txt" in text
        assert "not an exclusion" in text.lower() or "no science exclusions" in text.lower()

    def test_expected_baselines_role_ci_gate_only(self):
        for rel in (
            "benchmarks/datasets/astex_diverse.yaml",
            "python/flexaidds/dataset_runner/datasets/astex_diverse.yaml",
            "benchmarks/astex_diverse/manifest.yaml",
        ):
            text = _read(rel)
            assert "expected_baselines_role: ci_gate_only" in text, rel
            assert "0.70" in text

    def test_claim_receipt_placeholders_honest(self):
        text = _read("docs/ASTEX_CLAIM_RECEIPT.md")
        assert "UNSET" in text
        assert "no invented" in text.lower() or "not a measured" in text.lower()
        assert "N = 85" in text or "N=85" in text
        assert "2HR7" in text
        assert "ci_gate_only" in text


class TestZhaoDoiAndSupportMatrix:
    def test_ribosome_does_not_live_cite_werner_as_zhao(self):
        cpp = _read("LIB/NATURaL/RibosomeElongation.cpp")
        for line in cpp.splitlines():
            if "10.1021/jp109255g" not in line:
                continue
            low = line.lower()
            documented = (
                "doi-check: documented" in low
                or "no doi was invented" in low
                or "werner" in low
                or "do not attach" in low
            )
            assert documented, f"live Zhao attribution of Werner DOI: {line}"

    def test_openmp_fallback_clamps_num_threads(self):
        text = _read("LIB/tENCoM/tencm_cpu_fallback.cpp")
        assert text.count("num_threads(n_threads)") >= 2

    def test_support_matrix_matches_cmake_floors(self):
        text = _read("docs/SUPPORT_MATRIX.md")
        assert "GCC >= 14" in text
        assert "Clang >= 18" in text
        assert "Apple Clang >= 16" in text or "AppleClang >= 16" in text


def test_reproduce_astex85_dry_run_documented():
    sh = _read("scripts/reproduce_astex85.sh")
    assert "--dry-run" in sh
    assert "ASTEX_CLAIM_RECEIPT.md" in sh
    ci = _read(".github/workflows/ci.yml")
    assert "reproduce_astex85.sh --dry-run" in ci
