"""Polar-desolv lever: default OFF is bit-stable; high weight de-inverts 1G9V.

Drives the shipped engine via tools/probe_cf (SCORE_NATIVE + NATIVE_ONLY),
not a reimplementation of vcfunction. Skips when binary or Astex 1G9V inputs
are missing (CI without docking data).

Requires:
  - build/FlexAIDdS + build/probe_cf (or PATH)
  - ~/.flexaidds/benchmarks/astex_diverse/1G9V/{1G9V_apo.pdb,1G9V_ligand.sdf}
  - optional elected decoy: $FLEXAIDDS_LOCAL_ROOT/canary_clashfix2_*/out/1G9V/elected_pose.pdb
    or env FLEXAIDDS_1G9V_DECOY_PDB
"""
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
BENCH = Path.home() / ".flexaidds" / "benchmarks" / "astex_diverse" / "1G9V"
LOCAL = Path(os.environ.get("FLEXAIDDS_LOCAL_ROOT", Path.home() / "flexaidds_results"))


def _binaries():
    probe = REPO / "build" / "probe_cf"
    engine = REPO / "build" / "FlexAIDdS"
    if not probe.is_file() or not engine.is_file():
        return None, None
    return probe, engine


def _decoy_pdb() -> Path | None:
    env = os.environ.get("FLEXAIDDS_1G9V_DECOY_PDB")
    if env and Path(env).is_file():
        return Path(env)
    # Prefer the clashfix2 canary elected pose used in DIAG_1G9V
    candidates = sorted(LOCAL.glob("canary_clashfix2_*/out/1G9V/elected_pose.pdb"))
    if candidates:
        return candidates[-1]
    return None


def _config_path() -> Path | None:
    env = os.environ.get("FLEXAIDDS_1G9V_DOCK_CONFIG")
    if env and Path(env).is_file():
        return Path(env)
    cands = sorted(LOCAL.glob("canary_clashfix2_*/out/1G9V/dock_config.json"))
    return cands[-1] if cands else None


def _probe(probe: Path, engine: Path, pose: Path, weight: float) -> dict:
    rec = BENCH / "1G9V_apo.pdb"
    lig = BENCH / "1G9V_ligand.sdf"
    # Clean scorer-related env so concurrent shell experiments cannot leak.
    env = os.environ.copy()
    for k in list(env):
        if k.startswith("FLEXAIDDS_") and k not in (
            "FLEXAIDDS_LOCAL_ROOT",
            "FLEXAIDDS_1G9V_DECOY_PDB",
            "FLEXAIDDS_1G9V_DOCK_CONFIG",
            "FLEXAIDDS_POLAR_DESOLV_TEST_WEIGHT",
        ):
            env.pop(k, None)
    if weight > 0:
        env["FLEXAIDDS_POLAR_DESOLV_WEIGHT"] = str(weight)
    cmd = [
        str(probe),
        "--receptor",
        str(rec),
        "--ligand",
        str(lig),
        "--pose",
        str(pose),
        "--binary",
        str(engine),
        "--data-dir",
        str(engine.parent),
        "--pdb",
        "1G9V",
    ]
    cfg = _config_path()
    if cfg is not None:
        cmd.extend(["--config", str(cfg)])
    r = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=120)
    assert r.returncode == 0, f"probe_cf failed: {r.stderr[-500:]}\n{r.stdout[-500:]}"
    for line in reversed(r.stdout.strip().splitlines()):
        line = line.strip()
        if line.startswith("{") and "cf_total" in line:
            return json.loads(line)
    raise AssertionError(f"no JSON in probe_cf stdout: {r.stdout[-300:]}")


@pytest.mark.skipif(
    not (BENCH / "1G9V_apo.pdb").is_file() or not (BENCH / "1G9V_ligand.sdf").is_file(),
    reason="Astex 1G9V docking data not installed under ~/.flexaidds/benchmarks",
)
@pytest.mark.skipif(_binaries()[0] is None, reason="build/probe_cf + FlexAIDdS missing")
def test_polar_desolv_default_matches_native_cf():
    """w=0 must reproduce SCORE_NATIVE ~ -69.42 on crystal pose (no polar term)."""
    probe, engine = _binaries()
    nat = _probe(probe, engine, BENCH / "1G9V_ligand.sdf", 0.0)
    assert nat["cf_total"] == pytest.approx(-69.420582, abs=0.05)
    assert nat["cf_com"] == pytest.approx(-129.2074, abs=0.1)


@pytest.mark.skipif(
    not (BENCH / "1G9V_apo.pdb").is_file() or _decoy_pdb() is None or _binaries()[0] is None,
    reason="need 1G9V data + elected decoy PDB + built probe_cf",
)
def test_polar_desolv_high_weight_de_inverts_1g9v():
    """High polar-desolv weight must flip CF(native) < CF(elected decoy)."""
    probe, engine = _binaries()
    decoy = _decoy_pdb()
    assert decoy is not None
    nat0 = _probe(probe, engine, BENCH / "1G9V_ligand.sdf", 0.0)
    dec0 = _probe(probe, engine, decoy, 0.0)
    # Baseline invert: decoy preferred (more negative)
    assert dec0["cf_total"] < nat0["cf_total"], "baseline expected invert on this decoy"

    w = float(os.environ.get("FLEXAIDDS_POLAR_DESOLV_TEST_WEIGHT", "120"))
    nat = _probe(probe, engine, BENCH / "1G9V_ligand.sdf", w)
    dec = _probe(probe, engine, decoy, w)
    # Term must move sas (and thus total) for both when w>0
    assert nat["cf_sas"] > nat0["cf_sas"] + 1.0
    assert dec["cf_sas"] > dec0["cf_sas"] + 1.0
    # De-invert: native preferred
    assert nat["cf_total"] < dec["cf_total"], (
        f"w={w} did not flip: native={nat['cf_total']} elected={dec['cf_total']}"
    )
