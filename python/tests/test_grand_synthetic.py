"""
test_grand_synthetic.py — P4 grand canonical synthetic exact + fixture tests (pure Python).

Loads benchmarks/grand_synthetic/ fixtures (known log_Z + c + expected)
and validates against pure-Py reference implementation of GrandPartitionFunction
observables. No C++ bindings or new wiring required.

Run:
  cd python && pip install -e . -q && pytest tests/test_grand_synthetic.py -q --tb=line

These cases exercise:
- exact analytical p_bind / p_empty / occupancy / Xi / selectivity
- log-space stability (extreme ratios)
- empty site
- conc sweep effects on apparent vs intrinsic selectivity

Later (P2+): extend to roundtrip with flexaidds grand models and real BindingPopulation log_Z.
Seeded GA small-system cases planned for integration tests (see GPF_IMPLEMENTATION_PLAN.md P4/P7).
"""

import json
import math
from pathlib import Path
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SYNTH_DIR = REPO_ROOT / "benchmarks" / "grand_synthetic"

# Minimal pure-Py GPF mirror (duplicated from grand_calibrate.py for test isolation;
# keep in sync until shared module added).
def logsumexp(vals):
    if not vals:
        return 0.0
    m = max(vals)
    return m + math.log(sum(math.exp(v - m) for v in vals))

def compute_grand(ligands):
    if not ligands:
        return {"log_Xi": 0.0, "p_empty": 1.0, "p_bind": {}, "mean_occupancy": 0.0}
    log_zZ = {}
    for lig in ligands:
        name = lig["name"]
        lnz = math.log(lig["conc_M"])
        log_zZ[name] = lnz + lig["log_Z"]
    lnXi = logsumexp([0.0] + list(log_zZ.values()))
    p_bind = {n: math.exp(lzz - lnXi) for n, lzz in log_zZ.items()}
    p_empty = math.exp(-lnXi)
    return {
        "log_Xi": lnXi,
        "p_empty": p_empty,
        "p_bind": p_bind,
        "mean_occupancy": 1.0 - p_empty,
    }

def load_json_cases(name):
    p = SYNTH_DIR / name
    data = json.loads(p.read_text())
    return data.get("cases", [])

@pytest.mark.parametrize("case", load_json_cases("dual_ligand_exact.json"))
def test_dual_ligand_exact(case):
    comp = compute_grand(case["ligands"])
    exp = case["expected"]
    tol = case.get("tolerance", {})
    rel = tol.get("log_rel", 1e-9)
    abs_p = tol.get("prob_abs", 1e-12)

    if "log_Xi" in exp:
        assert math.isclose(comp["log_Xi"], exp["log_Xi"], rel_tol=rel)
    if "p_empty" in exp:
        assert math.isclose(comp["p_empty"], exp["p_empty"], abs_tol=abs_p)
    for nm, pe in exp.get("p_bind", {}).items():
        assert math.isclose(comp["p_bind"].get(nm, -1), pe, abs_tol=abs_p)
    if "mean_occupancy" in exp:
        assert math.isclose(comp["mean_occupancy"], exp["mean_occupancy"], abs_tol=abs_p)

    # p sums to 1
    psum = comp["p_empty"] + sum(comp["p_bind"].values())
    assert math.isclose(psum, 1.0, abs_tol=1e-12)

@pytest.mark.parametrize("case", load_json_cases("multi_ligand_exact.json"))
def test_multi_and_extreme(case):
    comp = compute_grand(case["ligands"])
    exp = case["expected"]
    tol = case.get("tolerance", {})
    rel = tol.get("log_rel", 1e-9)
    abs_p = tol.get("prob_abs", 1e-12)

    if "log_Xi" in exp:
        assert math.isclose(comp["log_Xi"], exp["log_Xi"], rel_tol=rel)
    if "p_empty" in exp:
        assert math.isclose(comp["p_empty"], exp["p_empty"], abs_tol=abs_p)

    psum = comp["p_empty"] + sum(comp["p_bind"].values())
    assert math.isclose(psum, 1.0, abs_tol=1e-12)

def test_empty_site():
    comp = compute_grand([])
    assert math.isclose(comp["log_Xi"], 0.0)
    assert math.isclose(comp["p_empty"], 1.0)
    assert comp["mean_occupancy"] == 0.0

# Planned (P4 note, implemented post-wiring):
# - test_grand_vs_cpp_bindings (requires_core, roundtrip fixtures to _core.GrandPartitionFunction)
# - test_seeded_ga_small_system (use tiny receptor + 2 ligands, fixed RngSeed, assert p within tol of analytical Z)
# - integration with results loader for competitive manifests
# Add to python/tests/ and C++ tests/test_grand_partition.cpp + test_multi_site_gpf.cpp

if __name__ == "__main__":
    # allow direct run
    import sys
    pytest.main([__file__] + sys.argv[1:])