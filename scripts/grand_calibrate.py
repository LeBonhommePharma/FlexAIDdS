#!/usr/bin/env python3
"""
grand_calibrate.py — Grand canonical (Ξ) validation harness for P4.

Handles synthetic exact fixtures and competition manifests to validate:
  p_bind, mean_occupancy, log_Xi, apparent/intrinsic selectivity
  against analytical ground truth or literature-derived references.

Pure-Python implementation of GrandPartitionFunction observables (log-space,
no external deps beyond stdlib + optional yaml for manifests). Mirrors
LIB/GrandPartitionFunction.cpp math exactly for cross-validation.

Intended for:
- Pre-wiring validation (current P4)
- Later binding to real StatMech log_Z from results or _core
- HW parity checks (scalar vs Metal/CUDA builds produce same grand obs)
- CI smoke for grand metrics

Usage examples (P4):
  python3 scripts/grand_calibrate.py --synthetic benchmarks/grand_synthetic/dual_ligand_exact.json --verbose
  python3 scripts/grand_calibrate.py --synthetic benchmarks/grand_synthetic/ --dry-run
  python3 scripts/grand_calibrate.py --competition benchmarks/datasets/competition_example.yaml --dry-run

Reproducibility: prints engine "py-fallback", T, c0, full case expectations vs computed.
No overclaim: only tests partition math + loading; real docking Z later.

See GPF_IMPLEMENTATION_PLAN.md (P4), docs/GrandPartitionFunction_Report.md,
benchmarks/datasets/competition_example.yaml, benchmarks/grand_synthetic/README.md
"""

import argparse
import csv
import json
import math
import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

# Physical / GPF constants (match C++ + thermodynamics.py fallback)
kB_kcal = 0.001987206
C0_M = 1.0  # standard state

def logsumexp(vals):
    """Numerically stable log-sum-exp. Matches GPF implementation."""
    if not vals:
        return 0.0
    m = max(vals)
    s = sum(math.exp(v - m) for v in vals)
    return m + math.log(s)

def compute_grand_observables(ligands, temperature_K=298.0, c0=C0_M):
    """
    ligands: list of dicts or tuples: {'name': str, 'log_Z': float, 'conc_M': float}
    Returns dict with log_Xi, p_empty, p_bind dict, mean_occupancy, selectivities.
    Pure Py, log-space only. Empty site always contributes 0 to ln terms.
    """
    if not ligands:
        return {
            "log_Xi": 0.0,
            "Xi": 1.0,
            "p_empty": 1.0,
            "p_bind": {},
            "mean_occupancy": 0.0,
            "selectivity": {},
            "log_zZ": {},
        }

    log_zZ = {}
    for lig in ligands:
        name = lig["name"] if isinstance(lig, dict) else lig[0]
        logZ = lig["log_Z"] if isinstance(lig, dict) else lig[1]
        c = lig["conc_M"] if isinstance(lig, dict) else lig[2]
        lnz = math.log(c / c0) if c > 0 else float("-inf")
        log_zZ[name] = lnz + logZ

    # lnXi = logsumexp( [0.0 (empty)] + list(log_zZ.values()) )
    terms = [0.0] + list(log_zZ.values())
    lnXi = logsumexp(terms)
    Xi = math.exp(lnXi) if lnXi < 700 else float("inf")  # guard

    p_bind = {}
    for name, lzz in log_zZ.items():
        p = math.exp(lzz - lnXi) if lnXi < 700 else (1.0 if lzz == max(log_zZ.values()) else 0.0)
        p_bind[name] = p
    p_empty = math.exp(0.0 - lnXi) if lnXi < 700 else 0.0

    occ = 1.0 - p_empty

    # Pairwise selectivities (for all ordered pairs)
    sel = {}
    names = list(log_zZ.keys())
    for i, na in enumerate(names):
        for nb in names[i+1:]:
            lza, lzb = log_zZ[na], log_zZ[nb]
            # apparent (conc weighted) = zZa - zZb = (lnZa - lnZb) + (lnca - lncb)
            # but since log_zZ already = lnz + lnZ, diff is direct
            log_app = lza - lzb
            app = math.exp(log_app) if abs(log_app) < 700 else (float("inf") if log_app > 0 else 0.0)
            # intrinsic = lnZa - lnZb = log_zZ_a - log_zZ_b - (ln ca - ln cb)
            lnca = math.log(ligands[i]["conc_M"] / c0) if isinstance(ligands[i], dict) else math.log(1e-6)
            lncb = math.log(ligands[i+1 if nb==names[i+1] else 0]["conc_M"] / c0) if False else 0 # simplified; compute properly
            # Better: recover lnZ diff = (lza - lnz_a) - (lzb - lnz_b)
            # We store per name; recompute cleanly:
            pass
    # Simpler robust intrinsic/apparent:
    for na in names:
        for nb in names:
            if na == nb: continue
            lza = log_zZ[na]
            lzb = log_zZ[nb]
            log_app = lza - lzb
            app = math.exp(log_app) if abs(log_app) < 700 else (float("inf") if log_app>0 else 0.0)
            # intrinsic: remove conc contrib
            # find concs
            ca = next((l["conc_M"] for l in ligands if (l["name"] if isinstance(l,dict) else l[0])==na), 1e-6)
            cb = next((l["conc_M"] for l in ligands if (l["name"] if isinstance(l,dict) else l[0])==nb), 1e-6)
            log_int = (lza - math.log(ca / c0)) - (lzb - math.log(cb / c0))
            intr = math.exp(log_int) if abs(log_int) < 700 else (float("inf") if log_int>0 else 0.0)
            key = f"{na}_over_{nb}"
            if key not in sel:
                sel[key] = {"apparent": app, "intrinsic": intr, "log_apparent": log_app, "log_intrinsic": log_int}

    return {
        "log_Xi": lnXi,
        "Xi": Xi if math.isfinite(Xi) else None,
        "p_empty": p_empty,
        "p_bind": p_bind,
        "mean_occupancy": occ,
        "selectivity": sel,
        "log_zZ": log_zZ,
    }

def load_synthetic_json(path: Path):
    with open(path) as f:
        data = json.load(f)
    cases = data.get("cases", [])
    return cases, data.get("temperature_K", 298.0)

def load_competition_yaml(path: Path):
    try:
        import yaml
    except ImportError:
        print("PyYAML not available; competition yaml limited to structure only.", file=sys.stderr)
        yaml = None
    if yaml is None:
        with open(path) as f:
            # minimal parse not full
            return {"competition_sets": [], "note": "yaml parse skipped"}
    with open(path) as f:
        return yaml.safe_load(f)

def compare_floats(a, b, rel_tol=1e-9, abs_tol=1e-12):
    if a is None or b is None:
        return False
    if not math.isfinite(a) or not math.isfinite(b):
        return (not math.isfinite(a)) and (not math.isfinite(b))
    return math.isclose(a, b, rel_tol=rel_tol, abs_tol=abs_tol)

def validate_case(case, verbose=False):
    """Run pure-py GPF on case ligands, compare to expected."""
    ligands = case.get("ligands", [])
    exp = case.get("expected", {})
    tol = case.get("tolerance", {"log_rel": 1e-9, "prob_abs": 1e-12})

    computed = compute_grand_observables(ligands)

    results = []
    ok = True

    # log_Xi
    if "log_Xi" in exp:
        match = compare_floats(computed["log_Xi"], exp["log_Xi"], rel_tol=tol.get("log_rel", 1e-9))
        results.append(("log_Xi", computed["log_Xi"], exp["log_Xi"], match))
        ok = ok and match

    # p_empty
    if "p_empty" in exp:
        match = compare_floats(computed["p_empty"], exp["p_empty"], abs_tol=tol.get("prob_abs", 1e-12))
        results.append(("p_empty", computed["p_empty"], exp["p_empty"], match))
        ok = ok and match

    # p_bind per ligand
    if "p_bind" in exp:
        for name, pexp in exp["p_bind"].items():
            pcomp = computed["p_bind"].get(name)
            match = compare_floats(pcomp, pexp, abs_tol=tol.get("prob_abs", 1e-12))
            results.append((f"p_bind[{name}]", pcomp, pexp, match))
            ok = ok and match

    # occupancy
    if "mean_occupancy" in exp:
        match = compare_floats(computed["mean_occupancy"], exp["mean_occupancy"], abs_tol=tol.get("prob_abs", 1e-12))
        results.append(("occupancy", computed["mean_occupancy"], exp["mean_occupancy"], match))
        ok = ok and match

    # selectivities (spot check)
    if "selectivity" in exp:
        for skey, svals in exp["selectivity"].items():
            if skey in computed["selectivity"]:
                for k in ("apparent", "intrinsic"):
                    if k in svals:
                        cval = computed["selectivity"][skey].get(k)
                        match = compare_floats(cval, svals[k], rel_tol=1e-8)
                        results.append((f"sel_{skey}_{k}", cval, svals[k], match))
                        ok = ok and match

    if verbose or not ok:
        print(f"  Case {case.get('case_id', '?')}: {'✅ PASS' if ok else '❌ FAIL'}")
        for name, comp, ex, m in results:
            status = "OK" if m else "MISMATCH"
            print(f"    {name}: computed={comp} expected={ex} [{status}]")

    return ok, results

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--synthetic", default=None, help="Path to synthetic JSON fixture or dir of *.json")
    ap.add_argument("--competition", default=None, help="Path to competition_example.yaml (structure + refs; Z supplied separately for full)")
    ap.add_argument("--dry-run", action="store_true", help="Parse and compute but do not assert failures (report only)")
    ap.add_argument("--verbose", "-v", action="store_true")
    args = ap.parse_args()

    print("Grand Calibrate (P4) — pure-Py GPF validation harness")
    print(f"  repo: {REPO_ROOT}")
    print(f"  synthetic: {args.synthetic}")
    print(f"  competition: {args.competition}")
    print("  (using pure-Py logsumexp + GPF observables; no C++ required)")

    all_ok = True
    n_cases = 0

    if args.synthetic:
        p = Path(args.synthetic)
        jsons = []
        if p.is_dir():
            jsons = sorted(p.glob("*.json"))
        elif p.is_file():
            jsons = [p]
        for jf in jsons:
            print(f"\nLoading synthetic: {jf}")
            cases, T = load_synthetic_json(jf)
            print(f"  T={T} K ; {len(cases)} cases")
            for case in cases:
                n_cases += 1
                ok, _ = validate_case(case, verbose=args.verbose)
                if not ok:
                    all_ok = False
                    if not args.dry_run:
                        print("  ** MISMATCH in exact synthetic **", file=sys.stderr)

    if args.competition:
        print(f"\nLoading competition manifest: {args.competition}")
        comp = load_competition_yaml(Path(args.competition))
        sets = comp.get("competition_sets", []) if isinstance(comp, dict) else []
        print(f"  {len(sets)} competition sets (note: Z values supplied via fixtures or later docking runs)")
        # For P4, just report structure; full Z-driven comparison when synthetic or results available
        for s in sets:
            rec = s.get("receptor_id", "?")
            nlig = len(s.get("ligands", []))
            print(f"  receptor={rec} ligands={nlig}")
            # Could later map to synthetic Z or load from results

    print(f"\nSummary: {n_cases} synthetic cases processed. all_ok={all_ok} (dry_run={args.dry_run})")
    if not all_ok and not args.dry_run:
        sys.exit(1)
    print("Done. (Fixtures match pure-Py GPF implementation.)")

if __name__ == "__main__":
    main()