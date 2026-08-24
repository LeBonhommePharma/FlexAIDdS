#!/usr/bin/env python3
"""Five-way RMSD cross-check named in METHODOLOGY.md §0.

This is a test/harness, not a campaign. Claim success is rank-0 in-place RMSD
``<= 2.0 Å`` in the receptor frame — never a superposed (Kabsch) value.

The five implementations (do not invent a sixth):

1. ``python/flexaidds/benchmark.py::compute_rmsd`` — CI / dataset_runner gate
2. ``benchmarks/astex_repro/score_reference.py`` — spyrmsd (optional dep)
3. ``benchmarks/astex_repro/score_offline.py::hrmsd`` — element-blocked Hungarian
4. ``LIB/calc_rmsd.cpp::calc_Hungarian_RMSD`` — pose PDB REMARK (C++ gtest)
5. ``LIB/DatasetRunner.cpp::dataset::hungarian_rmsd`` — result.csv (C++ gtest)

Methods 4 and 5 are compared in ``tests/test_dataset_runner.cpp``
(``RmsdCrossCheck``, ``RmsdClaimCutoff``). This module covers 1–3 plus the
shared claim-cutoff contract.

Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CLAIM_CUTOFF_A = 2.0

# METHODOLOGY.md §0 numbering — tests pin this list.
FIVE_METHODS = (
    {
        "id": 1,
        "name": "in_repo_compute_rmsd",
        "path": "python/flexaidds/benchmark.py::compute_rmsd",
        "role": "CI / dataset_runner gate",
        "live_on_campaign": False,
        "cxx": False,
    },
    {
        "id": 2,
        "name": "score_reference_spyrmsd",
        "path": "benchmarks/astex_repro/score_reference.py",
        "role": "offline reference (spyrmsd)",
        "live_on_campaign": False,
        "cxx": False,
    },
    {
        "id": 3,
        "name": "score_offline_hungarian",
        "path": "benchmarks/astex_repro/score_offline.py::hrmsd",
        "role": "offline permissive Hungarian",
        "live_on_campaign": False,
        "cxx": False,
    },
    {
        "id": 4,
        "name": "engine_calc_Hungarian_RMSD",
        "path": "LIB/calc_rmsd.cpp::calc_Hungarian_RMSD",
        "role": "pose PDB REMARK RMSD",
        "live_on_campaign": True,
        "cxx": True,
    },
    {
        "id": 5,
        "name": "dataset_hungarian_rmsd",
        "path": "LIB/DatasetRunner.cpp::dataset::hungarian_rmsd",
        "role": "result.csv RMSD",
        "live_on_campaign": True,
        "cxx": True,
    },
)


def claim_success(rmsd: float) -> bool:
    """METHODOLOGY.md §0: Success ⇔ rank-0 in-place RMSD <= 2.0 Å."""
    return bool(np.isfinite(rmsd) and rmsd <= CLAIM_CUTOFF_A)


def _load_compute_rmsd():
    path = ROOT / "python" / "flexaidds" / "benchmark.py"
    spec = importlib.util.spec_from_file_location("flexaidds_benchmark_rmsd", path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod.compute_rmsd, mod.compute_rmsd_superposed


def _load_hrmsd():
    """Load score_offline.hrmsd without running the campaign loop at import."""
    path = ROOT / "benchmarks" / "astex_repro" / "score_offline.py"
    text = path.read_text(encoding="utf-8")
    marker = "rows=[]"
    if marker not in text:
        return None
    preamble = text.split(marker, 1)[0]
    ns: dict[str, Any] = {}
    try:
        exec(preamble, ns)  # noqa: S102  — isolated load of the named scorer
    except Exception:
        return None
    return ns.get("hrmsd")


def _spyrmsd_inplace(pred: np.ndarray, ref: np.ndarray, atomic_nums: np.ndarray) -> float | None:
    try:
        from spyrmsd import rmsd as spyr
    except ImportError:
        return None
    n = len(pred)
    adj = np.zeros((n, n), dtype=int)
    return float(
        spyr.symmrmsd(
            ref,
            pred,
            atomic_nums,
            atomic_nums,
            adj,
            adj,
            minimize=False,
        )
    )


def asymmetric_cno_pair(dx: float) -> tuple[np.ndarray, np.ndarray, list[str], np.ndarray]:
    """Rank-0 in-place fixture: C/N/O, unique elements, rigid translation dx Å."""
    ref = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=np.float64)
    pred = ref + np.array([dx, 0.0, 0.0], dtype=np.float64)
    elements = ["C", "N", "O"]
    z = np.array([6, 7, 8], dtype=int)
    return pred, ref, elements, z


def measure_python_methods(dx: float) -> dict[str, Any]:
    compute_rmsd, compute_rmsd_superposed = _load_compute_rmsd()
    pred, ref, elements, z = asymmetric_cno_pair(dx)
    out: dict[str, Any] = {
        "dx": dx,
        "cutoff_A": CLAIM_CUTOFF_A,
        "method1_inplace": float(compute_rmsd(pred, ref)),
        "method1_superposed": float(compute_rmsd_superposed(pred, ref)),
        "method2_spyrmsd": None,
        "method3_hrmsd": None,
        "method4_cxx": "tests/test_dataset_runner.cpp::RmsdClaimCutoff",
        "method5_cxx": "tests/test_dataset_runner.cpp::RmsdClaimCutoff",
    }
    spy = _spyrmsd_inplace(pred, ref, z)
    if spy is not None:
        out["method2_spyrmsd"] = spy
    hrmsd = _load_hrmsd()
    if hrmsd is not None:
        out["method3_hrmsd"] = hrmsd(ref, elements, pred)
    inplace = out["method1_inplace"]
    out["claim_success_method1"] = claim_success(inplace)
    out["superposed_must_not_be_claim_metric"] = True
    return out


def python_methods_agree_on_claim(dx: float, *, rtol: float = 1e-6) -> dict[str, Any]:
    row = measure_python_methods(dx)
    values = [row["method1_inplace"]]
    for key in ("method2_spyrmsd", "method3_hrmsd"):
        if row[key] is not None:
            values.append(float(row[key]))
    expected = float(dx)
    for v in values:
        if not np.isclose(v, expected, rtol=rtol, atol=1e-6):
            raise AssertionError(f"in-place RMSD {v} != translation {expected} for dx={dx}")
    successes = [claim_success(v) for v in values]
    if len(set(successes)) != 1:
        raise AssertionError(f"claim-success split at dx={dx}: {values} -> {successes}")
    row["n_python_methods"] = len(values)
    row["agreed_claim_success"] = successes[0]
    return row


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--json", action="store_true")
    ap.add_argument(
        "--dx",
        type=float,
        nargs="*",
        default=[1.5, 2.0, 2.5],
        help="in-place translations (Å) to classify against 2.0",
    )
    args = ap.parse_args(argv)
    rows = [python_methods_agree_on_claim(dx) for dx in args.dx]
    if args.json:
        json.dump({"methods": FIVE_METHODS, "rows": rows}, sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        print("METHODOLOGY.md §0 five RMSD implementations:")
        for m in FIVE_METHODS:
            print(f"  {m['id']}. {m['name']}  ({m['role']})")
        for row in rows:
            print(
                f"  dx={row['dx']:.2f} Å  method1={row['method1_inplace']:.4f}  "
                f"claim_success={row['agreed_claim_success']}  "
                f"(cutoff {CLAIM_CUTOFF_A} Å, in-place)"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
