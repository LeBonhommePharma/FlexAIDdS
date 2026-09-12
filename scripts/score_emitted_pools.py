#!/usr/bin/env python3
"""Score saved cluster-head files with an explicit per-restart pose budget.

The supplied scoring module owns RMSD, ligand matching, sentinel rejection and
deduplication. This adapter changes only the enumeration budget. It never runs
an engine and does not infer the budget from the current shell environment.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import offline_elector as oe


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_scoring_code(path):
    """Load an explicitly chosen trusted scorer, retaining its RMSD convention."""
    path = Path(path).resolve(strict=True)
    digest = sha256_file(path)
    spec = importlib.util.spec_from_file_location("_emitted_pool_scorer_" + digest[:16], path)
    if spec is None or spec.loader is None:
        raise ValueError("cannot import scoring code: " + str(path))
    module = importlib.util.module_from_spec(spec)
    # No writes, including __pycache__, in a completed campaign's code directory.
    previous = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        spec.loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = previous
    if not callable(getattr(module, "Reference", None)) or not callable(getattr(module, "pool_tiers", None)):
        raise ValueError("scoring code must expose Reference and pool_tiers")
    return module, digest


def score_saved_pool(prefixes, reference, scoring_code, pose_budget=50, top_n=10):
    """Keep the explicit budget through cell -> pool -> every restart enumerator."""
    oe.validate_pose_budget(pose_budget)
    if type(top_n) is not int or top_n < 1:
        raise ValueError("top_n must be a positive integer")
    prefixes = [str(Path(p).resolve()) for p in prefixes]
    if not prefixes or len(set(prefixes)) != len(prefixes):
        raise ValueError("provide at least one unique restart prefix")
    reference = Path(reference).resolve(strict=True)
    scorer, score_digest = load_scoring_code(scoring_code)
    poses, stats = oe.build_pool(prefixes, budget=pose_budget)
    ref = scorer.Reference(str(reference))
    tiers = scorer.pool_tiers(poses, ref, lambda pose: pose["cf"], top_n=top_n)
    truncation = any(s["cf_truncated"] or s["fo_truncated"] for s in stats)
    parse_failures = sum(len(s["parse_errors"]) for s in stats)
    empty_prefixes = [s["prefix"] for s in stats if s["enumerated"] == 0]
    # pool_tiers may intentionally filter/dedup and may silently omit bad poses.
    # A count is observable; the reason for each downstream drop is not inferred.
    n_scored = tiers["n_pool"]
    if type(n_scored) is not int or not 0 <= n_scored <= len(poses):
        raise ValueError("scorer returned an invalid pool size")
    raw = sum(s["enumerated"] for s in stats)
    assert raw == len(poses) + parse_failures
    raw_paths = [p["path"] for p in poses]
    raw_paths.extend(e["path"] for s in stats for e in s["parse_errors"])
    input_receipts = []
    for path in raw_paths:
        sidecar = Path(path).with_suffix(".mcf")
        input_receipts.append({
            "path": path, "sha256": sha256_file(path),
            "mcf_sha256": sha256_file(sidecar) if sidecar.exists() else None,
        })
    return {
        "schema": "flexaidds-offline-pool-budget/1",
        "pose_budget": pose_budget,
        "pose_budget_source": "explicit_argument",
        "prefixes": prefixes,
        "enum_stats": stats,
        "n_enumerated": raw,
        "n_parsed": len(poses),
        "n_parse_failed": parse_failures,
        "n_pool_scored": n_scored,
        "n_not_in_scored_pool": len(poses) - n_scored,
        "empty_prefixes": empty_prefixes,
        "offline_enumerator_truncated": truncation,
        "status": ("missing_or_empty_prefixes" if empty_prefixes else
                   "parse_failures" if parse_failures else
                   "scored" if n_scored else "empty_pool"),
        "top_n": top_n,
        "oracle_spyrmsd_A": tiers["oracle"],
        "oracle_naive_A": tiers["oracle_naive"],
        "topn_spyrmsd_A": tiers["top_n"],
        "topn_naive_A": tiers["top_n_naive"],
        "score_code": str(Path(scoring_code).resolve()),
        "score_code_sha256": score_digest,
        "reference": str(reference),
        "reference_sha256": sha256_file(reference),
        "enumerator_sha256": sha256_file(oe.__file__),
        "adapter_sha256": sha256_file(__file__),
        "pose_inputs": input_receipts,
        "interpretation": (
            "Saved cluster-head pool only. An offline truncation flag cannot "
            "prove absence of upstream write/cluster caps. Downstream filtering "
            "and failure reasons remain the supplied scorer's responsibility."
        ),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", action="append", required=True,
                        help="Emitted restart prefix, repeat in original restart order")
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--scoring-code", type=Path, required=True,
                        help="Trusted frozen score_electors.py (Reference + pool_tiers)")
    parser.add_argument("--pose-budget", type=int, default=oe.POSE_BUDGET,
                        help="Per-restart CF/DP+FO budget, default 50; range 1..5000")
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--out", type=Path, required=True,
                        help="New JSON output path; existing files are never overwritten")
    args = parser.parse_args(argv)
    try:
        oe.validate_pose_budget(args.pose_budget)
        if args.out.exists():
            raise ValueError("output already exists: " + str(args.out))
        if not args.out.parent.is_dir():
            raise ValueError("output directory does not exist: " + str(args.out.parent))
        result = score_saved_pool(args.prefix, args.reference, args.scoring_code,
                                  pose_budget=args.pose_budget, top_n=args.top_n)
        encoded = json.dumps(result, indent=2, allow_nan=False) + "\n"
        with args.out.open("x") as handle:
            handle.write(encoded)
    except (ValueError, OSError, ImportError) as exc:
        parser.error(str(exc))
    print(f"budget={result['pose_budget']} enumerated={result['n_enumerated']} "
          f"parsed={result['n_parsed']} scored={result['n_pool_scored']} out={args.out}")
    return 0 if result["status"] == "scored" else 2


if __name__ == "__main__":
    raise SystemExit(main())
