# v131_safe_common.py — shared path resolution for Lane A (v131_safe) launchers
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

from __future__ import annotations

import os
import subprocess
from typing import Iterable


REPO = "/Users/lp.more/Projects/FlexAIDdS"


def git_root(script_dir: str) -> str:
    return subprocess.check_output(
        ["git", "-C", script_dir, "rev-parse", "--show-toplevel"], text=True
    ).strip()


def resolve_worktree() -> str:
    candidates = [
        os.environ.get("FLEXAIDDS_V131_WORKTREE", ""),
        f"{REPO}/../FlexAIDdS_v131_safe",
        "/Users/lp.more/.grok/worktrees/projects-flexaidds/FlexAIDdS_v131_safe",
    ]
    for path in candidates:
        if path and os.path.isdir(os.path.join(path, "build_lto")):
            return path
    for path in candidates:
        if path and os.path.isdir(path):
            return path
    return f"{REPO}/../FlexAIDdS_v131_safe"


def asset_roots(worktree: str, git_root_path: str) -> list[str]:
    roots: list[str] = []
    for root in (git_root_path, worktree, REPO):
        if root and root not in roots:
            roots.append(root)
    return roots


def resolve_oracle_dir(worktree: str, git_root_path: str) -> str:
    rel = "benchmarks/astex_diverse/astex_diverse"
    for root in asset_roots(worktree, git_root_path):
        path = os.path.join(root, rel)
        if os.path.isdir(path):
            return path
    return os.path.join(REPO, rel)


def resolve_oracle_asset(worktree: str, git_root_path: str, *parts: str) -> str:
    for root in asset_roots(worktree, git_root_path):
        path = os.path.join(root, "benchmarks", "astex_diverse", "astex_diverse", *parts)
        if os.path.exists(path):
            return path
    return os.path.join(REPO, "benchmarks", "astex_diverse", "astex_diverse", *parts)


def patch_pair_paths(pair: dict, worktree: str, git_root_path: str) -> dict:
    """Rewrite manifest paths for bf8cf1d2 holo/data targets when Projects tree is stale."""
    rid = pair.get("receptor_id", "")
    out = dict(pair)
    astex = resolve_oracle_dir(worktree, git_root_path)

    def astex_path(*parts: str) -> str:
        resolved = resolve_oracle_asset(worktree, git_root_path, *parts)
        if os.path.exists(resolved):
            return resolved
        return os.path.join(astex, *parts)

    if rid == "1G9V":
        out["receptor_pdb"] = astex_path("1G9V", "1G9V_apo.pdb")
        out["ligand_sdf"] = astex_path("1G9V", "1G9V_ligand.sdf")
        out["oracle_site_pdb"] = astex_path("1G9V", "1G9V_binding_site.pdb")
    elif rid == "1TW6":
        out["receptor_pdb"] = astex_path("1TW6", "1TW6_holo.pdb")
        out["ligand_sdf"] = astex_path("1TW6", "1TW6_ligand.sdf")
        out["oracle_site_pdb"] = astex_path("1TW6", "1TW6_binding_site.pdb")
    elif rid == "1HNN":
        out["receptor_pdb"] = astex_path("1HNN", "1HNN_apo.pdb")
        out["ligand_sdf"] = astex_path("1HNN", "1HNN_ligand.sdf")
        out["oracle_site_pdb"] = astex_path("1HNN", "1HNN_ligand_centered_site.pdb")
    return out


def patch_manifest(manifest: dict, worktree: str, git_root_path: str) -> dict:
    out = dict(manifest)
    out["pairs"] = [
        patch_pair_paths(p, worktree, git_root_path) for p in manifest["pairs"]
    ]
    out["astex_diverse_dir"] = resolve_oracle_dir(worktree, git_root_path)
    return out


def validate_lane_a_assets(worktree: str, git_root_path: str) -> None:
    required = (
        ("1G9V", "1G9V_apo.pdb"),
        ("1TW6", "1TW6_holo.pdb"),
        ("1HNN", "1HNN_ligand_centered_site.pdb"),
    )
    missing = []
    for parts in required:
        path = resolve_oracle_asset(worktree, git_root_path, *parts)
        if not os.path.isfile(path):
            missing.append(path)
    if missing:
        raise FileNotFoundError(
            "Missing v131 Lane A assets (run build_v131_safe.sh / sync bf8cf1d2 data):\n"
            + "\n".join(f"  - {p}" for p in missing)
        )


def v127_protocol_env(
    binary: str,
    build: str,
    cache: str,
    oracle_dir: str,
) -> dict[str, str]:
    return {
        "FLEXAIDDS_BINARY":                binary,
        "FLEXAIDDS_BUILD":                 build,
        "FLEXAIDDS_REPO":                  REPO,
        "FLEXAIDDS_ORACLE_SITE_DIR":       oracle_dir,
        "FLEXAIDDS_RESTARTS":              "5",
        "FLEXAIDDS_PARALLEL_RESTARTS":     "1",
        "FLEXAIDDS_EVAL_SCALE_DIHEDRAL":   "1",
        "FLEXAIDDS_CONSENSUS_SCORER":      "1",
        "FLEXAIDDS_SEED_ELITISM":          "1",
        "FLEXAIDDS_N_ELITE":               "1",
        "FLEXAIDDS_BUDGET_SCALE":          "1",
        "FLEXAIDDS_SOFTCORE_WAL":          "1",
        "FLEXAIDDS_SOFTCORE_FLOOR":        "0.5",
        "FLEXAIDDS_T_HOT":                 "500",
        "FLEXAIDDS_NATIVE_SEED_FRAC":      "0.90",
        "FLEXAIDDS_VCT_R0":                "4",
        "FLEXAIDDS_RECEPTOR_ROTAMER_PREP": "0",
        "FLEXAIDDS_DATA_DIR":              build,
        "FLEXAIDDS_ALLOW_CONCURRENT":      "1",
        "FLEXAIDDS_BENCH_CACHE":           cache,
        "OMP_WAIT_POLICY":                 "passive",
        "OMP_PLACES":                      "cores",
        "OMP_PROC_BIND":                   "spread",
    }


ENV_POP_KEYS: tuple[str, ...] = (
    "FLEXAIDDS_USE_DP", "FLEXAIDDS_FINE_GRID",
    "FLEXAIDDS_FORCE_RIGID", "FLEXAIDDS_USE_SHANNON",
    "FLEXAIDDS_VCT_NORM",
    "FLEXAIDDS_SHARING_ALPHA", "FLEXAIDDS_BOOM_FRAC",
    "FLEXAIDDS_RING_FLEX",
    "FLEXAIDDS_THERMO", "FLEXAIDDS_HVIB",
    "FLEXAIDDS_PRIORITY_TARGETS", "FLEXAIDDS_FREQSEL",
)


def scrub_env(env: dict) -> dict:
    out = dict(env)
    for k in ENV_POP_KEYS:
        out.pop(k, None)
    return out