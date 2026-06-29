# v132_common.py — v132 manifest validation + v127 protocol env
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

from __future__ import annotations

import json
import os


REPO = "/Users/lp.more/Projects/FlexAIDdS"
ASTEX = f"{REPO}/benchmarks/astex_diverse/astex_diverse"

LANE_A_REQUIRED = (
    ("1G9V", "1G9V_apo.pdb"),
    ("1G9V", "1G9V_holo.pdb"),
    ("1TW6", "1TW6_holo.pdb"),
    ("1HNN", "1HNN_ligand_centered_site.pdb"),
)


def validate_lane_a_assets() -> None:
    missing = []
    for parts in LANE_A_REQUIRED:
        path = os.path.join(ASTEX, *parts)
        if not os.path.isfile(path):
            missing.append(path)
    if missing:
        raise FileNotFoundError(
            "Missing v132 Lane A assets:\n" + "\n".join(f"  - {p}" for p in missing)
        )


def validate_manifest(manifest_path: str) -> None:
    manifest = json.load(open(manifest_path))
    missing = []
    for pair in manifest.get("pairs", []):
        rid = pair.get("receptor_id", "?")
        for key in ("receptor_pdb", "ligand_sdf", "oracle_site_pdb"):
            path = pair.get(key, "")
            if path and not os.path.isfile(path):
                missing.append(f"{rid}.{key}={path}")
    if missing:
        raise FileNotFoundError(
            f"Manifest preflight failed ({manifest_path}):\n"
            + "\n".join(f"  - {m}" for m in missing)
        )


def v132_protocol_env(binary: str, build: str, cache: str, oracle_dir: str) -> dict[str, str]:
    """v130-class recipe: r0=4, Fix-B ON, CRG ON, sulfo binary @ HEAD."""
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
        "FLEXAIDDS_CRG":                   "1",
        "FLEXAIDDS_CRG_RMSD_MAX":          "2.5",
        "FLEXAIDDS_CRG_CF_WINDOW":         "15",
        "FLEXAIDDS_DATA_DIR":              build,
        "FLEXAIDDS_ALLOW_CONCURRENT":      "1",
        "FLEXAIDDS_BENCH_CACHE":           cache,
        "OMP_WAIT_POLICY":                 "passive",
        "OMP_PLACES":                      "cores",
        "OMP_PROC_BIND":                   "spread",
    }