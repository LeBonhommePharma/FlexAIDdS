#!/usr/bin/env python3
# lib_v132_ablation.py — shared step definitions for v132 oracle ablation ladder
#
# Audit prescription (2026-06-29):
#   Restore consensus ON for oracle-ceiling campaigns, then ablate one knob per
#   run (consensus already isolated vs v131; then binary safe vs HEAD, logsumexp
#   binary, hbond weight) — not another combined knob turn like v131.
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

REPO = Path(os.environ.get("FLEXAIDDS_REPO", "/Users/lp.more/Projects/FlexAIDdS"))
GIT_ROOT = Path(
    os.environ.get(
        "FLEXAIDDS_GIT_ROOT",
        "/Users/lp.more/.grok/worktrees/projects-flexaidds/perfornance-swarm",
    )
)
RESULTS = Path(
    os.environ.get(
        "FLEXAIDDS_RESULTS_ROOT",
        "/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results",
    )
)
ORACLE_DIR = REPO / "benchmarks/astex_diverse/astex_diverse"
JSON_V131 = GIT_ROOT / "benchmarks/datasets/benchmark_astex_native_85_v131.json"
WORKTREE_V131_SAFE = Path(
    os.environ.get(
        "FLEXAIDDS_V131_WORKTREE",
        str(REPO / "../FlexAIDdS_v131_safe"),
    )
)
WORKTREE_V127B = Path(os.environ.get("FLEXAIDDS_V127B_WORKTREE", str(REPO / "../FlexAIDdS_v127b_repro")))

REF_V130 = "v130_20260629_0548_sulfo_expB_full85"
REF_V131 = "v131_20260629_0835_r07_nofixb_full85"
REF_V127 = "v127_20260629_0139_optB_full85"
REF_V124 = "v124_full85_20260626_0413_consensus_guard"

ENV_SNAPSHOT_KEYS = (
    "FLEXAIDDS_BINARY",
    "FLEXAIDDS_ORACLE_SITE_DIR",
    "FLEXAIDDS_RESTARTS",
    "FLEXAIDDS_PARALLEL_RESTARTS",
    "FLEXAIDDS_EVAL_SCALE_DIHEDRAL",
    "FLEXAIDDS_CONSENSUS_SCORER",
    "FLEXAIDDS_SEED_ELITISM",
    "FLEXAIDDS_N_ELITE",
    "FLEXAIDDS_BUDGET_SCALE",
    "FLEXAIDDS_SOFTCORE_WAL",
    "FLEXAIDDS_SOFTCORE_FLOOR",
    "FLEXAIDDS_T_HOT",
    "FLEXAIDDS_NATIVE_SEED_FRAC",
    "FLEXAIDDS_VCT_R0",
    "FLEXAIDDS_RECEPTOR_ROTAMER_PREP",
    "FLEXAIDDS_DATA_DIR",
    "FLEXAIDDS_BENCH_CACHE",
    "FLEXAIDDS_ALLOW_CONCURRENT",
    "FLEXAIDDS_HBOND_WEIGHT",
)

STRIP_ENV_KEYS = (
    "FLEXAIDDS_USE_DP",
    "FLEXAIDDS_FINE_GRID",
    "FLEXAIDDS_FORCE_RIGID",
    "FLEXAIDDS_USE_SHANNON",
    "FLEXAIDDS_VCT_NORM",
    "FLEXAIDDS_SHARING_ALPHA",
    "FLEXAIDDS_BOOM_FRAC",
    "FLEXAIDDS_RING_FLEX",
    "FLEXAIDDS_THERMO",
    "FLEXAIDDS_HVIB",
    "FLEXAIDDS_PRIORITY_TARGETS",
    "FLEXAIDDS_FREQSEL",
)


@dataclass(frozen=True)
class AblationStep:
    step_id: str
    label: str
    description: str
    binary_path: str
    binary_src: str
    runner_path: str
    data_dir: str
    build_dir: str
    git_cwd: str
    json_pairs: str
    cache_suffix: str
    ablation_knob: str
    ablation_delta_vs_consensus_on: str
    reference_runs: dict
    audit_rationale: str
    env_overrides: dict = field(default_factory=dict)
    build_script: str | None = None
    min_commit: str | None = None
    success_gate: str = "compare_vs_v130_v131"


def resolve_v131_safe_worktree() -> Path:
    for path in (
        WORKTREE_V131_SAFE,
        Path("/Users/lp.more/.grok/worktrees/projects-flexaidds/FlexAIDdS_v131_safe"),
        REPO / "../FlexAIDdS_v131_safe",
    ):
        if (path / "build_lto").is_dir() or path.is_dir():
            return path.resolve()
    return WORKTREE_V131_SAFE.resolve()


def resolve_v127b_worktree() -> Path:
    for path in (
        WORKTREE_V127B,
        Path("/Users/lp.more/.grok/worktrees/projects-flexaidds/FlexAIDdS_v127b_repro"),
        REPO / "../FlexAIDdS_v127b_repro",
    ):
        if (path / "build_lto").is_dir() or path.is_dir():
            return path.resolve()
    return WORKTREE_V127B.resolve()


def base_oracle_env(*, binary: str, build: str, data_dir: str, cache: str) -> dict:
    return {
        "FLEXAIDDS_BINARY": binary,
        "FLEXAIDDS_BUILD": build,
        "FLEXAIDDS_REPO": str(REPO),
        "FLEXAIDDS_ORACLE_SITE_DIR": str(ORACLE_DIR),
        "FLEXAIDDS_RESTARTS": "5",
        "FLEXAIDDS_PARALLEL_RESTARTS": "1",
        "FLEXAIDDS_EVAL_SCALE_DIHEDRAL": "1",
        "FLEXAIDDS_CONSENSUS_SCORER": "1",
        "FLEXAIDDS_SEED_ELITISM": "1",
        "FLEXAIDDS_N_ELITE": "1",
        "FLEXAIDDS_BUDGET_SCALE": "1",
        "FLEXAIDDS_SOFTCORE_WAL": "1",
        "FLEXAIDDS_SOFTCORE_FLOOR": "0.5",
        "FLEXAIDDS_T_HOT": "500",
        "FLEXAIDDS_NATIVE_SEED_FRAC": "0.90",
        "FLEXAIDDS_VCT_R0": "4",
        "FLEXAIDDS_RECEPTOR_ROTAMER_PREP": "0",
        "FLEXAIDDS_DATA_DIR": data_dir,
        "FLEXAIDDS_ALLOW_CONCURRENT": "1",
        "FLEXAIDDS_BENCH_CACHE": cache,
        "OMP_WAIT_POLICY": "passive",
        "OMP_PLACES": "cores",
        "OMP_PROC_BIND": "spread",
    }


def build_steps() -> list[AblationStep]:
    head_build = str(REPO / "build_lto")
    head_binary_src = f"{head_build}/FlexAIDdS"
    head_runner = f"{head_build}/benchmark_datasets"
    safe_root = resolve_v131_safe_worktree()
    safe_build = str(safe_root / "build_lto")
    v127b_root = resolve_v127b_worktree()
    v127b_build = str(v127b_root / "build_lto")

    consensus_on_refs = {
        "v130_consensus_on": REF_V130,
        "v131_consensus_off": REF_V131,
        "v127_optB": REF_V127,
        "v124_guard": REF_V124,
    }

    return [
        AblationStep(
            step_id="consensus_on",
            label="v132a_consensus_on_full85",
            description=(
                "HEAD binary + v131 data pipeline (holo 1TW6, 1HNN expB, sulfo remap) "
                "with v130 oracle protocol: consensus ON, r0=4. Primary audit candidate "
                "for highest honest oracle-ceiling number."
            ),
            binary_path="/tmp/FlexAIDdS_v132a",
            binary_src=head_binary_src,
            runner_path=head_runner,
            data_dir=head_build,
            build_dir=head_build,
            git_cwd=str(REPO),
            json_pairs=str(JSON_V131),
            cache_suffix="v132a_consensus_on",
            ablation_knob="consensus",
            ablation_delta_vs_consensus_on="baseline step (consensus=1, r0=4)",
            reference_runs=consensus_on_refs,
            audit_rationale=(
                "Scientific robustness audit: restore CONSENSUS_SCORER=1 for oracle-ceiling "
                "campaigns. Isolates v131's consensus=OFF + r0=7 combined turn."
            ),
            min_commit="04ff1735",
        ),
        AblationStep(
            step_id="safe_binary",
            label="v132b_safe_binary_full85",
            description=(
                "v131_safe binary (82ad51f4 + sulfo + holo data cherry-picks, pre-P0 gaboom) "
                "with identical consensus_on env. Isolates HEAD scoring/H-bond/VCT vs safe tree."
            ),
            binary_path="/tmp/FlexAIDdS_v132b",
            binary_src=f"{safe_build}/FlexAIDdS",
            runner_path=f"{safe_build}/benchmark_datasets",
            data_dir=safe_build,
            build_dir=safe_build,
            git_cwd=str(safe_root),
            json_pairs=str(JSON_V131),
            cache_suffix="v132b_safe_binary",
            ablation_knob="binary_safe_vs_head",
            ablation_delta_vs_consensus_on="binary → v131_safe worktree only",
            reference_runs={
                "v132a_consensus_on": "prior_ladder_step",
                **consensus_on_refs,
            },
            audit_rationale=(
                "Ablation ladder: binary safe (pre-27e68e51 / pre-ba5364d3) vs HEAD while "
                "holding consensus ON and v131 data fixes constant."
            ),
            build_script=str(GIT_ROOT / "scripts/build_v131_safe.sh"),
        ),
        AblationStep(
            step_id="logsumexp_only",
            label="v132c_logsumexp_only_full85",
            description=(
                "a4056163 logsumexp-only binary (no ba5364d3 H-bond/VCT patch) with "
                "consensus_on env + v131 JSON. Isolates H-bond/VCT scoring kernel delta."
            ),
            binary_path="/tmp/FlexAIDdS_v132c",
            binary_src=f"{v127b_build}/FlexAIDdS",
            runner_path=f"{v127b_build}/benchmark_datasets",
            data_dir=v127b_build,
            build_dir=v127b_build,
            git_cwd=str(v127b_root),
            json_pairs=str(JSON_V131),
            cache_suffix="v132c_logsumexp_only",
            ablation_knob="logsumexp_vs_hbond_vct",
            ablation_delta_vs_consensus_on="binary → a4056163 logsumexp-only worktree",
            reference_runs={
                "v132a_consensus_on": "prior_ladder_step",
                "v127b_smoke": "v127b 3-target smoke (if run)",
                **consensus_on_refs,
            },
            audit_rationale=(
                "Reproduction-gap audit step 7: ablate H-bond/VCT (ba5364d3) while keeping "
                "logsumexp selector; mirrors launch_v127b_smoke at full-85 scale."
            ),
            build_script=str(GIT_ROOT / "scripts/build_v127b_logsumexp.sh"),
        ),
        AblationStep(
            step_id="hbond_zero",
            label="v132d_hbond_weight_zero_full85",
            description=(
                "HEAD binary + consensus_on env + FLEXAIDDS_HBOND_WEIGHT=0. Approximate "
                "hbond-off ablation without dock_config surgery (v122-style probe)."
            ),
            binary_path="/tmp/FlexAIDdS_v132d",
            binary_src=head_binary_src,
            runner_path=head_runner,
            data_dir=head_build,
            build_dir=head_build,
            git_cwd=str(REPO),
            json_pairs=str(JSON_V131),
            cache_suffix="v132d_hbond_zero",
            ablation_knob="hbond",
            ablation_delta_vs_consensus_on="FLEXAIDDS_HBOND_WEIGHT=0 only",
            reference_runs={
                "v132a_consensus_on": "prior_ladder_step",
                "v122_hbond_off": "historical v122 protocol",
            },
            audit_rationale=(
                "Scientific robustness audit: hbond-off ablation if oracle rate regresses. "
                "Uses HBOND_WEIGHT=0 env lever (search path still enabled; document caveat)."
            ),
            env_overrides={"FLEXAIDDS_HBOND_WEIGHT": "0"},
            min_commit="04ff1735",
        ),
    ]


def step_by_id(step_id: str) -> AblationStep:
    for step in build_steps():
        if step.step_id == step_id:
            return step
    known = ", ".join(s.step_id for s in build_steps())
    raise SystemExit(f"ERROR: unknown step '{step_id}'. Choose from: {known}")


def ladder_step_ids() -> list[str]:
    return [s.step_id for s in build_steps()]