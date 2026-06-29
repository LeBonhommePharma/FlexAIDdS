#!/usr/bin/env python3
# lib_worker_orders.py — Grok Build worker campaign guards (shared by launch scripts)
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO = SCRIPT_DIR.parent
ORDERS_REPO = REPO / "docs/dev/grok_build_worker_orders.json"
ORDERS_RESULTS = Path(
    os.environ.get(
        "FLEXAIDDS_RESULTS_ROOT",
        "/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results",
    )
) / "grok_build_worker_orders.json"

CAMPAIGN_SCRIPT_MAP = {
    "v132_ablation_ladder": "scripts/queue_v132_ablation_ladder.py",
    "v132_ablation_step": "scripts/launch_v132_ablation.py",
    "v131_head_full85": "scripts/launch_v131_full85.py",
}


def load_orders() -> dict:
    for path in (ORDERS_RESULTS, ORDERS_REPO):
        if path.is_file():
            return json.loads(path.read_text())
    return {}


def bisect_complete() -> bool:
    summary = Path(
        os.environ.get(
            "FLEXAIDDS_RESULTS_ROOT",
            "/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results",
        )
    ) / "vcontacts_bisect_summary.json"
    return summary.is_file()


def assert_campaign_allowed(campaign: str, *, script_name: str = "") -> None:
    """Exit non-zero if worker orders block this campaign."""
    orders = load_orders()
    if not orders:
        return

    blocked = set(orders.get("blocked_campaigns") or [])
    if campaign not in blocked:
        return

    phase = orders.get("current_phase", "unknown")
    msg = (
        f"BLOCKED by Grok Build worker orders: campaign={campaign!r} "
        f"(current_phase={phase!r}).\n"
        f"  Read: docs/dev/GROK_BUILD_WORKER_DISPATCH_20260629.md\n"
        f"  NOW: bisect Vcontacts smoke-12×3. Do NOT run v132 ablation or HEAD full-85.\n"
        f"  v131_safe binary is the escape hatch to 78/85 while bisect runs."
    )
    if script_name:
        msg += f"\n  Script: {script_name}"
    if not bisect_complete():
        msg += "\n  Gate: vcontacts_bisect_summary.json not yet written."
    print(msg, file=sys.stderr)
    raise SystemExit(2)


def assert_script_allowed(script_path: str | Path) -> None:
    """Block launch if script is on the blocked list and bisect incomplete."""
    orders = load_orders()
    if not orders or bisect_complete():
        return

    script = Path(script_path).name
    blocked_scripts = orders.get("blocked_scripts") or []
    for entry in blocked_scripts:
        if Path(entry).name == script:
            for camp, mapped in CAMPAIGN_SCRIPT_MAP.items():
                if mapped.endswith(script):
                    assert_campaign_allowed(camp, script_name=str(script_path))
            assert_campaign_allowed("combined_knob_turn", script_name=str(script_path))