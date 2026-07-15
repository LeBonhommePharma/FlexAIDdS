#!/usr/bin/env python3
"""Deck-style red-pair barplots from LIVE 3dsig_r10 docks (not archived medians).

Scans three-engine arm OUTs for per-target result.csv / pose PDBs, computes
S_top10 + 10k bootstrap medians via bootstrap_3dsig_s_top10 helpers, and writes
pure-SVG bars (FlexAID vs FlexAIDdS). Competitors are optional and only drawn
when a competitors JSON supplies values.

Works with partial arms:
  - zero complete cases → status JSON + empty-aware SVG (exit 0)
  - A only / B0 only / B missing → plot available red bars, flag incomplete
  - --demo → synthetic pilot8-like rates for layout smoke (not a claim)

Usage (after A/B0/B finish on pilot8 or full85):

  export FLEXAIDDS_LOCAL_ROOT="${FLEXAIDDS_LOCAL_ROOT:-$HOME/flexaidds_results}"
  ROOT="${FLEXAIDDS_ROOT:-$HOME/Projects/FlexAIDdS}"
  TE="$FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine"

  python3 "$ROOT/scripts/plot_3dsig_live_bars.py" \\
    --campaign-root "$TE" \\
    --campaign 3dsig_r10 \\
    --out-dir "$ROOT/docs/figures/3dsig_live" \\
    --bootstraps 10000

  # Layout smoke without waiting for docks:
  python3 scripts/plot_3dsig_live_bars.py --demo --out-dir /tmp/3dsig_live_demo

Copyright 2026 Le Bonhomme Pharma
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple
from xml.sax.saxutils import escape

# ---------------------------------------------------------------------------
# Import bootstrap helpers without requiring package install
# ---------------------------------------------------------------------------
_SCRIPT_DIR = Path(__file__).resolve().parent
_BOOT_PATH = _SCRIPT_DIR / "bootstrap_3dsig_s_top10.py"


def _load_bootstrap_mod():
    spec = importlib.util.spec_from_file_location("bootstrap_3dsig_s_top10", _BOOT_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {_BOOT_PATH}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# Colors match archived deck: red pair first
COLOR_FLEXAID = "#E53935"
COLOR_FLEXAIDDS = "#E53935"
COLOR_VINA = "#1E88E5"
COLOR_FLEXX = "#00ACC1"
COLOR_RDOCK = "#8E24AA"
COLOR_OTHER = "#607D8B"

COMPETITOR_COLORS = {
    "vina": COLOR_VINA,
    "flexx": COLOR_FLEXX,
    "rdock": COLOR_RDOCK,
    "rDock": COLOR_RDOCK,
}


def default_campaign_root() -> Path:
    local = Path.home() / "flexaidds_results"
    env_local = Path(
        __import__("os").environ.get("FLEXAIDDS_LOCAL_ROOT", str(local))
    ).expanduser()
    # Prefer local-first three_engine tree
    for base in (env_local, local):
        te = base / "campaigns" / "three_engine"
        if te.is_dir():
            return te
    return env_local / "campaigns" / "three_engine"


def inventory_arm(arm_dir: Path) -> Dict[str, Any]:
    """Shallow inventory of arm_dir (no recursive CloudDocs walk)."""
    info: Dict[str, Any] = {
        "arm_dir": str(arm_dir),
        "exists": arm_dir.is_dir(),
        "case_dirs": [],
        "result_csv": [],
        "result_csv_with_rmsd": [],
        "result_csv_empty": [],
        "pose_pdb_cases": [],
        "receipt": None,
    }
    if not arm_dir.is_dir():
        return info
    receipt = arm_dir / "RUN_RECEIPT.json"
    if receipt.is_file():
        try:
            info["receipt"] = json.loads(receipt.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            info["receipt"] = {"error": str(exc)}
    for child in sorted(arm_dir.iterdir()):
        if not child.is_dir():
            continue
        name = child.name
        if len(name) < 4:
            continue
        info["case_dirs"].append(name)
        rc = child / "result.csv"
        if rc.is_file():
            info["result_csv"].append(name)
            # quick check: any numeric rmsd-like field
            try:
                text = rc.read_text(errors="replace")
                has = False
                for line in text.splitlines()[1:]:
                    parts = line.split(",")
                    # rmsd_top1 is col index 4 in current schema; also mode_rmsd_*
                    if "mode_rmsd_" in text.splitlines()[0]:
                        hdr = text.splitlines()[0].split(",")
                        for i, h in enumerate(hdr):
                            if h.startswith("mode_rmsd_") and i < len(parts):
                                if parts[i].strip() not in ("", "NA", "None"):
                                    has = True
                                    break
                    # fallback: non-empty rmsd_top1 / rmsd_bcr
                    if not has and len(parts) > 5:
                        # arm,engine,matrix,pdb,rmsd_top1,rmsd_bcr,...
                        if parts[4].strip() not in ("", "NA") or parts[5].strip() not in (
                            "",
                            "NA",
                        ):
                            has = True
                if has:
                    info["result_csv_with_rmsd"].append(name)
                else:
                    info["result_csv_empty"].append(name)
            except OSError:
                info["result_csv_empty"].append(name)
        # pose PDBs (ranked, not INI)
        poses = [
            p
            for p in child.glob(f"{name}_r*_*.pdb")
            if not p.name.endswith("_INI.pdb")
        ]
        if not poses:
            poses = [
                p
                for p in child.glob("*_[0-9]*.pdb")
                if not p.name.endswith("_INI.pdb") and "_r" in p.name
            ]
        if poses:
            info["pose_pdb_cases"].append(name)
    return info


def arm_stats(
    boot,
    arm_dir: Path,
    n_boot: int,
    seed: int,
    thresh: float,
) -> Dict[str, Any]:
    cases = boot.load_arm_dir(arm_dir) if arm_dir.is_dir() else {}
    success = [boot.s_top10(v, thresh) for v in cases.values()]
    stats = boot.bootstrap_median(success, n_boot, seed)
    stats["metric"] = "S_top10"
    stats["thresh_A"] = thresh
    stats["n_pdbs_loaded"] = len(cases)
    stats["pdb_ids"] = sorted(cases.keys())
    stats["case_success"] = {
        pid: bool(boot.s_top10(rms, thresh)) for pid, rms in sorted(cases.items())
    }
    stats["arm_dir"] = str(arm_dir)
    stats["complete"] = bool(stats["n_cases"] and stats["median"] is not None)
    return stats


def svg_bars(
    title: str,
    methods: Sequence[str],
    vals: Sequence[Optional[float]],
    colors: Sequence[str],
    footer: str,
    subtitle: str = "",
    p05: Optional[Sequence[Optional[float]]] = None,
    p95: Optional[Sequence[Optional[float]]] = None,
) -> str:
    W, H = 720, 440
    L, R, T, B = 70, 30, 58 if subtitle else 50, 100
    plot_w = W - L - R
    plot_h = H - T - B
    n = len(methods)
    gap = 12
    bar_w = (plot_w - gap * (n + 1)) / max(n, 1)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
        f'viewBox="0 0 {W} {H}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{W/2}" y="28" text-anchor="middle" font-family="Helvetica,Arial,sans-serif" '
        f'font-size="15" font-weight="bold">{escape(title)}</text>',
    ]
    if subtitle:
        parts.append(
            f'<text x="{W/2}" y="46" text-anchor="middle" font-size="11" '
            f'fill="#555" font-family="Helvetica,Arial,sans-serif">{escape(subtitle)}</text>'
        )
    parts += [
        f'<line x1="{L}" y1="{T}" x2="{L}" y2="{T+plot_h}" stroke="#222" stroke-width="1.5"/>',
        f'<line x1="{L}" y1="{T+plot_h}" x2="{L+plot_w}" y2="{T+plot_h}" stroke="#222" stroke-width="1.5"/>',
    ]
    for yv in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
        y = T + plot_h * (1.0 - yv)
        parts.append(
            f'<line x1="{L}" y1="{y:.1f}" x2="{L+plot_w}" y2="{y:.1f}" '
            f'stroke="#eee" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{L-8}" y="{y+4:.1f}" text-anchor="end" font-size="11" '
            f'font-family="Helvetica,Arial,sans-serif" fill="#444">{yv:.1f}</text>'
        )
    parts.append(
        f'<text x="18" y="{T+plot_h/2}" transform="rotate(-90 18 {T+plot_h/2})" '
        f'text-anchor="middle" font-size="12" font-family="Helvetica,Arial,sans-serif">'
        f"Success rate (bootstrap median)</text>"
    )

    plotted_vals: List[float] = []
    for i, (m, v, c) in enumerate(zip(methods, vals, colors)):
        x = L + gap + i * (bar_w + gap)
        if v is None:
            # empty bar placeholder
            parts.append(
                f'<rect x="{x:.1f}" y="{T:.1f}" width="{bar_w:.1f}" height="{plot_h:.1f}" '
                f'fill="#f5f5f5" stroke="#bbb" stroke-width="0.8" stroke-dasharray="4 3"/>'
            )
            parts.append(
                f'<text x="{x+bar_w/2:.1f}" y="{T+plot_h/2:.1f}" text-anchor="middle" '
                f'font-size="11" fill="#888" font-family="Helvetica,Arial,sans-serif">n/a</text>'
            )
        else:
            plotted_vals.append(v)
            bh = plot_h * max(0.0, min(1.0, v))
            y = T + plot_h - bh
            parts.append(
                f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{bh:.1f}" '
                f'fill="{c}" stroke="#111" stroke-width="0.8"/>'
            )
            parts.append(
                f'<text x="{x+bar_w/2:.1f}" y="{y-6:.1f}" text-anchor="middle" '
                f'font-size="12" font-weight="bold" font-family="Helvetica,Arial,sans-serif">'
                f"{v:.2f}</text>"
            )
            # optional CI whiskers
            if p05 and p95 and i < len(p05) and i < len(p95):
                lo, hi = p05[i], p95[i]
                if lo is not None and hi is not None:
                    y_lo = T + plot_h * (1.0 - max(0.0, min(1.0, lo)))
                    y_hi = T + plot_h * (1.0 - max(0.0, min(1.0, hi)))
                    cx = x + bar_w / 2
                    parts.append(
                        f'<line x1="{cx:.1f}" y1="{y_hi:.1f}" x2="{cx:.1f}" y2="{y_lo:.1f}" '
                        f'stroke="#111" stroke-width="1.2"/>'
                    )
                    parts.append(
                        f'<line x1="{cx-6:.1f}" y1="{y_hi:.1f}" x2="{cx+6:.1f}" y2="{y_hi:.1f}" '
                        f'stroke="#111" stroke-width="1.2"/>'
                    )
                    parts.append(
                        f'<line x1="{cx-6:.1f}" y1="{y_lo:.1f}" x2="{cx+6:.1f}" y2="{y_lo:.1f}" '
                        f'stroke="#111" stroke-width="1.2"/>'
                    )
        parts.append(
            f'<text x="{x+bar_w/2:.1f}" y="{T+plot_h+16}" text-anchor="middle" '
            f'font-size="11" font-family="Helvetica,Arial,sans-serif" '
            f'transform="rotate(20 {x+bar_w/2:.1f} {T+plot_h+16})">{escape(m)}</text>'
        )

    # Δ line when first two red-pair values present
    if len(vals) >= 2 and vals[0] is not None and vals[1] is not None:
        delta = vals[1] - vals[0]
        parts.append(
            f'<text x="{W/2}" y="{H-36}" text-anchor="middle" font-size="12" '
            f'fill="#B71C1C" font-family="Helvetica,Arial,sans-serif">'
            f"Δ FlexAIDdS−FlexAID = {delta:+.2f}  ·  red = FlexAID pair</text>"
        )
    elif not plotted_vals:
        parts.append(
            f'<text x="{W/2}" y="{H-36}" text-anchor="middle" font-size="12" '
            f'fill="#666" font-family="Helvetica,Arial,sans-serif">'
            f"No complete S_top10 cases yet — waiting on docks</text>"
        )

    parts.append(
        f'<text x="8" y="{H-10}" font-size="9" fill="#555" font-family="Helvetica,Arial,sans-serif">'
        f"{escape(footer)}</text>"
    )
    parts.append("</svg>")
    return "\n".join(parts)


def load_competitors(path: Optional[Path]) -> List[Tuple[str, float, str]]:
    """Optional competitors JSON: [{name, median, color?}, ...] or {name: median}."""
    if path is None or not path.is_file():
        return []
    data = json.loads(path.read_text())
    out: List[Tuple[str, float, str]] = []
    if isinstance(data, dict) and "methods" in data:
        data = data["methods"]
    if isinstance(data, dict):
        for name, val in data.items():
            if name.startswith("_"):
                continue
            if isinstance(val, dict):
                med = float(val["median"])
                col = val.get("color") or COMPETITOR_COLORS.get(name, COLOR_OTHER)
            else:
                med = float(val)
                col = COMPETITOR_COLORS.get(name, COLOR_OTHER)
            out.append((str(name), med, col))
    elif isinstance(data, list):
        for row in data:
            name = str(row.get("name") or row.get("method"))
            med = float(row["median"])
            col = row.get("color") or COMPETITOR_COLORS.get(name, COLOR_OTHER)
            out.append((name, med, col))
    return out


def choose_flexaid_arm(
    stats_a: Dict[str, Any],
    stats_b0: Dict[str, Any],
    prefer: str,
) -> Tuple[str, Dict[str, Any]]:
    """Pick FlexAID control arm: prefer A (JCIM) unless --flexaid-arm B0 or A empty."""
    prefer = prefer.upper()
    if prefer == "B0":
        return "B0", stats_b0
    if prefer == "A":
        if stats_a.get("complete") or not stats_b0.get("complete"):
            return "A", stats_a
        return "B0", stats_b0
    # auto
    if stats_a.get("complete"):
        return "A", stats_a
    if stats_b0.get("complete"):
        return "B0", stats_b0
    # partial: prefer whichever has more loaded cases
    if (stats_a.get("n_cases") or 0) >= (stats_b0.get("n_cases") or 0):
        return "A", stats_a
    return "B0", stats_b0


def demo_stats(seed: int = 20170715) -> Dict[str, Dict[str, Any]]:
    """Synthetic pilot8-like rates for layout dry-run (NOT a claim)."""
    # Fixed demo values near archived Astex pair for visual check
    return {
        "A": {
            "n_cases": 8,
            "n_success": 5,
            "point": 0.625,
            "median": 0.625,
            "p05": 0.375,
            "p95": 0.875,
            "n_boot": 10000,
            "complete": True,
            "pdb_ids": ["DEMO1", "DEMO2", "DEMO3", "DEMO4", "DEMO5", "DEMO6", "DEMO7", "DEMO8"],
            "note": "synthetic --demo",
        },
        "B0": {
            "n_cases": 0,
            "n_success": 0,
            "point": None,
            "median": None,
            "p05": None,
            "p95": None,
            "n_boot": 10000,
            "complete": False,
            "pdb_ids": [],
            "note": "synthetic --demo unused",
        },
        "B": {
            "n_cases": 8,
            "n_success": 6,
            "point": 0.750,
            "median": 0.750,
            "p05": 0.500,
            "p95": 1.000,
            "n_boot": 10000,
            "complete": True,
            "pdb_ids": ["DEMO1", "DEMO2", "DEMO3", "DEMO4", "DEMO5", "DEMO6", "DEMO7", "DEMO8"],
            "note": "synthetic --demo",
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--campaign-root",
        type=Path,
        default=None,
        help="three_engine root containing A/B0/B (default: $FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine)",
    )
    ap.add_argument("--campaign", default="3dsig_r10", help="Campaign subdir name (default 3dsig_r10)")
    ap.add_argument("--arm-a", type=Path, default=None, help="Override path to arm A OUT")
    ap.add_argument("--arm-b0", type=Path, default=None, help="Override path to arm B0 OUT")
    ap.add_argument("--arm-b", type=Path, default=None, help="Override path to arm B OUT")
    ap.add_argument(
        "--flexaid-arm",
        choices=("auto", "A", "B0"),
        default="auto",
        help="Which arm labels as FlexAID red bar (default auto: A if complete else B0)",
    )
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument("--bootstraps", type=int, default=10_000)
    ap.add_argument("--seed", type=int, default=20170715)
    ap.add_argument("--thresh", type=float, default=2.0)
    ap.add_argument(
        "--competitors",
        type=Path,
        default=None,
        help="Optional JSON of competitor bootstrap medians (only plotted if present)",
    )
    ap.add_argument(
        "--demo",
        action="store_true",
        help="Synthetic pilot8-like rates for layout smoke (not a live claim)",
    )
    ap.add_argument(
        "--inventory-only",
        action="store_true",
        help="Print inventory JSON and exit (no plot)",
    )
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[1]
    out_dir = args.out_dir or (root / "docs" / "figures" / "3dsig_live")
    out_dir.mkdir(parents=True, exist_ok=True)

    camp_root = args.campaign_root or default_campaign_root()
    arm_a = args.arm_a or (camp_root / "A" / args.campaign)
    arm_b0 = args.arm_b0 or (camp_root / "B0" / args.campaign)
    arm_b = args.arm_b or (camp_root / "B" / args.campaign)

    inv = {
        "campaign_root": str(camp_root),
        "campaign": args.campaign,
        "A": inventory_arm(arm_a),
        "B0": inventory_arm(arm_b0),
        "B": inventory_arm(arm_b),
    }
    inv_path = out_dir / "3dsig_live_inventory.json"
    inv_path.write_text(json.dumps(inv, indent=2) + "\n")
    print("wrote", inv_path)

    if args.inventory_only:
        print(json.dumps(inv, indent=2))
        return 0

    boot = _load_bootstrap_mod()

    if args.demo:
        all_stats = demo_stats(args.seed)
        source_note = "SYNTHETIC --demo (layout only; not a live re-dock claim)"
    else:
        all_stats = {
            "A": arm_stats(boot, arm_a, args.bootstraps, args.seed, args.thresh),
            "B0": arm_stats(boot, arm_b0, args.bootstraps, args.seed, args.thresh),
            "B": arm_stats(boot, arm_b, args.bootstraps, args.seed, args.thresh),
        }
        source_note = f"LIVE {args.campaign} · S_top10 RMSD<{args.thresh:g}Å · {args.bootstraps} bootstrap"

    fa_label, fa_stats = choose_flexaid_arm(
        all_stats["A"], all_stats["B0"], args.flexaid_arm
    )
    ds_stats = all_stats["B"]

    methods: List[str] = [
        f"FlexAID ({fa_label})",
        "FlexAIDdS (B)",
    ]
    medians: List[Optional[float]] = [
        fa_stats.get("median"),
        ds_stats.get("median"),
    ]
    p05s: List[Optional[float]] = [fa_stats.get("p05"), ds_stats.get("p05")]
    p95s: List[Optional[float]] = [fa_stats.get("p95"), ds_stats.get("p95")]
    colors: List[str] = [COLOR_FLEXAID, COLOR_FLEXAIDDS]

    for name, med, col in load_competitors(args.competitors):
        methods.append(name)
        medians.append(med)
        p05s.append(None)
        p95s.append(None)
        colors.append(col)

    n_fa = fa_stats.get("n_cases") or 0
    n_ds = ds_stats.get("n_cases") or 0
    n_success_note = (
        f"n={n_fa}/{n_ds} cases (FlexAID/FlexAIDdS) · "
        f"success={fa_stats.get('n_success')}/{ds_stats.get('n_success')}"
    )
    title = f"3Dsig red pair — LIVE {args.campaign}"
    if args.demo:
        title += " [DEMO]"
    elif not (fa_stats.get("complete") and ds_stats.get("complete")):
        title += " [PARTIAL]"

    subtitle = n_success_note
    footer = (
        f"{source_note} · campaign_root={camp_root} · "
        f"not archived medians · regenerate after A→B0→B complete"
    )

    svg = svg_bars(
        title=title,
        methods=methods,
        vals=medians,
        colors=colors,
        footer=footer,
        subtitle=subtitle,
        p05=p05s,
        p95=p95s,
    )
    svg_path = out_dir / f"3dsig_live_{args.campaign}_red_pair.svg"
    svg_path.write_text(svg + "\n")
    print("wrote", svg_path)

    # Per-arm bootstrap detail
    meta = {
        "note": (
            "Live 3Dsig red-pair S_top10 bootstrap medians from new docks"
            if not args.demo
            else "DEMO synthetic layout only"
        ),
        "protocol": "docs/implementation/3dsig_red_pair_protocol.md",
        "metric": "S_top10",
        "thresh_A": args.thresh,
        "n_boot": args.bootstraps,
        "seed": args.seed,
        "campaign": args.campaign,
        "campaign_root": str(camp_root),
        "flexaid_arm_used": fa_label,
        "demo": bool(args.demo),
        "arms": all_stats,
        "inventory": inv if not args.demo else None,
        "methods": methods,
        "medians": medians,
        "svg": str(svg_path),
        "blocked_on": _blocked_on(inv if not args.demo else None, all_stats, args.demo),
        "regen_command": _regen_command(camp_root, args.campaign, out_dir),
    }
    meta_path = out_dir / f"3dsig_live_{args.campaign}_bootstrap.json"
    meta_path.write_text(json.dumps(meta, indent=2, default=str) + "\n")
    print("wrote", meta_path)

    # Human status
    print("--- status ---")
    for arm in ("A", "B0", "B"):
        st = all_stats[arm]
        print(
            f"  {arm}: complete={st.get('complete')} n={st.get('n_cases')} "
            f"median={st.get('median')} dir={st.get('arm_dir', arm_a if arm=='A' else arm_b0 if arm=='B0' else arm_b)}"
        )
    print("blocked_on:", meta["blocked_on"])
    print("regen:", meta["regen_command"])

    # Exit 0 even if partial — pipeline must ship before docks finish
    return 0


def _blocked_on(
    inv: Optional[Dict[str, Any]],
    stats: Dict[str, Dict[str, Any]],
    demo: bool,
) -> List[str]:
    if demo:
        return []
    blocked: List[str] = []
    if not stats["A"].get("complete") and not stats["B0"].get("complete"):
        blocked.append("FlexAID arm (A or B0): need complete S_top10 cases (pose PDBs with REMARK RMSD or mode_rmsd_* columns)")
    if not stats["B"].get("complete"):
        blocked.append("FlexAIDdS arm B: not started or no complete S_top10 cases yet")
    if inv:
        for arm in ("A", "B0", "B"):
            empty = inv[arm].get("result_csv_empty") or []
            poses = inv[arm].get("pose_pdb_cases") or []
            if empty and not poses:
                blocked.append(
                    f"{arm}: result.csv present for {empty} but no RMSD / no ranked pose PDBs "
                    f"(engine wrote .cad/.par.res only — packaging gap)"
                )
    return blocked


def _regen_command(camp_root: Path, campaign: str, out_dir: Path) -> str:
    return (
        f'python3 scripts/plot_3dsig_live_bars.py '
        f'--campaign-root "{camp_root}" '
        f'--campaign {campaign} '
        f'--out-dir "{out_dir}" '
        f'--bootstraps 10000'
    )


if __name__ == "__main__":
    raise SystemExit(main())
