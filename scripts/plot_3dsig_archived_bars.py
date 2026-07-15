#!/usr/bin/env python3
"""Rebuild Morency 3Dsig 2017 red-pair barplots from archived medians (pure SVG).

No matplotlib/numpy dependency. Numbers from Morency_LP_3Dsig_2017.pdf and
R2 notes_05-05017.dat. Figure fast path only — not a live re-dock claim.

Usage:
  python3 scripts/plot_3dsig_archived_bars.py
  python3 scripts/plot_3dsig_archived_bars.py --out-dir docs/figures/3dsig_2017

Copyright 2026 Le Bonhomme Pharma
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import List
from xml.sax.saxutils import escape

ARCHIVED = {
    "astex_diverse": {
        "title": "Astex Diverse (N=85) — self-dock",
        "source": "Morency_LP_3Dsig_2017.pdf bar labels",
        "methods": ["FlexAID", "FlexAIDdS", "Vina", "FlexX", "rDock"],
        "median": [0.66, 0.69, 0.82, 0.78, 0.88],
        "colors": ["#E53935", "#E53935", "#1E88E5", "#00ACC1", "#8E24AA"],
    },
    "astex_nn_flrp": {
        "title": "Astex Non-Native FLRP (N=1112)",
        "source": "PDF + notes_05-05017.dat",
        "methods": ["FlexAID", "FlexAIDdS", "Vina", "FlexX", "rDock"],
        "median": [0.48, 0.52, 0.44, 0.50, 0.64],
        "colors": ["#E53935", "#E53935", "#1E88E5", "#00ACC1", "#8E24AA"],
    },
    "astex_nn_flfp": {
        "title": "Astex Non-Native FLFP (N=1112)",
        "source": "PDF + notes_05-05017.dat",
        "methods": ["FlexAID", "FlexAIDdS", "Vina", "rDock"],
        "median": [0.50, 0.54, 0.50, 0.74],
        "colors": ["#E53935", "#E53935", "#1E88E5", "#8E24AA"],
    },
    "hap2_flrp": {
        "title": "HAP2 FLRP (N=76)",
        "source": "PDF + notes_05-05017.dat",
        "methods": ["FlexAID", "FlexAIDdS", "Vina", "FlexX", "rDock"],
        "median": [0.18, 0.22, 0.13, 0.11, 0.22],
        "colors": ["#E53935", "#E53935", "#1E88E5", "#00ACC1", "#8E24AA"],
    },
    "hap2_flfp": {
        "title": "HAP2 FLFP (N=76)",
        "source": "PDF + notes_05-05017.dat",
        "methods": ["FlexAID", "FlexAIDdS", "Vina", "rDock"],
        "median": [0.29, 0.36, 0.33, 0.20],
        "colors": ["#E53935", "#E53935", "#1E88E5", "#8E24AA"],
    },
}


def svg_bars(
    title: str,
    methods: List[str],
    vals: List[float],
    colors: List[str],
    source: str,
) -> str:
    W, H = 720, 420
    L, R, T, B = 70, 30, 50, 90
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
        # axes
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
    for i, (m, v, c) in enumerate(zip(methods, vals, colors)):
        x = L + gap + i * (bar_w + gap)
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
        parts.append(
            f'<text x="{x+bar_w/2:.1f}" y="{T+plot_h+16}" text-anchor="middle" '
            f'font-size="11" font-family="Helvetica,Arial,sans-serif" '
            f'transform="rotate(20 {x+bar_w/2:.1f} {T+plot_h+16})">{escape(m)}</text>'
        )
    if len(vals) >= 2:
        delta = vals[1] - vals[0]
        parts.append(
            f'<text x="{W/2}" y="{H-28}" text-anchor="middle" font-size="12" '
            f'fill="#B71C1C" font-family="Helvetica,Arial,sans-serif">'
            f"Δ FlexAIDdS−FlexAID = {delta:+.2f}  ·  red = FlexAID pair</text>"
        )
    parts.append(
        f'<text x="8" y="{H-8}" font-size="9" fill="#555" font-family="Helvetica,Arial,sans-serif">'
        f"Archived · {escape(source)} · S_top10 RMSD&lt;2Å · 10×2e6 · 10k bootstrap · not a live re-dock"
        f"</text>"
    )
    parts.append("</svg>")
    return "\n".join(parts)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=None)
    args = ap.parse_args()
    root = Path(__file__).resolve().parents[1]
    out = args.out_dir or (root / "docs" / "figures" / "3dsig_2017")
    out.mkdir(parents=True, exist_ok=True)

    written = []
    for key, panel in ARCHIVED.items():
        svg = svg_bars(
            panel["title"],
            panel["methods"],
            panel["median"],
            panel["colors"],
            panel["source"],
        )
        path = out / f"3dsig_archived_{key}.svg"
        path.write_text(svg + "\n")
        written.append(str(path))
        print("wrote", path)

    meta = {
        "note": "Archived 3Dsig 2017 medians for figure rebuild only",
        "protocol": "docs/implementation/3dsig_red_pair_protocol.md",
        "deck": "Morency_LP_3Dsig_2017.pdf",
        "panels": ARCHIVED,
        "files": written,
    }
    meta_path = out / "3dsig_archived_medians.json"
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")
    print("wrote", meta_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
