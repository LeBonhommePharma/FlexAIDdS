#!/usr/bin/env python3
"""Toy μVT validation: MOR (fentanyl+naloxone) and 5-HT2A (5-MeO-DMT+5-HT).

Uses synthetic log_Z values only — not experimental affinity claims.
Prints ⟨N⟩, apparent K_i trends, and NRGsuite-style set_concentration usage.

Usage:
  python3 scripts/muvt_toy_validation.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from flexaidds.grand_canonical import (  # noqa: E402
    CompetitiveSite,
    kB_kcal,
    plot_occupancy_curve,
    set_concentration,
)


def nvt_F(log_Z: float, T: float) -> float:
    return -kB_kcal * T * log_Z


def run_mor() -> None:
    T = 310.0
    log_Z_fen, log_Z_nal = 15.0, 12.0
    print("=== MOR toy: fentanyl + naloxone (synthetic log_Z) ===")
    print(f"  NVT F_fentanyl = {nvt_F(log_Z_fen, T):.3f} kcal/mol")
    print(f"  NVT F_naloxone = {nvt_F(log_Z_nal, T):.3f} kcal/mol")

    site = CompetitiveSite(T)
    site.add("fentanyl", log_Z_fen, 1e-9)
    site.add("naloxone", log_Z_nal, 1e-6)
    print(f"  p_f={site.binding_probability('fentanyl'):.4f}  "
          f"p_n={site.binding_probability('naloxone'):.4f}  "
          f"⟨N⟩={site.mean_N():.4f}")
    print(f"  collapse={site.ligand_entropy_collapse():.3f}  "
          f"S_mix={site.mixing_entropy():.6f} kcal/mol/K")

    # NRGsuite command: set_concentration([L1, L2, ...])
    set_concentration(site, [1e-9, 1e-4], names=["fentanyl", "naloxone"])
    print("  after set_concentration([1e-9, 1e-4]):")
    print(f"  p_f={site.binding_probability('fentanyl'):.4f}  "
          f"p_n={site.binding_probability('naloxone'):.4f}  "
          f"⟨N⟩={site.mean_N():.4f}")

    curve = site.occupancy_vs_concentration(
        "naloxone", [10 ** x for x in range(-9, -2)]
    )
    print("  naloxone titration (c_M, p_species, mean_N):")
    for pt in curve:
        print(f"    {pt.concentration_M:.1e}  {pt.p_species:.4f}  {pt.mean_N:.4f}")
    ax = plot_occupancy_curve(curve, title="MOR: naloxone titration (toy)")
    if ax is not None:
        out = ROOT / "results" / "muvt_mor_occupancy_toy.png"
        out.parent.mkdir(parents=True, exist_ok=True)
        ax.figure.savefig(out, dpi=120, bbox_inches="tight")
        print(f"  wrote {out}")


def run_5ht2a() -> None:
    T = 310.0
    log_Z_dmt, log_Z_ht = 14.0, 11.0
    print("\n=== 5-HT2A toy: 5-MeO-DMT + 5-HT endogenous (synthetic log_Z) ===")
    site = CompetitiveSite(T)
    site.add("5-MeO-DMT", log_Z_dmt, 1e-8)
    site.add("5-HT", log_Z_ht, 1e-6)
    print(f"  p_dmt={site.binding_probability('5-MeO-DMT'):.4f}  "
          f"p_ht={site.binding_probability('5-HT'):.4f}  "
          f"⟨N⟩={site.mean_N():.4f}")
    set_concentration(site, {"5-HT": 1e-3})
    print(f"  after set_concentration(5-HT=1 mM): "
          f"p_ht={site.binding_probability('5-HT'):.4f}  "
          f"collapse={site.ligand_entropy_collapse():.3f}")


def main() -> int:
    run_mor()
    run_5ht2a()
    print("\nNRGsuite integration:")
    print("  set_concentration([L1, L2, ...])  # Molar units")
    print("  See docs/theory.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
