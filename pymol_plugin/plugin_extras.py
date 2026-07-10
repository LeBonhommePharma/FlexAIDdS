"""Plugin polish: status, publication view, CF contacts, RMSD overlay, T-scan.

Display-only post-hoc helpers. Never change ranking or scoring.
"""

from __future__ import annotations

import math
import os
from pathlib import Path
from typing import Dict, List, Optional

try:
    from pymol import cmd
except ImportError as exc:
    raise ImportError("PyMOL not available") from exc

from .session import SESSION, object_name


def flexaids_status(*_args, **_kwargs) -> None:
    """Dump current plugin session state."""
    print("FlexAID∆S session status")
    print(f"  prefix:      {SESSION.prefix}")
    print(f"  temperature: {SESSION.temperature_K} K")
    print(f"  output_dir:  {SESSION.output_dir}")
    print(f"  n_modes:     {SESSION.n_modes}")
    if SESSION.result is None:
        print("  result:      (none loaded)")
        return
    for mode in SESSION.result.binding_modes:
        src = (mode.metadata or {}).get("ledger_source", "?")
        f = mode.free_energy
        f_str = f"{f:.2f}" if f is not None else "N/A"
        cf = mode.best_cf
        cf_str = f"{cf:.2f}" if cf is not None else "N/A"
        objs = len(SESSION.object_names_for(mode.mode_id))
        print(
            f"  mode {mode.mode_id}: F={f_str} CF={cf_str} "
            f"n_poses={mode.n_poses} objs={objs} ledger={src}"
        )
    recent = os.environ.get("FLEXAIDDS_RESULTS", "")
    if recent:
        print(f"  FLEXAIDDS_RESULTS: {recent}")


def publication_view(
    mode_id: int = 1,
    ray: int = 0,
    output_png: str = "",
) -> None:
    """One-click publication view for the best pose of a mode."""
    mid = int(mode_id)
    mode = SESSION.get_mode(mid)
    if mode is None:
        print("ERROR: No loaded mode. Use flexaids_load first.")
        return
    names = SESSION.object_names_for(mid)
    if not names:
        print(f"ERROR: No PyMOL objects for mode {mid}.")
        return
    best = mode.best_pose()
    if best is not None:
        best_name = object_name(SESSION.prefix, mid, best.pose_rank)
    else:
        best_name = names[0]

    for n in names:
        cmd.disable(n)
    cmd.enable(best_name)
    cmd.hide("everything", best_name)
    cmd.show("cartoon", f"{best_name} and polymer")
    cmd.show("sticks", f"{best_name} and organic")
    cmd.color("grey80", f"{best_name} and polymer")
    cmd.util.cnc(f"{best_name} and organic")
    cmd.bg_color("white")
    cmd.set("ray_shadows", 0)
    cmd.set("ambient", 0.4)
    cmd.set("specular", 0.0)
    cmd.orient(best_name)
    cmd.zoom(best_name, buffer=2.0)
    print(f"Publication view: {best_name}")

    out = str(output_png).strip()
    if out:
        if int(ray):
            cmd.ray(1200, 900)
        cmd.png(out, dpi=300)
        print(f"  wrote {out}")


def show_cf_contacts(mode_id: int = 1) -> None:
    """Print CF-component summary from REMARKs and label optres on best pose.

    Does not invent PLIP interactions — only surfaces existing CF.* REMARK data.
    """
    mid = int(mode_id)
    mode = SESSION.get_mode(mid)
    if mode is None:
        print("ERROR: No loaded mode. Use flexaids_load first.")
        return
    best = mode.best_pose()
    if best is None:
        print("ERROR: No poses in mode.")
        return
    remarks = best.remarks or {}
    components = {
        k: v for k, v in remarks.items()
        if k.startswith("cf_") or k in ("cf", "cf_com", "cf_hbond", "cf_wal", "cf_sas", "cf_con", "cf_gist", "cf_metal", "cf_elec")
    }
    # Normalise dotted keys that became cf_com etc.
    print(f"CF component map for mode {mid} pose {best.pose_rank} (CF proxy terms):")
    for key in sorted(remarks):
        if key.startswith("cf") and key != "cf_app":
            print(f"  {key:12s} = {remarks[key]}")
    if not any(k.startswith("cf") for k in remarks):
        print("  (no CF.* REMARK lines found)")

    # Label optimizable residues if present
    optres_lines = [
        k for k in remarks
        if "optimizable" in k or "residue" in k
    ]
    names = SESSION.object_names_for(mid)
    if names:
        best_name = object_name(SESSION.prefix, mid, best.pose_rank)
        try:
            cmd.distance(
                f"{SESSION.prefix}_contacts_m{mid}",
                f"({best_name} and organic)",
                f"({best_name} and polymer and name CA)",
                mode=2,
                cutoff=4.0,
            )
            print(f"  dashed contacts object: {SESSION.prefix}_contacts_m{mid}")
        except Exception as exc:
            print(f"  (could not create contact dashes: {exc})")


def show_rmsd_overlay(mode_id: int = 1) -> None:
    """Print dual RMSD (raw / symmetry-corrected) and success badge."""
    mid = int(mode_id)
    mode = SESSION.get_mode(mid)
    if mode is None:
        print("ERROR: No loaded mode. Use flexaids_load first.")
        return
    best = mode.best_pose()
    if best is None:
        print("ERROR: No poses.")
        return
    raw = best.rmsd_raw
    sym = best.rmsd_sym
    print(f"RMSD overlay — mode {mid} pose {best.pose_rank}:")
    print(f"  rmsd_raw = {raw if raw is not None else 'N/A'} Å")
    print(f"  rmsd_sym = {sym if sym is not None else 'N/A'} Å")
    metric = sym if sym is not None else raw
    if metric is not None and metric <= 2.0:
        print("  badge: RMSD ≤ 2.0 Å (pose geometry success)")
        print("  note: full benchmark success also requires PoseBusters when available")
    elif metric is not None:
        print(f"  badge: RMSD = {metric:.2f} Å (> 2.0 — not a geometric success)")
    else:
        print("  badge: no reference RMSD in REMARKs")


def temperature_scan(
    mode_id: int = 1,
    temperatures: str = "298,310",
) -> None:
    """Recompute pure-Python StatMech F/H/S at multiple T from pose CF ensemble.

    Display-only; does not change stored ranking. Labels results as
    ensemble estimate from CF proxy.
    """
    from flexaidds.thermodynamics import StatMechEngine

    mid = int(mode_id)
    mode = SESSION.get_mode(mid)
    if mode is None:
        print("ERROR: No loaded mode. Use flexaids_load first.")
        return
    energies = []
    for p in mode.poses:
        if p.cf is not None:
            energies.append(float(p.cf))
        elif p.cf_app is not None:
            energies.append(float(p.cf_app))
    if not energies:
        print("ERROR: No CF values for temperature scan.")
        return

    temps = []
    for tok in str(temperatures).split(","):
        tok = tok.strip()
        if tok:
            temps.append(float(tok))
    if not temps:
        temps = [298.0, 310.0]

    print(
        f"T-scan mode {mid} (ensemble estimate from CF proxy; "
        f"n_samples={len(energies)}):"
    )
    print(f"  {'T(K)':>8}  {'F':>10}  {'H':>10}  {'S':>12}  {'-TS':>10}")
    rows = []
    for T in temps:
        eng = StatMechEngine(T)
        for e in energies:
            eng.add_sample(e)
        td = eng.compute()
        rows.append((T, td.free_energy, td.mean_energy, td.entropy, -T * td.entropy))
        print(
            f"  {T:8.1f}  {td.free_energy:10.3f}  {td.mean_energy:10.3f}  "
            f"{td.entropy:12.6f}  {-T * td.entropy:10.3f}"
        )
    if len(rows) >= 2:
        # Ranking shift: only one mode here; report ΔF across T
        df = rows[-1][1] - rows[0][1]
        print(f"  ΔF ({rows[0][0]:.0f}→{rows[-1][0]:.0f} K) = {df:+.3f} kcal/mol")


def show_cleft(sphere_pdb: str) -> None:
    """Load a FlexAID cleft sphere PDB into PyMOL."""
    path = Path(str(sphere_pdb).strip())
    if not path.is_file():
        print(f"ERROR: Sphere file not found: {sphere_pdb}")
        return
    obj = f"{SESSION.prefix}_cleft_spheres"
    try:
        cmd.delete(obj)
    except Exception:
        pass
    cmd.load(str(path), obj)
    cmd.show("spheres", obj)
    cmd.color("yellow", obj)
    cmd.set("sphere_transparency", 0.5, obj)
    print(f"Loaded cleft spheres as '{obj}'")


def colorbar_legend(metric: str = "cf") -> None:
    """Print the burgundy↔purple color legend for the active metric."""
    metric = str(metric).strip().lower()
    unit = "kcal/mol (CF proxy)" if metric == "cf" else "kcal/mol (ensemble F)"
    print(f"Color legend — metric={metric} [{unit}]")
    print("  burgundy  = lower / more favourable")
    print("  purple    = higher / less favourable")
    print("  (spectrum implemented via set_color on pose objects)")


def recent_results(*_args, **_kwargs) -> None:
    """List last few result dirs under FLEXAIDDS_RESULTS (if set)."""
    base = os.environ.get("FLEXAIDDS_RESULTS")
    if not base:
        print("FLEXAIDDS_RESULTS is not set.")
        return
    root = Path(base)
    if not root.is_dir():
        print(f"FLEXAIDDS_RESULTS does not exist: {base}")
        return
    # Prefer recent directories by mtime
    dirs = [p for p in root.iterdir() if p.is_dir()]
    dirs.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    print(f"Recent results under {base}:")
    for d in dirs[:5]:
        print(f"  {d}")
