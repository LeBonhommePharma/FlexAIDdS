#!/usr/bin/env python3
"""Build a wall-oracle panel with saturating soft-wall contacts.

Falsemin poses often have cf_wal well below WAL_CONTACT_CAP (50), so
FLEXAIDDS_WAL_COERCIVE (remove cap) is a no-op. This script:

1. Loads production dock configs (ops/gates/configs).
2. For each clean probe, translates the falsemin ligand toward the receptor
   COM by a ladder of offsets (Å) and scores with probe_cf.
3. Keeps the offset that maximizes cf_wal (or first that reaches ≥ min_wal).
4. Writes a panel TSV + poses under out-dir for wall_coercive_oracle re-run.

Usage:
  python3 scripts/wall_saturating_panel.py \\
    --probe-cf build/probe_cf --binary build/FlexAIDdS --data-dir . \\
    --out-dir ~/flexaidds_results/workorders/wall_sat_panel
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, List, Optional, Sequence, Tuple

CLEAN = ("1J3J", "1K3U", "1L7F", "1N1M", "1M2Z")


def repo_root() -> Path:
    p = Path(__file__).resolve().parents[1]
    return p if (p / "AGENTS.md").exists() else Path.cwd()


def com_pdb(path: Path) -> Optional[Tuple[float, float, float]]:
    xs: List[float] = []
    ys: List[float] = []
    zs: List[float] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        try:
            xs.append(float(line[30:38]))
            ys.append(float(line[38:46]))
            zs.append(float(line[46:54]))
        except ValueError:
            continue
    if not xs:
        return None
    n = float(len(xs))
    return (sum(xs) / n, sum(ys) / n, sum(zs) / n)


def translate_pdb(src: Path, dst: Path, dx: float, dy: float, dz: float) -> None:
    lines_out: List[str] = []
    for line in src.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                x = float(line[30:38]) + dx
                y = float(line[38:46]) + dy
                z = float(line[46:54]) + dz
                line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
            except ValueError:
                pass
        lines_out.append(line)
    dst.write_text("\n".join(lines_out) + "\n", encoding="utf-8")


def probe(
    probe_cf: Path,
    *,
    receptor: Path,
    pose: Path,
    ligand: Path,
    binary: Path,
    data_dir: Path,
    config: Path,
    label: str,
) -> Optional[dict]:
    env = os.environ.copy()
    env["FLEXAIDDS_WAL_COERCIVE"] = "0"
    env.pop("FLEXAIDDS_COM_BURIAL_CAP", None)
    cmd = [
        str(probe_cf),
        "--receptor",
        str(receptor),
        "--pose",
        str(pose),
        "--ligand",
        str(ligand),
        "--binary",
        str(binary),
        "--data-dir",
        str(data_dir),
        "--config",
        str(config),
        "--pdb",
        label,
    ]
    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=600, env=env, check=False
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    for line in reversed(proc.stdout.splitlines()):
        line = line.strip()
        if not line.startswith("{"):
            continue
        try:
            obj = json.loads(line)
            if "cf_total" in obj:
                return obj
        except (json.JSONDecodeError, TypeError, ValueError):
            continue
    return None


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    repo = repo_root()
    ap.add_argument("--repo", type=Path, default=repo)
    ap.add_argument("--probe-cf", type=Path, default=None)
    ap.add_argument("--binary", type=Path, default=None)
    ap.add_argument("--data-dir", type=Path, default=None)
    ap.add_argument("--panel", nargs="*", default=list(CLEAN))
    ap.add_argument(
        "--offsets",
        type=float,
        nargs="*",
        default=[0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
        help="Å translation of decoy toward receptor COM",
    )
    ap.add_argument("--min-wal", type=float, default=45.0, help="target cf_wal for saturation")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(argv)
    repo = args.repo.resolve()
    probe_cf = args.probe_cf or (repo / "build" / "probe_cf")
    binary = args.binary or (repo / "build" / "FlexAIDdS")
    data_dir = args.data_dir or repo
    out = args.out_dir.expanduser()
    poses = out / "poses"
    poses.mkdir(parents=True, exist_ok=True)

    rows: List[dict] = []
    for pdb in args.panel:
        rec = repo / f"benchmarks/astex_diverse/astex_diverse/{pdb}/{pdb}_apo.pdb"
        lig = repo / f"benchmarks/astex_diverse/astex_diverse/{pdb}/{pdb}_ligand.sdf"
        dec0 = repo / f"diagnostic/refs/{pdb}/falsemin_armA.pdb"
        cfg = repo / f"ops/gates/configs/{pdb}_dock_config.json"
        if not all(p.is_file() for p in (rec, lig, dec0, cfg, probe_cf, binary)):
            rows.append({"pdb": pdb, "status": "skip_missing"})
            continue
        rcom = com_pdb(rec)
        lcom = com_pdb(dec0)
        if not rcom or not lcom:
            rows.append({"pdb": pdb, "status": "skip_com"})
            continue
        # Unit vector decoy COM → receptor COM (burial)
        vx, vy, vz = rcom[0] - lcom[0], rcom[1] - lcom[1], rcom[2] - lcom[2]
        norm = (vx * vx + vy * vy + vz * vz) ** 0.5
        if norm < 1e-6:
            rows.append({"pdb": pdb, "status": "skip_zero_vector"})
            continue
        ux, uy, uz = vx / norm, vy / norm, vz / norm

        best: Optional[dict] = None
        ladder: List[dict] = []
        for off in args.offsets:
            pose_p = poses / f"{pdb}_burial_{off:.1f}.pdb"
            translate_pdb(dec0, pose_p, ux * off, uy * off, uz * off)
            obj = probe(
                probe_cf,
                receptor=rec,
                pose=pose_p,
                ligand=lig,
                binary=binary,
                data_dir=data_dir,
                config=cfg,
                label=f"{pdb}_b{off:.1f}",
            )
            if not obj:
                ladder.append({"offset": off, "status": "probe_fail"})
                continue
            entry = {
                "offset": off,
                "cf_total": obj.get("cf_total"),
                "cf_wal": obj.get("cf_wal"),
                "pose": str(pose_p),
            }
            ladder.append(entry)
            wal = float(obj.get("cf_wal") or 0.0)
            if best is None or wal > float(best.get("cf_wal") or -1e9):
                best = entry
            if wal >= args.min_wal:
                best = entry
                break
        rows.append(
            {
                "pdb": pdb,
                "status": "ok",
                "best_offset": best.get("offset") if best else None,
                "best_cf_wal": best.get("cf_wal") if best else None,
                "best_cf_total": best.get("cf_total") if best else None,
                "best_pose": best.get("pose") if best else None,
                "ladder": ladder,
            }
        )
        print(
            f"[{pdb}] best_wal={best.get('cf_wal') if best else None} "
            f"offset={best.get('offset') if best else None}",
            flush=True,
        )

    # Write saturating manifest for wall_coercive_oracle (decoy = best pose)
    man = out / "saturating_panel_manifest.tsv"
    lines = [
        "# pdb\treceptor\tnative_sdf\tdecoy_pdb\tdock_config",
    ]
    n_sat = 0
    for r in rows:
        if r.get("status") != "ok" or not r.get("best_pose"):
            continue
        pdb = r["pdb"]
        wal = float(r.get("best_cf_wal") or 0)
        if wal >= args.min_wal:
            n_sat += 1
        rec = f"benchmarks/astex_diverse/astex_diverse/{pdb}/{pdb}_apo.pdb"
        nat = f"benchmarks/astex_diverse/astex_diverse/{pdb}/{pdb}_ligand.sdf"
        cfg = f"ops/gates/configs/{pdb}_dock_config.json"
        # store absolute decoy path
        lines.append(f"{pdb}\t{rec}\t{nat}\t{r['best_pose']}\t{cfg}")
    man.write_text("\n".join(lines) + "\n", encoding="utf-8")

    summary = {
        "n_ok": sum(1 for r in rows if r.get("status") == "ok"),
        "n_saturating": n_sat,
        "min_wal": args.min_wal,
        "manifest": str(man),
        "rows": [
            {
                k: v
                for k, v in r.items()
                if k != "ladder"
            }
            for r in rows
        ],
    }
    (out / "saturating_panel_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    # CSV flat
    with (out / "saturating_panel.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=[
                "pdb",
                "status",
                "best_offset",
                "best_cf_wal",
                "best_cf_total",
                "best_pose",
            ],
        )
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in w.fieldnames})

    print(json.dumps({k: summary[k] for k in ("n_ok", "n_saturating", "min_wal", "manifest")}, indent=2))
    return 0 if summary["n_ok"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
