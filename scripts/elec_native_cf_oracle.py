#!/usr/bin/env python3
"""W1.2 score-only: FLEXAIDDS_USE_ELEC=0 vs 1 native-vs-decoy CF via probe_cf.

Must not mass-invert clean probes when elec is ON (default remains OFF).
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Optional

PANEL = ("1G9V", "1M2Z", "1N1M", "1J3J", "1K3U", "1L7F")


def resolve_receptor(repo: Path, pdb: str) -> Optional[Path]:
    for c in (
        Path.home() / f".flexaidds/benchmarks/astex_diverse/{pdb}/{pdb}_apo.pdb",
        repo / f"benchmarks/astex_diverse/astex_diverse/{pdb}/{pdb}_apo.pdb",
        repo / f"benchmarks/astex_diverse/{pdb}/{pdb}_apo.pdb",
    ):
        if c.is_file():
            return c
    return None


def probe(
    probe_cf: Path,
    *,
    receptor: Path,
    pose: Path,
    ligand: Path,
    binary: Path,
    data_dir: Path,
    config: Path,
    use_elec: bool,
    label: str,
) -> Optional[float]:
    env = os.environ.copy()
    env["FLEXAIDDS_USE_ELEC"] = "1" if use_elec else "0"
    env.pop("FLEXAIDDS_WAL_COERCIVE", None)
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
        "--pdb",
        label,
    ]
    if config.is_file():
        cmd.extend(["--config", str(config)])
    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=600, env=env, check=False
        )
    except (OSError, subprocess.TimeoutExpired) as e:
        print(f"error {label}: {e}", file=sys.stderr)
        return None
    for line in reversed(proc.stdout.splitlines()):
        line = line.strip()
        if line.startswith("{"):
            try:
                o = json.loads(line)
                if "cf_total" in o:
                    return float(o["cf_total"])
            except (json.JSONDecodeError, TypeError, ValueError):
                pass
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--probe-cf", type=Path, required=True)
    ap.add_argument("--binary", type=Path, required=True)
    ap.add_argument("--data-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--panel", nargs="*", default=list(PANEL))
    args = ap.parse_args()
    repo = args.repo.resolve()
    cfg = repo / "diagnostic" / "probe_config.json"
    out = args.out_dir.expanduser()
    out.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    for pdb in args.panel:
        d = repo / "diagnostic" / "refs" / pdb
        nat, dec = d / "native.sdf", d / "falsemin_armA.pdb"
        rec = resolve_receptor(repo, pdb)
        if not (rec and nat.is_file() and dec.is_file()):
            rows.append({"pdb": pdb, "status": "skip"})
            continue
        print(f"[{pdb}] elec OFF/ON …")
        n0 = probe(
            args.probe_cf,
            receptor=rec,
            pose=nat,
            ligand=nat,
            binary=args.binary,
            data_dir=args.data_dir,
            config=cfg,
            use_elec=False,
            label=f"{pdb}_n0",
        )
        d0 = probe(
            args.probe_cf,
            receptor=rec,
            pose=dec,
            ligand=nat,
            binary=args.binary,
            data_dir=args.data_dir,
            config=cfg,
            use_elec=False,
            label=f"{pdb}_d0",
        )
        n1 = probe(
            args.probe_cf,
            receptor=rec,
            pose=nat,
            ligand=nat,
            binary=args.binary,
            data_dir=args.data_dir,
            config=cfg,
            use_elec=True,
            label=f"{pdb}_n1",
        )
        d1 = probe(
            args.probe_cf,
            receptor=rec,
            pose=dec,
            ligand=nat,
            binary=args.binary,
            data_dir=args.data_dir,
            config=cfg,
            use_elec=True,
            label=f"{pdb}_d1",
        )

        def win(n, d):
            return n is not None and d is not None and (n - d) <= 0

        w0, w1 = win(n0, d0), win(n1, d1)
        rows.append(
            {
                "pdb": pdb,
                "status": "ok",
                "cf_n_off": n0,
                "cf_d_off": d0,
                "dCF_off": None if n0 is None or d0 is None else n0 - d0,
                "native_wins_off": w0,
                "cf_n_on": n1,
                "cf_d_on": d1,
                "dCF_on": None if n1 is None or d1 is None else n1 - d1,
                "native_wins_on": w1,
                "inverted_by_elec": w0 and not w1,
            }
        )
    ok = [r for r in rows if r.get("status") == "ok"]
    n_inv = sum(1 for r in ok if r.get("inverted_by_elec"))
    summary = {
        "n_scored": len(ok),
        "n_inverted_by_elec": n_inv,
        "mass_invert": n_inv > max(1, len(ok) // 2),
        "pass_no_mass_invert": n_inv <= max(1, len(ok) // 2),
        "rows": rows,
    }
    (out / "elec_oracle.json").write_text(json.dumps(summary, indent=2) + "\n")
    with (out / "elec_oracle.csv").open("w", newline="") as fh:
        keys = sorted({k for r in rows for k in r})
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    md = [
        "# Elec OFF vs ON native-CF oracle",
        "",
        f"Scored: {len(ok)}; inverted_by_elec: **{n_inv}**; "
        f"mass_invert: **{summary['mass_invert']}**; "
        f"PASS (no mass invert): **{summary['pass_no_mass_invert']}**",
        "",
        "| pdb | dCF_off | win_off | dCF_on | win_on | inverted |",
        "|-----|---------|---------|--------|--------|----------|",
    ]
    for r in rows:
        if r.get("status") != "ok":
            md.append(f"| {r['pdb']} | skip | | | | |")
            continue
        md.append(
            f"| {r['pdb']} | {r.get('dCF_off')} | {r.get('native_wins_off')} | "
            f"{r.get('dCF_on')} | {r.get('native_wins_on')} | {r.get('inverted_by_elec')} |"
        )
    (out / "elec_oracle.md").write_text("\n".join(md) + "\n")
    print(json.dumps({k: summary[k] for k in summary if k != "rows"}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
