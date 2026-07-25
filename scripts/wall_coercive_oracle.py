#!/usr/bin/env python3
"""W2 score-only wall oracle: FLEXAIDDS_WAL_COERCIVE OFF vs ON via probe_cf.

Uses **production** per-target dock configs under ops/gates/configs/{PDB}_dock_config.json
(and optionally ops/gates/panel_manifest.tsv) so LOCCLF / optres pocket pruning
matches claim CF. Do **not** default to diagnostic/probe_config.json (1G9V ligand
path + wrong sas_weight/vct_r0 → ~200× inflated CF).

Panel default: methodology clean probes 1J3J 1K3U 1L7F 1N1M 1M2Z (1G9V excluded).

Usage:
  python3 scripts/wall_coercive_oracle.py \\
    --probe-cf build/probe_cf --binary build/FlexAIDdS --data-dir . \\
    --manifest ops/gates/panel_manifest.tsv \\
    --out-dir ~/flexaidds_results/workorders/wall_oracle_prod

Does not rebuild binaries. No GA.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

# OPS campaign methodology STEP 2 clean probes (1G9V chain-control only).
CLEAN_PANEL = ("1J3J", "1K3U", "1L7F", "1N1M", "1M2Z")


def find_repo_root() -> Path:
    p = Path(__file__).resolve().parents[1]
    if (p / "AGENTS.md").exists():
        return p
    return Path.cwd()


def resolve_dock_config(repo: Path, pdb: str) -> Optional[Path]:
    """Production dock_config for LOCCLF-equivalent CF (never diagnostic/probe_config)."""
    candidates = [
        repo / "ops" / "gates" / "configs" / f"{pdb}_dock_config.json",
        repo / "ops" / "gates" / "configs" / f"{pdb}.json",
    ]
    for c in candidates:
        if c.is_file():
            return c
    return None


def resolve_receptor(repo: Path, pdb: str) -> Optional[Path]:
    candidates = [
        repo / f"benchmarks/astex_diverse/astex_diverse/{pdb}/{pdb}_apo.pdb",
        repo / f"benchmarks/astex_diverse/{pdb}/{pdb}_apo.pdb",
        Path.home() / f".flexaidds/benchmarks/astex_diverse/{pdb}/{pdb}_apo.pdb",
    ]
    for c in candidates:
        if c.is_file():
            return c
    return None


def resolve_native_ligand(repo: Path, pdb: str) -> Optional[Path]:
    """Crystal ligand SDF (topology for PDB decoys + native pose when no diagnostic)."""
    candidates = [
        repo / f"benchmarks/astex_diverse/astex_diverse/{pdb}/{pdb}_ligand.sdf",
        repo / f"benchmarks/astex_diverse/{pdb}/{pdb}_ligand.sdf",
        repo / "diagnostic" / "refs" / pdb / "native.sdf",
    ]
    for c in candidates:
        if c.is_file():
            return c
    return None


def resolve_decoy(repo: Path, pdb: str) -> Optional[Path]:
    d = repo / "diagnostic" / "refs" / pdb / "falsemin_armA.pdb"
    return d if d.is_file() else None


def load_manifest(repo: Path, path: Path) -> List[Dict[str, Path]]:
    """Parse panel_manifest.tsv → list of {pdb, receptor, native, decoy, config}."""
    rows: List[Dict[str, Path]] = []
    text = path.read_text(encoding="utf-8")
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        pdb, rec, nat, dec, cfg = parts[0], parts[1], parts[2], parts[3], parts[4]
        if "PLACEHOLDER" in dec:
            continue

        def _abs(p: str) -> Path:
            pp = Path(p)
            return pp if pp.is_absolute() else (repo / pp)

        rows.append(
            {
                "pdb": pdb,
                "receptor": _abs(rec),
                "native": _abs(nat),
                "decoy": _abs(dec),
                "config": _abs(cfg),
            }
        )
    return rows


def probe_cf_total(
    probe: Path,
    *,
    receptor: Path,
    pose: Path,
    ligand: Path,
    binary: Path,
    data_dir: Path,
    config: Path,
    wal_coercive: bool,
    label: str,
) -> Optional[float]:
    """Score one pose; requires production config (fail-closed)."""
    if not config.is_file():
        print(f"  missing config for {label}: {config}", file=sys.stderr)
        return None
    env = os.environ.copy()
    env["FLEXAIDDS_WAL_COERCIVE"] = "1" if wal_coercive else "0"
    env.pop("FLEXAIDDS_COM_BURIAL_CAP", None)
    cmd = [
        str(probe),
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
            cmd,
            capture_output=True,
            text=True,
            timeout=600,
            env=env,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as e:
        print(f"  probe_cf error {label}: {e}", file=sys.stderr)
        return None
    for line in reversed(proc.stdout.splitlines()):
        line = line.strip()
        if not line.startswith("{"):
            continue
        try:
            obj = json.loads(line)
            if "cf_total" in obj:
                return float(obj["cf_total"])
        except (json.JSONDecodeError, TypeError, ValueError):
            continue
    if proc.returncode != 0:
        print(
            f"  probe_cf fail {label} rc={proc.returncode}: {proc.stderr[-300:]}",
            file=sys.stderr,
        )
    return None


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    repo = find_repo_root()
    ap.add_argument("--repo", type=Path, default=repo)
    ap.add_argument("--probe-cf", type=Path, default=None)
    ap.add_argument("--binary", type=Path, default=None)
    ap.add_argument("--data-dir", type=Path, default=None)
    ap.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Override single config for all targets (discouraged; prefer per-target)",
    )
    ap.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help="panel_manifest.tsv (default: ops/gates/panel_manifest.tsv if present)",
    )
    ap.add_argument(
        "--panel",
        nargs="*",
        default=None,
        help=f"PDB codes (default: clean {CLEAN_PANEL})",
    )
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(argv)

    repo = args.repo.expanduser().resolve()
    probe = args.probe_cf
    if probe is None:
        for c in (repo / "build" / "probe_cf",):
            if c.is_file():
                probe = c
                break
    if probe is None or not Path(probe).is_file():
        print("error: probe_cf not found", file=sys.stderr)
        return 2

    binary = args.binary
    if binary is None:
        for c in (
            repo / "build" / "FlexAIDdS",
            Path.home() / "flexaidds_results/baseline_engine/FlexAIDdS",
        ):
            if c.is_file():
                binary = c
                break
    if binary is None or not Path(binary).is_file():
        print("error: FlexAIDdS binary not found", file=sys.stderr)
        return 2

    data_dir = args.data_dir or repo
    out_dir = args.out_dir.expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Prefer manifest (production paths)
    manifest_path = args.manifest
    if manifest_path is None:
        cand = repo / "ops" / "gates" / "panel_manifest.tsv"
        if cand.is_file():
            manifest_path = cand

    work: List[Dict[str, Any]] = []
    if manifest_path is not None and Path(manifest_path).is_file():
        want = set(args.panel) if args.panel else set(CLEAN_PANEL)
        for row in load_manifest(repo, Path(manifest_path)):
            if row["pdb"] not in want:
                continue
            work.append(row)
    else:
        for pdb in args.panel or list(CLEAN_PANEL):
            cfg = args.config or resolve_dock_config(repo, pdb)
            if cfg is None:
                work.append({"pdb": pdb, "status": "skip_missing_config"})
                continue
            work.append(
                {
                    "pdb": pdb,
                    "receptor": resolve_receptor(repo, pdb),
                    "native": resolve_native_ligand(repo, pdb),
                    "decoy": resolve_decoy(repo, pdb),
                    "config": cfg,
                }
            )

    rows: list[dict[str, Any]] = []
    for item in work:
        pdb = item["pdb"]
        if item.get("status") == "skip_missing_config":
            rows.append({"pdb": pdb, "status": "skip_missing_config"})
            continue
        rec = item.get("receptor")
        nat = item.get("native")
        dec = item.get("decoy")
        cfg = item.get("config")
        if isinstance(rec, str):
            rec = Path(rec)
        if isinstance(nat, str):
            nat = Path(nat)
        if isinstance(dec, str):
            dec = Path(dec)
        if isinstance(cfg, str):
            cfg = Path(cfg)

        if not rec or not Path(rec).is_file():
            rows.append({"pdb": pdb, "status": "skip_missing_receptor"})
            continue
        if not nat or not Path(nat).is_file():
            rows.append({"pdb": pdb, "status": "skip_missing_native"})
            continue
        if not dec or not Path(dec).is_file():
            rows.append({"pdb": pdb, "status": "skip_missing_decoy"})
            continue
        if not cfg or not Path(cfg).is_file():
            rows.append(
                {
                    "pdb": pdb,
                    "status": "skip_missing_config",
                    "config": str(cfg) if cfg else "",
                }
            )
            continue

        # Refuse diagnostic/probe_config as production (wrong LOCCLF)
        cfg_s = str(cfg)
        if "diagnostic/probe_config" in cfg_s.replace("\\", "/"):
            rows.append(
                {
                    "pdb": pdb,
                    "status": "fail_nonproduction_config",
                    "config": cfg_s,
                }
            )
            continue

        print(f"[{pdb}] production config={cfg} …")
        cf_n0 = probe_cf_total(
            Path(probe),
            receptor=Path(rec),
            pose=Path(nat),
            ligand=Path(nat),
            binary=Path(binary),
            data_dir=Path(data_dir),
            config=Path(cfg),
            wal_coercive=False,
            label=f"{pdb}_nat_off",
        )
        cf_d0 = probe_cf_total(
            Path(probe),
            receptor=Path(rec),
            pose=Path(dec),
            ligand=Path(nat),
            binary=Path(binary),
            data_dir=Path(data_dir),
            config=Path(cfg),
            wal_coercive=False,
            label=f"{pdb}_dec_off",
        )
        cf_n1 = probe_cf_total(
            Path(probe),
            receptor=Path(rec),
            pose=Path(nat),
            ligand=Path(nat),
            binary=Path(binary),
            data_dir=Path(data_dir),
            config=Path(cfg),
            wal_coercive=True,
            label=f"{pdb}_nat_on",
        )
        cf_d1 = probe_cf_total(
            Path(probe),
            receptor=Path(rec),
            pose=Path(dec),
            ligand=Path(nat),
            binary=Path(binary),
            data_dir=Path(data_dir),
            config=Path(cfg),
            wal_coercive=True,
            label=f"{pdb}_dec_on",
        )

        def delta(n: Optional[float], d: Optional[float]) -> Optional[float]:
            if n is None or d is None:
                return None
            return n - d

        d0 = delta(cf_n0, cf_d0)
        d1 = delta(cf_n1, cf_d1)
        rows.append(
            {
                "pdb": pdb,
                "status": "ok",
                "config": str(cfg),
                "cf_native_off": cf_n0,
                "cf_decoy_off": cf_d0,
                "dCF_off": d0,
                "native_wins_off": d0 is not None and d0 <= 0,
                "cf_native_on": cf_n1,
                "cf_decoy_on": cf_d1,
                "dCF_on": d1,
                "native_wins_on": d1 is not None and d1 <= 0,
                "rescued": (
                    d0 is not None
                    and d0 > 0
                    and d1 is not None
                    and d1 <= 0
                ),
                "identical_off_on": (
                    cf_n0 is not None
                    and cf_n1 is not None
                    and abs(cf_n0 - cf_n1) < 1e-6
                    and cf_d0 is not None
                    and cf_d1 is not None
                    and abs(cf_d0 - cf_d1) < 1e-6
                ),
            }
        )

    ok = [r for r in rows if r.get("status") == "ok"]
    n = len(ok)
    n_off = sum(1 for r in ok if r.get("native_wins_off"))
    n_on = sum(1 for r in ok if r.get("native_wins_on"))
    n_rescued = sum(1 for r in ok if r.get("rescued"))
    need = max(1, (n * 7 + 7) // 8) if n else 0
    regressed = [
        r
        for r in ok
        if r.get("native_wins_off") and not r.get("native_wins_on")
    ]
    wall_pass = n >= 1 and n_on >= need and not regressed

    csv_path = out_dir / "wall_oracle.csv"
    keys = sorted({k for r in rows for k in r})
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    md_lines = [
        "# Wall coercive oracle (score-only, production configs)",
        "",
        f"Binary: `{binary}`",
        f"probe_cf: `{probe}`",
        f"Manifest/config policy: **ops/gates/configs/{{PDB}}_dock_config.json** "
        f"(not diagnostic/probe_config.json)",
        f"Panel scored ok: **{n}**",
        f"Native wins OFF: **{n_off}/{n}**",
        f"Native wins ON: **{n_on}/{n}**",
        f"Rescued (fail OFF → pass ON): **{n_rescued}**",
        f"Need ≥ ceil(7n/8) = **{need}**",
        f"Regressed already-min: **{len(regressed)}**",
        f"**VERDICT: {'PASS' if wall_pass else 'FAIL'}**",
        "",
        "| PDB | config | cf_nat_off | cf_dec_off | dCF_off | dCF_on | win_off | win_on |",
        "|-----|--------|-----------:|-----------:|-------:|------:|:-------:|:------:|",
    ]
    for r in ok:
        cfg_name = Path(str(r.get("config", "") or "")).name
        md_lines.append(
            f"| {r['pdb']} | `{cfg_name}` | "
            f"{r.get('cf_native_off')} | {r.get('cf_decoy_off')} | "
            f"{r.get('dCF_off')} | {r.get('dCF_on')} | "
            f"{r.get('native_wins_off')} | {r.get('native_wins_on')} |"
        )
    skips = [r for r in rows if r.get("status") != "ok"]
    if skips:
        md_lines += ["", "## Skips / errors", ""]
        for r in skips:
            md_lines.append(f"- {r.get('pdb')}: {r.get('status')} {r.get('config','')}")

    md_lines += [
        "",
        "## Cadence",
        "",
        "- Phase: STEP 2 wall oracle",
        "- One variable: FLEXAIDDS_WAL_COERCIVE",
        f"- PASS/FAIL: **{'PASS' if wall_pass else 'FAIL'}**",
        "",
    ]
    if not wall_pass:
        md_lines.append(
            "**STOP before memetic / WALL_PILOT_PASS.** Re-diagnose wall / panel."
        )

    (out_dir / "wall_oracle.md").write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    summary = {
        "n_scored": n,
        "native_wins_off": n_off,
        "native_wins_on": n_on,
        "n_rescued": n_rescued,
        "need": need,
        "wall_pilot_pass": wall_pass,
        "binary": str(binary),
        "production_config_policy": "ops/gates/configs/{PDB}_dock_config.json",
    }
    (out_dir / "wall_oracle_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    print(f"Wrote {csv_path}")
    print(f"Wrote {out_dir / 'wall_oracle.md'}")
    return 0 if wall_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
