#!/usr/bin/env python3
"""W1.3 / W2 score-only wall oracle: FLEXAIDDS_WAL_COERCIVE OFF vs ON via probe_cf.

No GA. For each panel target with diagnostic native.sdf + falsemin_armA.pdb,
scores CF with the existing FlexAIDdS binary under WAL_COERCIVE=0 and =1.

Acceptance (plan): native CF <= decoy CF (dCF = native - decoy <= 0) on as many
panel targets as possible when ON, especially where OFF fails.

Usage:
  python3 scripts/wall_coercive_oracle.py \\
    --probe-cf /path/to/probe_cf \\
    --binary /path/to/FlexAIDdS \\
    --data-dir /path/to/data \\
    --out-dir ~/flexaidds_results/workorders/wall_oracle

Does not rebuild binaries. Safe to run while other campaigns use different OUT.
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

DEFAULT_PANEL = ("1G9V", "1M2Z", "1N1M", "1J3J", "1K3U", "1L7F", "1HNN", "1HP0")


def find_repo_root() -> Path:
    p = Path(__file__).resolve().parents[1]
    if (p / "AGENTS.md").exists():
        return p
    return Path.cwd()


def resolve_receptor(repo: Path, pdb: str) -> Optional[Path]:
    candidates = [
        Path.home() / f".flexaidds/benchmarks/astex_diverse/{pdb}/{pdb}_apo.pdb",
        repo / f"benchmarks/astex_diverse/astex_diverse/{pdb}/{pdb}_apo.pdb",
        repo / f"benchmarks/astex_diverse/{pdb}/{pdb}_apo.pdb",
    ]
    for c in candidates:
        if c.is_file():
            return c
    return None


def resolve_pair(repo: Path, pdb: str) -> tuple[Optional[Path], Optional[Path]]:
    d = repo / "diagnostic" / "refs" / pdb
    nat = d / "native.sdf"
    false = d / "falsemin_armA.pdb"
    return (nat if nat.is_file() else None, false if false.is_file() else None)


def probe_cf_total(
    probe: Path,
    *,
    receptor: Path,
    pose: Path,
    ligand: Path,
    binary: Path,
    data_dir: Path,
    config: Optional[Path],
    wal_coercive: bool,
    label: str,
) -> Optional[float]:
    env = os.environ.copy()
    env["FLEXAIDDS_WAL_COERCIVE"] = "1" if wal_coercive else "0"
    # Ensure CAP does not confound
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
        "--pdb",
        label,
    ]
    if config and config.is_file():
        cmd.extend(["--config", str(config)])
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
    # last JSON line with cf_total
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
    ap.add_argument(
        "--probe-cf",
        type=Path,
        default=None,
        help="probe_cf binary (default: repo/build/probe_cf or PATH)",
    )
    ap.add_argument("--binary", type=Path, default=None)
    ap.add_argument("--data-dir", type=Path, default=None)
    ap.add_argument("--config", type=Path, default=None)
    ap.add_argument("--panel", nargs="*", default=list(DEFAULT_PANEL))
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(argv)

    repo = args.repo.expanduser().resolve()
    probe = args.probe_cf
    if probe is None:
        for c in (repo / "build" / "probe_cf", Path("/Users/lp.more/Projects/FlexAIDdS/build/probe_cf")):
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
            Path("/Users/lp.more/Projects/FlexAIDdS/build/FlexAIDdS"),
            Path.home() / "flexaidds_results/baseline_engine/FlexAIDdS",
        ):
            if c.is_file():
                binary = c
                break
    if binary is None or not Path(binary).is_file():
        print("error: FlexAIDdS binary not found", file=sys.stderr)
        return 2

    data_dir = args.data_dir or Path(binary).parent
    config = args.config or (repo / "diagnostic" / "probe_config.json")
    out_dir = args.out_dir.expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    for pdb in args.panel:
        rec = resolve_receptor(repo, pdb)
        nat, false = resolve_pair(repo, pdb)
        if not rec or not nat or not false:
            rows.append(
                {
                    "pdb": pdb,
                    "status": "skip_missing_inputs",
                    "receptor": str(rec) if rec else "",
                    "native": str(nat) if nat else "",
                    "decoy": str(false) if false else "",
                }
            )
            continue
        print(f"[{pdb}] scoring OFF/ON …")
        cf_n0 = probe_cf_total(
            Path(probe),
            receptor=rec,
            pose=nat,
            ligand=nat,
            binary=Path(binary),
            data_dir=Path(data_dir),
            config=config if config.is_file() else None,
            wal_coercive=False,
            label=f"{pdb}_nat_off",
        )
        cf_d0 = probe_cf_total(
            Path(probe),
            receptor=rec,
            pose=false,
            ligand=nat,
            binary=Path(binary),
            data_dir=Path(data_dir),
            config=config if config.is_file() else None,
            wal_coercive=False,
            label=f"{pdb}_dec_off",
        )
        cf_n1 = probe_cf_total(
            Path(probe),
            receptor=rec,
            pose=nat,
            ligand=nat,
            binary=Path(binary),
            data_dir=Path(data_dir),
            config=config if config.is_file() else None,
            wal_coercive=True,
            label=f"{pdb}_nat_on",
        )
        cf_d1 = probe_cf_total(
            Path(probe),
            receptor=rec,
            pose=false,
            ligand=nat,
            binary=Path(binary),
            data_dir=Path(data_dir),
            config=config if config.is_file() else None,
            wal_coercive=True,
            label=f"{pdb}_dec_on",
        )

        def delta(n, d):
            if n is None or d is None:
                return None
            return n - d  # <=0 means native better (lower CF)

        d0 = delta(cf_n0, cf_d0)
        d1 = delta(cf_n1, cf_d1)
        native_wins_off = d0 is not None and d0 <= 0
        native_wins_on = d1 is not None and d1 <= 0
        rows.append(
            {
                "pdb": pdb,
                "status": "ok",
                "cf_native_off": cf_n0,
                "cf_decoy_off": cf_d0,
                "dCF_off": d0,
                "native_wins_off": native_wins_off,
                "cf_native_on": cf_n1,
                "cf_decoy_on": cf_d1,
                "dCF_on": d1,
                "native_wins_on": native_wins_on,
                "rescued": (not native_wins_off) and native_wins_on,
            }
        )

    ok = [r for r in rows if r.get("status") == "ok"]
    n_off = sum(1 for r in ok if r.get("native_wins_off"))
    n_on = sum(1 for r in ok if r.get("native_wins_on"))
    n_rescued = sum(1 for r in ok if r.get("rescued"))
    n = len(ok)
    bar_78 = n_on >= 7 if n >= 8 else (n_on >= max(1, (n * 7 + 7) // 8))
    wall_pass = n_on > n_off and (n_rescued > 0 or n_on == n)

    csv_path = out_dir / "wall_oracle.csv"
    fields = list(rows[0].keys()) if rows else ["pdb"]
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}))
        w.writeheader()
        for r in rows:
            w.writerow(r)

    md = [
        "# Wall coercive oracle (score-only)",
        "",
        f"Binary: `{binary}`",
        f"probe_cf: `{probe}`",
        f"Panel scored: {n}",
        f"Native wins OFF: **{n_off}/{n}**",
        f"Native wins ON (WAL_COERCIVE=1): **{n_on}/{n}**",
        f"Rescued (OFF fail → ON pass): **{n_rescued}**",
        f"≥7/8 style bar: **{'PASS' if bar_78 else 'FAIL'}** (scored n={n})",
        f"Wall pilot PASS recommendation: **{'YES' if wall_pass else 'NO'}** "
        f"(set FLEXAIDDS_WALL_PILOT_PASS only if YES)",
        "",
        "| pdb | dCF_off | win_off | dCF_on | win_on | rescued |",
        "|-----|---------|---------|--------|--------|---------|",
    ]
    for r in rows:
        if r.get("status") != "ok":
            md.append(f"| {r['pdb']} | — | — | — | — | skip |")
            continue
        md.append(
            f"| {r['pdb']} | {r.get('dCF_off')} | {r.get('native_wins_off')} | "
            f"{r.get('dCF_on')} | {r.get('native_wins_on')} | {r.get('rescued')} |"
        )
    md.append("")
    md_path = out_dir / "wall_oracle.md"
    md_path.write_text("\n".join(md) + "\n", encoding="utf-8")

    summary = {
        "n_scored": n,
        "native_wins_off": n_off,
        "native_wins_on": n_on,
        "n_rescued": n_rescued,
        "bar_7_of_8_style": bar_78,
        "wall_pilot_pass": wall_pass,
        "binary": str(binary),
        "rows": rows,
    }
    (out_dir / "wall_oracle.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: summary[k] for k in summary if k != "rows"}, indent=2))
    print(f"Wrote {csv_path}")
    print(f"Wrote {md_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
