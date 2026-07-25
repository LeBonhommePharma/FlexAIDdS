#!/usr/bin/env python3
"""
dock_any.py — Fast, strict, any-target / any-ligand docking entry for the flexaidds skill.

Flexible inputs:
  • Local receptor PDB + ligand MOL2/SDF/PDB
  • Optional RCSB PDB ID download (self-docking redock of a HET residue)

Strict by default:
  • resolve_build.py --check (hard-fail; sets FLEXAIDDS_REQUIRE_BUILD=1)
  • ensure_docking_data.py --check (then ensure if missing)
  • Softβ election OFF
  • Local-first OUT under $FLEXAIDDS_LOCAL_ROOT
  • Delegates to scripts/run_flexaidds.sh for provenance sidecars

Usage:
  python3 .grok/skills/flexaidds/scripts/dock_any.py \\
      --receptor target.pdb --ligand ligand.mol2 --temperature 298.15

  python3 .grok/skills/flexaidds/scripts/dock_any.py \\
      --pdb 1STP --ligand-res BTN --temperature 298.15

  # Preflight only (no dock):
  python3 .grok/skills/flexaidds/scripts/dock_any.py --receptor r.pdb --ligand l.mol2 --dry-run
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


def repo_root() -> Path:
    env = os.environ.get("FLEXAIDDS_ROOT", "").strip()
    if env:
        return Path(env).expanduser().resolve()
    here = Path(__file__).resolve()
    for p in [here.parents[i] for i in range(2, 8)]:
        if (p / "AGENTS.md").is_file() or (p / ".git").exists():
            return p
    return here.parents[4]


def run(cmd: list[str], *, cwd: Path | None = None, check: bool = True) -> subprocess.CompletedProcess[str]:
    print("+", " ".join(cmd), flush=True)
    return subprocess.run(cmd, cwd=cwd, text=True, check=check)


def preflight(root: Path, *, skip_ensure: bool) -> None:
    os.environ["FLEXAIDDS_REQUIRE_BUILD"] = "1"
    # Softβ OFF unless user already opted in explicitly to 1
    if os.environ.get("FLEXAIDDS_SOFTBETA_ELECTION", "").strip() not in ("1", "true", "yes"):
        os.environ["FLEXAIDDS_SOFTBETA_ELECTION"] = "0"
    if os.environ.get("FLEXAIDDS_ELECTION_SHANNON_F", "").strip() not in ("1", "true", "yes"):
        os.environ["FLEXAIDDS_ELECTION_SHANNON_F"] = "0"

    skill = root / ".grok" / "skills" / "flexaidds" / "scripts"
    run([sys.executable, str(skill / "resolve_build.py"), "--check", "--repo-root", str(root)])
    ensure = skill / "ensure_docking_data.py"
    if not skip_ensure:
        # check first; if fails, ensure then re-check
        chk = subprocess.run(
            [sys.executable, str(ensure), "--check"],
            cwd=root,
            text=True,
            check=False,
        )
        if chk.returncode != 0:
            run([sys.executable, str(ensure)], cwd=root)
            run([sys.executable, str(ensure), "--check"], cwd=root)


def local_out_dir(root: Path, label: str) -> Path:
    base = Path(
        os.environ.get("FLEXAIDDS_LOCAL_ROOT")
        or (Path.home() / "flexaidds_results")
    ).expanduser()
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = base / "dock_any" / f"{label}_{ts}"
    out.mkdir(parents=True, exist_ok=True)
    return out


def download_pdb(pdb_id: str, dest: Path) -> Path:
    pdb_id = pdb_id.strip().upper()
    dest.parent.mkdir(parents=True, exist_ok=True)
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"Downloading {url} → {dest}", flush=True)
    with urllib.request.urlopen(url, timeout=60) as resp:  # noqa: S310 — fixed RCSB host
        dest.write_bytes(resp.read())
    if dest.stat().st_size < 200:
        raise SystemExit(f"Download too small for {pdb_id}: {dest}")
    return dest


def split_pdb_self_dock(pdb_path: Path, ligand_res: str, work: Path) -> tuple[Path, Path]:
    """Naive ATOM/HETATM split for self-docking redock. Ligand written as PDB (run_flexaidds accepts)."""
    resn = ligand_res.strip().upper()
    rec_lines: list[str] = []
    lig_lines: list[str] = []
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            r = line[17:20].strip().upper()
            if r == resn and line.startswith("HETATM"):
                lig_lines.append(line)
            elif r == resn and line.startswith("ATOM  "):
                # rare: ligand as ATOM
                lig_lines.append(line)
            else:
                # strip waters/ions often helpful for redock TARGET
                if r in {"HOH", "WAT", "DOD", "NA", "CL", "MG", "ZN", "CA", "K", "SO4", "PO4"}:
                    continue
                rec_lines.append(line)
        elif line.startswith(("CRYST1", "TER", "END", "MODEL", "ENDMDL")):
            rec_lines.append(line)
    if not lig_lines:
        raise SystemExit(
            f"No HETATM/ATOM residue {resn!r} found in {pdb_path}. "
            "Pass --ligand-res matching the co-crystal ligand three-letter code."
        )
    if not rec_lines:
        raise SystemExit(f"No receptor atoms left after split for {pdb_path}")
    rec = work / "receptor.pdb"
    lig = work / "ligand.pdb"
    rec.write_text("\n".join(rec_lines) + "\n", encoding="utf-8")
    lig.write_text("\n".join(lig_lines) + "\nEND\n", encoding="utf-8")
    # Prefer MOL2 via openbabel when available
    lig_mol2 = work / "ligand.mol2"
    obabel = shutil.which("obabel")
    if obabel:
        subprocess.run(
            [obabel, str(lig), "-O", str(lig_mol2), "-h", "--partialcharge", "gasteiger"],
            check=False,
            capture_output=True,
            text=True,
        )
        if lig_mol2.is_file() and lig_mol2.stat().st_size > 50:
            return rec, lig_mol2
    return rec, lig


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--receptor", "-r", help="Receptor PDB path")
    p.add_argument("--ligand", "-l", help="Ligand MOL2/SDF/PDB path")
    p.add_argument("--pdb", help="RCSB PDB ID for self-docking download + split")
    p.add_argument("--ligand-res", default="", help="Three-letter ligand residue for --pdb split (e.g. BTN)")
    p.add_argument("--temperature", type=float, default=298.15)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--outdir", "-o", help="Output directory (default: $FLEXAIDDS_LOCAL_ROOT/dock_any/...)")
    p.add_argument("--dry-run", action="store_true", help="Preflight + prepare inputs only")
    p.add_argument("--skip-ensure", action="store_true", help="Skip ensure_docking_data (not recommended)")
    p.add_argument("--visualize", action="store_true", help="Pass --visualize to run_flexaidds.sh")
    p.add_argument("--create-bundle", action="store_true", help="CREATE_BUNDLE=1 for run_flexaidds.sh")
    args = p.parse_args()

    root = repo_root()
    preflight(root, skip_ensure=args.skip_ensure)

    work_label = "custom"
    receptor: Path | None = None
    ligand: Path | None = None

    if args.pdb:
        if not args.ligand_res:
            raise SystemExit("--pdb requires --ligand-res (HET three-letter code)")
        work_label = args.pdb.strip().upper()
        out = Path(args.outdir).expanduser() if args.outdir else local_out_dir(root, work_label)
        work = out / "inputs"
        work.mkdir(parents=True, exist_ok=True)
        raw = download_pdb(args.pdb, work / f"{work_label}.pdb")
        receptor, ligand = split_pdb_self_dock(raw, args.ligand_res, work)
        print(f"Self-docking mode: {work_label} / {args.ligand_res.upper()}")
    else:
        if not args.receptor or not args.ligand:
            raise SystemExit("Provide --receptor and --ligand, or --pdb with --ligand-res")
        receptor = Path(args.receptor).expanduser().resolve()
        ligand = Path(args.ligand).expanduser().resolve()
        if not receptor.is_file() or not ligand.is_file():
            raise SystemExit(f"Missing input files: {receptor} / {ligand}")
        work_label = f"{receptor.stem}_{ligand.stem}"
        out = Path(args.outdir).expanduser() if args.outdir else local_out_dir(root, work_label)
        out.mkdir(parents=True, exist_ok=True)

    assert receptor is not None and ligand is not None
    print(f"OUT: {out}")
    print(f"receptor: {receptor}")
    print(f"ligand:   {ligand}")
    print("docking_mode: self_docking (file/PDB path) — for cross-docking, pass apo receptor + cognate ligand files explicitly")

    if args.dry_run:
        print("DRY-RUN: preflight OK; not launching engine")
        (out / "dock_any_dry_run.txt").write_text(
            f"receptor={receptor}\nligand={ligand}\ntemperature={args.temperature}\n",
            encoding="utf-8",
        )
        return 0

    launcher = root / "scripts" / "run_flexaidds.sh"
    if not launcher.is_file():
        raise SystemExit(f"Missing {launcher}")

    env = os.environ.copy()
    env["FLEXAIDDS_REQUIRE_BUILD"] = "1"
    env["FLEXAIDDS_SOURCE"] = str(root)
    env["SKIP_REBUILD"] = env.get("SKIP_REBUILD", "1")
    if args.create_bundle:
        env["CREATE_BUNDLE"] = "1"
    if args.visualize:
        env["VISUALIZE"] = "1"

    cmd = [
        "bash",
        str(launcher),
        str(receptor),
        str(ligand),
        "--outdir",
        str(out),
        "--temperature",
        str(args.temperature),
        "--seed",
        str(args.seed),
    ]
    if args.visualize:
        cmd.append("--visualize")

    print("+", " ".join(cmd), flush=True)
    proc = subprocess.run(cmd, cwd=root, env=env)
    if proc.returncode != 0:
        print(
            "Engine/wrapper failed. Refuse docking-success language without result artifacts.",
            file=sys.stderr,
        )
        return proc.returncode

    # Deception-proof tip
    print(
        "\nDone. Claim success only if on-disk outputs under OUT pass your gates "
        "(RMSD≤2.0 Å + PoseBusters for pose claims; CF is a scoring proxy, not ΔG)."
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
