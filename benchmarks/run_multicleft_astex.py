#!/usr/bin/env python3
"""Restore the old Astex multi-cleft GA orchestration.

This launcher is intentionally narrow:
- Astex Diverse only.
- Get_Cleft generates ranked *_sph_<rank>.pdb files.
- benchmark_datasets runs one normal FlexAIDdS child per cleft variant.
- Each child receives FLEXAIDDS_CLEFT_SPHERE_FILE for a real cleft grid.
- In known-site mode, each child keeps the original oracle_site_pdb so the
  ligand IC frame remains anchored to the crystal pose while the search grid is
  still the ranked cleft sphere.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ASTEX_JSON = REPO / "benchmarks" / "datasets" / "benchmark_astex_native_85.json"

GET_CLEFT = Path(os.path.expanduser(
    os.environ.get("FLEXAIDDS_GET_CLEFT", "/Users/lp.more/Projects/Get_Cleft/Get_Cleft")
))
BENCHMARK_BIN = Path(os.path.expanduser(
    os.environ.get("FLEXAIDDS_BENCHMARK_BIN", str(REPO / "build_lto" / "benchmark_datasets"))
))
FLEXAIDDS_BIN = Path(os.path.expanduser(
    os.environ.get("FLEXAIDDS_BINARY", str(REPO / "build_lto" / "FlexAIDdS"))
))

SMOKE_CODES = ["1G9V", "1GM8", "1GPK"]


def die(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(2)


def run_checked(cmd: list[str], cwd: Path | None = None) -> None:
    proc = subprocess.run(cmd, cwd=str(cwd) if cwd else None, text=True)
    if proc.returncode != 0:
        die(f"command failed ({proc.returncode}): {' '.join(cmd)}")


def load_astex_entries() -> list[dict]:
    if not ASTEX_JSON.exists():
        die(f"missing Astex JSON: {ASTEX_JSON}")
    data = json.loads(ASTEX_JSON.read_text())
    pairs = data.get("pairs", [])
    if not pairs:
        die(f"Astex JSON has no pairs: {ASTEX_JSON}")
    return pairs


def select_entries(entries: list[dict], args: argparse.Namespace) -> list[dict]:
    if args.targets:
        wanted = [x.strip().upper() for x in args.targets.replace(",", " ").split() if x.strip()]
    elif args.scope == "smoke":
        wanted = SMOKE_CODES
    else:
        wanted = []

    if wanted:
        wanted_set = set(wanted)
        entries = [e for e in entries if e.get("receptor_id", "").upper() in wanted_set]

    if args.max_targets > 0:
        entries = entries[: args.max_targets]

    if not entries:
        die("no Astex entries selected")
    return entries


def generate_clefts(entry: dict, cleft_dir: Path, top_clefts: int, force: bool) -> list[Path]:
    code = entry["receptor_id"]
    receptor = Path(entry["receptor_pdb"])
    if not receptor.exists():
        die(f"missing receptor for {code}: {receptor}")

    cleft_dir.mkdir(parents=True, exist_ok=True)
    prefix = cleft_dir / code
    expected = [cleft_dir / f"{code}_sph_{rank}.pdb" for rank in range(1, top_clefts + 1)]

    def cleft_rank(path: Path) -> int:
        try:
            return int(path.stem.rsplit("_sph_", 1)[1])
        except (IndexError, ValueError):
            return top_clefts + 1

    existing = sorted(
        (p for p in cleft_dir.glob(f"{code}_sph_*.pdb") if cleft_rank(p) <= top_clefts),
        key=cleft_rank,
    )
    if not force and existing:
        return existing

    if force or not all(p.exists() for p in expected):
        for stale in cleft_dir.glob(f"{code}_*"):
            stale.unlink()
        run_checked([
            str(GET_CLEFT),
            "-p", str(receptor),
            "-t", str(top_clefts),
            "-s",
            "-o", str(prefix),
        ])

    sph_files = [p for p in expected if p.exists()]
    if not sph_files:
        die(f"Get_Cleft produced no sphere files for {code} in {cleft_dir}")
    return sph_files


def write_manifest(entries: list[dict], args: argparse.Namespace, out_dir: Path) -> Path:
    cleft_root = out_dir / "clefts"
    pairs: list[dict] = []

    for entry in entries:
        code = entry["receptor_id"]
        oracle_site = entry.get("oracle_site_pdb") or entry.get("binding_site_path") or ""
        if args.mode != "autonomous" and not oracle_site:
            die(f"{code}: known-site mode requested but source entry has no oracle_site_pdb")
        if oracle_site and not Path(oracle_site).exists():
            die(f"{code}: oracle site file is missing: {oracle_site}")

        sph_files = generate_clefts(entry, cleft_root / code, args.top_clefts, args.force_clefts)
        for rank, sph in enumerate(sph_files, 1):
            pair = {
                "receptor_id": f"{code}__clf{rank}",
                "ligand_id": entry.get("ligand_id", code),
                "receptor_pdb": entry["receptor_pdb"],
                "ligand_sdf": entry["ligand_sdf"],
                "rmsd_ref_sdf": entry.get("rmsd_ref_sdf", entry["ligand_sdf"]),
                "cleft_sphere_file": str(sph),
            }
            if args.mode != "autonomous":
                pair["oracle_site_pdb"] = oracle_site
            pairs.append(pair)

    manifest = {
        "schema_version": 1,
        "name": "astex_multicleft_restoration",
        "description": "Astex Diverse expanded to one independent GA per ranked Get_Cleft cleft.",
        "n_pairs": len(pairs),
        "pairs": pairs,
    }
    manifest_path = out_dir / "astex_multicleft_restoration.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest_path


def sanitized_env(omp_threads: int, restarts: int) -> dict[str, str]:
    env = os.environ.copy()
    for key in [
        "FLEXAIDDS_MULTI_CLEFT",
        "FLEXAIDDS_ORACLE_SITE",
        "FLEXAIDDS_ORACLE_SITE_DIR",
        "FLEXAIDDS_THERMO",
        "FLEXAIDDS_T_EFF",
        "FLEXAIDDS_TENCOM_SCALE",
        "FLEXAIDDS_GRID_CACHE_DIR",
    ]:
        env.pop(key, None)
    env["OMP_NUM_THREADS"] = str(omp_threads)
    env["FLEXAIDDS_PARALLEL_RESTARTS"] = "0"
    env["FLEXAIDDS_BINARY"] = str(FLEXAIDDS_BIN)
    env["FLEXAIDDS_RESTARTS"] = str(restarts)
    return env


def build_command(manifest: Path, args: argparse.Namespace, out_dir: Path) -> list[str]:
    results_dir = out_dir / "results"
    cmd = [
        str(BENCHMARK_BIN),
        "--benchmark", f"crossdock_json:{manifest}",
        "--output", str(results_dir),
        "--threads", str(args.workers),
        "--omp-threads", str(args.omp_threads),
        "--ga-generations", str(args.ga_generations),
        "--ga-population", str(args.ga_population),
        "--job-timeout-seconds", str(args.job_timeout_seconds),
        "--mode", args.mode,
    ]
    if args.force_results:
        cmd.append("--force")
    return cmd


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scope", choices=["smoke", "full"], default="smoke")
    parser.add_argument("--targets", default="", help="Comma/space separated Astex codes")
    parser.add_argument("--max-targets", type=int, default=0)
    parser.add_argument("--top-clefts", type=int, default=3)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--omp-threads", type=int, default=1)
    parser.add_argument("--ga-generations", type=int, default=500)
    parser.add_argument("--ga-population", type=int, default=1000)
    parser.add_argument("--restarts", type=int, default=1)
    parser.add_argument("--job-timeout-seconds", type=int, default=3600)
    parser.add_argument("--mode", choices=["autonomous", "oracle-ceiling"], default="oracle-ceiling")
    parser.add_argument("--output-dir", default="")
    parser.add_argument("--force-clefts", action="store_true")
    parser.add_argument("--force-results", action="store_true")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--detach", action="store_true")
    args = parser.parse_args()

    if args.workers * args.omp_threads > 4:
        die("RAM/CPU guard: workers * omp_threads must be <= 4 on the 18 GB MacBook profile")
    if args.top_clefts < 1:
        die("--top-clefts must be >= 1")
    if args.restarts < 1:
        die("--restarts must be >= 1")
    if not GET_CLEFT.exists():
        die(f"missing Get_Cleft binary: {GET_CLEFT}")
    if not BENCHMARK_BIN.exists():
        die(f"missing benchmark_datasets binary: {BENCHMARK_BIN}")
    if not FLEXAIDDS_BIN.exists():
        die(f"missing FlexAIDdS binary: {FLEXAIDDS_BIN}")

    out_dir = Path(args.output_dir) if args.output_dir else REPO / "results" / f"multicleft_{args.scope}_latest"
    out_dir.mkdir(parents=True, exist_ok=True)

    entries = select_entries(load_astex_entries(), args)
    manifest = write_manifest(entries, args, out_dir)
    cmd = build_command(manifest, args, out_dir)
    (out_dir / "run_command.txt").write_text(" ".join(cmd) + "\n")

    print("mission=astex_multicleft_restoration")
    print(f"entries={len(entries)}")
    print(f"top_clefts={args.top_clefts}")
    print(f"manifest={manifest}")
    print(f"results={out_dir / 'results'}")

    if args.prepare_only:
        return 0

    env = sanitized_env(args.omp_threads, args.restarts)
    log_path = out_dir / "benchmark_runner.log"
    run_cmd = cmd
    if shutil.which("caffeinate"):
        run_cmd = ["caffeinate", "-dimsu"] + cmd

    if args.detach:
        log = open(log_path, "a")
        proc = subprocess.Popen(
            run_cmd,
            cwd=str(REPO),
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        (out_dir / "launcher.pid").write_text(str(proc.pid) + "\n")
        print(f"detached_pid={proc.pid}")
        print(f"log={log_path}")
        return 0

    with open(log_path, "a") as log:
        proc = subprocess.run(run_cmd, cwd=str(REPO), env=env, stdout=log, stderr=subprocess.STDOUT)
    return proc.returncode


if __name__ == "__main__":
    raise SystemExit(main())
