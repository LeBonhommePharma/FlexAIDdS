#!/usr/bin/env python3
"""
Publishable blind-docking benchmark — autonomous mode (hardened).

Protocol (claim-safe defaults):
  - Mode:    autonomous (no native seeds; seed_fraction=0)
  - Site:    optional oracle site dir ONLY if FLEXAIDDS_ORACLE_SITE_DIR is set
             by the operator (not forced by this script)
  - Binary:  must be an existing file (FLEXAIDDS_BINARY)
  - Seed elitism: OFF unless FLEXAIDDS_SEED_ELITISM=1 is set by the operator
  - Softcore WAL: not forced (SOFTWA soft-core is already engine default)

Historical baselines (v123 session notes):
  v50b autonomous : 69/85 = 81.2%
  v127 oracle-ceiling : 78/85 = 91.8%  (NOT publishable)
  rDock (literature) : 88.2%

Paths resolve from FLEXAIDDS_* env vars / repo root — no machine-absolute hardcoding.
"""
from __future__ import annotations

import datetime
import json
import os
import subprocess
import sys
from pathlib import Path


def _repo_root() -> Path:
    env = os.environ.get("FLEXAIDDS_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    return Path(__file__).resolve().parents[1]


def main() -> int:
    repo = _repo_root()
    local_root = Path(
        os.environ.get("FLEXAIDDS_LOCAL_ROOT", str(Path.home() / "flexaidds_results"))
    ).expanduser()

    bin_path = Path(
        os.environ.get(
            "FLEXAIDDS_BINARY",
            str(local_root / "three_engine_entropy_q1" / "bin" / "C"),
        )
    ).expanduser()
    if not bin_path.is_file():
        print(
            f"[publishable_blind] ERROR: FLEXAIDDS_BINARY is not a file: {bin_path}\n"
            "  Build/pin a claim binary and export FLEXAIDDS_BINARY=/path/to/FlexAID",
            file=sys.stderr,
        )
        return 2

    data_dir = Path(
        os.environ.get("FLEXAIDDS_DATA_DIR", str(repo / "build"))
    ).expanduser()
    ds = Path(
        os.environ.get(
            "FLEXAIDDS_BENCHMARK_DATASETS",
            str(data_dir / "benchmark_datasets"),
        )
    ).expanduser()
    if not ds.exists():
        print(
            f"[publishable_blind] ERROR: benchmark runner not found: {ds}\n"
            "  Set FLEXAIDDS_BENCHMARK_DATASETS or build ENABLE_BENCHMARK_DATASETS.",
            file=sys.stderr,
        )
        return 2

    tag = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    out = Path(
        os.environ.get(
            "FLEXAIDDS_PUBLISHABLE_OUT",
            str(local_root / "campaigns" / f"publishable_{tag}_blind"),
        )
    ).expanduser()
    out.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env["FLEXAIDDS_BINARY"] = str(bin_path)
    env["FLEXAIDDS_DATA_DIR"] = str(data_dir)
    # Do NOT force ORACLE_SITE_DIR — publishable blind should not silently
    # attach crystal sites unless the operator set the env.
    # Do NOT force SEED_ELITISM / SOFTCORE_WAL — claim knobs stay operator-owned.
    env.setdefault("FLEXAIDDS_RESTARTS", "5")
    env.setdefault("FLEXAIDDS_PARALLEL_RESTARTS", "1")
    env.setdefault("FLEXAIDDS_EVAL_SCALE_DIHEDRAL", "1")
    env.setdefault("FLEXAIDDS_SEED_ELITISM", "0")  # autonomous publishable default
    env.setdefault("FLEXAIDDS_BUDGET_SCALE", "1")
    env.setdefault("FLEXAIDDS_ALLOW_CONCURRENT", "1")

    # Refuse claim OUT contamination
    out_s = str(out)
    if "claim" in out_s.lower() and os.environ.get("FLEXAIDDS_ALLOW_CLAIM_OUT") != "1":
        print(
            f"[publishable_blind] ERROR: refusing OUT path that looks like a claim dir: {out}\n"
            "  Set FLEXAIDDS_ALLOW_CLAIM_OUT=1 to override.",
            file=sys.stderr,
        )
        return 2

    cmd = [
        str(ds),
        "--benchmark",
        "astex",
        "--mode",
        "autonomous",
        "--restarts",
        env["FLEXAIDDS_RESTARTS"],
        "--threads",
        env.get("FLEXAIDDS_THREADS", "2"),
        "--omp-threads",
        env.get("FLEXAIDDS_OMP_THREADS", "5"),
        "--output",
        str(out),
    ]

    # Provenance before launch
    try:
        import hashlib

        h = hashlib.sha256()
        with open(bin_path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        bin_sha = h.hexdigest()
    except OSError:
        bin_sha = ""

    provenance = {
        "launcher": "scripts/launch_publishable_blind.py",
        "mode": "autonomous",
        "binary": str(bin_path),
        "binary_sha256": bin_sha,
        "data_dir": str(data_dir),
        "runner": str(ds),
        "output": str(out),
        "oracle_site_dir_set": bool(env.get("FLEXAIDDS_ORACLE_SITE_DIR")),
        "seed_elitism": env.get("FLEXAIDDS_SEED_ELITISM", "0"),
        "cmd": cmd,
        "started_utc": datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
    }
    (out / "PROVENANCE.json").write_text(json.dumps(provenance, indent=2) + "\n")

    foreground = os.environ.get("FLEXAIDDS_PUBLISHABLE_FOREGROUND", "0") == "1"
    log_path = out / "stdout.log"

    if foreground or not hasattr(os, "fork"):
        print(f"[publishable_blind] Running foreground → {out}")
        print(f"[publishable_blind] Binary: {bin_path}")
        with open(log_path, "w") as log:
            proc = subprocess.run(cmd, env=env, cwd=str(repo), stdout=log, stderr=log)
        return int(proc.returncode)

    # Double-fork daemon (POSIX), with pid file
    log = open(log_path, "w", buffering=1)
    pid = os.fork()
    if pid != 0:
        print(f"[publishable_blind] Launched parent PID {pid}")
        print(f"[publishable_blind] Binary: {bin_path}")
        print(f"[publishable_blind] Mode: autonomous (seed_elitism={env.get('FLEXAIDDS_SEED_ELITISM')})")
        print(f"[publishable_blind] Output: {out}")
        print(f"[publishable_blind] Provenance: {out / 'PROVENANCE.json'}")
        return 0

    os.setsid()
    pid2 = os.fork()
    if pid2 != 0:
        sys.exit(0)

    (out / "worker.pid").write_text(f"{os.getpid()}\n")
    os.chdir(repo)
    proc = subprocess.Popen(
        cmd,
        env=env,
        stdout=log,
        stderr=log,
        start_new_session=True,
    )
    log.write(f"[publishable_blind] Worker PID {proc.pid}\n")
    log.flush()
    (out / "worker.pid").write_text(f"{proc.pid}\n")
    rc = proc.wait()
    log.close()
    sys.exit(rc)


if __name__ == "__main__":
    sys.exit(main())
