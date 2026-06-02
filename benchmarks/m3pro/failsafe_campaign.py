#!/usr/bin/env python3
"""Failsafe local runner for FlexAIDdS M3 Pro (and portable) benchmark repetitions.

Execution policy (when remote sync is enabled):
  - execute binaries + hot paths from fast local storage (never from slow remote);
  - sync every completed run, status, logs, and final analysis to the configured remote
    (iCloud on M3 Pro, or any other durable location via --remote-base).

New in this version (production-grade portability):
  - --remote-base + --no-remote-sync for general use beyond the original M3 Pro rig.
  - Full FLEXAIDDS_* env integration (source ~/.flexaidds_env and it just works).
  - Auto repo detection, configurable lock, $TMPDIR defaults, etc.

Use --preflight-only for safe dry-runs on any machine / by Codex agents.

NOTE (P1.10 / canonical protocol): For the production full first runs at exact 298 K / 310 K (Astex Diverse + Non-Native) that deliver the best BindingMode with full thermo ledger, Metal accel, iCloud-only artifacts, and validity gating, ALWAYS use the flexaid-docking skill canonical launcher:
  bash .grok/skills/flexaid-docking/scripts/launch_full_benchmark.sh <dataset> <298|310> <name>
See .grok/skills/flexaid-docking/SKILL.md "M3 Pro iCloud Canonical Best-BindingMode Protocol" (and the 4 launched full-*-TS dirs). This failsafe script is for repetition debugging / custom OMP tests only.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import selectors
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


DEFAULT_DATASETS = ("astex", "astex_nonnative", "hap2")
DEFAULT_REPO = Path("/Users/lp.more/Projects/FlexAIDdS")
DEFAULT_BUILD = Path("/private/tmp/flexaidds-benchmark-omp-build")
DEFAULT_ICLOUD_BASE = Path.home() / "Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS"
EXPECTED_MIN_SYSTEMS = {
    "astex": 85,
    "astex_diverse": 85,
    "astex_nonnative": 100,
    "hap2": 59,
}
CACHE_SUBDIRS = {
    "astex": ("astex_diverse",),
    "astex_diverse": ("astex_diverse",),
    "astex_nonnative": ("astex_nonnative", "astex_diverse"),
    "hap2": ("hap2",),
}
FATAL_PATTERNS = (
    "Segmentation fault",
    "segmentation fault",
    "Bus error",
    "Trace/BPT trap",
    "command not found",
    "dyld:",
    "Abort trap",
)


class CampaignError(RuntimeError):
    pass


@dataclass
class RunMetrics:
    total_systems: int
    successful: int
    success_rate: float
    results_csv: Path
    summary_csv: Path | None


CURRENT_PGID: int | None = None
STOP_REQUESTED = False


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def local_now_tag() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def is_icloud_path(path: Path) -> bool:
    s = str(path)
    return "Mobile Documents" in s or "com~apple~CloudDocs" in s


def require_local_path(name: str, path: Path, strict: bool = True) -> None:
    if is_icloud_path(path):
        msg = f"{name} is on iCloud: {path}"
        if strict:
            raise CampaignError(f"{name} must be local APFS, not iCloud: {path}")
        else:
            # On pure iCloud-only machines (M3 Pro 18GB policy), the installed build can live on iCloud
            # as long as actual GA execution uses a fast local --local-base.
            print(f"[failsafe] WARNING: {msg} — acceptable on iCloud-only rigs. Execution will still use local hot paths.")


def mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def write_json(path: Path, payload: dict) -> None:
    mkdir(path.parent)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def signal_handler(signum: int, _frame) -> None:
    global STOP_REQUESTED
    STOP_REQUESTED = True
    if CURRENT_PGID is not None:
        try:
            os.killpg(CURRENT_PGID, signal.SIGTERM)
        except ProcessLookupError:
            pass


def run_logged(
    cmd: list[str],
    *,
    cwd: Path,
    env: dict[str, str],
    log_path: Path,
    timeout_s: int | None,
) -> int:
    """Run a subprocess in its own process group, streaming output to a log."""
    global CURRENT_PGID
    mkdir(log_path.parent)
    started = time.monotonic()
    with log_path.open("a", encoding="utf-8", errors="replace") as log:
        log.write(f"\n[{utc_now()}] RUN {' '.join(cmd)}\n")
        log.flush()
        proc = subprocess.Popen(
            cmd,
            cwd=str(cwd),
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            preexec_fn=os.setsid,
        )
        CURRENT_PGID = os.getpgid(proc.pid)
        sel = selectors.DefaultSelector()
        assert proc.stdout is not None
        sel.register(proc.stdout, selectors.EVENT_READ)

        try:
            while True:
                if STOP_REQUESTED:
                    raise CampaignError("stop requested")
                if timeout_s and time.monotonic() - started > timeout_s:
                    try:
                        os.killpg(CURRENT_PGID, signal.SIGTERM)
                    except ProcessLookupError:
                        pass
                    time.sleep(5)
                    if proc.poll() is None:
                        try:
                            os.killpg(CURRENT_PGID, signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                    raise CampaignError(f"timeout after {timeout_s}s: {' '.join(cmd)}")

                events = sel.select(timeout=1.0)
                for key, _ in events:
                    line = key.fileobj.readline()
                    if line:
                        log.write(line)
                        log.flush()

                if proc.poll() is not None:
                    remaining = proc.stdout.read()
                    if remaining:
                        log.write(remaining)
                    rc = proc.wait()
                    log.write(f"[{utc_now()}] EXIT {rc}\n")
                    log.flush()
                    return rc
        finally:
            CURRENT_PGID = None
            try:
                sel.unregister(proc.stdout)
            except Exception:
                pass


def rsync_copy(src: Path, dst: Path | None, *, delete: bool = False) -> None:
    """Production-grade rsync helper.

    When dst is None (--no-remote-sync or no remote configured), this is a silent no-op.
    This makes the entire remote sync subsystem optional in a bulletproof way.
    """
    if dst is None:
        return
    if not src.exists():
        return
    mkdir(dst)
    rsync = shutil.which("rsync")
    if rsync:
        cmd = [rsync, "-a"]
        if delete:
            cmd.append("--delete")
        cmd.extend([str(src) + "/", str(dst) + "/"])
        subprocess.run(cmd, check=True)
    else:
        shutil.copytree(src, dst, dirs_exist_ok=True)


def stage_runtime_data(repo: Path, build: Path) -> dict:
    """Stage all runtime .dat/.def/.lst files onto local APFS hot paths."""
    source = repo / "WRK"
    if not source.is_dir():
        raise CampaignError(f"missing WRK runtime data directory: {source}")

    files = sorted(
        p for p in source.iterdir()
        if p.is_file() and p.suffix.lower() in {".dat", ".def", ".lst"}
    )
    if not files:
        raise CampaignError(f"no .dat/.def/.lst runtime files found in {source}")

    staged_by_name = {p.name: p for p in files}
    if "rotobs.lst" not in staged_by_name:
        for candidate in (
            repo / "WRK" / "build_iee" / "rotobs.lst",
            repo / "build-ci-fix" / "rotobs.lst",
            repo / "build-xcode" / "Release" / "rotobs.lst",
            repo / ".grok" / "skills" / "flexaid-docking" / "data" / "rotobs.lst",
        ):
            if candidate.is_file():
                staged_by_name["rotobs.lst"] = candidate
                break
    if "rotobs.lst" not in staged_by_name:
        raise CampaignError("missing required runtime data file rotobs.lst")

    files = [staged_by_name[name] for name in sorted(staged_by_name)]
    destinations = []
    for dst in (build, build / "WRK", Path("/private/tmp") / "WRK", Path(tempfile.gettempdir()) / "WRK"):
        if dst not in destinations:
            destinations.append(dst)
    for dst in destinations:
        mkdir(dst)
        for src in files:
            shutil.copy2(src, dst / src.name)

    return {
        "source": str(source),
        "destinations": [str(p) for p in destinations],
        "file_count": len(files),
        "files": [p.name for p in files],
    }


def seed_cache(remote_cache: Path | None, local_cache: Path, datasets: list[str], *, include_smoke: bool) -> None:
    """Safe cache seeder — no-op if no remote configured."""
    if remote_cache is None:
        mkdir(local_cache)
        return
    mkdir(local_cache)
    wanted: set[str] = set()
    if include_smoke:
        wanted.add("astex_diverse")
    for dataset in datasets:
        wanted.update(CACHE_SUBDIRS.get(dataset, (dataset,)))
    for subdir in sorted(wanted):
        src = remote_cache / subdir
        dst = local_cache / subdir
        if src.exists():
            print(f"[failsafe] seeding cache subdir {subdir}")
            rsync_copy(src, dst, delete=False)
        else:
            print(f"[failsafe] cache subdir absent, runner may download: {subdir}")


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f))


def find_single(patterns: Iterable[str], directory: Path) -> Path | None:
    matches: list[Path] = []
    for pattern in patterns:
        matches.extend(sorted(directory.glob(pattern)))
    return matches[0] if matches else None


def parse_run_metrics(run_dir: Path) -> RunMetrics:
    results_csv = find_single(("*_results.csv",), run_dir)
    if results_csv is None:
        raise CampaignError(f"missing aggregate results CSV in {run_dir}")
    rows = read_csv_rows(results_csv)
    if not rows:
        raise CampaignError(f"empty aggregate results CSV: {results_csv}")

    total = len(rows)
    successes = 0
    for row in rows:
        try:
            successes += int(float(row.get("success", "0") or "0"))
        except ValueError:
            pass

    summary_csv = find_single(("*_summary.csv",), run_dir)
    if summary_csv is not None:
        summary_rows = read_csv_rows(summary_csv)
        if summary_rows:
            sr = summary_rows[0]
            try:
                total = int(float(sr.get("total_systems", total)))
                successes = int(float(sr.get("successful", successes)))
            except ValueError:
                pass

    success_rate = successes / total if total else 0.0
    return RunMetrics(total, successes, success_rate, results_csv, summary_csv)


def scan_fatal_logs(run_dir: Path) -> list[str]:
    findings: list[str] = []
    for path in run_dir.rglob("*.log"):
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for pattern in FATAL_PATTERNS:
            if pattern in text:
                findings.append(f"{path}: {pattern}")
                break
    return findings


def validate_run(
    *,
    dataset: str,
    run_dir: Path,
    require_nonzero_success: bool,
    expected_min_systems: int | None,
) -> RunMetrics:
    metrics = parse_run_metrics(run_dir)
    if expected_min_systems and metrics.total_systems < expected_min_systems:
        raise CampaignError(
            f"{dataset}: expected at least {expected_min_systems} systems, got {metrics.total_systems}"
        )
    fatal = scan_fatal_logs(run_dir)
    if fatal:
        raise CampaignError("fatal log signatures:\n" + "\n".join(fatal[:20]))
    if require_nonzero_success and metrics.successful == 0:
        raise CampaignError(f"{dataset}: zero successful targets after completed run")
    marker = {
        "dataset": dataset,
        "validated_at": utc_now(),
        "total_systems": metrics.total_systems,
        "successful": metrics.successful,
        "success_rate": metrics.success_rate,
        "results_csv": str(metrics.results_csv),
        "summary_csv": str(metrics.summary_csv) if metrics.summary_csv else None,
    }
    write_json(run_dir / "RUN_OK.json", marker)
    return metrics


def write_status(status_file: Path, **updates) -> None:
    payload = {
        "updated_at": utc_now(),
        **updates,
    }
    write_json(status_file, payload)


def acquire_lock(lock_dir: Path) -> None:
    try:
        lock_dir.mkdir()
        (lock_dir / "pid").write_text(str(os.getpid()) + "\n")
        return
    except FileExistsError:
        pid_file = lock_dir / "pid"
        try:
            pid = int(pid_file.read_text().strip())
        except Exception as exc:
            raise CampaignError(f"campaign lock exists without readable pid: {lock_dir}") from exc
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            shutil.rmtree(lock_dir, ignore_errors=True)
            lock_dir.mkdir()
            pid_file.write_text(str(os.getpid()) + "\n")
            return
        except PermissionError as exc:
            raise CampaignError(f"campaign lock held by pid {pid}: {lock_dir}") from exc
        raise CampaignError(f"campaign lock held by live pid {pid}: {lock_dir}")


def run_analysis(repo: Path, local_results: Path, local_analysis: Path, n_bootstrap: int, log: Path) -> int:
    analyzer = repo / "benchmarks/m3pro/analyze_repetitions.py"
    if not analyzer.exists():
        return 127
    cmd = [
        sys.executable,
        str(analyzer),
        "--results-dir",
        str(local_results / "tier2"),
        "--n-bootstrap",
        str(n_bootstrap),
        "--output-dir",
        str(local_analysis),
    ]
    env = os.environ.copy()
    return run_logged(cmd, cwd=repo, env=env, log_path=log, timeout_s=3600)


def _resolve_repo(default: Path | None) -> Path:
    """Best-effort portable repo root detection for Codex / multi-user / CI use.

    Order: explicit --repo > FLEXAIDDS_REPO env > git rev-parse --show-toplevel > cwd > original M3 default.
    """
    if default:
        return default
    env = os.environ.get("FLEXAIDDS_REPO")
    if env:
        return Path(env)
    try:
        root = subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"], text=True, stderr=subprocess.DEVNULL
        ).strip()
        if root:
            return Path(root)
    except Exception:
        pass
    cwd = Path.cwd()
    if (cwd / "LIB" / "flexaid.h").exists() or (cwd / "CMakeLists.txt").exists():
        return cwd
    # Last resort: keep the (user-specific) default so existing M3 setups don't break
    return default or DEFAULT_REPO


def _load_flexaidds_env() -> dict[str, str]:
    """Bulletproof loader for FLEXAIDDS_* variables.

    Priority:
      1. Already-present os.environ values (user did `source ~/.flexaidds_env` before python)
      2. Safe parse of ~/.flexaidds_env (or $FLEXAIDDS_ENV_FILE) without executing shell

    Never raises. Returns only FLEXAIDDS_* keys. Handles common quoting.
    This matches the set -a / source pattern used by the companion .sh scripts.
    """


def _warn_non_portable_remote(remote_base: Path | None) -> None:
    """Gentle production warning for the M3-specific iCloud policy on other platforms."""
    if remote_base is None:
        return
    if is_icloud_path(remote_base) and platform.system() != "Darwin":
        print(
            "[failsafe] WARNING: iCloud-style remote path detected on non-macOS. "
            "The 'local APFS only + rsync to remote' policy and path checks are tuned for "
            "the original M3 Pro + iCloud 2TB environment. Things may still work, but "
            "you are in less-tested territory.",
            file=sys.stderr,
        )
    env: dict[str, str] = {}

    # 1. Respect already-exported environment (highest priority)
    for k, v in os.environ.items():
        if k.startswith("FLEXAIDDS_"):
            env[k] = v

    if env:
        return env  # User already sourced — use it as-is

    # 2. Safe file parse (no shell, no eval)
    env_file = os.environ.get("FLEXAIDDS_ENV_FILE")
    if not env_file:
        env_file = str(Path.home() / ".flexaidds_env")

    p = Path(env_file)
    if not p.exists():
        return env

    try:
        for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("export "):
                line = line[7:].strip()
            if "=" not in line:
                continue
            k, v = line.split("=", 1)
            k = k.strip()
            if not k.startswith("FLEXAIDDS_"):
                continue
            v = v.strip()
            # Strip common quoting
            if (v.startswith('"') and v.endswith('"')) or (v.startswith("'") and v.endswith("'")):
                v = v[1:-1]
            env[k] = v
    except Exception as exc:
        # Production-grade: never let a bad env file kill the script
        print(f"[failsafe] WARNING: failed to parse {env_file}: {exc}", file=sys.stderr)

    return env


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__ + "\n\n"
        "Portability note: Developed for a specific M3 Pro + iCloud-only setup. "
        "All critical paths are overridable via CLI flags or FLEXAIDDS_* environment variables "
        "(see the companion .sh scripts). Use --preflight-only for safe testing on other machines/Codex agents."
    )
    p.add_argument("--repo", type=Path, default=None,
                   help="Path to FlexAIDdS repo root (auto-detected via FLEXAIDDS_REPO, git, or cwd if omitted)")
    p.add_argument("--build", type=Path, default=DEFAULT_BUILD)
    p.add_argument("--icloud-base", type=Path, default=DEFAULT_ICLOUD_BASE,
                   help="Legacy alias for --remote-base (M3 Pro iCloud users)")
    p.add_argument("--remote-base", type=Path, default=None,
                   help="Remote durable storage base (generalization of iCloud). All results/logs/data are rsynced here after each step.")
    p.add_argument("--no-remote-sync", action="store_true",
                   help="Disable all remote (iCloud / other) synchronization. Useful for pure local or CI runs.")
    p.add_argument("--run-id", default=f"m3pro_failsafe_10rep_{local_now_tag()}")
    p.add_argument("--local-base", type=Path, default=None)
    p.add_argument("--local-cache", type=Path, default=Path(tempfile.gettempdir()) / "flexaidds_benchmark_cache")
    p.add_argument("--datasets", nargs="+", default=list(DEFAULT_DATASETS))
    p.add_argument("--runs", type=int, default=10)
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--omp-threads", type=int, default=6)
    # GA generations: benchmark plan §6.1 canonical spec = 500 gen × 1000 chrom = 510k evals.
    # Shannon HSC early-stop (H < 1.3863 nats) terminates converged runs before 500.
    # Use --ga-generations 2000 explicitly only for exploratory / hard-landscape runs.
    p.add_argument("--ga-generations", type=int, default=500)
    p.add_argument("--ga-population", type=int, default=1000)
    p.add_argument("--temperature", type=float, default=300.0)
    p.add_argument("--clustering", default="FO")
    p.add_argument("--bootstrap", type=int, default=10000)
    p.add_argument("--run-timeout-hours", type=float, default=168.0)
    p.add_argument("--dock-timeout-hours", type=float, default=12.0)
    p.add_argument("--resume", action="store_true")
    p.add_argument("--skip-smoke", action="store_true")
    p.add_argument("--seed-cache", action="store_true")
    p.add_argument("--mirror-cache", action="store_true")
    p.add_argument("--preflight-only", action="store_true")
    p.add_argument("--lock-dir", type=Path, default=None,
                   help="Global campaign lock dir (defaults to $TMPDIR/flexaidds_campaign.lock for portability)")
    return p.parse_args()


def main() -> int:
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)
    args = parse_args()

    flex_env = _load_flexaidds_env() or {}

    # Apply FLEXAIDDS_* env as high-priority fallbacks (matches .sh campaign behavior)
    # This makes the Python script "just work" after the user has sourced ~/.flexaidds_env
    if not args.repo and "FLEXAIDDS_REPO" in flex_env:
        args.repo = Path(flex_env["FLEXAIDDS_REPO"])
    if str(args.build) == str(DEFAULT_BUILD) and "FLEXAIDDS_BUILD" in flex_env:
        args.build = Path(flex_env["FLEXAIDDS_BUILD"])
    if str(args.icloud_base) == str(DEFAULT_ICLOUD_BASE) and "FLEXAIDDS_ICLOUD" in flex_env:
        args.icloud_base = Path(flex_env["FLEXAIDDS_ICLOUD"])
    if str(args.local_cache) == str(Path(tempfile.gettempdir()) / "flexaidds_benchmark_cache") and "FLEXAIDDS_BENCHMARK_DATA" in flex_env:
        args.local_cache = Path(flex_env["FLEXAIDDS_BENCHMARK_DATA"])

    repo = _resolve_repo(args.repo).resolve()
    build = args.build.resolve()
    binary = build / "benchmark_datasets"
    docking_binary = build / "FlexAID"
    local_base = (args.local_base or Path("/private/tmp") / "flexaidds_campaigns" / args.run_id).resolve()
    local_cache = args.local_cache.resolve()
    local_results = local_base / "results"
    local_logs = local_base / "logs"
    local_analysis = local_results / "analysis"
    # Remote (durable) base resolution — production-grade generalization
    # Priority: --remote-base > --icloud-base (legacy) > FLEXAIDDS_ICLOUD env > hard default
    remote_base = args.remote_base or args.icloud_base
    if args.no_remote_sync:
        remote_base = None

    if remote_base:
        remote_base = remote_base.resolve()
        remote_results = remote_base / "results" / args.run_id
        remote_logs = remote_base / "logs" / args.run_id
        remote_cache = remote_base / "benchmark_data"
    else:
        remote_results = remote_logs = remote_cache = None

    # Keep old variable names in this scope for minimal diff in the rest of the (large) function
    # while supporting the new flags. This preserves exact behavior for existing M3 Pro callers.
    icloud_base = remote_base or args.icloud_base   # for any legacy direct uses
    icloud_results = remote_results
    icloud_logs = remote_logs
    icloud_cache = remote_cache

    # status + master log are ALWAYS local-first (robust even with --no-remote-sync)
    status_file = local_logs / "campaign_status.json"
    master_log = local_logs / "campaign.log"
    lock_dir = (args.lock_dir or Path("/private/tmp") / "flexaidds_campaign.lock").resolve()
    run_timeout_s = int(args.run_timeout_hours * 3600) if args.run_timeout_hours > 0 else None
    dock_timeout_s = int(args.dock_timeout_hours * 3600) if args.dock_timeout_hours > 0 else 0

    _warn_non_portable_remote(remote_base)

    try:
        require_local_path("repo", repo, strict=False)
        require_local_path("build", build, strict=False)
        require_local_path("binary", binary, strict=False)
        require_local_path("local_base", local_base, strict=True)
        require_local_path("local_cache", local_cache, strict=True)
        if not binary.is_file() or not os.access(binary, os.X_OK):
            raise CampaignError(f"benchmark_datasets not executable: {binary}")
        if not docking_binary.is_file() or not os.access(docking_binary, os.X_OK):
            raise CampaignError(f"FlexAID not executable: {docking_binary}")
        acquire_lock(lock_dir)

        for p in (local_cache, local_results, local_logs, icloud_results, icloud_logs):
            if p is not None:
                mkdir(p)

        runtime_data = stage_runtime_data(repo, build)

        env = os.environ.copy()
        env.update(
            {
                "FLEXAIDDS_REPO": str(repo),
                "FLEXAIDDS_BUILD": str(build),
                "FLEXAIDDS_BINARY": str(docking_binary),
                "FLEXAIDDS_BENCHMARK_DATA": str(local_cache),
                "FLEXAIDDS_RESULTS": str(local_results),
                "FLEXAIDDS_LOGS": str(local_logs),
                "HOME": str(Path("/private/tmp/flexaidds-campaign-home")),
                "OMP_NUM_THREADS": str(args.omp_threads),
                "OMP_PLACES": "cores",
                "OMP_PROC_BIND": "spread",
                "OMP_WAIT_POLICY": "passive",
            }
        )
        mkdir(Path(env["HOME"]))

        manifest = {
            "run_id": args.run_id,
            "pid": os.getpid(),
            "created_at": utc_now(),
            "repo": str(repo),
            "git_head": subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=repo, text=True).strip(),
            "git_branch": subprocess.check_output(["git", "branch", "--show-current"], cwd=repo, text=True).strip(),
            "build": str(build),
            "binary": str(binary),
            "docking_binary": str(docking_binary),
            "local_base": str(local_base),
            "local_cache": str(local_cache),
            "local_results": str(local_results),
            "icloud_results": str(icloud_results),
            "icloud_logs": str(icloud_logs),
            "datasets": args.datasets,
            "runs": args.runs,
            "workers": args.workers,
            "omp_threads": args.omp_threads,
            "ga_generations": args.ga_generations,
            "ga_population": args.ga_population,
            "temperature": args.temperature,
            "clustering": args.clustering,
            "bootstrap": args.bootstrap,
            "run_timeout_seconds": run_timeout_s,
            "dock_timeout_seconds": dock_timeout_s,
            "runtime_data": runtime_data,
        }
        write_json(local_logs / "campaign_manifest.json", manifest)
        if icloud_logs is not None:
            write_json(icloud_logs / "campaign_manifest.json", manifest)
        write_status(status_file, state="preflight", run_id=args.run_id, manifest=manifest)

        print(f"[failsafe] run_id={args.run_id}")
        print(f"[failsafe] local_base={local_base}")
        if remote_base:
            label = "iCloud" if is_icloud_path(remote_base) else "remote"
            print(f"[failsafe] {label} results={remote_results}")
        else:
            print("[failsafe] remote sync disabled (--no-remote-sync)")
        if args.seed_cache:
            print("[failsafe] seeding selected local benchmark cache from iCloud")
            seed_cache(icloud_cache, local_cache, list(args.datasets), include_smoke=not args.skip_smoke)
        else:
            mkdir(local_cache)
            print(f"[failsafe] using local APFS cache without iCloud seed: {local_cache}")

        help_rc = run_logged([str(binary), "--help"], cwd=repo, env=env, log_path=master_log, timeout_s=120)
        if help_rc not in (0, 1):
            raise CampaignError(f"benchmark_datasets --help returned {help_rc}")

        if not args.skip_smoke:
            smoke_dir = local_results / "preflight_smoke"
            if smoke_dir.exists():
                shutil.rmtree(smoke_dir)
            mkdir(smoke_dir)
            smoke_list = local_base / "smoke_targets.txt"
            smoke_list.write_text("1G9V\n")
            print("[failsafe] running 1-target smoke test")
            rc = run_logged(
                [
                    str(binary),
                    "--benchmark",
                    f"pdb_list:{smoke_list}",
                    "--threads",
                    "1",
                    "--cache",
                    str(local_cache),
                    "--output",
                    str(smoke_dir),
                    "--ga-generations",
                    "1",
                    "--ga-population",
                    "20",
                    "--omp-threads",
                    str(args.omp_threads),
                    "--job-timeout-seconds",
                    str(min(dock_timeout_s, 1800) if dock_timeout_s else 1800),
                    "--temperature",
                    str(args.temperature),
                    "--clustering",
                    args.clustering,
                ],
                cwd=repo,
                env=env,
                log_path=master_log,
                timeout_s=1800,
            )
            if rc != 0:
                raise CampaignError(f"smoke test returned {rc}")
            validate_run(
                dataset="preflight_smoke",
                run_dir=smoke_dir,
                require_nonzero_success=True,
                expected_min_systems=1,
            )
            if icloud_results is not None:
                rsync_copy(smoke_dir, icloud_results / "preflight_smoke", delete=True)
                rsync_copy(local_logs, icloud_logs, delete=False)

        if args.preflight_only:
            if icloud_logs is not None:
                rsync_copy(local_logs, icloud_logs, delete=False)
            write_status(status_file, state="preflight_complete", run_id=args.run_id, manifest=manifest)
            print("[failsafe] preflight complete")
            return 0

        completed = 0
        failed = 0
        total_runs = len(args.datasets) * args.runs
        write_status(status_file, state="running", run_id=args.run_id, completed=0, failed=0, total=total_runs)

        for dataset in args.datasets:
            for run_num in range(1, args.runs + 1):
                if STOP_REQUESTED:
                    raise CampaignError("stop requested")
                run_name = f"run{run_num:02d}"
                local_run = local_results / "tier2" / dataset / run_name
                remote_run = (icloud_results / "tier2" / dataset / run_name) if icloud_results is not None else None

                if args.resume and remote_run is not None and (remote_run / "RUN_OK.json").exists():
                    print(f"[failsafe] {dataset} {run_name}: already validated in remote")
                    completed += 1
                    continue

                if local_run.exists():
                    shutil.rmtree(local_run)
                if args.resume and remote_run is not None and remote_run.exists():
                    print(f"[failsafe] {dataset} {run_name}: restoring partial remote state")
                    rsync_copy(remote_run, local_run, delete=False)
                else:
                    mkdir(local_run)
                stale_failed = local_run / "RUN_FAILED.json"
                if stale_failed.exists():
                    stale_failed.replace(local_run / f"RUN_FAILED.interrupted_{local_now_tag()}.json")
                write_status(
                    status_file,
                    state="running",
                    run_id=args.run_id,
                    active_dataset=dataset,
                    active_run=run_name,
                    completed=completed,
                    failed=failed,
                    total=total_runs,
                )

                print(f"[failsafe] {dataset} {run_name}: starting")
                rc = run_logged(
                    [
                        str(binary),
                        "--benchmark",
                        dataset,
                        "--threads",
                        str(args.workers),
                        "--cache",
                        str(local_cache),
                        "--output",
                        str(local_run),
                        "--ga-generations",
                        str(args.ga_generations),
                        "--ga-population",
                        str(args.ga_population),
                        "--omp-threads",
                        str(args.omp_threads),
                        "--job-timeout-seconds",
                        str(dock_timeout_s),
                        "--temperature",
                        str(args.temperature),
                        "--clustering",
                        args.clustering,
                    ],
                    cwd=repo,
                    env=env,
                    log_path=master_log,
                    timeout_s=run_timeout_s,
                )
                if rc != 0:
                    failed += 1
                    write_json(local_run / "RUN_FAILED.json", {"dataset": dataset, "run": run_name, "exit_code": rc})
                    rsync_copy(local_run, icloud_run, delete=True)
                    rsync_copy(local_logs, icloud_logs, delete=False)
                    raise CampaignError(f"{dataset} {run_name}: benchmark_datasets returned {rc}")

                metrics = validate_run(
                    dataset=dataset,
                    run_dir=local_run,
                    require_nonzero_success=True,
                    expected_min_systems=EXPECTED_MIN_SYSTEMS.get(dataset),
                )
                print(
                    f"[failsafe] {dataset} {run_name}: "
                    f"{metrics.successful}/{metrics.total_systems} success "
                    f"({metrics.success_rate * 100:.2f}%)"
                )
                completed += 1
                rsync_copy(local_run, icloud_run, delete=True)
                rsync_copy(local_logs, icloud_logs, delete=False)
                if args.mirror_cache:
                    rsync_copy(local_cache, icloud_cache, delete=False)
                write_status(
                    status_file,
                    state="running",
                    run_id=args.run_id,
                    active_dataset=dataset,
                    active_run=run_name,
                    completed=completed,
                    failed=failed,
                    total=total_runs,
                    last_success_rate=metrics.success_rate,
                )

        print("[failsafe] running bootstrap analysis")
        analysis_rc = run_analysis(repo, local_results, local_analysis, args.bootstrap, master_log)
        if analysis_rc != 0:
            raise CampaignError(f"bootstrap analysis failed with exit code {analysis_rc}")
        rsync_copy(local_results, icloud_results, delete=False)
        rsync_copy(local_logs, icloud_logs, delete=False)
        write_status(status_file, state="complete", run_id=args.run_id, completed=completed, failed=failed, total=total_runs)
        print("[failsafe] complete")
        return 0
    except Exception as exc:
        try:
            write_status(status_file, state="failed", run_id=args.run_id, error=str(exc))
            if local_results.exists():
                rsync_copy(local_results, icloud_results, delete=False)
            if local_logs.exists():
                rsync_copy(local_logs, icloud_logs, delete=False)
        except Exception:
            pass
        print(f"[failsafe] ERROR: {exc}", file=sys.stderr)
        return 1
    finally:
        try:
            lock_dir.rmdir()
        except OSError:
            pass


if __name__ == "__main__":
    raise SystemExit(main())
