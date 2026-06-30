"""Launch queue planning for large-N benchmark campaigns (astex_nonnative + posex_cd)."""

from __future__ import annotations

import json
import logging
import os
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

LARGE_N_ENTRIES = {
    "astex_nonnative": 1113,
    "posex_cd": 1312,
}


def _repo_root(explicit: Optional[Path] = None) -> Path:
    if explicit is not None:
        return explicit
    return Path(__file__).resolve().parents[3]


def build_launch_plan(repo_root: Optional[Path] = None) -> Dict[str, Any]:
    """Return full-N launch commands for astex_nonnative (1113) and posex_cd (1312)."""
    root = _repo_root(repo_root)
    launcher = root / ".grok/skills/flexaid-docking/scripts/launch_full_benchmark.sh"
    posex_json = root / "benchmarks/datasets/posex_cd_1312.json"
    return {
        "astex_nonnative": {
            "n_entries": LARGE_N_ENTRIES["astex_nonnative"],
            "launcher": str(launcher),
            "command": [
                "bash", str(launcher), "astex_nonnative", "298", "astex_nonnative_298K_vht",
            ],
        },
        "posex_cd": {
            "n_pairs": LARGE_N_ENTRIES["posex_cd"],
            "posex_json": str(posex_json),
            "command": [
                "benchmark_datasets",
                f"--benchmark=crossdock_json:{posex_json}",
                "--threads", "8",
            ],
        },
    }


def build_run_status(
    sibling_count: int,
    *,
    repo_root: Optional[Path] = None,
    scratch: Optional[Path] = None,
) -> Dict[str, Any]:
    """Build run_status.json payload for conditioned or post-wait launch."""
    root = _repo_root(repo_root)
    plan = build_launch_plan(root)
    checked_at = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

    if sibling_count > 0:
        return {
            "status": "waiting_for_astex_diverse_siblings",
            "sibling_process_count": sibling_count,
            "datasets_pending": ["astex_nonnative_1113", "posex_cd_1312"],
            "checked_at": checked_at,
            "next_action": "benchmarks/m3pro/queue_large_benchmarks.sh",
            "launch_plan": plan,
        }

    return {
        "status": "launched_full",
        "astex_nonnative": plan["astex_nonnative"],
        "posex_cd": plan["posex_cd"],
        "launched_at": checked_at,
        "launch_plan": plan,
    }


def count_astex_diverse_siblings() -> int:
    """Best-effort count of active Astex Diverse / astex benchmark processes."""
    try:
        result = subprocess.run(
            ["pgrep", "-fl", "astex_diverse|benchmark_datasets.*astex"],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            return 0
        return len([ln for ln in result.stdout.strip().splitlines() if ln.strip()])
    except Exception:
        return 0


def astex_nonnative_running() -> bool:
    """True when an astex_nonnative full benchmark is already in flight."""
    try:
        result = subprocess.run(
            ["pgrep", "-fl", "astex_nonnative_298K|benchmark astex_nonnative"],
            capture_output=True,
            text=True,
            check=False,
        )
        return result.returncode == 0 and bool(result.stdout.strip())
    except Exception:
        return False


def wait_for_astex_nonnative_clear(max_wait_s: int = 0, interval_s: int = 120) -> bool:
    """Poll until no astex_nonnative benchmark process remains. Returns True when clear."""
    if max_wait_s <= 0:
        return not astex_nonnative_running()
    elapsed = 0
    while elapsed < max_wait_s:
        if not astex_nonnative_running():
            return True
        time.sleep(interval_s)
        elapsed += interval_s
    return not astex_nonnative_running()


def write_run_status(path: Path, status: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(status, indent=2))


def _benchmark_datasets_available() -> bool:
    return shutil.which("benchmark_datasets") is not None


def execute_launches(
    repo_root: Optional[Path] = None,
    scratch: Optional[Path] = None,
    *,
    dry_run: bool = False,
) -> Dict[str, Any]:
    """Execute full-N launch commands; return metadata for run_status.json."""
    root = _repo_root(repo_root)
    scratch = scratch or Path(os.environ.get("SCRATCH", "/tmp"))
    plan = build_launch_plan(root)
    results: Dict[str, Any] = {"astex_nonnative": {}, "posex_cd": {}}

    skip_nn = os.environ.get("SKIP_ASTEX_NONNATIVE", "").strip() in {"1", "true", "yes"}
    nn_cmd = plan["astex_nonnative"]["command"]
    launcher = Path(nn_cmd[1]) if len(nn_cmd) > 1 else None
    if dry_run:
        results["astex_nonnative"] = {"command": nn_cmd, "dry_run": True}
    elif skip_nn or astex_nonnative_running():
        results["astex_nonnative"] = {
            "command": nn_cmd,
            "skipped": "astex_nonnative already running or SKIP_ASTEX_NONNATIVE set",
        }
    elif launcher is None or not launcher.is_file():
        results["astex_nonnative"] = {
            "command": nn_cmd,
            "skipped": f"launcher missing: {launcher}",
        }
    else:
        try:
            proc = subprocess.run(nn_cmd, cwd=str(root), check=False)
            results["astex_nonnative"] = {"command": nn_cmd, "returncode": proc.returncode}
        except Exception as exc:
            logger.warning("astex_nonnative launch failed: %s", exc)
            results["astex_nonnative"] = {"command": nn_cmd, "error": str(exc)}

    posex_json = Path(plan["posex_cd"]["posex_json"])
    posex_out = scratch / f"posex_cd_298K_{int(time.time())}"
    posex_out.mkdir(parents=True, exist_ok=True)
    results["posex_cd"]["output_dir"] = str(posex_out)

    if dry_run:
        results["posex_cd"]["command"] = plan["posex_cd"]["command"]
        results["posex_cd"]["dry_run"] = True
    elif not posex_json.is_file():
        results["posex_cd"]["skipped"] = f"missing posex json: {posex_json}"
    elif not _benchmark_datasets_available():
        results["posex_cd"]["skipped"] = "benchmark_datasets not on PATH"
    else:
        posex_cmd = [
            "benchmark_datasets",
            f"--benchmark=crossdock_json:{posex_json}",
            "--output", str(posex_out),
            "--threads", "8",
        ]
        log_out = posex_out / "binary.log"
        err_out = posex_out / "stderr.log"
        try:
            with open(log_out, "a") as log_f, open(err_out, "a") as err_f:
                proc = subprocess.Popen(
                    posex_cmd,
                    cwd=str(root),
                    stdout=log_f,
                    stderr=err_f,
                    stdin=subprocess.DEVNULL,
                )
            results["posex_cd"].update({"command": posex_cmd, "pid": proc.pid})
            write_run_status(posex_out / "run_status.json", {
                "status": "launched",
                "dataset": "crossdock_json:posex_cd_1312",
                "n_pairs": LARGE_N_ENTRIES["posex_cd"],
                "output_dir": str(posex_out),
                "pid": proc.pid,
            })
        except Exception as exc:
            logger.warning("posex_cd launch failed: %s", exc)
            results["posex_cd"].update({"command": posex_cmd, "error": str(exc)})

    return results


def main(argv: Optional[List[str]] = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description="Large-N benchmark launch queue planner")
    parser.add_argument("--scratch", default=os.environ.get("SCRATCH", "/tmp"))
    parser.add_argument("--repo-root", default=None)
    parser.add_argument("--sibling-count", type=int, default=None)
    parser.add_argument("--write-status", type=Path, default=None)
    parser.add_argument("--execute", action="store_true", help="Run launch commands after planning")
    parser.add_argument("--dry-run", action="store_true", help="Plan only, do not invoke binaries")
    args = parser.parse_args(argv)

    repo_root = Path(args.repo_root) if args.repo_root else None
    scratch = Path(args.scratch)
    sibling_count = (
        args.sibling_count
        if args.sibling_count is not None
        else count_astex_diverse_siblings()
    )

    status = build_run_status(sibling_count, repo_root=repo_root, scratch=scratch)

    # Write status before execution so verification artifacts survive launch failures.
    if args.write_status:
        write_run_status(args.write_status, status)

    if args.execute and status["status"] == "launched_full":
        try:
            launch_meta = execute_launches(repo_root, scratch, dry_run=args.dry_run)
        except Exception as exc:
            logger.warning("execute_launches failed: %s", exc)
            launch_meta = {"error": str(exc)}
        status["launch_results"] = launch_meta
        if args.write_status:
            write_run_status(args.write_status, status)

    print(json.dumps(status, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())