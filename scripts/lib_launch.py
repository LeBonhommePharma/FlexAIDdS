#!/usr/bin/env python3
# lib_launch.py — shared helpers for FlexAIDdS benchmark launch scripts
#
# Import from any launch_vN.py:
#   from lib_launch import priority_from_prev_run, priority_for_fix
#   from lib_launch import dock_preflight  # Sol #9 fail-closed gate
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

import csv
import glob
import json
import os
import signal
import subprocess
import sys

# Fail-closed multi-session dock gate (benchmark_coord.py) — hard-required.
# Soft-missing the module would let every launch_v*.py skip Sol #9.
_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
if _SCRIPTS not in sys.path:
    sys.path.insert(0, _SCRIPTS)
try:
    from benchmark_coord import (  # noqa: E402
        MAX_WORKERS,
        preflight as dock_preflight,
        release_lock as dock_release_lock,
        status as dock_status,
    )
except ImportError as _coord_imp_err:  # pragma: no cover
    raise ImportError(
        "Sol #9 dock coordination requires scripts/benchmark_coord.py — "
        "fail-closed (cannot soft-disable). Fix the import or PYTHONPATH."
    ) from _coord_imp_err


def priority_from_prev_run(prev_result_dir, lo=1.8, hi=2.5):
    """Return comma-separated PDB IDs with best_cluster_rmsd in (lo, hi) from a
    previous run — i.e., near-misses that almost made it under the 2 Å cutoff.

    Pass the result to FLEXAIDDS_PRIORITY_TARGETS so those targets run first in
    the next campaign and get fresh eyes early, before the worker pool disperses
    onto easier wins.

    Args:
        prev_result_dir: root output dir of a prior run (e.g.
            ~/flexaidds_results/v41_20260613_zshannonselect).
        lo: RMSD lower bound (exclusive). Default 1.8 Å.
        hi: RMSD upper bound (exclusive). Default 2.5 Å.

    Returns:
        Comma-separated string of PDB IDs, e.g. "1HP0,1Q4G,1Q41".
        Empty string if no near-misses found or dir doesn't exist.

    Example::
        prio = priority_from_prev_run("~/flexaidds_results/v41_...", lo=1.8, hi=2.5)
        if prio:
            env["FLEXAIDDS_PRIORITY_TARGETS"] = prio
    """
    prev_result_dir = os.path.expanduser(prev_result_dir)
    if not os.path.isdir(prev_result_dir):
        return ""

    near_misses = []
    for csv_path in glob.glob(f"{prev_result_dir}/*/result.csv"):
        pdb_id = os.path.basename(os.path.dirname(csv_path))
        try:
            with open(csv_path, newline="") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    try:
                        # rmsd_hungarian = top-1 selected pose (benchmark metric).
                        # best_cluster_rmsd = oracle minimum (not what we report).
                        # Near-miss detection should use the benchmark metric.
                        rmsd = float(
                            row.get("rmsd_hungarian")
                            or row.get("best_cluster_rmsd")
                            or 999
                        )
                        if lo < rmsd < hi:
                            near_misses.append(pdb_id)
                            break  # one hit per target is enough
                    except (ValueError, TypeError):
                        pass
        except OSError:
            pass

    return ",".join(near_misses)


def priority_for_fix(result_dir, fix_types, max_targets=10):
    """Return comma-separated PDB IDs whose failure mode matches any of fix_types.

    Reads the failure_modes.json produced by failure_classify.py.  If the JSON
    does not exist yet, runs the classifier on-the-fly (requires failure_classify
    to be importable from the same scripts/ directory).

    Args:
        result_dir: root output dir of a completed benchmark run.
        fix_types:  list of failure-mode strings the next fix addresses, e.g.
                    ['selection_miss', 'CF_false_minimum'].
                    Prefix matching is supported: 'search_failure' matches
                    'search_failure', 'search_failure:rigid', etc.
        max_targets: maximum number of PDB IDs to return (default: 10).

    Returns:
        Comma-separated string of PDB IDs ordered by ascending RMSD
        (closest to the 2 Å success cutoff first — hardest near-misses first).
        Empty string if no matches or dir doesn't exist.

    Column-semantics note:
        This function reads RMSD from rmsd_hungarian (top-1 selected pose),
        which is the benchmark metric.  See failure_classify.py docstring.

    Example::
        prio = priority_for_fix(
            "~/flexaidds_results/v40_...",
            fix_types=["selection_miss", "CF_false_minimum"],
        )
        if prio:
            env["FLEXAIDDS_PRIORITY_TARGETS"] = prio
    """
    result_dir = os.path.expanduser(result_dir)
    if not os.path.isdir(result_dir):
        return ""

    json_path = os.path.join(result_dir, "failure_modes.json")

    # Auto-generate the JSON if it doesn't exist yet
    if not os.path.isfile(json_path):
        try:
            import sys
            scripts_dir = os.path.dirname(os.path.abspath(__file__))
            if scripts_dir not in sys.path:
                sys.path.insert(0, scripts_dir)
            from failure_classify import classify_run  # noqa: PLC0415
            classifications, _ = classify_run(result_dir)
            with open(json_path, "w") as fh:
                json.dump(classifications, fh, indent=2)
                fh.write("\n")
        except Exception:
            # If classify fails (e.g. no result CSVs yet), return empty
            return ""

    try:
        with open(json_path) as fh:
            classifications = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return ""

    def _mode_matches(mode, fix_types):
        if mode is None:
            return False
        for ft in fix_types:
            if mode == ft or mode.startswith(ft + ":") or ft.startswith(mode + ":"):
                return True
        return False

    matched = []
    for pdb_id, entry in classifications.items():
        if entry.get("success"):
            continue
        if _mode_matches(entry.get("failure_mode"), fix_types):
            matched.append((entry.get("rmsd", 999.0), pdb_id))

    # Sort ascending by RMSD (nearest to success threshold first)
    matched.sort(key=lambda x: x[0])
    top = [pid for _, pid in matched[:max_targets]]
    return ",".join(top)


def _infer_workers(cmd, env, workers):
    """Resolve worker count for Sol #9 preflight (default 1 if unknown)."""
    if workers is not None:
        return int(workers)
    if env:
        for key in ("FLEXAIDDS_BENCH_THREADS", "BENCH_THREADS"):
            raw = env.get(key)
            if raw:
                try:
                    return int(raw)
                except (TypeError, ValueError):
                    pass
    if cmd:
        for i, arg in enumerate(cmd):
            if arg in ("--threads", "--workers") and i + 1 < len(cmd):
                try:
                    return int(cmd[i + 1])
                except (TypeError, ValueError):
                    pass
    return 1


def _infer_binary(cmd, env, binary):
    """Resolve engine binary path for stamp + preflight."""
    if binary:
        return binary
    if env and env.get("FLEXAIDDS_BINARY"):
        return env["FLEXAIDDS_BINARY"]
    if cmd:
        for arg in cmd:
            if not arg or arg.startswith("-"):
                continue
            # skip known wrappers
            base = os.path.basename(arg)
            if base in ("caffeinate", "env", "nice", "nohup", "time"):
                continue
            if os.path.isfile(arg) and os.access(arg, os.X_OK):
                return arg
    raise RuntimeError(
        "launch_session_isolated: cannot resolve FLEXAIDDS_BINARY for Sol #9 "
        "preflight (pass binary= or set env FLEXAIDDS_BINARY)"
    )


def launch_session_isolated(
    cmd,
    env,
    output_dir,
    *,
    stdout_log=None,
    stderr_log=None,
    cwd=None,
    workers=None,
    binary=None,
    owner=None,
    skip_dock_preflight=False,
    purpose="dock",
):
    """Launch *cmd* as a fully detached daemon — SIGHUP-immune, isolated session.

    Implements a POSIX double-fork so the grandchild process outlives the
    calling terminal.  The grandchild inherits *env*, runs in *cwd* (defaults
    to *output_dir*), and logs to *stdout_log* / *stderr_log* (defaults:
    ``{output_dir}/stdout.log`` and ``{output_dir}/stderr.log``).

    **Sol #9:** Before any fork, runs fail-closed ``dock_preflight`` (hold,
    mkdir lock, disk floor, WORKERS≤4, binary stamp into *output_dir*). On
    success, ``env['FLEXAIDDS_BINARY']`` is rewritten to the stamped path when
    a stamp is produced. Lock token is written to
    ``{output_dir}/BENCHMARK_DOCK_LOCK_TOKEN`` for later ``dock_release_lock``.
    Pass ``skip_dock_preflight=True`` only for unit tests that mock the gate.

    Fix B: subprocess.Popen is called with ``start_new_session=True`` so the
    worker is placed in its own process group *and* POSIX session — it cannot
    receive SIGHUP from a terminal close regardless of setsid() in the parent.

    The PID of the spawned subprocess is written to
    ``{output_dir}/benchmark.pid`` and returned to the caller so it can be
    embedded in a provenance JSON.

    Args:
        cmd:         Command list passed to subprocess.Popen.
        env:         Environment dict for the subprocess.
        output_dir:  Root result directory (created if absent).
        stdout_log:  Path for stdout capture (default: output_dir/stdout.log).
        stderr_log:  Path for stderr capture (default: output_dir/stderr.log).
        cwd:         Working directory for the daemon (default: output_dir).
        workers:     Optional explicit worker count (else inferred).
        binary:      Optional engine binary path (else FLEXAIDDS_BINARY / cmd).
        owner:       Lock owner label (default: launch_session_isolated).
        skip_dock_preflight: Emergency/test-only bypass (default False).
        purpose:     Lock purpose metadata.

    Returns:
        int — PID of the spawned benchmark_datasets subprocess.

    Raises:
        RuntimeError if Sol #9 preflight refuses or the pipe IPC handshake fails.
    """
    from pathlib import Path

    output_dir = os.path.expanduser(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # ── Sol #9 fail-closed gate (before any fork / dock child) ──────────────
    if not skip_dock_preflight:
        if dock_preflight is None:  # pragma: no cover — import hard-fails above
            raise RuntimeError(
                "launch_session_isolated: dock_preflight unavailable (Sol #9)"
            )
        n_workers = _infer_workers(cmd, env, workers)
        bin_path = _infer_binary(cmd, env, binary)
        pf = dock_preflight(
            out=Path(output_dir),
            workers=n_workers,
            binary=Path(bin_path),
            owner=owner or "launch_session_isolated",
            purpose=purpose,
            acquire=True,
            stamp=True,
        )
        if not pf.ok:
            raise RuntimeError(
                f"Sol #9 dock preflight refused: {pf.reason}"
            )
        # Persist token for orderly release; point engine at stamped binary.
        token_path = os.path.join(output_dir, "BENCHMARK_DOCK_LOCK_TOKEN")
        if pf.lock_token:
            with open(token_path, "w", encoding="utf-8") as tf:
                tf.write(pf.lock_token + "\n")
        if env is not None and pf.stamped_binary:
            env = dict(env)
            env["FLEXAIDDS_BINARY"] = pf.stamped_binary
        receipt = {
            "ok": True,
            "reason": pf.reason,
            "lock_token": pf.lock_token,
            "stamped_binary": pf.stamped_binary,
            "workers": n_workers,
            "owner": owner or "launch_session_isolated",
        }
        with open(
            os.path.join(output_dir, "DOCK_PREFLIGHT.json"),
            "w",
            encoding="utf-8",
        ) as rf:
            json.dump(receipt, rf, indent=2)
            rf.write("\n")

    pid_file   = os.path.join(output_dir, "benchmark.pid")
    stdout_log = stdout_log or os.path.join(output_dir, "stdout.log")
    stderr_log = stderr_log or os.path.join(output_dir, "stderr.log")
    cwd        = cwd or output_dir

    # Pipe carries the grandchild's subprocess PID back to the grandparent.
    r_fd, w_fd = os.pipe()

    if os.fork() > 0:
        # ── Grandparent: read PID from pipe, write pid_file, return ────────
        os.close(w_fd)
        data = b""
        while True:
            chunk = os.read(r_fd, 64)
            if not chunk:
                break
            data += chunk
        os.close(r_fd)
        if not data.strip():
            raise RuntimeError(
                f"launch_session_isolated: no PID received from daemon — "
                f"check {stdout_log} and {stderr_log}"
            )
        child_pid = int(data.strip())
        with open(pid_file, "w") as f:
            f.write(str(child_pid) + "\n")
        return child_pid

    # ── First child ─────────────────────────────────────────────────────────
    os.close(r_fd)
    os.setsid()                            # new session, detach from terminal

    if os.fork() > 0:
        os._exit(0)                        # first child exits; grandchild lives

    # ── Grandchild (daemon) — fully detached ────────────────────────────────
    signal.signal(signal.SIGHUP, signal.SIG_IGN)   # belt-and-suspenders
    os.chdir(cwd)

    with open(stdout_log, "w") as out, open(stderr_log, "w") as err:
        p = subprocess.Popen(
            cmd,
            stdout=out,
            stderr=err,
            env=env,
            start_new_session=True,        # Fix B: own session → immune to SIGHUP
        )

    # Hand PID back to grandparent via the pipe, then close.
    os.write(w_fd, str(p.pid).encode())
    os.close(w_fd)
    os._exit(0)
