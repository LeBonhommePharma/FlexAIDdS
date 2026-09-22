#!/usr/bin/env python3
"""Start a command in its own session, fully detached. A portable setsid(1).

WHY THIS EXISTS -- measured 2026-09-21
    Three campaign relaunches were dispatched with `nohup setsid bash ...` and all
    three died instantly, leaving 41-byte logs reading

        nohup: setsid: No such file or directory

    setsid(1) is a util-linux program and does not ship on macOS. The syscall it
    wraps, setsid(2), is POSIX and is exposed by Python as os.setsid(), so the
    capability was always available -- only the binary was missing.

    The reason a plain `nohup cmd &` was not enough: nohup only arranges for
    SIGHUP to be ignored. The child remains in the parent's session and process
    group, so anything that signals the GROUP -- which is what happened when this
    session's daemon restarted -- takes the child with it. Detaching requires a
    new session, which requires setsid(2) from a process that is not a process
    group leader. Hence the double fork.

    launchd would also solve it and would additionally survive reboot, but writing
    a LaunchAgent is out of scope here: that path is protected.

WHAT IT GUARANTEES
    * the command runs in a NEW session with no controlling terminal, so neither
      the parent's death nor a signal to the parent's process group reaches it
    * stdin is /dev/null and stdout/stderr are redirected to the log path given,
      opened in append mode so a restart cannot truncate a previous run's log
    * the child's pid is written to <log>.pid BEFORE exec, so "did it start" can
      be answered from disk rather than inferred
    * exit status reports whether the SPAWN succeeded, which is not the same
      thing as the command succeeding -- verify the log content separately
      (scripts/verify_launch.sh exists for that and encodes the known failure
      signatures)

WHAT IT DOES NOT SURVIVE
    A reboot. Nothing short of a launchd job does. The mitigation for that is a
    working resume gate, which this project now has: a relaunch after the reboot
    of 2026-09-21 skipped 83 of 85 completed cells instead of re-docking them.

Usage:
    spawn_detached.py --log <path> [--cwd <dir>] [--env K=V ...] -- <cmd> [args...]
"""

from __future__ import annotations

import argparse
import os
import sys


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--log", required=True, help="stdout+stderr destination (append mode)")
    ap.add_argument("--cwd")
    ap.add_argument("--env", action="append", default=[], metavar="K=V")
    ap.add_argument("cmd", nargs=argparse.REMAINDER)
    a = ap.parse_args()

    cmd = a.cmd[1:] if a.cmd and a.cmd[0] == "--" else a.cmd
    if not cmd:
        print("no command given", file=sys.stderr)
        return 2

    # Fail BEFORE forking if the command does not exist. This is the exact defect
    # being fixed: `nohup setsid ...` reported success from the shell's point of
    # view and left the failure buried in a log nobody read for hours.
    exe = cmd[0]
    if "/" in exe:
        if not (os.path.isfile(exe) and os.access(exe, os.X_OK)):
            print(f"not executable: {exe}", file=sys.stderr)
            return 127
    else:
        from shutil import which
        if which(exe) is None:
            print(f"not on PATH: {exe}", file=sys.stderr)
            return 127

    logpath = os.path.abspath(a.log)
    os.makedirs(os.path.dirname(logpath), exist_ok=True)
    pidpath = logpath + ".pid"

    env = dict(os.environ)
    for kv in a.env:
        k, _, v = kv.partition("=")
        env[k] = v

    # First fork: the parent returns to the caller immediately. The child is not a
    # process-group leader, which is the precondition setsid(2) requires.
    r, w = os.pipe()
    pid = os.fork()
    if pid > 0:
        os.close(w)
        got = os.read(r, 64).decode().strip()
        os.close(r)
        os.waitpid(pid, 0)              # reap the intermediate; it exits at once
        if not got.isdigit():
            print(f"spawn failed: {got or 'no pid reported'}", file=sys.stderr)
            return 1
        print(got)
        return 0

    os.close(r)
    try:
        os.setsid()                     # NEW SESSION: the whole point
        # Second fork: the session leader exits, so the final process can never
        # reacquire a controlling terminal.
        pid2 = os.fork()
        if pid2 > 0:
            os.write(w, str(pid2).encode())
            os.close(w)
            os._exit(0)

        os.close(w)
        if a.cwd:
            os.chdir(a.cwd)
        fd_in = os.open(os.devnull, os.O_RDONLY)
        fd_out = os.open(logpath, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o644)
        os.dup2(fd_in, 0)
        os.dup2(fd_out, 1)
        os.dup2(fd_out, 2)
        with open(pidpath, "w") as fh:
            fh.write(f"{os.getpid()}\n")
        os.execvpe(cmd[0], cmd, env)
    except Exception as exc:            # noqa: BLE001 - must reach the parent
        try:
            os.write(w, f"error: {exc}".encode())
            os.close(w)
        except OSError:
            pass
        os._exit(127)
    os._exit(127)


if __name__ == "__main__":
    sys.exit(main())
