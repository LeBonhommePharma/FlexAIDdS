#!/usr/bin/env python3
"""v130 live progress watcher — run anytime: python3 v130_watch.py

Scans full_v130/ for completed result.csv files, bins by RMSD tier,
shows currently-running target from the log tail, and estimates ETA.
Usage: python3 v130_watch.py [--watch N]   (--watch N: refresh every N sec)
"""

import csv, glob, os, sys, time, signal
from datetime import datetime, timedelta

BASE  = "/Users/lp.more/Projects/FlexAIDdS/benchmarks/astex_repro"
DIR   = f"{BASE}/full_v130"
LOG   = f"{BASE}/v130.log"
LOG2  = f"{BASE}/v130_1hq2.log"
TOTAL = 85   # 84 + 1HQ2 queued separately

V130_PID  = 98837
QUEUE_PID = 6618


def pid_alive(pid):
    try: os.kill(pid, 0); return True
    except OSError: return False

def read_results():
    ok, nm, fail, sent = [], [], [], []
    for r in sorted(glob.glob(f"{DIR}/*/result.csv")):
        pdb = r.split(os.sep)[-2]
        try:
            with open(r) as f:
                row = next(csv.DictReader(f))
            # column may be rmsd_to_crystal or rmsd
            rmsd = float(row.get("rmsd_to_crystal") or row.get("rmsd") or -1)
            wt   = float(row.get("wall_time_s") or 0)
        except Exception:
            sent.append((pdb, -1, 0)); continue
        if rmsd < 0:        sent.append((pdb, rmsd, wt))
        elif rmsd < 2.0:    ok.append((pdb, rmsd, wt))
        elif rmsd < 2.5:    nm.append((pdb, rmsd, wt))
        else:               fail.append((pdb, rmsd, wt))
    return ok, nm, fail, sent

def current_target():
    """Read last few KB of log to find the active PDB code."""
    for logfile in [LOG, LOG2]:
        if not os.path.exists(logfile): continue
        try:
            with open(logfile, "rb") as f:
                f.seek(0, 2); size = f.tell()
                f.seek(max(0, size - 8192))
                tail = f.read().decode(errors="replace")
            for line in reversed(tail.splitlines()):
                if "EVAL-SCALE" in line or "EVAL-BUDGET" in line:
                    # e.g. [EVAL-SCALE] 1G9V: ...
                    parts = line.split("]")
                    if len(parts) > 1:
                        code = parts[1].strip().split(":")[0].strip()
                        if len(code) == 4: return code, logfile
        except Exception:
            pass
    return "?", LOG

def report():
    ok, nm, fail, sent = read_results()
    done = len(ok) + len(nm) + len(fail) + len(sent)
    cur, _ = current_target()

    v130_up  = pid_alive(V130_PID)
    queue_up = pid_alive(QUEUE_PID)

    avg_wt = None
    all_done = ok + nm + fail + sent
    times = [wt for _, _, wt in all_done if wt > 0]
    if times:
        avg_wt = sum(times) / len(times)
        remaining = TOTAL - done
        eta_s = avg_wt * remaining
        eta_str = str(timedelta(seconds=int(eta_s)))
    else:
        eta_str = "n/a"

    now = datetime.now().strftime("%H:%M:%S")
    print(f"\n{'='*58}")
    print(f"  v130 progress  [{now}]")
    print(f"{'='*58}")
    print(f"  PIDs   v130={V130_PID} {'✓' if v130_up else '✗ DEAD'}   "
          f"1HQ2-queue={QUEUE_PID} {'✓' if queue_up else '✗ DEAD'}")
    print(f"  Done   {done}/{TOTAL}   running: {cur}")
    print(f"  ✓ success  (<2.0Å): {len(ok):3d}  —  "
          + ", ".join(p for p, *_ in ok) if ok else f"  ✓ success  (<2.0Å):   0")
    if nm:
        print(f"  ~ near-miss (2-2.5): {len(nm):2d}  —  "
              + ", ".join(f"{p}({r:.2f})" for p, r, _ in nm))
    print(f"  ✗ fail   (≥2.5Å): {len(fail):3d}")
    if fail:
        for p, r, _ in sorted(fail, key=lambda x: x[1]):
            print(f"       {p}  RMSD={r:.3f}Å")
    if sent:
        print(f"  ⚠ sentinel (no pose): {len(sent):2d}  — "
              + ", ".join(p for p, *_ in sent))
    if avg_wt:
        print(f"  avg wall/target: {avg_wt/60:.1f} min   ETA: ~{eta_str}")
    print(f"{'='*58}\n")

def main():
    interval = None
    args = sys.argv[1:]
    if "--watch" in args:
        idx = args.index("--watch")
        try: interval = int(args[idx + 1])
        except (IndexError, ValueError): interval = 60

    if interval:
        signal.signal(signal.SIGINT, lambda *_: sys.exit(0))
        while True:
            os.system("clear")
            report()
            time.sleep(interval)
    else:
        report()

if __name__ == "__main__":
    main()
