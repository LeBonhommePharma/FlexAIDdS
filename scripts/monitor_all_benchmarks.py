#!/usr/bin/env python3
"""monitor_all_benchmarks.py — production multi-campaign benchmark monitor.

Discovers concurrent docking / orchestrator processes and reports progress for
known FlexAIDdS campaign namespaces (oracle-ceiling, C0 full85, FlexAID A/B0/B).

Never kills processes. Never dual-launches campaigns.

Usage:
  python3 scripts/monitor_all_benchmarks.py
  python3 scripts/monitor_all_benchmarks.py --json-out /tmp/bench_status.json
  python3 scripts/monitor_all_benchmarks.py --scratch "$SCRATCH"

Environment (optional):
  FLEXAIDDS_ROOT              repo root (default: git toplevel or cwd)
  FLEXAIDDS_ICLOUD            iCloud benchmarks root
  FLEXAIDDS_RESULTS           results root (default: $FLEXAIDDS_ICLOUD/results)
  FLEXAIDDS_QUEUE_ROOT        three-engine queue root
  FLEXAIDDS_ORACLE_CAMPAIGN   oracle-ceiling campaign dir
  FLEXAIDDS_MONITOR_SCRATCH   default snapshot directory

Copyright 2026 Le Bonhomme Pharma
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Optional


# ─── path resolution ─────────────────────────────────────────────────────────

def repo_root() -> Path:
    env = os.environ.get("FLEXAIDDS_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
        return Path(out)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return Path.cwd().resolve()


def _first_existing(candidates: Iterable[Path]) -> Optional[Path]:
    for p in candidates:
        try:
            if p.is_dir():
                return p.resolve()
        except OSError:
            continue
    return None


def icloud_root() -> Path:
    env = os.environ.get("FLEXAIDDS_ICLOUD")
    if env:
        p = Path(env).expanduser()
        if p.is_dir():
            return p.resolve()
    home = Path.home()
    found = _first_existing(
        [
            home
            / "Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks",
            home / "Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS",
        ]
    )
    if found:
        return found
    return (
        home
        / "Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks"
    )


def results_root() -> Path:
    env = os.environ.get("FLEXAIDDS_RESULTS")
    if env:
        p = Path(env).expanduser()
        # Prefer env only if campaigns/ (or parent campaigns) look real
        if p.is_dir() and (
            (p / "campaigns").is_dir()
            or any(p.glob("**/C0_full85*"))
            or any(p.glob("**/three_engine/**"))
        ):
            return p.resolve()
    ic = icloud_root()
    for cand in (ic / "results", ic):
        if (cand / "campaigns").is_dir() or cand.is_dir():
            if (cand / "campaigns").is_dir():
                return cand.resolve()
    return (ic / "results").resolve()


def queue_root() -> Path:
    env = os.environ.get("FLEXAIDDS_QUEUE_ROOT")
    if env:
        p = Path(env).expanduser()
        if p.is_dir():
            return p.resolve()
    ic = icloud_root()
    found = _first_existing(
        [
            ic / "queues/three_engine_entropy_q1",
            ic / "queues/three_engine_entropy_q1".replace("q1", "q1"),
        ]
    )
    return found or (ic / "queues/three_engine_entropy_q1")


def oracle_campaign() -> Path:
    env = os.environ.get("FLEXAIDDS_ORACLE_CAMPAIGN")
    if env:
        return Path(env).expanduser().resolve()
    return Path.home() / "flexaidds_results/oracle_ceiling_restore_v43proto_r3"


# ─── system helpers ───────────────────────────────────────────────────────────

def mem_snapshot() -> dict[str, float]:
    try:
        ps = int(subprocess.check_output(["pagesize"], text=True).strip())
    except (subprocess.CalledProcessError, FileNotFoundError):
        ps = 4096
    try:
        raw = subprocess.check_output(["vm_stat"], text=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return {}
    d: dict[str, int] = {}
    for line in raw.splitlines()[1:]:
        if ":" not in line:
            continue
        k, v = line.split(":", 1)
        try:
            d[k.strip()] = int(v.strip().rstrip("."))
        except ValueError:
            pass

    def gb(pages: float) -> float:
        return pages * ps / 1e9

    free = d.get("Pages free", 0)
    pur = d.get("Pages purgeable", 0)
    ina = d.get("Pages inactive", 0)
    return {
        "free_GB": gb(free),
        "purgeable_GB": gb(pur),
        "inactive_GB": gb(ina),
        "available_est_GB": gb(free + pur + 0.3 * ina),
    }


PROCESS_PATTERNS = [
    re.compile(r"bin/[ABC]/FlexAID(?:dS)?\b"),
    re.compile(r"/FlexAIDdS(\s|$)"),
    re.compile(r"/FlexAID(\s|$)"),
    re.compile(r"benchmark_datasets\b"),
    re.compile(r"serial_AB_after_C0"),
    re.compile(r"throughput_maximizer"),
    re.compile(r"run_flexaid_arm_pilot8"),
    re.compile(r"run_[AB0]+_pilot8\.sh"),
    re.compile(r"ram_guard_resume_pilot"),
    re.compile(r"generate_flexaid_inp\.py"),
    re.compile(r"--mode\s+oracle-ceiling"),
    re.compile(r"C0_full85"),
    re.compile(r"defined-cleft-redock"),
]

# Ignore search/noise PIDs that only *mention* campaign paths
NOISE_CMD = re.compile(
    r"^(rg|grep|ugrep|find|bfs|mdfind|spotlight)\b|"
    r"\brg\s+-|"
    r"monitor_all_benchmarks"
)


@dataclass
class ProcInfo:
    pid: int
    rss_mb: float
    pcpu: float
    state: str
    etime: str
    command: str
    campaign_hint: str = ""


def classify_command(cmd: str) -> str:
    if "oracle_ceiling" in cmd or "oracle-ceiling" in cmd:
        return "oracle_ceiling"
    if "C0_full85" in cmd or "defined_cleft_nativeseed_forbidden" in cmd:
        return "C0_full85"
    if "three_engine/A" in cmd or re.search(r"\barm[= ]A\b|/bin/A/FlexAID", cmd):
        return "arm_A"
    if "three_engine/B0" in cmd or "/bin/B/FlexAID" in cmd and "B0" in cmd:
        return "arm_B0"
    if "three_engine/B/" in cmd or re.search(r"run_B_pilot8|arm[= ]B\b", cmd):
        return "arm_B"
    if "throughput_maximizer" in cmd or "serial_AB" in cmd:
        return "ab_watcher"
    if "benchmark_datasets" in cmd:
        return "benchmark_datasets"
    if "FlexAIDdS" in cmd:
        return "FlexAIDdS"
    if "FlexAID" in cmd:
        return "FlexAID"
    return "other"


def list_benchmark_processes() -> list[ProcInfo]:
    try:
        raw = subprocess.check_output(
            ["ps", "-axo", "pid=,rss=,pcpu=,state=,etime=,command="],
            text=True,
            stderr=subprocess.DEVNULL,
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        return []
    out: list[ProcInfo] = []
    self_markers = ("monitor_all_benchmarks",)
    for line in raw.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split(None, 5)
        if len(parts) < 6:
            continue
        pid_s, rss_s, pcpu_s, state, etime, cmd = parts
        if any(m in cmd for m in self_markers):
            continue
        if NOISE_CMD.search(cmd):
            continue
        # zsh -c wrappers that only run rg over campaign paths
        if "rg -l" in cmd or "rg -n" in cmd:
            continue
        if not any(p.search(cmd) for p in PROCESS_PATTERNS):
            continue
        try:
            pid = int(pid_s)
            rss_mb = int(rss_s) / 1024.0
            pcpu = float(pcpu_s)
        except ValueError:
            continue
        out.append(
            ProcInfo(
                pid=pid,
                rss_mb=round(rss_mb, 1),
                pcpu=pcpu,
                state=state,
                etime=etime,
                command=cmd[:500],
                campaign_hint=classify_command(cmd),
            )
        )
    return out


# ─── campaign progress ────────────────────────────────────────────────────────

def count_result_csv(root: Path) -> int:
    if not root.is_dir():
        return 0
    try:
        return sum(1 for _ in root.rglob("result.csv") if _.is_file())
    except OSError:
        return 0


def list_result_targets(root: Path) -> list[str]:
    if not root.is_dir():
        return []
    names: list[str] = []
    try:
        for p in root.rglob("result.csv"):
            if p.is_file():
                names.append(p.parent.name)
    except OSError:
        pass
    return sorted(set(names))


def read_bcr_stats(root: Path) -> dict[str, Any]:
    """best_cluster_rmsd stats from result.csv files."""
    ok = 0
    bad_neg = 0
    le2 = 0
    n = 0
    fails: list[str] = []
    if not root.is_dir():
        return {"n": 0, "ok_bcr": 0, "bad_neg1": 0, "bcr_le2": 0, "fails": []}
    try:
        paths = list(root.rglob("result.csv"))
    except OSError:
        paths = []
    for rc in paths:
        try:
            with open(rc, newline="") as f:
                rows = list(csv.DictReader(f))
            if not rows:
                continue
            r = rows[0]
            bcr = r.get("best_cluster_rmsd") or r.get("BCR") or r.get("best_cluster_RMSD")
            n += 1
            if bcr is None:
                continue
            v = float(bcr)
            if v < 0:
                bad_neg += 1
            else:
                ok += 1
                if v <= 2.0:
                    le2 += 1
                else:
                    fails.append(rc.parent.name)
        except (OSError, ValueError, csv.Error):
            continue
    return {
        "n": n,
        "ok_bcr": ok,
        "bad_neg1": bad_neg,
        "bcr_le2": le2,
        "bcr_rate": (le2 / n) if n else None,
        "fails": fails[:20],
    }


def read_claim_partial(root: Path) -> dict[str, Any]:
    """Lightweight S1/S2-style counts from common CSV columns."""
    n = 0
    s1 = 0
    s2 = 0
    seed_echo = 0
    native_seed = 0
    if not root.is_dir():
        return {}
    try:
        paths = list(root.rglob("result.csv"))
    except OSError:
        return {}
    for rc in paths:
        try:
            with open(rc, newline="") as f:
                rows = list(csv.DictReader(f))
            if not rows:
                continue
            r = rows[0]
            n += 1
            # S1: success_rmsd or hungarian/rmsd fields
            sr = r.get("success_rmsd") or r.get("S1")
            if sr is not None and str(sr).strip() in ("1", "True", "true"):
                s1 += 1
            else:
                for key in ("rmsd_hungarian", "rmsd", "elected_rmsd", "rmsd_to_crystal"):
                    if key in r and r[key] not in (None, "", "-1"):
                        try:
                            if float(r[key]) <= 2.0 and float(r[key]) >= 0:
                                s1 += 1
                                break
                        except ValueError:
                            pass
            pb = r.get("success_pb") or r.get("pb_pass")
            if pb is not None and str(pb).strip() in ("1", "True", "true"):
                s2 += 1
            for sk, acc in (
                ("seed_echo", "seed_echo"),
                ("native_pose_seeded", "native_seed"),
            ):
                if sk in r:
                    try:
                        if float(r[sk]) != 0:
                            if sk == "seed_echo":
                                seed_echo += 1
                            else:
                                native_seed += 1
                    except ValueError:
                        if str(r[sk]).strip() not in ("0", "False", "false", ""):
                            if sk == "seed_echo":
                                seed_echo += 1
                            else:
                                native_seed += 1
        except (OSError, csv.Error):
            continue
    if n == 0:
        return {"n": 0}
    return {
        "n": n,
        "S1_n": s1,
        "S1_rate": s1 / n,
        "S2_n": s2,
        "S2_rate": s2 / n,
        "seed_echo_nonzero": seed_echo,
        "native_pose_seeded_nonzero": native_seed,
    }


def pid_file_status(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {"path": str(path), "exists": False, "live": False}
    try:
        pid = int(path.read_text().strip().split()[0])
    except (OSError, ValueError, IndexError):
        return {"path": str(path), "exists": True, "live": False, "pid": None}
    live = False
    try:
        os.kill(pid, 0)
        live = True
    except OSError:
        live = False
    return {"path": str(path), "exists": True, "pid": pid, "live": live}


@dataclass
class CampaignStatus:
    name: str
    path: str
    exists: bool
    result_csv_n: int
    target_total: Optional[int] = None
    metrics: dict[str, Any] = field(default_factory=dict)
    bcr: dict[str, Any] = field(default_factory=dict)
    pid_files: list[dict[str, Any]] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)


def resolve_c0_out() -> Path:
    """Locate C0 OUT across results roots (env may lag iCloud layout)."""
    name = "C0_full85_defined_cleft_nativeseed_forbidden"
    candidates = [
        results_root() / "campaigns" / name,
        icloud_root() / "results/campaigns" / name,
        Path.home()
        / "Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks"
        / "results/campaigns"
        / name,
    ]
    for c in candidates:
        if c.is_dir():
            return c
    return candidates[0]


def campaign_c0() -> CampaignStatus:
    out = resolve_c0_out()
    q = queue_root()
    st = CampaignStatus(
        name="C0_full85",
        path=str(out),
        exists=out.is_dir(),
        result_csv_n=count_result_csv(out),
        target_total=85,
        metrics=read_claim_partial(out),
        bcr=read_bcr_stats(out),
        pid_files=[
            pid_file_status(q / "logs/C0_full85.pid"),
            pid_file_status(q / "logs/C0_full85.lock"),
        ],
    )
    if st.result_csv_n >= 85:
        st.notes.append("N>=85 — ready for full claim aggregate")
    elif st.result_csv_n >= 5:
        st.notes.append("partial N>=5")
    return st


def campaign_arm(arm: str) -> CampaignStatus:
    base_candidates = [
        results_root() / "campaigns/three_engine",
        icloud_root() / "results/campaigns/three_engine",
        Path.home()
        / "Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks"
        / "results/campaigns/three_engine",
    ]
    te = next((b for b in base_candidates if b.is_dir()), base_candidates[0])
    out = te / arm / "pilot8"
    smoke = te / arm / "smoke"
    q = queue_root()
    n = count_result_csv(out)
    n_smoke = count_result_csv(smoke)
    st = CampaignStatus(
        name=f"FlexAID_{arm}_pilot8",
        path=str(out),
        exists=out.is_dir() or smoke.is_dir(),
        result_csv_n=n,
        target_total=8,
        metrics=read_claim_partial(out),
        bcr=read_bcr_stats(out),
        pid_files=[pid_file_status(q / f"logs/run_{arm}_pilot8.lock")],
        notes=[f"smoke_result_csv={n_smoke}"],
    )
    return st


def campaign_oracle() -> CampaignStatus:
    out = oracle_campaign()
    st = CampaignStatus(
        name="oracle_ceiling_v43proto_r3",
        path=str(out),
        exists=out.is_dir(),
        result_csv_n=count_result_csv(out),
        target_total=85,
        metrics=read_claim_partial(out),
        bcr=read_bcr_stats(out),
    )
    if st.bcr.get("bad_neg1", 0):
        st.notes.append("best_cluster_rmsd=-1 present — consider patch_bcr_from_poses")
    if st.result_csv_n >= 85:
        st.notes.append("COMPLETE candidate — aggregate_oracle_ceiling")
    return st


def campaign_watchers() -> dict[str, Any]:
    q = queue_root()
    return {
        "ab_chain": pid_file_status(q / "logs/run_AB_pilot8_chain.pid"),
        "throughput_maximizer": pid_file_status(q / "logs/throughput_maximizer.pid"),
        "c0": pid_file_status(q / "logs/C0_full85.pid"),
    }


def heavy_dock_count(procs: Iterable[ProcInfo]) -> int:
    """Count true dock workers (not caffeinate / sh / benchmark_datasets parent)."""
    n = 0
    for p in procs:
        cmd = p.command
        if cmd.strip().startswith("caffeinate") or cmd.strip().startswith("sh -c"):
            continue
        if "benchmark_datasets" in cmd and "FlexAID" not in cmd:
            continue
        is_binary = bool(
            re.search(r"bin/[ABC]/FlexAID(?:dS)?\b|/FlexAIDdS(\s|$)|/FlexAID(\s|$)", cmd)
        )
        if is_binary and p.rss_mb >= 100:
            n += 1
    return n


def maybe_aggregate_oracle(root: Path, scratch: Path) -> Optional[Path]:
    """If N>=85 and aggregate script exists, write rates JSON."""
    n = count_result_csv(root)
    if n < 85:
        if n >= 5:
            partial = scratch / "oracle_ceiling_rates_partial.json"
            bcr = read_bcr_stats(root)
            partial.write_text(json.dumps({"N": n, "bcr": bcr, "partial": True}, indent=2) + "\n")
            return partial
        return None
    script = repo_root() / "scripts/aggregate_oracle_ceiling.py"
    out = scratch / "oracle_ceiling_rates.json"
    if script.is_file():
        try:
            subprocess.check_call(
                [sys.executable, str(script), str(root), "--json", str(out)],
                cwd=str(repo_root()),
            )
            return out
        except subprocess.CalledProcessError:
            pass
    # fallback from BCR columns
    bcr = read_bcr_stats(root)
    payload = {
        "N": n,
        "ceiling_n": bcr.get("bcr_le2"),
        "ceiling_rate": bcr.get("bcr_rate"),
        "failed_BCR_ids": bcr.get("fails"),
        "source": "monitor_fallback",
    }
    out.write_text(json.dumps(payload, indent=2) + "\n")
    return out


# ─── report ───────────────────────────────────────────────────────────────────

def build_report() -> dict[str, Any]:
    procs = list_benchmark_processes()
    campaigns = [
        asdict(campaign_oracle()),
        asdict(campaign_c0()),
        asdict(campaign_arm("A")),
        asdict(campaign_arm("B0")),
        asdict(campaign_arm("B")),
    ]
    return {
        "ts_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "memory": mem_snapshot(),
        "processes": [asdict(p) for p in procs],
        "process_count": len(procs),
        "heavy_dock_count": heavy_dock_count(procs),
        "pid_files": campaign_watchers(),
        "campaigns": campaigns,
        "policy": {
            "never_kill_healthy": True,
            "never_dual_launch_same_out": True,
            "max_heavy_docks_recommended": 1,
        },
        "paths": {
            "repo": str(repo_root()),
            "results": str(results_root()),
            "queue": str(queue_root()),
            "oracle": str(oracle_campaign()),
        },
    }


def format_text(report: dict[str, Any]) -> str:
    lines: list[str] = []
    lines.append(f"=== benchmark multi-monitor {report['ts_utc']} ===")
    mem = report.get("memory") or {}
    if mem:
        lines.append(
            f"RAM free={mem.get('free_GB', 0):.2f}G "
            f"avail_est={mem.get('available_est_GB', 0):.2f}G "
            f"purgeable={mem.get('purgeable_GB', 0):.2f}G"
        )
    lines.append(
        f"processes={report['process_count']} heavy_docks≈{report['heavy_dock_count']}"
    )
    if report["heavy_dock_count"] > 1:
        lines.append("WARNING: >1 heavy dock process — RAM thrash risk on 18GB hosts")
    lines.append("--- LIVE PROCESSES ---")
    if not report["processes"]:
        lines.append("(none matched)")
    for p in report["processes"]:
        lines.append(
            f"  pid={p['pid']} rss={p['rss_mb']}MB cpu={p['pcpu']}% "
            f"state={p['state']} etime={p['etime']} hint={p['campaign_hint']}"
        )
        lines.append(f"    {p['command'][:160]}")
    lines.append("--- PID FILES ---")
    for k, v in (report.get("pid_files") or {}).items():
        lines.append(f"  {k}: live={v.get('live')} pid={v.get('pid')}")
    lines.append("--- CAMPAIGNS ---")
    for c in report["campaigns"]:
        tot = c.get("target_total") or "?"
        lines.append(
            f"  [{c['name']}] N={c['result_csv_n']}/{tot} exists={c['exists']}"
        )
        lines.append(f"    path={c['path']}")
        bcr = c.get("bcr") or {}
        if bcr.get("n"):
            rate = bcr.get("bcr_rate")
            rate_s = f"{rate:.3f}" if isinstance(rate, float) else "n/a"
            lines.append(
                f"    BCR: n={bcr.get('n')} le2={bcr.get('bcr_le2')} "
                f"rate={rate_s} neg1={bcr.get('bad_neg1')}"
            )
        m = c.get("metrics") or {}
        if m.get("n"):
            lines.append(
                f"    claim-ish: S1={m.get('S1_n')}/{m.get('n')} "
                f"S2={m.get('S2_n')}/{m.get('n')} "
                f"seed_echo_nz={m.get('seed_echo_nonzero')} "
                f"native_seed_nz={m.get('native_pose_seeded_nonzero')}"
            )
        for note in c.get("notes") or []:
            lines.append(f"    note: {note}")
    lines.append("--- POLICY ---")
    lines.append("never kill healthy · never dual-launch same OUT · one heavy GA preferred")
    return "\n".join(lines) + "\n"


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--json-out", type=Path, default=None, help="Write full JSON report")
    ap.add_argument(
        "--scratch",
        type=Path,
        default=None,
        help="Scratch dir for monitor_HHMM.txt and aggregates",
    )
    ap.add_argument("--quiet", action="store_true", help="JSON only to stdout if --json-out -")
    args = ap.parse_args(argv)

    scratch = args.scratch
    if scratch is None:
        env = os.environ.get("FLEXAIDDS_MONITOR_SCRATCH")
        scratch = Path(env).expanduser() if env else None
    if scratch is not None:
        scratch.mkdir(parents=True, exist_ok=True)

    report = build_report()
    text = format_text(report)

    if scratch is not None:
        hhmm = time.strftime("%H%M")
        snap = scratch / f"monitor_{hhmm}.txt"
        snap.write_text(text)
        (scratch / "monitor_latest.txt").write_text(text)
        jpath = scratch / "monitor_latest.json"
        jpath.write_text(json.dumps(report, indent=2) + "\n")
        report["snapshot_txt"] = str(snap)
        report["snapshot_json"] = str(jpath)

        # Oracle aggregate side effects (read-only re: processes)
        oracle = oracle_campaign()
        if oracle.is_dir():
            agg = maybe_aggregate_oracle(oracle, scratch)
            if agg:
                report["oracle_aggregate"] = str(agg)
                text += f"oracle_aggregate={agg}\n"

        # C0 partial rates
        c0 = resolve_c0_out()
        if c0.is_dir():
            n = count_result_csv(c0)
            if n >= 5:
                partial = {
                    "campaign": "C0_full85",
                    "N": n,
                    "metrics": read_claim_partial(c0),
                    "bcr": read_bcr_stats(c0),
                    "ts_utc": report["ts_utc"],
                }
                ppath = scratch / "c0_partial_rates.json"
                ppath.write_text(json.dumps(partial, indent=2) + "\n")
                report["c0_partial"] = str(ppath)

        snap.write_text(format_text(report) if "oracle_aggregate" in report else text)

    if args.json_out:
        if str(args.json_out) == "-":
            print(json.dumps(report, indent=2))
        else:
            args.json_out.parent.mkdir(parents=True, exist_ok=True)
            args.json_out.write_text(json.dumps(report, indent=2) + "\n")

    if not args.quiet:
        print(text if scratch is None else format_text(report), end="")

    # Exit codes: 0 ok; 2 warning multi-heavy; never fail hard on incomplete N
    if report["heavy_dock_count"] > 1:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
