#!/usr/bin/env python3
"""benchmark_ops_monitor.py — automated ops + finished-run deep analysis.

Single entry point for:
  1) Live computing status (PIDs, RAM, campaign N/totals)
  2) Deep forensic analysis of newly finished result.csv targets

Designed for scheduler / cron / agent loops. Never kills dock processes.
Never dual-launches campaigns.

**Anti-hang (production):** prefer ``~/flexaidds_results`` campaign trees.
Never ``rglob`` / recursive ``find`` under CloudDocs. Optional iCloud status
writes are best-effort with short timeouts via ``icloud_safe_io``.

Usage:
  python3 scripts/benchmark_ops_monitor.py
  python3 scripts/benchmark_ops_monitor.py --scratch /path/to/scratch
  python3 scripts/benchmark_ops_monitor.py --json-out status.json

Environment:
  FLEXAIDDS_ROOT, FLEXAIDDS_ICLOUD, FLEXAIDDS_RESULTS, FLEXAIDDS_QUEUE_ROOT
  FLEXAIDDS_LOCAL_ROOT, FLEXAIDDS_MONITOR_SCRATCH

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
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# scripts/ on sys.path for icloud_safe_io
_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))
try:
    import icloud_safe_io as _icio
except ImportError:  # pragma: no cover
    _icio = None  # type: ignore


# ─── paths ────────────────────────────────────────────────────────────────────

def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


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


def local_root() -> Path:
    if _icio is not None:
        return _icio.local_root()
    env = os.environ.get("FLEXAIDDS_LOCAL_ROOT", "").strip()
    return Path(env).expanduser() if env else Path.home() / "flexaidds_results"


def icloud_root() -> Path:
    env = os.environ.get("FLEXAIDDS_ICLOUD")
    if env:
        return Path(env).expanduser()
    return (
        Path.home()
        / "Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks"
    )


def results_root() -> Path:
    """iCloud durable results (may hang if scanned deeply — prefer local)."""
    env = os.environ.get("FLEXAIDDS_RESULTS")
    if env:
        return Path(env).expanduser()
    return icloud_root() / "results"


def local_campaigns_root() -> Path:
    return local_root() / "campaigns"


def queue_root() -> Path:
    env = os.environ.get("FLEXAIDDS_QUEUE_ROOT")
    if env:
        return Path(env).expanduser()
    # Prefer local queue staging for pid files / logs
    local_q = local_root() / "three_engine_entropy_q1"
    if (local_q / "logs").is_dir() or (local_q / "bin").is_dir():
        return local_q
    return icloud_root() / "queues/three_engine_entropy_q1"


def default_scratch() -> Path:
    env = os.environ.get("FLEXAIDDS_MONITOR_SCRATCH")
    if env:
        return Path(env).expanduser()
    # Prefer local ops dir over /tmp that vanishes
    return local_root() / "logs" / "ops_monitor"

# ─── ops: memory + processes ──────────────────────────────────────────────────

def mem_snapshot() -> Dict[str, float]:
    try:
        ps = int(subprocess.check_output(["pagesize"], text=True).strip())
    except Exception:
        ps = 4096
    try:
        raw = subprocess.check_output(["vm_stat"], text=True)
    except Exception:
        return {}
    d: Dict[str, int] = {}
    for line in raw.splitlines()[1:]:
        if ":" not in line:
            continue
        k, v = line.split(":", 1)
        try:
            d[k.strip()] = int(v.strip().rstrip("."))
        except ValueError:
            pass

    def gb(n: float) -> float:
        return n * ps / 1e9

    free = d.get("Pages free", 0)
    pur = d.get("Pages purgeable", 0)
    ina = d.get("Pages inactive", 0)
    return {
        "free_GB": round(gb(free), 3),
        "purgeable_GB": round(gb(pur), 3),
        "available_est_GB": round(gb(free + pur + 0.3 * ina), 3),
    }


PROC_RE = re.compile(
    r"bin/[ABC]/FlexAID(?:dS)?\b|/FlexAIDdS(\s|$)|benchmark_datasets\b|"
    r"serial_AB_after_C0|throughput_maximizer|run_flexaid_arm_pilot8|"
    r"C0_claim_clean|C0_full85"
)


def list_procs() -> List[Dict[str, Any]]:
    try:
        raw = subprocess.check_output(
            ["ps", "-axo", "pid=,rss=,pcpu=,state=,etime=,command="],
            text=True,
            stderr=subprocess.DEVNULL,
        )
    except Exception:
        return []
    out: List[Dict[str, Any]] = []
    for line in raw.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split(None, 5)
        if len(parts) < 6:
            continue
        pid_s, rss_s, pcpu_s, state, etime, cmd = parts
        if "benchmark_ops_monitor" in cmd or "monitor_all_benchmarks" in cmd:
            continue
        if re.search(r"\brg\s|^\s*(rg|grep|ugrep)\b", cmd):
            continue
        if not PROC_RE.search(cmd):
            continue
        try:
            out.append(
                {
                    "pid": int(pid_s),
                    "rss_mb": round(int(rss_s) / 1024.0, 1),
                    "pcpu": float(pcpu_s),
                    "state": state,
                    "etime": etime,
                    "command": cmd[:400],
                    "hint": _hint(cmd),
                }
            )
        except ValueError:
            continue
    return out


def _hint(cmd: str) -> str:
    if "claim_g2000" in cmd or "C0_claim" in cmd:
        return "C0_claim"
    if "C0_full85_defined_cleft" in cmd or "defined_cleft_nativeseed" in cmd:
        return "C0_legacy"
    if "three_engine/A" in cmd or "/bin/A/FlexAID" in cmd:
        return "arm_A"
    if "three_engine/B0" in cmd:
        return "arm_B0"
    if "three_engine/B/" in cmd or "TEMPER 21" in cmd:
        return "arm_B"
    if "FlexAIDdS" in cmd:
        return "FlexAIDdS"
    if "benchmark_datasets" in cmd:
        return "runner"
    return "other"


def current_dock_target(procs: List[Dict[str, Any]]) -> Optional[str]:
    for p in procs:
        cmd = p["command"]
        if "FlexAID" not in cmd and "FlexAIDdS" not in cmd:
            continue
        m = re.search(r"astex_diverse/([0-9A-Za-z]{4})/", cmd)
        if m:
            return m.group(1)
        # local OUT: .../campaigns/.../1IA1/r4/1IA1 or .../1IA1/1IA1_dockin
        m = re.search(r"/([0-9A-Za-z]{4})/(?:r\d+/)?(?:\1_|dock_config)", cmd)
        if m:
            return m.group(1)
        m = re.search(r"/([0-9A-Za-z]{4})/[0-9A-Za-z]{4}_dockin", cmd)
        if m:
            return m.group(1)
    return None


def pid_file(path: Path) -> Dict[str, Any]:
    if not path.is_file():
        return {"exists": False, "live": False, "pid": None}
    try:
        pid = int(path.read_text().strip().split()[0])
    except Exception:
        return {"exists": True, "live": False, "pid": None}
    live = False
    try:
        os.kill(pid, 0)
        live = True
    except OSError:
        live = False
    return {"exists": True, "live": live, "pid": pid, "path": str(path)}


# ─── campaign scan ────────────────────────────────────────────────────────────

CAMPAIGN_SPECS = [
    {
        "id": "C0_claim",
        "rel": "campaigns/C0_full85_claim_g2000_popmod_20260715",
        "total": 85,
        "claim": True,
    },
    {
        "id": "C0_legacy",
        "rel": "campaigns/C0_full85_defined_cleft_nativeseed_forbidden",
        "total": 85,
        "claim": False,
    },
    {
        "id": "A_pilot8",
        "rel": "campaigns/three_engine/A/pilot8",
        "total": 8,
        "claim": True,
    },
    {
        "id": "B_pilot8",
        "rel": "campaigns/three_engine/B/pilot8",
        "total": 8,
        "claim": True,
    },
    {
        "id": "B0_pilot8",
        "rel": "campaigns/three_engine/B0/pilot8",
        "total": 8,
        "claim": False,
    },
]


def _ffloat(x: Any) -> Optional[float]:
    try:
        if x is None or x == "":
            return None
        return float(x)
    except (TypeError, ValueError):
        return None


def parse_result_csv(path: Path) -> Dict[str, Any]:
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return {"empty": True, "path": str(path)}
    r = rows[0]
    rh = _ffloat(r.get("rmsd_hungarian"))
    rtc = _ffloat(r.get("rmsd_to_crystal"))
    bcr = _ffloat(r.get("best_cluster_rmsd"))
    sr_raw = r.get("success_rmsd")
    s1 = str(sr_raw).strip() in ("1", "True", "true")
    if not s1 and rh is not None and 0 <= rh <= 2.0:
        s1 = True
    pb = str(r.get("success_pb") or r.get("pb_pass") or "").strip() in (
        "1",
        "True",
        "true",
    )
    bcr_ok = bcr is not None and 0 <= bcr <= 2.0
    bcr_neg = bcr is not None and bcr < 0
    nposes = r.get("num_poses")
    try:
        nposes_i = int(float(nposes)) if nposes not in (None, "") else None
    except ValueError:
        nposes_i = None
    cf = _ffloat(r.get("best_score") or r.get("elected_cf"))
    return {
        "path": str(path),
        "pdb_id": r.get("pdb_id") or path.parent.name,
        "rmsd_hungarian": rh,
        "rmsd_to_crystal": rtc,
        "best_cluster_rmsd": bcr,
        "success_rmsd": sr_raw,
        "s1": int(s1),
        "s2": int(s1 and pb),
        "pb_pass": pb,
        "bcr_le2": int(bcr_ok),
        "bcr_neg1": int(bcr_neg),
        "election_gap": int(bcr_ok and not s1),
        "packaging_bug": int(bcr_neg or (nposes_i == 0 and path.parent.exists())),
        "best_score": cf,
        "num_poses": nposes_i,
        "wall_time_s": _ffloat(r.get("wall_time_s")),
        "seed_echo": r.get("seed_echo"),
        "native_pose_seeded": r.get("native_pose_seeded"),
        "pb_backend": r.get("pb_backend"),
        "elected_pose_path": r.get("elected_pose_path"),
        "mtime": path.stat().st_mtime,
    }


def classify(rec: Dict[str, Any]) -> List[str]:
    tags: List[str] = []
    if rec.get("empty"):
        return ["empty_csv"]
    if rec["s1"]:
        tags.append("S1_HIT")
    else:
        tags.append("S1_MISS")
    if rec["s2"]:
        tags.append("S2_HIT")
    else:
        tags.append("S2_MISS")
    if rec["bcr_le2"]:
        tags.append("BCR_OK")
    elif rec["bcr_neg1"]:
        tags.append("BCR_SENTINEL")
    else:
        tags.append("BCR_MISS")
    if rec["election_gap"]:
        tags.append("ELECTION_GAP")
    if rec["packaging_bug"]:
        tags.append("PACKAGING_BUG")
    cf = rec.get("best_score")
    if cf is not None and abs(cf) > 400 and not rec["s1"]:
        tags.append("PATHOLOGICAL_CF")
    se = str(rec.get("seed_echo") or "0")
    ns = str(rec.get("native_pose_seeded") or "0")
    if se not in ("0", "0.0", "False", "false", ""):
        tags.append("SEED_ECHO")
    if ns not in ("0", "0.0", "False", "false", ""):
        tags.append("NATIVE_SEEDED")
    return tags


def scan_campaign(camp_id: str, root: Path, total: int) -> Dict[str, Any]:
    """Scan one campaign. Uses one-level glob only (never rglob). CloudDocs-safe."""
    empty = {
        "id": camp_id,
        "path": str(root),
        "exists": False,
        "N": 0,
        "total": total,
        "results": [],
        "storage": "unknown",
    }
    try:
        exists = root.is_dir()
    except OSError:
        return empty
    if not exists:
        return empty

    storage = "clouddocs" if (_icio and _icio.is_clouddocs(root)) else "local"
    if _icio is not None:
        csv_paths = _icio.safe_glob_result_csvs(root, timeout_s=20.0)
    else:
        try:
            csv_paths = [
                p
                for p in sorted(root.glob("*/result.csv"))
                if "incomplete" not in p.parent.name
            ]
        except OSError:
            csv_paths = []

    results: List[Dict[str, Any]] = []
    for rc in csv_paths:
        try:
            results.append(parse_result_csv(rc))
        except Exception as e:
            results.append(
                {"pdb_id": rc.parent.name, "error": str(e), "path": str(rc)}
            )
    n = len(results)
    s1 = sum(1 for r in results if r.get("s1"))
    s2 = sum(1 for r in results if r.get("s2"))
    bcr = sum(1 for r in results if r.get("bcr_le2"))
    gap = sum(1 for r in results if r.get("election_gap"))
    neg = sum(1 for r in results if r.get("bcr_neg1"))
    return {
        "id": camp_id,
        "path": str(root),
        "exists": True,
        "N": n,
        "total": total,
        "S1": s1,
        "S2": s2,
        "BCR_le2": bcr,
        "election_gap": gap,
        "bcr_neg1": neg,
        "S1_rate": (s1 / n) if n else None,
        "S2_rate": (s2 / n) if n else None,
        "BCR_rate": (bcr / n) if n else None,
        "results": results,
        "storage": storage,
    }


# ─── deep analysis state ──────────────────────────────────────────────────────

def load_state(path: Path) -> Dict[str, Any]:
    if path.is_file():
        try:
            return json.loads(path.read_text())
        except Exception:
            pass
    return {"analyzed": {}, "last_run_utc": None}


def save_state(path: Path, state: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(state, indent=2) + "\n")


def dock_config_budget(target_dir: Path) -> Dict[str, Any]:
    """Pull pop/gen from dock_config.json (shallow only — no **/ under CloudDocs)."""
    # Prefer top-level then one restart subdir; never rglob.
    candidates = [target_dir / "dock_config.json"]
    try:
        for sub in sorted(target_dir.glob("r*/dock_config.json"))[:5]:
            candidates.append(sub)
    except OSError:
        pass
    for p in candidates:
        try:
            if not p.is_file():
                continue
            d = json.loads(p.read_text())
            ga = d.get("ga") or d.get("genetic_algorithm") or {}
            return {
                "num_chromosomes": ga.get("num_chromosomes") or ga.get("population"),
                "num_generations": ga.get("num_generations") or ga.get("generations"),
                "temperature": d.get("temperature")
                or (d.get("thermodynamics") or {}).get("temperature"),
                "clustering": (d.get("clustering_algorithm") or d.get("clustering")),
                "path": str(p),
            }
        except Exception:
            continue
    return {}


def analyze_new(
    campaigns: List[Dict[str, Any]], state: Dict[str, Any]
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    analyzed = state.setdefault("analyzed", {})
    new_items: List[Dict[str, Any]] = []
    for camp in campaigns:
        if not camp.get("exists"):
            continue
        for rec in camp.get("results") or []:
            if rec.get("error") or rec.get("empty"):
                continue
            key = f"{camp['id']}/{rec.get('pdb_id')}"
            mtime = rec.get("mtime") or 0
            prev = analyzed.get(key)
            if prev and prev.get("mtime") == mtime:
                continue
            tags = classify(rec)
            budget = dock_config_budget(Path(rec["path"]).parent)
            item = {
                **rec,
                "campaign": camp["id"],
                "tags": tags,
                "budget": budget,
                "analyzed_utc": utc_now(),
            }
            new_items.append(item)
            analyzed[key] = {
                "mtime": mtime,
                "s1": rec.get("s1"),
                "s2": rec.get("s2"),
                "bcr": rec.get("best_cluster_rmsd"),
                "tags": tags,
                "ts": utc_now(),
            }
    state["last_run_utc"] = utc_now()
    return new_items, state


# ─── reporting ────────────────────────────────────────────────────────────────

def format_ops_brief(
    mem: Dict[str, float],
    procs: List[Dict[str, Any]],
    campaigns: List[Dict[str, Any]],
    pid_files: Dict[str, Any],
    new_n: int,
    cur: Optional[str],
) -> str:
    lines = [
        f"# Benchmark ops + results — {utc_now()}",
        "",
        f"RAM free={mem.get('free_GB', 0):.2f}G avail_est={mem.get('available_est_GB', 0):.2f}G",
        f"Live processes={len(procs)} current_target={cur or 'none'} new_finishes={new_n}",
        "",
        "## Campaigns",
    ]
    for c in campaigns:
        if not c.get("exists") and c["N"] == 0:
            lines.append(f"- **{c['id']}**: missing/empty")
            continue
        rate = c.get("S1_rate")
        rate_s = f"{100*rate:.1f}%" if rate is not None else "n/a"
        lines.append(
            f"- **{c['id']}**: N={c['N']}/{c['total']}  "
            f"S1={c.get('S1',0)}/{c['N']} ({rate_s})  "
            f"S2={c.get('S2',0)}  BCR≤2={c.get('BCR_le2',0)}  "
            f"gap={c.get('election_gap',0)}  neg1={c.get('bcr_neg1',0)}"
        )
    lines.append("")
    lines.append("## PID files")
    for k, v in pid_files.items():
        lines.append(f"- {k}: live={v.get('live')} pid={v.get('pid')}")
    lines.append("")
    lines.append("## Live workers")
    if not procs:
        lines.append("- (none)")
    for p in procs[:8]:
        lines.append(
            f"- pid={p['pid']} rss={p['rss_mb']}MB cpu={p['pcpu']}% "
            f"hint={p['hint']} etime={p['etime']}"
        )
    lines.append("")
    # completion
    claim = next((c for c in campaigns if c["id"] == "C0_claim"), None)
    legacy = next((c for c in campaigns if c["id"] == "C0_legacy"), None)
    a = next((c for c in campaigns if c["id"] == "A_pilot8"), None)
    b = next((c for c in campaigns if c["id"] == "B_pilot8"), None)
    claim_done = claim and claim["N"] >= claim["total"]
    if claim_done and a and a["N"] >= 8 and b and b["N"] >= 8:
        lines.append("**COMPLETE** claim path C0+A+B pilots.")
    elif claim_done:
        lines.append("**C0 claim COMPLETE** — next: A pilot8 → B@TEMPER21 (serial, iCloud).")
    elif claim and claim["N"] == 0 and legacy and legacy["N"] > 0:
        lines.append(
            "**NOT COMPLETE** — dirty legacy C0 still has results; "
            "clean claim OUT empty (launch run_C0_claim_clean.sh if reboot intended)."
        )
    else:
        lines.append("**NOT COMPLETE** — campaigns still running or not started.")
    lines.append("")
    lines.append(
        "CF = scoring proxy. S1 = primary. BCR = diagnostic. No dual-launch. One heavy GA."
    )
    return "\n".join(lines) + "\n"


def format_new_targets_md(new_items: List[Dict[str, Any]]) -> str:
    if not new_items:
        return ""
    lines = [f"## New finishes ({len(new_items)})", ""]
    for it in new_items:
        tags = ",".join(it.get("tags") or [])
        lines.append(
            f"### {it.get('campaign')}/{it.get('pdb_id')}  `{tags}`"
        )
        lines.append(
            f"- HUNG={it.get('rmsd_hungarian')} XTAL={it.get('rmsd_to_crystal')} "
            f"BCR={it.get('best_cluster_rmsd')} CF={it.get('best_score')} "
            f"poses={it.get('num_poses')} wall_s={it.get('wall_time_s')}"
        )
        lines.append(
            f"- S1={it.get('s1')} S2={it.get('s2')} PB={it.get('pb_pass')} "
            f"seed_echo={it.get('seed_echo')} native_seeded={it.get('native_pose_seeded')} "
            f"backend={it.get('pb_backend')}"
        )
        b = it.get("budget") or {}
        if b:
            lines.append(
                f"- dock_config: chroms={b.get('num_chromosomes')} "
                f"gen={b.get('num_generations')} T={b.get('temperature')} "
                f"clust={b.get('clustering')}"
            )
        lines.append("")
    return "\n".join(lines)


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scratch", type=Path, default=None)
    ap.add_argument("--json-out", type=Path, default=None)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    scratch = args.scratch or default_scratch()
    scratch.mkdir(parents=True, exist_ok=True)

    # Prefer local campaign trees; iCloud is durable mirror only.
    res = results_root()
    local_camps = local_campaigns_root()
    q = queue_root()
    local_log = local_root() / "logs" / "C0_claim"

    mem = mem_snapshot()
    procs = list_procs()
    cur = current_dock_target(procs)
    # Pid files: local first, then queue (may be iCloud — only open single small files)
    def _pid_prefer(*paths: Path) -> Dict[str, Any]:
        for p in paths:
            info = pid_file(p)
            if info.get("exists"):
                return info
        return {"exists": False, "live": False, "pid": None}

    pid_files = {
        "C0_claim": _pid_prefer(
            local_log / "C0_claim_clean.pid",
            q / "logs/C0_claim_clean.pid",
        ),
        "C0_legacy": _pid_prefer(q / "logs/C0_full85.pid"),
        "ab_chain": _pid_prefer(q / "logs/run_AB_pilot8_chain.pid"),
        "throughput": _pid_prefer(q / "logs/throughput_maximizer.pid"),
        "sync_loop": _pid_prefer(local_log / "sync_ops_loop.pid"),
    }

    campaigns: List[Dict[str, Any]] = []
    for spec in CAMPAIGN_SPECS:
        # Local campaign dir wins when present
        camp_name = Path(spec["rel"]).name
        local_root_c = local_camps / camp_name
        icloud_root_c = res / spec["rel"]
        try:
            local_ok = local_root_c.is_dir()
        except OSError:
            local_ok = False
        root = local_root_c if local_ok else icloud_root_c
        campaigns.append(scan_campaign(spec["id"], root, spec["total"]))
    # Discover extra local campaigns (e.g. C0_v137_*) without walking iCloud
    try:
        known = {Path(s["rel"]).name for s in CAMPAIGN_SPECS}
        if local_camps.is_dir():
            for d in sorted(local_camps.iterdir()):
                if not d.is_dir() or d.name in known:
                    continue
                if d.name.startswith(("C0_", "DPFO_", "three_engine")):
                    campaigns.append(scan_campaign(d.name, d, 85))
    except OSError:
        pass
    state_path = scratch / "finished_run_analysis_state.json"
    state = load_state(state_path)
    new_items, state = analyze_new(campaigns, state)
    save_state(state_path, state)

    # Write deep log append
    deep_log = scratch / "FINISHED_RUN_DEEP_LOG.md"
    with deep_log.open("a") as f:
        f.write(f"\n---\n# {utc_now()}\n\n")
        if new_items:
            f.write(format_new_targets_md(new_items))
        else:
            f.write(
                f"no new finished targets | current={cur} | "
                + " ".join(
                    f"{c['id']}={c['N']}/{c['total']}" for c in campaigns if c.get("exists") or c["N"]
                )
                + "\n"
            )

    brief = format_ops_brief(mem, procs, campaigns, pid_files, len(new_items), cur)
    if new_items:
        brief += "\n" + format_new_targets_md(new_items)

    (scratch / "finished_run_latest.md").write_text(brief)
    hhmm = time.strftime("%H%M")
    (scratch / f"monitor_{hhmm}.txt").write_text(brief)
    (scratch / "monitor_latest.txt").write_text(brief)

    # Mirror brief to *local* ops log first (never block on CloudDocs).
    try:
        llog = local_root() / "logs" / "ops"
        llog.mkdir(parents=True, exist_ok=True)
        (llog / "benchmark_ops_latest.md").write_text(brief)
    except OSError:
        pass
    # Optional thin iCloud write — best-effort only, no retry loops
    try:
        if _icio is None or not _icio.is_clouddocs(q):
            qlog = q / "logs"
            qlog.mkdir(parents=True, exist_ok=True)
            (qlog / "benchmark_ops_latest.md").write_text(brief)
    except OSError:
        pass

    payload = {
        "ts_utc": utc_now(),
        "memory": mem,
        "processes": procs,
        "current_target": cur,
        "pid_files": pid_files,
        "campaigns": [
            {k: v for k, v in c.items() if k != "results"} for c in campaigns
        ],
        "new_finishes": [
            {
                "campaign": i.get("campaign"),
                "pdb_id": i.get("pdb_id"),
                "tags": i.get("tags"),
                "s1": i.get("s1"),
                "s2": i.get("s2"),
                "bcr": i.get("best_cluster_rmsd"),
                "cf": i.get("best_score"),
            }
            for i in new_items
        ],
        "paths": {
            "results_icloud": str(res),
            "campaigns_local": str(local_camps),
            "queue": str(q),
            "scratch": str(scratch),
            "local_root": str(local_root()),
        },
        "anti_hang": "local_first_no_clouddocs_rglob",
    }
    (scratch / "monitor_latest.json").write_text(json.dumps(payload, indent=2) + "\n")
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(payload, indent=2) + "\n")

    if not args.quiet:
        print(brief)

    heavy = sum(
        1
        for p in procs
        if p["rss_mb"] >= 200 and "FlexAID" in p["command"]
    )
    return 2 if heavy > 1 else 0


if __name__ == "__main__":
    sys.exit(main())
