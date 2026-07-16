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
  FORCE_COLOR=1 python3 scripts/benchmark_ops_monitor.py   # ANSI even when piped

Environment:
  FLEXAIDDS_ROOT, FLEXAIDDS_ICLOUD, FLEXAIDDS_RESULTS, FLEXAIDDS_QUEUE_ROOT
  FLEXAIDDS_LOCAL_ROOT, FLEXAIDDS_MONITOR_SCRATCH
  FORCE_COLOR / NO_COLOR  — control ANSI coloring on stdout

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
    r"run_flexaid_arm_pilot8|run_3dsig_red_pair_serial|"
    r"throughput_maximizer"
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
    """Process labels match campaign path IDs — no ad-hoc nicknames."""
    if "run_3dsig_red_pair_full85" in cmd:
        return "three_engine/full85_serial"
    if "generate_flexaid_inp" in cmd:
        return "three_engine/prep"
    if "ProcessLigand" in cmd:
        return "three_engine/ProcessLigand"
    if "--full85" in cmd and "run_flexaid_arm" in cmd:
        if " A " in f" {cmd} " or cmd.rstrip().endswith(" A") or "arms A" in cmd:
            return "three_engine/A/full85"
        if "B0" in cmd:
            return "three_engine/B0/full85"
        if re.search(r"\bB\b", cmd) and "B0" not in cmd:
            return "three_engine/B/full85"
        return "three_engine/full85_arm"
    if "three_engine/A/" in cmd or "/work/A/" in cmd or "/bin/A/FlexAID" in cmd:
        return "three_engine/A"
    if "three_engine/B0/" in cmd or "/work/B0/" in cmd:
        return "three_engine/B0"
    if "three_engine/B/" in cmd or "/work/B/" in cmd or "/bin/B/FlexAID" in cmd:
        return "three_engine/B"
    if "run_3dsig_red_pair_serial" in cmd:
        return "three_engine/serial"
    if "run_flexaid_arm_pilot8" in cmd:
        return "three_engine/arm_launcher"
    if "FlexAIDdS" in cmd:
        return "FlexAIDdS"
    if "benchmark_datasets" in cmd:
        return "benchmark_datasets"
    return "other"


def current_dock_target(procs: List[Dict[str, Any]]) -> Optional[str]:
    for p in procs:
        cmd = p["command"]
        if "FlexAID" not in cmd and "FlexAIDdS" not in cmd and "generate_flexaid" not in cmd:
            continue
        for pat in (
            r"/3dsig_full85[^/]*/([0-9A-Za-z]{4})/",
            r"/3dsig_full85_r\d+/([0-9A-Za-z]{4})/",
            r"/work(?:_scratch_[^/]+)?/(?:A|B0|B|C)/([0-9A-Za-z]{4})/",
            r"astex_diverse/([0-9A-Za-z]{4})/",
            r"/([0-9A-Za-z]{4})/(?:r\d+/)?(?:\1_|dock_config)",
            r"/([0-9A-Za-z]{4})/[0-9A-Za-z]{4}_dockin",
        ):
            m = re.search(pat, cmd)
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
#
# Live science: full85 red-pair ONLY (A → B0 → B serial). Pilot8 is not tracked.
# Campaign id from FLEXAID_CAMPAIGN (default: latest scratch pin on disk).
# C0_claim / C0_legacy paths are REMOVED from ops — suspended; do not surface.


def _active_full85_campaign() -> str:
    env = os.environ.get("FLEXAIDDS_OPS_CAMPAIGN") or os.environ.get("FLEXAID_CAMPAIGN")
    if env and env.strip():
        return env.strip()
    # Prefer newest three_engine/A/3dsig_full85_* with local results root
    root = local_campaigns_root() / "three_engine" / "A"
    try:
        cands = sorted(
            (
                p
                for p in root.iterdir()
                if p.is_dir() and p.name.startswith("3dsig_full85")
            ),
            key=lambda p: p.stat().st_mtime,
            reverse=True,
        )
        if cands:
            return cands[0].name
    except OSError:
        pass
    return "3dsig_full85_scratch_3b2fa57cc"


def _full85_campaign_specs() -> List[Dict[str, Any]]:
    camp = _active_full85_campaign()
    return [
        {
            "id": f"three_engine/A/{camp}",
            "rel": f"three_engine/A/{camp}",
            "total": 85,
            "panel": "full85",
            "arm": "A",
            "description": f"FlexAID TEMPER0 CLUSTA CF · {camp} R=1",
        },
        {
            "id": f"three_engine/B0/{camp}",
            "rel": f"three_engine/B0/{camp}",
            "total": 85,
            "panel": "full85",
            "arm": "B0",
            "description": f"master TEMPER0 CLUSTA CF · {camp} R=1",
        },
        {
            "id": f"three_engine/B/{camp}",
            "rel": f"three_engine/B/{camp}",
            "total": 85,
            "panel": "full85",
            "arm": "B",
            "description": f"master TEMPER21 CLUSTA FO · {camp} R=1",
        },
    ]


# Evaluated at import for static tooling; main() rebinds via CAMPAIGN_SPECS =
CAMPAIGN_SPECS = _full85_campaign_specs()


def _ffloat(x: Any) -> Optional[float]:
    try:
        if x is None or x == "":
            return None
        return float(x)
    except (TypeError, ValueError):
        return None


def _truthy_int(val: Any) -> Optional[int]:
    """Parse 0/1 / true/false; None if absent."""
    if val is None or val == "":
        return None
    s = str(val).strip().lower()
    if s in ("1", "true", "yes"):
        return 1
    if s in ("0", "false", "no"):
        return 0
    try:
        return 1 if int(float(s)) != 0 else 0
    except (TypeError, ValueError):
        return None


def parse_result_csv(path: Path) -> Dict[str, Any]:
    """Parse DatasetRunner or classic FlexAID parse_flexaid_arm_results result.csv."""
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return {"empty": True, "path": str(path)}
    r = rows[0]
    # RMSD fields: DatasetRunner + classic arm parser
    rh = _ffloat(r.get("rmsd_hungarian") or r.get("rmsd_top1"))
    rtc = _ffloat(r.get("rmsd_to_crystal") or r.get("rmsd_top1"))
    bcr = _ffloat(r.get("best_cluster_rmsd") or r.get("rmsd_bcr"))
    # S1: prefer explicit success_s1 / success_rmsd; else top1 ≤ 2.0 Å
    s1_flag = _truthy_int(r.get("success_s1"))
    if s1_flag is None:
        s1_flag = _truthy_int(r.get("success_rmsd"))
    if s1_flag is None:
        s1_flag = int(rh is not None and 0 <= rh <= 2.0)
    s1 = int(s1_flag)
    # S_top10 (3Dsig red-pair primary when present)
    s_top10_flag = _truthy_int(r.get("success_s_top10"))
    if s_top10_flag is None:
        # Derive from mode_rmsd_0..9 if columns exist
        modes = [_ffloat(r.get(f"mode_rmsd_{i}")) for i in range(10)]
        if any(m is not None for m in modes):
            s_top10_flag = int(
                any(m is not None and 0 <= m <= 2.0 for m in modes)
            )
        else:
            s_top10_flag = 0
    s_top10 = int(s_top10_flag)
    s3_flag = _truthy_int(r.get("success_s3"))
    if s3_flag is None:
        s3_flag = int(bcr is not None and 0 <= bcr <= 2.0)
    pb = str(r.get("success_pb") or r.get("pb_pass") or "").strip() in (
        "1",
        "True",
        "true",
    )
    bcr_ok = bcr is not None and 0 <= bcr <= 2.0
    bcr_neg = bcr is not None and bcr < 0
    nposes = r.get("num_poses") or r.get("n_poses")
    try:
        nposes_i = int(float(nposes)) if nposes not in (None, "") else None
    except ValueError:
        nposes_i = None
    cf = _ffloat(r.get("best_score") or r.get("elected_cf") or r.get("score_top1"))
    wall = _ffloat(r.get("wall_time_s") or r.get("wall_s"))
    wall_src = "csv" if wall is not None and wall > 0 else None
    if wall is None or wall <= 0:
        # Launcher file or pose-mtime proxy for historical empty wall_s cells
        wall, wall_src = _wall_from_target_dir(path.parent, r.get("pdb_id") or path.parent.name)
    return {
        "path": str(path),
        "pdb_id": r.get("pdb_id") or path.parent.name,
        "rmsd_hungarian": rh,
        "rmsd_to_crystal": rtc,
        "rmsd_top1": _ffloat(r.get("rmsd_top1")) or rh,
        "best_cluster_rmsd": bcr,
        "success_rmsd": r.get("success_rmsd") or r.get("success_s1"),
        "s1": s1,
        "s_top10": s_top10,
        "s2": int(s1 and pb),
        "s3": int(s3_flag),
        "pb_pass": pb,
        "bcr_le2": int(bcr_ok),
        "bcr_neg1": int(bcr_neg),
        "election_gap": int(bcr_ok and not s1),
        "packaging_bug": int(bcr_neg or (nposes_i == 0 and path.parent.exists())),
        "best_score": cf,
        "num_poses": nposes_i,
        "wall_time_s": wall,
        "wall_src": wall_src,
        "seed_echo": r.get("seed_echo"),
        "native_pose_seeded": r.get("native_pose_seeded"),
        "pb_backend": r.get("pb_backend"),
        "elected_pose_path": r.get("elected_pose_path") or r.get("elected_path"),
        "mtime": path.stat().st_mtime,
    }


def _wall_from_target_dir(
    target_dir: Path, pdb_id: str
) -> Tuple[Optional[float], Optional[str]]:
    """Recover wall seconds from wall_s.txt / wall_timing.json / pose mtimes."""
    for name in ("wall_s.txt", "wall_timing.json"):
        wp = target_dir / name
        try:
            if not wp.is_file():
                continue
            text = wp.read_text().strip()
            if name.endswith(".json"):
                d = json.loads(text)
                v = _ffloat(d.get("wall_s"))
                if v is not None and v > 0:
                    return v, "launcher"
            else:
                v = _ffloat(text.split()[0] if text else None)
                if v is not None and v > 0:
                    return v, "launcher"
        except (OSError, ValueError, TypeError, json.JSONDecodeError):
            continue
    pdb = str(pdb_id or target_dir.name).upper()
    mt: List[float] = []
    try:
        for pat in (f"{pdb}_r*_*.pdb", f"{pdb}_*.pdb"):
            for p in target_dir.glob(pat):
                try:
                    if p.is_file():
                        mt.append(p.stat().st_mtime)
                except OSError:
                    continue
            if len(mt) >= 2:
                break
    except OSError:
        return None, None
    if len(mt) < 2:
        return None, None
    span = max(mt) - min(mt)
    # Reject absurd spans (bulk copy / iCloud re-touch); keep per-target docking-scale
    if span <= 1.0 or span > 48 * 3600:
        return None, None
    return float(span), "mtime_proxy"


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
    s_top10 = sum(1 for r in results if r.get("s_top10"))
    bcr = sum(1 for r in results if r.get("bcr_le2"))
    gap = sum(1 for r in results if r.get("election_gap"))
    neg = sum(1 for r in results if r.get("bcr_neg1"))
    pack = sum(1 for r in results if r.get("packaging_bug"))
    bcr_vals = sorted(
        r["best_cluster_rmsd"]
        for r in results
        if isinstance(r.get("best_cluster_rmsd"), (int, float))
        and r["best_cluster_rmsd"] is not None
        and r["best_cluster_rmsd"] >= 0
    )
    top1_vals = sorted(
        r["rmsd_hungarian"]
        for r in results
        if isinstance(r.get("rmsd_hungarian"), (int, float))
        and r["rmsd_hungarian"] is not None
        and r["rmsd_hungarian"] >= 0
    )

    def _median(xs: List[float]) -> Optional[float]:
        if not xs:
            return None
        m = len(xs) // 2
        return float(xs[m]) if len(xs) % 2 else 0.5 * (xs[m - 1] + xs[m])

    wall_vals = sorted(
        float(r["wall_time_s"])
        for r in results
        if isinstance(r.get("wall_time_s"), (int, float))
        and r["wall_time_s"] is not None
        and r["wall_time_s"] > 0
    )
    wall_srcs = {
        str(r.get("wall_src") or "?")
        for r in results
        if isinstance(r.get("wall_time_s"), (int, float)) and r.get("wall_time_s")
    }
    wall_sum = sum(wall_vals) if wall_vals else None
    wall_med = _median(wall_vals)
    wall_mean = (sum(wall_vals) / len(wall_vals)) if wall_vals else None
    # Project remaining wall at median for unfinished campaign
    wall_eta = None
    if wall_med is not None and total > n:
        wall_eta = wall_med * (total - n)

    return {
        "id": camp_id,
        "path": str(root),
        "exists": True,
        "N": n,
        "total": total,
        "S1": s1,
        "S2": s2,
        "S_top10": s_top10,
        "BCR_le2": bcr,
        "election_gap": gap,
        "bcr_neg1": neg,
        "packaging_bug": pack,
        "S1_rate": (s1 / n) if n else None,
        "S2_rate": (s2 / n) if n else None,
        "S_top10_rate": (s_top10 / n) if n else None,
        "BCR_rate": (bcr / n) if n else None,
        "BCR_median": _median(bcr_vals),
        "top1_median": _median(top1_vals),
        "BCR_best": bcr_vals[0] if bcr_vals else None,
        "wall_n": len(wall_vals),
        "wall_sum_s": wall_sum,
        "wall_median_s": wall_med,
        "wall_mean_s": wall_mean,
        "wall_min_s": wall_vals[0] if wall_vals else None,
        "wall_max_s": wall_vals[-1] if wall_vals else None,
        "wall_eta_remain_s": wall_eta,
        "wall_src": "+".join(sorted(wall_srcs)) if wall_srcs else None,
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


# ─── reporting (color + dashboard) ────────────────────────────────────────────

class _Ansi:
    """Terminal colors (TTY or FORCE_COLOR=1). Markdown files stay emoji-only."""

    RESET = "\033[0m"
    BOLD = "\033[1m"
    DIM = "\033[2m"
    UNDER = "\033[4m"
    RED = "\033[31m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"
    BLUE = "\033[34m"
    MAGENTA = "\033[35m"
    CYAN = "\033[36m"
    WHITE = "\033[37m"
    BRIGHT_RED = "\033[91m"
    BRIGHT_GREEN = "\033[92m"
    BRIGHT_YELLOW = "\033[93m"
    BRIGHT_BLUE = "\033[94m"
    BRIGHT_CYAN = "\033[96m"
    BRIGHT_WHITE = "\033[97m"
    BG_RED = "\033[41m"
    BG_GREEN = "\033[42m"
    BG_YELLOW = "\033[43m"
    BG_BLUE = "\033[44m"
    BG_MAGENTA = "\033[45m"
    BG_CYAN = "\033[46m"
    BG_DARK = "\033[100m"

    @classmethod
    def enabled(cls) -> bool:
        if os.environ.get("NO_COLOR", "").strip():
            return False
        if os.environ.get("FORCE_COLOR", "").strip() in ("1", "true", "yes"):
            return True
        return bool(getattr(sys.stdout, "isatty", lambda: False)())

    @classmethod
    def paint(cls, text: str, *codes: str) -> str:
        if not cls.enabled() or not codes:
            return text
        return "".join(codes) + text + cls.RESET


# ── badge / bar helpers ───────────────────────────────────────────────────────

def _progress_bar(n: int, total: int, width: int = 20) -> str:
    total = max(int(total or 0), 1)
    n = max(0, min(int(n or 0), total))
    filled = int(round(width * n / total))
    empty = width - filled
    return "█" * filled + "░" * empty


def _science_bar(n: int, hits: int, total: int, width: int = 20) -> str:
    """Dual-tone bar: green hit fraction inside docked fill, rest filled cyan/gray."""
    total = max(int(total or 0), 1)
    n = max(0, min(int(n or 0), total))
    hits = max(0, min(int(hits or 0), n))
    fill = int(round(width * n / total))
    hit_w = int(round(width * hits / total)) if n else 0
    hit_w = min(hit_w, fill)
    rest = fill - hit_w
    empty = width - fill
    return "▓" * hit_w + "█" * rest + "░" * empty


def _rate_badge(rate: Optional[float], *, good: float = 0.5, ok: float = 0.2) -> str:
    if rate is None:
        return "⬜ n/a"
    pct = 100.0 * rate
    if rate >= good:
        return f"🟢 {pct:.0f}%"
    if rate >= ok:
        return f"🟡 {pct:.0f}%"
    if rate > 0:
        return f"🟠 {pct:.0f}%"
    return f"🔴 {pct:.0f}%"


def _count_badge(hits: int, n: int, *, invert_zero_ok: bool = False) -> str:
    if n <= 0:
        return "⬜ —"
    if hits <= 0:
        return "🔴 0" if not invert_zero_ok else "🟢 0"
    if hits >= n:
        return f"🟢 {hits}/{n}"
    return f"🟡 {hits}/{n}"


def _rmsd_badge(x: Optional[float], *, unit: str = "Å") -> str:
    """Traffic-light RMSD: ≤2 green, ≤3 yellow, ≤5 orange, >5 red."""
    if x is None:
        return "—"
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "—"
    if v < 0:
        return f"⬛ {v:.2f}{unit}"
    s = f"{v:.2f}{unit}"
    if v <= 2.0:
        return f"🟢 {s}"
    if v <= 3.0:
        return f"🟡 {s}"
    if v <= 5.0:
        return f"🟠 {s}"
    return f"🔴 {s}"


def _cf_badge(x: Any) -> str:
    """CF scoring proxy badge — flag pathological magnitudes (not ΔG)."""
    cf = _ffloat(x)
    if cf is None:
        return "—"
    s = f"{cf:.1f}"
    if abs(cf) > 400:
        return f"🔴 {s} pathol"
    if abs(cf) > 200:
        return f"🟠 {s}"
    return f"⚪ {s}"


def _ram_badge(free_g: float, avail_g: float) -> str:
    if free_g < 0.5 or avail_g < 2.0:
        return f"🔴 CRITICAL free={free_g:.2f}G avail≈{avail_g:.2f}G"
    if free_g < 1.5 or avail_g < 4.0:
        return f"🟡 TIGHT free={free_g:.2f}G avail≈{avail_g:.2f}G"
    return f"🟢 OK free={free_g:.2f}G avail≈{avail_g:.2f}G"


def _bcr_buckets(results: List[Dict[str, Any]]) -> Dict[str, int]:
    """Histogram of BCR quality for denser science readout."""
    buckets = {"≤2": 0, "2–3": 0, "3–5": 0, ">5": 0, "neg/na": 0}
    for r in results or []:
        b = r.get("best_cluster_rmsd")
        if not isinstance(b, (int, float)):
            buckets["neg/na"] += 1
            continue
        if b < 0:
            buckets["neg/na"] += 1
        elif b <= 2.0:
            buckets["≤2"] += 1
        elif b <= 3.0:
            buckets["2–3"] += 1
        elif b <= 5.0:
            buckets["3–5"] += 1
        else:
            buckets[">5"] += 1
    return buckets


def _bucket_strip(buckets: Dict[str, int], n: int) -> str:
    if n <= 0:
        return "⬜ no data"
    parts = [
        f"🟢≤2:{buckets.get('≤2', 0)}",
        f"🟡2–3:{buckets.get('2–3', 0)}",
        f"🟠3–5:{buckets.get('3–5', 0)}",
        f"🔴>5:{buckets.get('>5', 0)}",
    ]
    if buckets.get("neg/na"):
        parts.append(f"⬛na:{buckets['neg/na']}")
    return " ".join(parts)


def _eta_str(
    results: List[Dict[str, Any]],
    n: int,
    total: int,
    wall_median_s: Optional[float] = None,
    wall_eta_remain_s: Optional[float] = None,
) -> str:
    """Rough remaining wall-time from finished target wall_time_s medians."""
    if n <= 0 or total <= 0:
        return "—"
    if n >= total:
        return "done"
    if wall_eta_remain_s is not None:
        return f"~{_fmt_wall(wall_eta_remain_s)}"
    med = wall_median_s
    if med is None:
        walls = sorted(
            float(r["wall_time_s"])
            for r in (results or [])
            if isinstance(r.get("wall_time_s"), (int, float)) and r["wall_time_s"] > 0
        )
        if len(walls) < 1:
            return "calc…"
        med = walls[len(walls) // 2]
    return f"~{_fmt_wall(med * (total - n))}"


def _proc_is_full85(p: Dict[str, Any]) -> bool:
    h = str(p.get("hint") or "")
    cmd = str(p.get("command") or "")
    return (
        "full85" in h
        or "full85" in cmd
        or "3dsig_full85" in cmd
        or "full85_serial" in h
    )


def _arm_hit(arm: str, p: Dict[str, Any]) -> bool:
    if not arm:
        return False
    h = str(p.get("hint") or "")
    cmd = str(p.get("command") or "")
    return (
        f"/{arm}/" in h
        or f"/{arm}/" in cmd
        or f" {arm} " in f" {cmd} "
        or f"arms {arm}" in cmd
        or f"--arms {arm}" in cmd
        or re.search(rf"(?:^|\s){re.escape(arm)}(?:\s|$)", cmd) is not None
        or h.endswith(f"/{arm}")
        or h.endswith(f"/{arm}/full85")
        or f"/{arm}/full85" in h
        or f"/bin/{arm}/" in cmd
    )


def _campaign_phase(c: Dict[str, Any], procs: List[Dict[str, Any]]) -> str:
    """RUNNING / PREP / DONE / EMPTY / FAIL-SCIENCE — panel-aware (no pilot false-RUN)."""
    cid = str(c.get("id") or "")
    n, tot = int(c.get("N") or 0), int(c.get("total") or 0)
    arm = str(c.get("arm") or (cid.split("/")[1] if "/" in cid else ""))
    panel = str(c.get("panel") or ("full85" if "full85" in cid else "pilot8"))
    is_full85 = "full85" in panel or "full85" in cid

    prep = False
    live = False
    full85_launcher_up = any(
        "full85" in str(p.get("hint") or "") for p in procs
    )
    for p in procs:
        h = str(p.get("hint") or "")
        cmd = str(p.get("command") or "")
        arm_hit = _arm_hit(arm, p)
        full85_cmd = _proc_is_full85(p)
        if is_full85 and "full85_serial" in h:
            live = True
        if arm_hit and is_full85 and full85_cmd:
            live = True
        if arm_hit and (not is_full85) and not full85_cmd:
            live = True
        if arm_hit and (
            "prep" in h or "generate_flexaid" in cmd or "ProcessLigand" in cmd
        ):
            # prep for full85 only attributes to full85 panel
            if is_full85 or not full85_launcher_up:
                prep = True
        if arm_hit and "FlexAID" in cmd:
            if is_full85:
                # attribute GA to full85 when launcher/work is full85, or ambiguous under arm
                if full85_cmd or full85_launcher_up or "3dsig_full85" in cmd:
                    live = True
            else:
                # pilot: ignore FlexAID owned by full85 campaign
                if not full85_cmd and not full85_launcher_up:
                    live = True

    # Completed dock count takes priority over mis-attributed live workers
    if tot > 0 and n >= tot and not (
        live and not is_full85 and any(
            _arm_hit(arm, p) and not _proc_is_full85(p) and "run_flexaid_arm" in str(p.get("command") or "")
            for p in procs
        )
    ):
        if live and is_full85:
            return "🔵 RUNNING"
        st = int(c.get("S_top10") or 0)
        pack = int(c.get("packaging_bug") or 0)
        if st == 0:
            if pack > 0 and pack >= n:
                return "🔴 DONE·PACK FAIL"
            return "🔴 DONE·SCI FAIL"
        return "🟢 DONE"

    if prep and n < tot:
        return "🟡 PREP"
    if live:
        return "🔵 RUNNING"
    if not c.get("exists") and n == 0:
        return "⬛ EMPTY"
    if n > 0:
        return "🟡 PARTIAL"
    return "⬜ PENDING"


def _phase_short(phase: str) -> str:
    """Compact phase token for pipeline strip."""
    if "RUNNING" in phase:
        return "RUN"
    if "PREP" in phase:
        return "PREP"
    if "PACK FAIL" in phase:
        return "PACK"
    if "SCI FAIL" in phase:
        return "FAIL"
    if "DONE" in phase:
        return "DONE"
    if "EMPTY" in phase:
        return "EMPTY"
    if "PARTIAL" in phase:
        return "PART"
    return "PEND"


def _fmt_A(x: Optional[float]) -> str:
    if x is None:
        return "—"
    return f"{x:.2f}Å"


def _fmt_wall(seconds: Optional[float]) -> str:
    """Human walltime: 90s / 12.5m / 3.2h."""
    if seconds is None:
        return "—"
    try:
        s = float(seconds)
    except (TypeError, ValueError):
        return "—"
    if s < 0:
        return "—"
    if s < 90:
        return f"{s:.0f}s"
    if s < 3600:
        return f"{s / 60:.1f}m"
    if s < 86400:
        return f"{s / 3600:.2f}h"
    return f"{s / 86400:.2f}d"


def _normalize_campaigns(
    campaigns: List[Dict[str, Any]],
) -> Dict[str, List[Dict[str, Any]]]:
    by_panel: Dict[str, List[Dict[str, Any]]] = {}
    for c in campaigns:
        panel = c.get("panel") or ("full85" if "full85" in c.get("id", "") else "pilot8")
        cc = dict(c)
        cc["panel"] = panel
        cc["arm"] = cc.get("arm") or cc["id"].split("/")[1]
        by_panel.setdefault(panel, []).append(cc)
    # full85 first (primary science), then pilot8, then anything else
    order = ["full85", "pilot8"]
    ordered: Dict[str, List[Dict[str, Any]]] = {}
    for p in order:
        if p in by_panel:
            ordered[p] = by_panel[p]
    for p, v in by_panel.items():
        if p not in ordered:
            ordered[p] = v
    return ordered


def _pipeline_strip(
    camps: List[Dict[str, Any]], procs: List[Dict[str, Any]], cur: Optional[str]
) -> List[str]:
    """A → B0 → B visual pipeline for one panel."""
    lines: List[str] = []
    bits: List[str] = []
    for c in camps:
        arm = c.get("arm") or "?"
        n, tot = int(c.get("N") or 0), int(c.get("total") or 1)
        phase = _campaign_phase(c, procs)
        short = _phase_short(phase)
        bar = _science_bar(n, int(c.get("S_top10") or 0), tot, width=12)
        pct = 100.0 * n / tot if tot else 0.0
        st = int(c.get("S_top10") or 0)
        emoji = phase.split()[0] if phase else "⬜"
        cur_mark = "◀" if (cur and "RUNNING" in phase) else ""
        bits.append(
            f"{emoji}**{arm}** `{bar}` {n}/{tot} ({pct:.0f}%) "
            f"S10={st} · {short}{cur_mark}"
        )
    lines.append("  " + "  →  ".join(bits))
    if cur:
        lines.append(f"  🎯 docking now: **`{cur}`**")
    return lines


def _headline_block(
    campaigns: List[Dict[str, Any]],
    procs: List[Dict[str, Any]],
    cur: Optional[str],
) -> List[str]:
    by_id = {c["id"]: c for c in campaigns}
    live_any = any(
        "three_engine" in p.get("hint", "")
        or "FlexAID" in p.get("command", "")
        or "full85" in p.get("command", "")
        or "generate_flexaid" in p.get("command", "")
        for p in procs
    )

    def _arm_docked(c: Optional[Dict[str, Any]]) -> bool:
        return bool(c and c.get("exists") and c.get("N", 0) >= c.get("total", 1))

    full85_active = any(_proc_is_full85(p) for p in procs)
    full85_partial = any(
        c.get("panel") == "full85" and int(c.get("N") or 0) > 0 for c in campaigns
    )
    # Resolve A/B0/B by panel+arm (campaign name is dynamic)
    def _arm(arm: str) -> Optional[Dict[str, Any]]:
        for c in campaigns:
            if c.get("panel") == "full85" and c.get("arm") == arm:
                return c
        return None

    fa, fb0, fb = _arm("A"), _arm("B0"), _arm("B")
    camp = _active_full85_campaign()
    lines: List[str] = []
    if full85_active or full85_partial or fa or fb0 or fb:
        na, nb0, nb = (
            int((fa or {}).get("N") or 0),
            int((fb0 or {}).get("N") or 0),
            int((fb or {}).get("N") or 0),
        )
        sta = int((fa or {}).get("S_top10") or 0)
        st0 = int((fb0 or {}).get("S_top10") or 0)
        st1 = int((fb or {}).get("S_top10") or 0)
        if live_any or full85_active:
            lines.append(
                f"> ### 🔵 NOT COMPLETE — full85 LIVE\n"
                f"> campaign=`{camp}` · **A** {na}/85 (S10={sta}) · **B0** {nb0}/85 "
                f"(S10={st0}) · **B** {nb}/85 (S10={st1}) · target=`{cur or 'prep/queue'}`"
            )
        elif _arm_docked(fb0) and _arm_docked(fb) and _arm_docked(fa):
            if st0 == 0 and st1 == 0 and sta == 0:
                lines.append(
                    f"> ### 🔴 DOCKING COMPLETE — SCIENCE GATE FAIL (full85)\n"
                    f"> `{camp}` S_top10 A={sta}/{na} B0={st0}/{nb0} B={st1}/{nb} · CF proxy ≠ ΔG"
                )
            else:
                lines.append(
                    f"> ### 🟢 DOCKING COMPLETE (full85)\n"
                    f"> `{camp}` S_top10 A={sta}/{na} B0={st0}/{nb0} B={st1}/{nb}"
                )
        else:
            lines.append(
                f"> ### 🟡 NOT COMPLETE — full85 partial\n"
                f"> `{camp}` A {na}/85 · B0 {nb0}/85 · B {nb}/85"
            )
    else:
        lines.append(
            f"> ### ⬛ IDLE — no full85 campaign activity\n"
            f"> expected campaign=`{camp}` · serial A→B0→B · Softβ OFF"
        )
    return lines


def format_ops_brief(
    mem: Dict[str, float],
    procs: List[Dict[str, Any]],
    campaigns: List[Dict[str, Any]],
    pid_files: Dict[str, Any],
    new_n: int,
    cur: Optional[str],
) -> str:
    free_g = float(mem.get("free_GB") or 0)
    avail_g = float(mem.get("available_est_GB") or 0)
    live_n = len(procs)
    lines: List[str] = []

    # ── hero header ───────────────────────────────────────────────────────
    lines.append("```")
    lines.append("╔══════════════════════════════════════════════════════════════════════╗")
    lines.append("║  🧬  FlexAIDdS  ·  BENCHMARK OPS MONITOR  ·  three_engine red-pair  ║")
    lines.append(f"║  ⏱  {utc_now():<62} ║")
    lines.append("╚══════════════════════════════════════════════════════════════════════╝")
    lines.append("```")
    lines.append("")

    # ── headline first (most important) ───────────────────────────────────
    lines.extend(_headline_block(campaigns, procs, cur))
    lines.append("")

    # ── resource strip ────────────────────────────────────────────────────
    heavy_flex = sum(
        1
        for p in procs
        if "FlexAID" in (p.get("command") or "")
        and "caffeinate" not in (p.get("command") or "")
        and float(p.get("rss_mb") or 0) >= 80
    )
    if live_n == 0:
        live_s = "⚪ idle · 0 workers"
    elif heavy_flex > 1:
        live_s = f"🔴 **DUAL-LAUNCH RISK** · {live_n} procs · **{heavy_flex} heavy FlexAID**"
    else:
        live_s = f"🔵 **{live_n}** proc(s) · heavy GA={heavy_flex} · serial OK"

    lines.append("### 📡 Systems")
    lines.append("")
    lines.append("| 🔋 RAM | ⚙️ Workers | 🎯 Target | 🆕 New finishes |")
    lines.append("|--------|------------|-----------|-----------------|")
    lines.append(
        f"| {_ram_badge(free_g, avail_g)} | {live_s} | "
        f"**`{cur or '—'}`** | "
        f"{'🆕 **' + str(new_n) + '**' if new_n else '0'} |"
    )
    lines.append("")

    if heavy_flex > 1:
        lines.append(
            "> ⚠️ **ALERT:** multiple heavy FlexAID processes — dual-launch is forbidden. "
            "Inspect live workers before starting anything."
        )
        lines.append("")

    by_panel = _normalize_campaigns(campaigns)

    # ── per-panel dashboards ──────────────────────────────────────────────
    for panel, camps in by_panel.items():
        if panel == "full85":
            title = f"🚀 Full Astex Diverse 85  ·  R=1  ·  **PRIMARY**  ·  `{_active_full85_campaign()}`"
            sub = "serial A → B0 → B · Softβ election OFF · CF scoring proxy ≠ ΔG · wall_s per target"
        else:
            title = f"📦 {panel}"
            sub = ""
        lines.append(f"## {title}")
        if sub:
            lines.append(f"_{sub}_")
        lines.append("")

        # pipeline strip
        lines.append("**Pipeline**")
        lines.append("")
        lines.extend(_pipeline_strip(camps, procs, cur if panel == "full85" else None))
        lines.append("")

        # main metrics table
        lines.append(
            "| Arm | Phase | Progress (▓=S10 hit) | S1 | **S_top10** | BCR≤2 | "
            "med BCR | best BCR | pack | ETA | stor |"
        )
        lines.append(
            "|:---:|:------|:---------------------|:--:|:-----------:|:-----:|"
            ":-------:|:--------:|:----:|:---:|:----:|"
        )
        for c in camps:
            n, tot = int(c.get("N") or 0), int(c.get("total") or 1)
            phase = _campaign_phase(c, procs)
            bar = _science_bar(n, int(c.get("S_top10") or 0), tot, width=18)
            pct = 100.0 * n / tot if tot else 0.0
            arm = c.get("arm") or "?"
            s1 = int(c.get("S1") or 0)
            st = int(c.get("S_top10") or 0)
            bcr = int(c.get("BCR_le2") or 0)
            pack = int(c.get("packaging_bug") or 0)
            stor = c.get("storage") or "?"
            eta = _eta_str(
                c.get("results") or [],
                n,
                tot,
                wall_median_s=c.get("wall_median_s"),
                wall_eta_remain_s=c.get("wall_eta_remain_s"),
            )
            lines.append(
                f"| **{arm}** | {phase} | `{bar}` **{n}/{tot}** ({pct:.0f}%) "
                f"| {_count_badge(s1, n)} "
                f"| {_count_badge(st, n)} "
                f"| {_count_badge(bcr, n)} "
                f"| {_rmsd_badge(c.get('BCR_median'))} "
                f"| {_rmsd_badge(c.get('BCR_best'))} "
                f"| {'🔴 **' + str(pack) + '**' if pack else '🟢 0'} "
                f"| {eta} "
                f"| {stor} |"
            )
        lines.append("")

        # computational walltime per method (primary ask for ops)
        lines.append(
            "**⏱ Compute walltime** (FlexAID CPU wall · not ΔG · "
            "src: `launcher` measured · `csv` · `mtime_proxy` pose span)"
        )
        lines.append("")
        lines.append(
            "| Arm | method | N timed | **Σ wall** | med / target | mean | min | max | "
            "proj remain | src |"
        )
        lines.append(
            "|:---:|--------|--------:|-----------:|-------------:|-----:|----:|----:|"
            "-----------:|-----|"
        )
        for c in camps:
            arm = c.get("arm") or "?"
            desc = (c.get("description") or "").split("·")[0].strip() or c.get("id", "")
            wn = int(c.get("wall_n") or 0)
            src = c.get("wall_src") or ("—" if wn == 0 else "?")
            lines.append(
                f"| **{arm}** | {desc} "
                f"| {wn}/{int(c.get('N') or 0)} "
                f"| **{_fmt_wall(c.get('wall_sum_s'))}** "
                f"| {_fmt_wall(c.get('wall_median_s'))} "
                f"| {_fmt_wall(c.get('wall_mean_s'))} "
                f"| {_fmt_wall(c.get('wall_min_s'))} "
                f"| {_fmt_wall(c.get('wall_max_s'))} "
                f"| {_fmt_wall(c.get('wall_eta_remain_s'))} "
                f"| `{src}` |"
            )
        lines.append("")

        # science KPI + BCR histogram (always visible — not collapsed)
        lines.append("**Science KPIs** (3Dsig: **S_top10** primary · S1 rank-0 · BCR cluster head)")
        lines.append("")
        for c in camps:
            arm = c.get("arm") or "?"
            n = int(c.get("N") or 0)
            desc = c.get("description") or ""
            if n <= 0:
                lines.append(
                    f"- **{arm}** `{c['id']}` — ⬛ no `result.csv` yet"
                    + (f" · _{desc}_" if desc else "")
                )
                continue
            buckets = _bcr_buckets(c.get("results") or [])
            lines.append(
                f"- **{arm}** `{c['id']}` · "
                f"S1 {_rate_badge(c.get('S1_rate'))} · "
                f"**S_top10 {_rate_badge(c.get('S_top10_rate'))}** · "
                f"BCR {_rate_badge(c.get('BCR_rate'))} · "
                f"Σwall {_fmt_wall(c.get('wall_sum_s'))} · "
                f"med {_fmt_wall(c.get('wall_median_s'))}/tgt · "
                f"gap={c.get('election_gap', 0)} neg1={c.get('bcr_neg1', 0)}"
            )
            lines.append(f"  - BCR hist: {_bucket_strip(buckets, n)}")
            if desc:
                lines.append(f"  - _{desc}_")
        lines.append("")

    # ── PID locks (live-first sort) ────────────────────────────────────────
    lines.append("## 🔒 PID / lock files")
    lines.append("")
    lines.append("| Lock | Status | PID |")
    lines.append("|------|:------:|-----|")
    # live locks first for scanability
    items = sorted(
        pid_files.items(),
        key=lambda kv: (0 if kv[1].get("live") else 1 if kv[1].get("exists") else 2, kv[0]),
    )
    for k, v in items:
        live = bool(v.get("live"))
        if live:
            badge = "🟢 **LIVE**"
        elif v.get("exists"):
            badge = "⚪ dead"
        else:
            badge = "⬛ none"
        pid = v.get("pid")
        lines.append(f"| `{k}` | {badge} | `{pid if pid is not None else '—'}` |")
    lines.append("")

    # ── live workers ──────────────────────────────────────────────────────
    lines.append("## ⚙️ Live workers")
    lines.append("")
    if not procs:
        lines.append("_No matching FlexAID / chain / prep processes._")
    else:
        lines.append("| PID | Role | CPU | RSS | State | Etime | Command |")
        lines.append("|-----|------|:---:|:---:|:-----:|------:|---------|")
        for p in procs[:14]:
            cpu = float(p.get("pcpu") or 0)
            rss = float(p.get("rss_mb") or 0)
            if cpu > 80:
                cpu_b = f"🔴 **{cpu:.0f}%**"
            elif cpu > 20:
                cpu_b = f"🟡 {cpu:.0f}%"
            else:
                cpu_b = f"🟢 {cpu:.0f}%"
            if rss > 1500:
                rss_b = f"🔴 **{rss:.0f}MB**"
            elif rss > 400:
                rss_b = f"🟡 {rss:.0f}MB"
            else:
                rss_b = f"🟢 {rss:.0f}MB"
            state = str(p.get("state") or "")
            state_b = f"🔵 {state}" if state.startswith("R") else f"⚪ {state}"
            hint = str(p.get("hint") or "")
            if "full85" in hint:
                role = f"🚀 `{hint}`"
            elif "prep" in hint or "ProcessLigand" in hint:
                role = f"🟡 `{hint}`"
            elif "FlexAID" in (p.get("command") or ""):
                role = f"⚙️ `{hint}`"
            else:
                role = f"`{hint}`"
            cmd = (p.get("command") or "")[:88].replace("|", "\\|")
            lines.append(
                f"| `{p['pid']}` | {role} | {cpu_b} | {rss_b} | {state_b} "
                f"| `{p.get('etime')}` | `{cmd}` |"
            )
        if len(procs) > 14:
            lines.append(f"| … | _{len(procs) - 14} more_ | | | | | |")
    lines.append("")

    # ── footer legend ─────────────────────────────────────────────────────
    lines.append("---")
    lines.append("")
    lines.append(
        "**Legend** · S1 = rank-0 ≤2 Å · **S_top10** = 3Dsig primary (ranks 0–9) · "
        "BCR = best cluster RMSD · CF = contact-function **scoring proxy** (not ΔG).  \n"
        "Progress bar: `▓` = S_top10 hits · `█` = docked without S10 · `░` = remaining.  \n"
        "RMSD: 🟢≤2 · 🟡≤3 · 🟠≤5 · 🔴>5 Å.  "
        "Phase: 🔵 running · 🟡 partial/prep · 🟢 done · 🔴 sci/pack fail · ⬛ empty.  \n"
        "Rules: **no dual-launch** · Softβ election default **OFF** · scope **three_engine only**."
    )
    lines.append("")
    return "\n".join(lines)


def format_new_targets_md(new_items: List[Dict[str, Any]]) -> str:
    """Compact color-coded table of newly finished targets (no per-target H2 spam)."""
    if not new_items:
        return ""
    lines: List[str] = []

    # summary counters
    n_hit = sum(1 for it in new_items if it.get("s1") or it.get("s_top10"))
    n_bcr = sum(1 for it in new_items if it.get("bcr_le2"))
    n_pack = sum(1 for it in new_items if "PACKAGING_BUG" in (it.get("tags") or []))
    n_path = sum(1 for it in new_items if "PATHOLOGICAL_CF" in (it.get("tags") or []))
    n_miss = len(new_items) - n_hit

    lines.append(f"## 🆕 New finishes · **{len(new_items)}** targets")
    lines.append("")
    lines.append(
        f"| 🟢 S1/S10 hits | 🔴 misses | 🟢 BCR≤2 | 🔴 pack bugs | 🔴 pathol CF |"
    )
    lines.append("|:--------------:|:---------:|:--------:|:------------:|:------------:|")
    lines.append(
        f"| **{n_hit}** | **{n_miss}** | **{n_bcr}** | **{n_pack}** | **{n_path}** |"
    )
    lines.append("")

    # sort: hits first, then by BCR ascending, then campaign/pdb
    def _sort_key(it: Dict[str, Any]) -> Tuple:
        hit = 0 if (it.get("s1") or it.get("s_top10") or it.get("bcr_le2")) else 1
        bcr = it.get("best_cluster_rmsd")
        bcr_k = float(bcr) if isinstance(bcr, (int, float)) and bcr >= 0 else 999.0
        return (hit, bcr_k, str(it.get("campaign") or ""), str(it.get("pdb_id") or ""))

    ordered = sorted(new_items, key=_sort_key)

    lines.append(
        "| Status | Target | Campaign | top1 | BCR | CF proxy | poses | wall | tags |"
    )
    lines.append(
        "|:------:|:------:|----------|-----:|----:|---------:|------:|-----:|------|"
    )
    for it in ordered:
        tags = it.get("tags") or []
        s1 = bool(it.get("s1"))
        s10 = bool(it.get("s_top10"))
        bcr_ok = bool(it.get("bcr_le2"))
        if s1 or s10 or bcr_ok:
            status = "🟢 HIT"
        elif "PACKAGING_BUG" in tags:
            status = "⬛ PACK"
        else:
            status = "🔴 MISS"

        bcr = it.get("best_cluster_rmsd")
        rh = it.get("rmsd_hungarian")
        tag_short = []
        for t in tags:
            if t in ("S1_HIT", "S2_HIT", "BCR_OK"):
                tag_short.append(f"🟢{t}")
            elif "MISS" in t or "BUG" in t or "PATHOL" in t or "SENTINEL" in t:
                tag_short.append(f"🔴{t}")
            elif t in ("ELECTION_GAP", "SEED_ECHO", "NATIVE_SEEDED"):
                tag_short.append(f"🟡{t}")
        # keep table tight — top 3 tags
        tag_s = " ".join(f"`{t}`" for t in tag_short[:3])
        if len(tag_short) > 3:
            tag_s += " …"

        wall = it.get("wall_time_s")
        if isinstance(wall, (int, float)) and wall > 0:
            wall_s = f"{wall:.0f}s" if wall < 3600 else f"{wall/3600:.1f}h"
        else:
            wall_s = "—"

        poses = it.get("num_poses")
        poses_s = str(poses) if poses is not None else "—"
        camp = str(it.get("campaign") or "")
        # shorten campaign path for readability
        camp_short = camp.replace("three_engine/", "te/")

        lines.append(
            f"| {status} | **`{it.get('pdb_id')}`** | `{camp_short}` "
            f"| {_rmsd_badge(rh if isinstance(rh, (int, float)) else None)} "
            f"| {_rmsd_badge(bcr if isinstance(bcr, (int, float)) else None)} "
            f"| {_cf_badge(it.get('best_score'))} "
            f"| {poses_s} | {wall_s} | {tag_s} |"
        )

    lines.append("")
    # detail footnotes only for hits or unusual tags (keep noise down)
    notables = [
        it
        for it in ordered
        if it.get("s1")
        or it.get("s_top10")
        or it.get("bcr_le2")
        or "ELECTION_GAP" in (it.get("tags") or [])
        or "SEED_ECHO" in (it.get("tags") or [])
    ]
    if notables:
        lines.append("<details><summary>Hit / notable detail (dock_config)</summary>")
        lines.append("")
        for it in notables:
            b = it.get("budget") or {}
            bud = ""
            if b:
                bud = (
                    f" · chroms={b.get('num_chromosomes')} gen={b.get('num_generations')} "
                    f"T={b.get('temperature')} clust={b.get('clustering')}"
                )
            lines.append(
                f"- **`{it.get('pdb_id')}`** `{it.get('campaign')}` · "
                f"S1={it.get('s1')} S10={it.get('s_top10')} PB={it.get('pb_pass')} "
                f"seed_echo={it.get('seed_echo')} native={it.get('native_pose_seeded')}"
                f"{bud}"
            )
        lines.append("")
        lines.append("</details>")
        lines.append("")
    return "\n".join(lines)


def colorize_brief_for_tty(text: str) -> str:
    """Apply rich ANSI highlighting for interactive terminals (files stay plain MD)."""
    if not _Ansi.enabled():
        return text
    A = _Ansi
    out = text
    replacements = [
        ("### 🔵 NOT COMPLETE — full85 LIVE", A.paint("### 🔵 NOT COMPLETE — full85 LIVE", A.BOLD, A.BRIGHT_CYAN)),
        ("### 🔵 NOT COMPLETE — three_engine LIVE", A.paint("### 🔵 NOT COMPLETE — three_engine LIVE", A.BOLD, A.BRIGHT_CYAN)),
        ("### 🟢 DOCKING COMPLETE (full85)", A.paint("### 🟢 DOCKING COMPLETE (full85)", A.BOLD, A.BRIGHT_GREEN)),
        ("### 🟢 DOCKING COMPLETE (pilot8)", A.paint("### 🟢 DOCKING COMPLETE (pilot8)", A.BOLD, A.BRIGHT_GREEN)),
        (
            "### 🔴 DOCKING COMPLETE — SCIENCE GATE FAIL (full85)",
            A.paint("### 🔴 DOCKING COMPLETE — SCIENCE GATE FAIL (full85)", A.BOLD, A.BRIGHT_RED),
        ),
        (
            "### 🔴 DOCKING COMPLETE — SCIENCE GATE FAIL (pilot8)",
            A.paint("### 🔴 DOCKING COMPLETE — SCIENCE GATE FAIL (pilot8)", A.BOLD, A.BRIGHT_RED),
        ),
        ("### 🟡 NOT COMPLETE — full85 partial", A.paint("### 🟡 NOT COMPLETE — full85 partial", A.BOLD, A.BRIGHT_YELLOW)),
        ("### 🟡 NOT COMPLETE — pilot8", A.paint("### 🟡 NOT COMPLETE — pilot8", A.BOLD, A.BRIGHT_YELLOW)),
        ("🔵 RUNNING", A.paint("🔵 RUNNING", A.BOLD, A.CYAN)),
        ("🔴 DONE·SCI FAIL", A.paint("🔴 DONE·SCI FAIL", A.BOLD, A.RED)),
        ("🔴 DONE·PACK FAIL", A.paint("🔴 DONE·PACK FAIL", A.BOLD, A.RED)),
        ("🟢 DONE", A.paint("🟢 DONE", A.BOLD, A.GREEN)),
        ("🟡 PREP", A.paint("🟡 PREP", A.YELLOW)),
        ("🟡 PARTIAL", A.paint("🟡 PARTIAL", A.YELLOW)),
        ("⬛ EMPTY", A.paint("⬛ EMPTY", A.DIM)),
        ("⬜ PENDING", A.paint("⬜ PENDING", A.DIM)),
        ("🟢 **LIVE**", A.paint("🟢 **LIVE**", A.BOLD, A.GREEN)),
        ("🔴 **DUAL-LAUNCH RISK**", A.paint("🔴 **DUAL-LAUNCH RISK**", A.BOLD, A.BG_RED, A.BRIGHT_WHITE)),
        ("🔴 CRITICAL", A.paint("🔴 CRITICAL", A.BOLD, A.RED)),
        ("🟡 TIGHT", A.paint("🟡 TIGHT", A.BOLD, A.YELLOW)),
        ("🟢 OK", A.paint("🟢 OK", A.GREEN)),
        ("🟢 HIT", A.paint("🟢 HIT", A.BOLD, A.GREEN)),
        ("🔴 MISS", A.paint("🔴 MISS", A.RED)),
        ("⬛ PACK", A.paint("⬛ PACK", A.DIM)),
        ("**PRIMARY**", A.paint("**PRIMARY**", A.BOLD, A.MAGENTA)),
    ]
    for old, new in replacements:
        out = out.replace(old, new)
    # legacy headline tokens (if any remain)
    out = out.replace(
        "**🔵 NOT COMPLETE",
        A.paint("**🔵 NOT COMPLETE", A.BOLD, A.CYAN),
    )
    out = out.replace(
        "**🟢 DOCKING COMPLETE",
        A.paint("**🟢 DOCKING COMPLETE", A.BOLD, A.GREEN),
    )
    out = out.replace(
        "**🔴 DOCKING COMPLETE — SCIENCE GATE FAIL",
        A.paint("**🔴 DOCKING COMPLETE — SCIENCE GATE FAIL", A.BOLD, A.RED),
    )
    out = out.replace(
        "**🟡 NOT COMPLETE",
        A.paint("**🟡 NOT COMPLETE", A.BOLD, A.YELLOW),
    )
    return out


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
    te_logs = local_root() / "logs" / "three_engine"

    mem = mem_snapshot()
    procs = list_procs()
    cur = current_dock_target(procs)

    def _pid_prefer(*paths: Path) -> Dict[str, Any]:
        for p in paths:
            info = pid_file(p)
            if info.get("exists"):
                return info
        return {"exists": False, "live": False, "pid": None}

    # PID keys use path-accurate campaign names (no C0_claim/C0_legacy)
    pid_files = {
        "three_engine/serial": _pid_prefer(
            te_logs / "run_3dsig_red_pair_serial.pid",
            q / "logs/run_3dsig_red_pair_serial.pid",
            q / "logs/run_AB_pilot8_chain.pid",
        ),
        "three_engine/full85_serial": _pid_prefer(
            te_logs / "run_3dsig_red_pair_full85.pid",
            q / "logs/run_3dsig_red_pair_full85.pid",
        ),
        "three_engine/A/pilot8": _pid_prefer(
            te_logs / "run_A_pilot8.lock", te_logs / "run_A_pilot8.pid"
        ),
        "three_engine/B0/pilot8": _pid_prefer(
            te_logs / "run_B0_pilot8.lock", te_logs / "run_B0_pilot8.pid"
        ),
        "three_engine/B/pilot8": _pid_prefer(
            te_logs / "run_B_pilot8.lock",
            te_logs / "run_B_pilot8_launcher.pid",
            te_logs / "run_B_pilot8.pid",
        ),
        "three_engine/A/full85": _pid_prefer(
            te_logs / "run_A_full85.lock", te_logs / "run_A_full85.pid"
        ),
        "three_engine/B0/full85": _pid_prefer(
            te_logs / "run_B0_full85.lock", te_logs / "run_B0_full85.pid"
        ),
        "three_engine/B/full85": _pid_prefer(
            te_logs / "run_B_full85.lock", te_logs / "run_B_full85.pid"
        ),
        "throughput_maximizer": _pid_prefer(q / "logs/throughput_maximizer.pid"),
    }

    # Rebind specs at runtime (FLEXAID_CAMPAIGN / newest full85 dir)
    campaign_specs = _full85_campaign_specs()
    campaigns: List[Dict[str, Any]] = []
    for spec in campaign_specs:
        # rel is under campaigns/, e.g. three_engine/A/3dsig_full85_scratch_*
        rel = Path(spec["rel"])
        local_root_c = local_camps / rel
        icloud_root_c = res / "campaigns" / rel
        # Also allow FLEXAIDDS_RESULTS already pointing at .../campaigns
        alt_icloud = res / rel
        try:
            if local_root_c.is_dir():
                root = local_root_c
            elif icloud_root_c.is_dir():
                root = icloud_root_c
            elif alt_icloud.is_dir():
                root = alt_icloud
            else:
                root = local_root_c  # missing → scan reports empty
        except OSError:
            root = local_root_c
        camp = scan_campaign(spec["id"], root, spec["total"])
        camp["description"] = spec.get("description", "")
        camp["panel"] = spec.get("panel", "")
        camp["arm"] = spec.get("arm", "")
        campaigns.append(camp)
    # Do NOT auto-discover C0_* or other trees — ops is three_engine red-pair only.
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
        # ANSI only on TTY/FORCE_COLOR; on-disk MD stays portable emoji markdown
        print(colorize_brief_for_tty(brief))

    heavy = sum(
        1
        for p in procs
        if p["rss_mb"] >= 200 and "FlexAID" in p["command"]
    )
    return 2 if heavy > 1 else 0


if __name__ == "__main__":
    sys.exit(main())
