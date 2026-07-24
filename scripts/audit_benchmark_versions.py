#!/usr/bin/env python3
"""Inventory benchmark/canary campaigns and report success rates + provenance.

Discovers docking result trees under local (and optional) roots, reads
RUN_RECEIPT / result.csv / canary FINAL or summary JSON, and emits a versioned
audit table (JSON + Markdown).

Success definitions (fail-closed, aligned with claim spirit):
  S1_rmsd        rmsd_to_crystal ≤ 2.0 Å when column present (diagnostic)
  genuine        seed_echo==0 AND pose_source not seed/ini AND
                 (rmsd_hungarian ≤ 2.0 if present else rmsd_to_crystal ≤ 2.0)
  claim_ready    claim_ready flag truthy when column present (optional)

Does NOT walk CloudDocs / iCloud Mobile Documents (AGENTS local-first rule).
Never uses Path.rglob under iCloud paths.

Usage:
  python3 scripts/audit_benchmark_versions.py
  python3 scripts/audit_benchmark_versions.py --roots ~/flexaidds_results \\
      --out-dir ~/flexaidds_results/workorders
  python3 scripts/audit_benchmark_versions.py --json-only

Exit: 0 always when inventory written; 2 on usage error.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Optional

RMSD_OK = 2.0
SEED_POSE_SOURCES = frozenset(
    {
        "ini_elitism",
        "seed_echo",
        "native_seed",
        "oracle_seed",
        "seed",
    }
)


def _default_roots() -> list[Path]:
    roots: list[Path] = []
    env = os.environ.get("FLEXAIDDS_LOCAL_ROOT", "").strip()
    if env:
        roots.append(Path(env).expanduser())
    else:
        roots.append(Path.home() / "flexaidds_results")
    # Common local campaign parents (never CloudDocs)
    extra = Path.home() / "flexaidds_results"
    if extra not in roots:
        roots.append(extra)
    return roots


def _f(row: dict[str, str], *keys: str) -> float:
    for k in keys:
        v = row.get(k)
        if v is None or v == "" or str(v).upper() == "NA":
            continue
        try:
            x = float(v)
            if math.isfinite(x):
                return x
        except (TypeError, ValueError):
            continue
    return float("nan")


def _flag0(row: dict[str, str], key: str) -> bool:
    return str(row.get(key, "")).strip() in ("0", "0.0", "False", "false", "")


def _truth(row: dict[str, str], key: str) -> bool:
    return str(row.get(key, "")).strip() in ("1", "True", "true", "YES", "yes")


def _is_seed_pose(row: dict[str, str]) -> bool:
    src = str(row.get("pose_source", "")).strip().lower()
    if src in SEED_POSE_SOURCES:
        return True
    if "seed" in src and "ga" not in src:
        return True
    return not _flag0(row, "seed_echo") and str(row.get("seed_echo", "")).strip() not in (
        "",
        "0",
        "0.0",
    )


def row_is_genuine(row: dict[str, str]) -> bool:
    """Genuine docking success for audit (seed_echo=0, non-seed source, RMSD≤2)."""
    se = str(row.get("seed_echo", "")).strip()
    if se not in ("0", "0.0", ""):
        # blank seed_echo: only accept if pose_source clearly ga
        if se != "" and se not in ("0", "0.0"):
            try:
                if float(se) != 0.0:
                    return False
            except ValueError:
                return False
    src = str(row.get("pose_source", "")).strip().lower()
    if src in SEED_POSE_SOURCES:
        return False
    # Prefer hungarian when present; fall back to ordered crystal RMSD
    rh = _f(row, "rmsd_hungarian")
    rc = _f(row, "rmsd_to_crystal")
    rmsd = rh if math.isfinite(rh) else rc
    if not math.isfinite(rmsd):
        return False
    return rmsd <= RMSD_OK


def row_s1(row: dict[str, str]) -> bool:
    """RMSD-only diagnostic (ordered crystal RMSD)."""
    rc = _f(row, "rmsd_to_crystal")
    if not math.isfinite(rc):
        rc = _f(row, "rmsd_hungarian")
    return math.isfinite(rc) and rc <= RMSD_OK


@dataclass
class TargetRow:
    pdb: str
    genuine: bool
    s1: bool
    seed_echo: str = ""
    pose_source: str = ""
    rmsd_hun: float = float("nan")
    rmsd_xtal: float = float("nan")
    best_cluster: float = float("nan")
    elected_cf: float = float("nan")
    cf_native: float = float("nan")
    claim_ready: bool = False


@dataclass
class CampaignAudit:
    name: str
    path: str
    kind: str  # campaign | canary | unknown
    n_targets: int = 0
    n_genuine: int = 0
    n_s1: int = 0
    genuine_rate: float = float("nan")
    s1_rate: float = float("nan")
    matrix_md5: str = ""
    git_commit: str = ""
    seed_elitism: str = ""
    mode: str = ""
    started_utc: str = ""
    binary_sha256: str = ""
    notes: str = ""
    targets: list[dict[str, Any]] = field(default_factory=list)
    provenance: dict[str, Any] = field(default_factory=dict)


def _read_receipt(camp: Path) -> dict[str, Any]:
    for name in ("RUN_RECEIPT.json", "out/RUN_RECEIPT.json", "provenance.json"):
        p = camp / name
        if p.is_file():
            try:
                return json.loads(p.read_text(encoding="utf-8", errors="replace"))
            except (OSError, ValueError):
                continue
    # canary provenance.txt
    pt = camp / "provenance.txt"
    if pt.is_file():
        d: dict[str, Any] = {}
        for line in pt.read_text(encoding="utf-8", errors="replace").splitlines():
            if " " in line:
                k, _, v = line.partition(" ")
                d[k.strip()] = v.strip()
        return d
    return {}


def _iter_result_csvs(camp: Path) -> Iterable[Path]:
    # Shallow: camp/*/result.csv or camp/out/*/result.csv
    for base in (camp, camp / "out"):
        if not base.is_dir():
            continue
        try:
            for child in sorted(base.iterdir()):
                if not child.is_dir():
                    continue
                r = child / "result.csv"
                if r.is_file():
                    yield r
        except OSError:
            continue
    # top-level summary CSVs
    for name in (
        "astex_diverse_results.csv",
        "astex_diverse_summary.csv",
        "summary.csv",
        "results.csv",
    ):
        p = camp / name
        if p.is_file():
            yield p
        p2 = camp / "out" / name
        if p2.is_file():
            yield p2


def _load_csv_rows(path: Path) -> list[dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8", errors="replace") as fh:
            return list(csv.DictReader(fh))
    except OSError:
        return []


def audit_campaign_dir(camp: Path) -> CampaignAudit:
    name = camp.name
    kind = "unknown"
    low = name.lower()
    if "canary" in low:
        kind = "canary"
    elif low.startswith("v_"):
        kind = "ab_arm"
    elif "_blind_" in low:
        kind = "probe"
    elif "campaign" in low or low.startswith("c0_") or "full85" in low:
        kind = "campaign"
    elif camp.parent.name in ("campaigns", "canary"):
        kind = "campaign" if camp.parent.name == "campaigns" else "canary"

    receipt = _read_receipt(camp)
    matrix = str(receipt.get("matrix_md5") or receipt.get("matrix") or "")
    commit = str(receipt.get("git_commit") or receipt.get("commit") or "")
    if len(matrix) > 12 and matrix[0:8].isalnum():
        matrix_short = matrix[:8]
    else:
        matrix_short = matrix

    # Prefer canary summary JSON when present
    for sj in (
        camp / "summary.json",
        camp / "summary_auto.json",
        camp / "canary_summary.json",
        camp / "base_summary.json",
    ):
        if sj.is_file():
            try:
                data = json.loads(sj.read_text(encoding="utf-8", errors="replace"))
            except (OSError, ValueError):
                data = {}
            rows_j = data.get("rows") or data.get("targets") or []
            if isinstance(rows_j, list) and rows_j:
                targets: list[TargetRow] = []
                for r in rows_j:
                    if not isinstance(r, dict):
                        continue
                    pdb = str(r.get("pdb") or r.get("pdb_id") or "")
                    # coerce to string map for helpers
                    smap = {k: str(v) if v is not None else "" for k, v in r.items()}
                    if "seed_echo" not in smap and r.get("genuine") is False:
                        smap.setdefault("seed_echo", "0")
                    if "rmsd_hungarian" not in smap and "rmsd_hun" in r:
                        smap["rmsd_hungarian"] = str(r.get("rmsd_hun", ""))
                    if "best_cluster_rmsd" not in smap and "best_cluster" in r:
                        smap["best_cluster_rmsd"] = str(r.get("best_cluster", ""))
                    g = bool(r.get("genuine")) if "genuine" in r else row_is_genuine(smap)
                    targets.append(
                        TargetRow(
                            pdb=pdb,
                            genuine=g,
                            s1=row_s1(smap),
                            seed_echo=smap.get("seed_echo", ""),
                            pose_source=smap.get("pose_source", ""),
                            rmsd_hun=_f(smap, "rmsd_hungarian", "rmsd_hun"),
                            rmsd_xtal=_f(smap, "rmsd_to_crystal"),
                            best_cluster=_f(smap, "best_cluster_rmsd", "best_cluster"),
                            elected_cf=_f(smap, "elected_cf"),
                            cf_native=_f(smap, "cf_native"),
                            claim_ready=_truth(smap, "claim_ready"),
                        )
                    )
                if targets:
                    ng = sum(1 for t in targets if t.genuine)
                    ns = sum(1 for t in targets if t.s1)
                    n = len(targets)
                    return CampaignAudit(
                        name=name,
                        path=str(camp),
                        kind=kind,
                        n_targets=n,
                        n_genuine=ng,
                        n_s1=ns,
                        genuine_rate=(ng / n) if n else float("nan"),
                        s1_rate=(ns / n) if n else float("nan"),
                        matrix_md5=matrix_short or matrix,
                        git_commit=commit[:12] if commit else "",
                        seed_elitism=str(receipt.get("seed_elitism", "")),
                        mode=str(receipt.get("mode", "")),
                        started_utc=str(receipt.get("started_utc", "")),
                        binary_sha256=str(receipt.get("binary_sha256", ""))[:16],
                        notes=f"from {sj.name}",
                        targets=[asdict(t) for t in targets],
                        provenance={
                            k: receipt.get(k)
                            for k in (
                                "matrix_md5",
                                "git_commit",
                                "seed_elitism",
                                "mode",
                                "pop",
                                "gen",
                                "restarts",
                            )
                            if k in receipt
                        },
                    )

    # CSV path
    seen_pdb: set[str] = set()
    targets2: list[TargetRow] = []
    for csv_path in _iter_result_csvs(camp):
        for row in _load_csv_rows(csv_path):
            pdb = str(row.get("pdb_id") or row.get("pdb") or row.get("code") or "").strip()
            if not pdb:
                # summary row without pdb — skip or treat as aggregate later
                continue
            if pdb in seen_pdb:
                continue
            seen_pdb.add(pdb)
            targets2.append(
                TargetRow(
                    pdb=pdb,
                    genuine=row_is_genuine(row),
                    s1=row_s1(row),
                    seed_echo=str(row.get("seed_echo", "")),
                    pose_source=str(row.get("pose_source", "")),
                    rmsd_hun=_f(row, "rmsd_hungarian"),
                    rmsd_xtal=_f(row, "rmsd_to_crystal"),
                    best_cluster=_f(row, "best_cluster_rmsd"),
                    elected_cf=_f(row, "elected_cf"),
                    cf_native=_f(row, "cf_native"),
                    claim_ready=_truth(row, "claim_ready"),
                )
            )

    # FINAL.md genuine line for canaries without CSV
    notes = ""
    final = camp / "FINAL.md"
    if final.is_file() and not targets2:
        text = final.read_text(encoding="utf-8", errors="replace")
        m = re.search(r"Genuine:\s*(\d+)\s*/\s*(\d+)", text, re.I)
        if m:
            ng, n = int(m.group(1)), int(m.group(2))
            return CampaignAudit(
                name=name,
                path=str(camp),
                kind=kind,
                n_targets=n,
                n_genuine=ng,
                n_s1=ng,
                genuine_rate=(ng / n) if n else float("nan"),
                s1_rate=(ng / n) if n else float("nan"),
                matrix_md5=matrix_short or matrix,
                git_commit=commit[:12] if commit else "",
                seed_elitism=str(receipt.get("seed_elitism", "")),
                mode=str(receipt.get("mode", "")),
                notes="parsed from FINAL.md only",
            )
        notes = "FINAL.md present but no parseable Genuine line / no result.csv"

    n = len(targets2)
    ng = sum(1 for t in targets2 if t.genuine)
    ns = sum(1 for t in targets2 if t.s1)
    return CampaignAudit(
        name=name,
        path=str(camp),
        kind=kind,
        n_targets=n,
        n_genuine=ng,
        n_s1=ns,
        genuine_rate=(ng / n) if n else float("nan"),
        s1_rate=(ns / n) if n else float("nan"),
        matrix_md5=matrix_short or matrix,
        git_commit=commit[:12] if commit else "",
        seed_elitism=str(receipt.get("seed_elitism", "")),
        mode=str(receipt.get("mode", "")),
        started_utc=str(receipt.get("started_utc", "")),
        binary_sha256=str(receipt.get("binary_sha256", ""))[:16],
        notes=notes,
        targets=[asdict(t) for t in targets2],
        provenance={
            k: receipt.get(k)
            for k in (
                "matrix_md5",
                "git_commit",
                "seed_elitism",
                "mode",
                "pop",
                "gen",
                "restarts",
            )
            if k in receipt
        },
    )


# Campaign / canary name patterns (local APFS only; never CloudDocs).
_CAMPAIGN_PREFIXES = (
    "canary_",
    "C0_",
    "c0_",
    "v_",  # A/B version arms e.g. v_comcap_softbeta_*
    "campaign_",
)
_CAMPAIGN_SUBSTR = (
    "full85",
    "_blind_",
    "astex",
    "claim_",
)


def _looks_like_campaign(name: str) -> bool:
    if name.startswith(".") or name in {
        "workorders",
        "logs",
        "data_9dc9",
        "matrix_pins",
        "pins",
        "probe_cache",
        "three_engine_entropy_q1",
        "baseline_engine",
        "hem_stripped_receptors",
    }:
        return False
    if any(name.startswith(p) for p in _CAMPAIGN_PREFIXES):
        return True
    low = name.lower()
    return any(s in low for s in _CAMPAIGN_SUBSTR)


def discover_campaigns(roots: list[Path]) -> list[Path]:
    found: list[Path] = []
    for root in roots:
        root = root.expanduser().resolve()
        if not root.is_dir():
            continue
        # Never walk under Mobile Documents / iCloud FileProvider
        root_s = str(root)
        if "Mobile Documents" in root_s or "CloudDocs" in root_s:
            continue
        try:
            for child in sorted(root.iterdir()):
                if not child.is_dir():
                    continue
                # Skip iCloud archive symlinks that land under CloudDocs
                try:
                    resolved = child.resolve()
                    rs = str(resolved)
                    if "Mobile Documents" in rs or "CloudDocs" in rs:
                        continue
                except OSError:
                    continue
                n = child.name
                if n == "campaigns":
                    try:
                        for c2 in sorted(child.iterdir()):
                            if c2.is_dir() and not c2.name.startswith("."):
                                try:
                                    r2 = str(c2.resolve())
                                    if "Mobile Documents" in r2 or "CloudDocs" in r2:
                                        continue
                                except OSError:
                                    continue
                                found.append(c2)
                    except OSError:
                        pass
                    continue
                if _looks_like_campaign(n):
                    found.append(child)
        except OSError:
            continue
    uniq: dict[str, Path] = {}
    for p in found:
        try:
            uniq[str(p.resolve())] = p
        except OSError:
            uniq[str(p)] = p
    return list(uniq.values())


def git_log_for_commit(repo: Path, commit: str, n: int = 5) -> str:
    if not commit or len(commit) < 7:
        return ""
    try:
        out = subprocess.run(
            ["git", "-C", str(repo), "log", "--oneline", f"-{n}", commit],
            capture_output=True,
            text=True,
            timeout=15,
            check=False,
        )
        return out.stdout.strip() if out.returncode == 0 else ""
    except (OSError, subprocess.TimeoutExpired):
        return ""


def git_log_range(repo: Path, older: str, newer: str, n: int = 40) -> list[str]:
    """Commits reachable from newer but not older (changelog between pins)."""
    if not older or not newer or len(older) < 7 or len(newer) < 7:
        return []
    if older.startswith(newer) or newer.startswith(older):
        return []
    try:
        out = subprocess.run(
            [
                "git",
                "-C",
                str(repo),
                "log",
                "--oneline",
                f"--max-count={n}",
                f"{older}..{newer}",
            ],
            capture_output=True,
            text=True,
            timeout=20,
            check=False,
        )
        if out.returncode != 0:
            return []
        return [ln for ln in out.stdout.splitlines() if ln.strip()]
    except (OSError, subprocess.TimeoutExpired):
        return []


def git_commit_timestamp(repo: Path, commit: str) -> str:
    if not commit or len(commit) < 7:
        return ""
    try:
        out = subprocess.run(
            [
                "git",
                "-C",
                str(repo),
                "show",
                "-s",
                "--format=%cI",
                commit,
            ],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
        return out.stdout.strip() if out.returncode == 0 else ""
    except (OSError, subprocess.TimeoutExpired):
        return ""


def git_tags(repo: Path) -> list[str]:
    try:
        out = subprocess.run(
            ["git", "-C", str(repo), "tag", "-l", "v*"],
            capture_output=True,
            text=True,
            timeout=15,
            check=False,
        )
        if out.returncode != 0:
            return []
        return [t for t in out.stdout.splitlines() if t.strip()]
    except (OSError, subprocess.TimeoutExpired):
        return []


def build_version_deltas(
    audits: list[CampaignAudit],
    repo: Path,
) -> list[dict[str, Any]]:
    """Order campaigns by git commit time and list code changes between pins.

    Groups by unique git_commit so A/B canaries on the same SHA share one delta
    edge; each edge carries the success-rate snapshot of every campaign on the
    newer pin for side-by-side comparison.
    """
    by_commit: dict[str, list[CampaignAudit]] = {}
    for a in audits:
        c = (a.git_commit or "").strip()
        if len(c) < 7:
            continue
        by_commit.setdefault(c, []).append(a)

    ordered: list[tuple[str, str]] = []  # (commit, timestamp)
    for c in by_commit:
        ts = git_commit_timestamp(repo, c)
        ordered.append((c, ts or "1970-01-01T00:00:00Z"))
    ordered.sort(key=lambda x: x[1])

    deltas: list[dict[str, Any]] = []
    for i in range(1, len(ordered)):
        older_c, older_ts = ordered[i - 1]
        newer_c, newer_ts = ordered[i]
        commits = git_log_range(repo, older_c, newer_c, n=50)
        older_camps = by_commit[older_c]
        newer_camps = by_commit[newer_c]

        def _snap(camps: list[CampaignAudit]) -> list[dict[str, Any]]:
            out: list[dict[str, Any]] = []
            for a in sorted(camps, key=lambda x: x.name):
                out.append(
                    {
                        "name": a.name,
                        "kind": a.kind,
                        "n_targets": a.n_targets,
                        "n_genuine": a.n_genuine,
                        "genuine_rate": a.genuine_rate
                        if math.isfinite(a.genuine_rate)
                        else None,
                        "s1_rate": a.s1_rate if math.isfinite(a.s1_rate) else None,
                        "matrix_md5": a.matrix_md5,
                    }
                )
            return out

        deltas.append(
            {
                "from_commit": older_c,
                "to_commit": newer_c,
                "from_ts": older_ts,
                "to_ts": newer_ts,
                "n_commits": len(commits),
                "commits": commits,
                "from_campaigns": _snap(older_camps),
                "to_campaigns": _snap(newer_camps),
            }
        )
    return deltas


def _fmt_rate(x: float | None) -> str:
    if x is None or not isinstance(x, (int, float)) or not math.isfinite(float(x)):
        return "n/a"
    return f"{float(x):.1%}"


def render_markdown(
    audits: list[CampaignAudit],
    *,
    repo: Path,
    generated: str,
    deltas: list[dict[str, Any]] | None = None,
) -> str:
    lines: list[str] = []
    lines.append("# Benchmark version audit")
    lines.append("")
    lines.append(f"Generated: {generated}")
    lines.append("")
    lines.append(
        "Metrics: **genuine** = seed_echo=0 + non-seed pose_source + RMSD≤2.0 Å; "
        "**S1** = rmsd_to_crystal≤2.0 (diagnostic only — not claim-ready)."
    )
    lines.append("")
    lines.append(
        "Local-first: only `$FLEXAIDDS_LOCAL_ROOT` / `~/flexaidds_results` trees "
        "(never CloudDocs `rglob`)."
    )
    lines.append("")
    lines.append("## Summary table")
    lines.append("")
    lines.append(
        "| Name | Kind | N | Genuine | Rate | S1 rate | Matrix | Commit | Seed elit. |"
    )
    lines.append("|------|------|---|---------|------|---------|--------|--------|------------|")
    for a in sorted(audits, key=lambda x: x.name):
        gr = _fmt_rate(a.genuine_rate if math.isfinite(a.genuine_rate) else None)
        sr = _fmt_rate(a.s1_rate if math.isfinite(a.s1_rate) else None)
        lines.append(
            f"| `{a.name}` | {a.kind} | {a.n_targets} | "
            f"{a.n_genuine}/{a.n_targets} | {gr} | {sr} | "
            f"`{a.matrix_md5 or '—'}` | `{a.git_commit or '—'}` | {a.seed_elitism or '—'} |"
        )
    lines.append("")

    # Success-rate leaderboard (by genuine rate, then N)
    scored = [a for a in audits if a.n_targets > 0 and math.isfinite(a.genuine_rate)]
    scored.sort(key=lambda x: (-x.genuine_rate, -x.n_targets, x.name))
    if scored:
        lines.append("## Success-rate ranking (genuine)")
        lines.append("")
        lines.append("| Rank | Name | Genuine | Rate | S1 | Matrix | Commit |")
        lines.append("|------|------|---------|------|----|--------|--------|")
        for i, a in enumerate(scored, 1):
            lines.append(
                f"| {i} | `{a.name}` | {a.n_genuine}/{a.n_targets} | "
                f"{_fmt_rate(a.genuine_rate)} | {_fmt_rate(a.s1_rate)} | "
                f"`{a.matrix_md5 or '—'}` | `{a.git_commit or '—'}` |"
            )
        lines.append("")

    lines.append("## Modifications between benchmark pins")
    lines.append("")
    lines.append(
        "Unique git commits ordered by committer date. Each step lists code changes "
        "(`git log older..newer`) and the success snapshots on each pin."
    )
    lines.append("")
    if deltas:
        for d in deltas:
            lines.append(
                f"### `{d['from_commit'][:12]}` → `{d['to_commit'][:12]}` "
                f"({d.get('n_commits', 0)} commits)"
            )
            lines.append("")
            lines.append(f"- From: `{d.get('from_ts', '')}`")
            lines.append(f"- To:   `{d.get('to_ts', '')}`")
            lines.append("- Campaigns on older pin:")
            for c in d.get("from_campaigns") or []:
                lines.append(
                    f"  - `{c['name']}`: genuine={c['n_genuine']}/{c['n_targets']} "
                    f"({_fmt_rate(c.get('genuine_rate'))}) matrix=`{c.get('matrix_md5') or '—'}`"
                )
            lines.append("- Campaigns on newer pin:")
            for c in d.get("to_campaigns") or []:
                lines.append(
                    f"  - `{c['name']}`: genuine={c['n_genuine']}/{c['n_targets']} "
                    f"({_fmt_rate(c.get('genuine_rate'))}) matrix=`{c.get('matrix_md5') or '—'}`"
                )
            commits = d.get("commits") or []
            if commits:
                lines.append("- Code changes:")
                lines.append("```")
                for ln in commits[:40]:
                    lines.append(ln)
                if len(commits) > 40:
                    lines.append(f"... ({len(commits) - 40} more)")
                lines.append("```")
            else:
                lines.append("- Code changes: _(empty range or unreachable commits)_")
            lines.append("")
    else:
        lines.append("_(fewer than two pinned commits — no inter-pin delta)_")
        lines.append("")

    lines.append("## Git tags in repo")
    lines.append("")
    tags = git_tags(repo)
    if tags:
        lines.append(", ".join(f"`{t}`" for t in tags[-20:]))
    else:
        lines.append("_(no v* tags found)_")
    lines.append("")
    lines.append("## Per-campaign notes + recent commits at pin")
    lines.append("")
    for a in sorted(audits, key=lambda x: x.name):
        lines.append(f"### `{a.name}`")
        lines.append("")
        lines.append(f"- Path: `{a.path}`")
        if a.notes:
            lines.append(f"- Notes: {a.notes}")
        if a.git_commit:
            log = git_log_for_commit(repo, a.git_commit, 8)
            if log:
                lines.append("- Git log at pin:")
                lines.append("```")
                lines.append(log)
                lines.append("```")
        fails = [
            t
            for t in a.targets
            if isinstance(t, dict) and not t.get("genuine") and t.get("pdb")
        ]
        if fails and len(fails) <= 40:
            lines.append(f"- Non-genuine targets ({len(fails)}):")
            for t in fails[:40]:
                lines.append(
                    f"  - {t.get('pdb')}: rmsd_hun={t.get('rmsd_hun')} "
                    f"bcr={t.get('best_cluster')} seed_echo={t.get('seed_echo')} "
                    f"src={t.get('pose_source')}"
                )
        elif fails:
            lines.append(f"- Non-genuine targets: {len(fails)} (list truncated in JSON)")
        lines.append("")
    lines.append("## Gate reminder")
    lines.append("")
    lines.append(
        "Astex-85 **accuracy claims** require genuine ≥1/3 on an agreed claim canary "
        "with seed_echo=0 (AGENTS / METHODOLOGY). Concurrent full-85 trees are not "
        "automatically claim-valid. S1 alone is never a claim headline."
    )
    lines.append("")
    return "\n".join(lines)


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--roots",
        nargs="*",
        default=None,
        help="Local result roots (default: $FLEXAIDDS_LOCAL_ROOT or ~/flexaidds_results)",
    )
    ap.add_argument(
        "--repo",
        type=Path,
        default=None,
        help="Git repo for tag/commit correlation (default: cwd or parents)",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Write audit_benchmark_versions.{json,md} here",
    )
    ap.add_argument("--json-only", action="store_true")
    args = ap.parse_args(argv)

    roots = [Path(r) for r in args.roots] if args.roots else _default_roots()
    repo = args.repo
    if repo is None:
        cwd = Path.cwd()
        for p in [cwd, *cwd.parents]:
            if (p / ".git").exists() and (p / "AGENTS.md").exists():
                repo = p
                break
        if repo is None:
            repo = cwd

    camps = discover_campaigns(roots)
    audits = [audit_campaign_dir(c) for c in camps]
    deltas = build_version_deltas(audits, repo)
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ")

    payload = {
        "generated_utc": generated,
        "roots": [str(r) for r in roots],
        "n_campaigns": len(audits),
        "campaigns": [asdict(a) for a in audits],
        "version_deltas": deltas,
        "git_tags": git_tags(repo),
        "metrics_contract": {
            "genuine": "seed_echo=0 + non-seed pose_source + (rmsd_hungarian|rmsd_to_crystal)≤2.0",
            "S1": "rmsd_to_crystal≤2.0 diagnostic only",
            "claim": "see benchmarks/protocols/admission_metrics_contract.md + aggregate_claim_metrics.py",
        },
    }

    out_dir = args.out_dir
    if out_dir is None:
        local = Path.home() / "flexaidds_results" / "workorders"
        out_dir = local if local.parent.is_dir() else Path.cwd()
    out_dir = out_dir.expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)

    jpath = out_dir / "audit_benchmark_versions.json"
    jpath.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")

    if not args.json_only:
        md = render_markdown(
            audits, repo=repo, generated=generated, deltas=deltas
        )
        mpath = out_dir / "audit_benchmark_versions.md"
        mpath.write_text(md, encoding="utf-8")
        print(f"Wrote {mpath}")
    print(f"Wrote {jpath}")
    print(f"Campaigns inventoried: {len(audits)}")
    print(f"Version deltas (unique commit edges): {len(deltas)}")
    for a in sorted(audits, key=lambda x: -x.n_targets)[:15]:
        gr = f"{a.genuine_rate:.0%}" if math.isfinite(a.genuine_rate) else "n/a"
        print(
            f"  {a.name}: n={a.n_targets} genuine={a.n_genuine} ({gr}) "
            f"commit={a.git_commit or '-'}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
