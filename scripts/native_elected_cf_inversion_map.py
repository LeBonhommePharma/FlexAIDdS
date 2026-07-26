#!/usr/bin/env python3
"""Score-only Native vs Elected CF inversion map (campaign Phase 1.5 / NEXT_CAMPAIGN_STEP).

Classifies each target SCORING-LOCKED | SEARCH-MISS | TIED under production LOCCLF.

Liveness: every score uses --config + --ligand; refuse diagnostic/probe_config.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

PANEL = ("1J3J", "1K3U", "1L7F", "1N1M", "1M2Z", "1OQ5", "1SQ5", "1YGC")
EPS_DEFAULT = 0.5


def repo_root() -> Path:
    p = Path(__file__).resolve().parents[1]
    return p if (p / "AGENTS.md").exists() else Path.cwd()


def resolve_config(repo: Path, pdb: str, pilot: Optional[Path] = None) -> Optional[Path]:
    candidates = [
        repo / "ops" / "gates" / "configs" / f"{pdb}_dock_config.json",
    ]
    if pilot is not None:
        candidates.append(pilot / pdb / "dock_config.json")
    for c in candidates:
        if c.is_file():
            return c
    return None


def resolve_native(repo: Path, pdb: str) -> Optional[Path]:
    for c in (
        Path.home() / f".flexaidds/benchmarks/astex_diverse/{pdb}/{pdb}_ligand.sdf",
        repo / f"benchmarks/astex_diverse/astex_diverse/{pdb}/{pdb}_ligand.sdf",
        repo / f"benchmarks/astex_diverse/{pdb}/{pdb}_ligand.sdf",
    ):
        if c.is_file():
            return c
    return None


def resolve_receptor(repo: Path, pdb: str) -> Optional[Path]:
    for c in (
        Path.home() / f".flexaidds/benchmarks/astex_diverse/{pdb}/{pdb}_apo.pdb",
        repo / f"benchmarks/astex_diverse/astex_diverse/{pdb}/{pdb}_apo.pdb",
        repo / f"benchmarks/astex_diverse/{pdb}/{pdb}_apo.pdb",
    ):
        if c.is_file():
            return c
    return None


def resolve_elected(pilot: Path, pdb: str) -> Optional[Path]:
    for c in (
        pilot / pdb / "elected_pose.pdb",
        pilot / pdb / f"{pdb}_elected.pdb",
    ):
        if c.is_file():
            return c
    # result.csv elected_pose_path
    rc = pilot / pdb / "result.csv"
    if rc.is_file():
        with rc.open() as f:
            row = next(csv.DictReader(f), None)
        if row and row.get("elected_pose_path"):
            p = Path(row["elected_pose_path"])
            if p.is_file():
                return p
    return None


def load_result_metrics(pilot: Path, pdb: str) -> Dict[str, Any]:
    rc = pilot / pdb / "result.csv"
    if not rc.is_file():
        return {}
    with rc.open() as f:
        row = next(csv.DictReader(f), None)
    if not row:
        return {}
    out: Dict[str, Any] = {}
    for k in (
        "rmsd_hungarian",
        "best_cluster_rmsd",
        "seed_echo",
        "elected_cf",
        "elected_pose_path",
    ):
        if k in row and row[k] not in (None, ""):
            try:
                out[k] = float(row[k]) if k != "elected_pose_path" else row[k]
            except ValueError:
                out[k] = row[k]
    return out


def probe_cf(
    probe: Path,
    *,
    receptor: Path,
    pose: Path,
    ligand: Path,
    config: Path,
    binary: Path,
    data_dir: Path,
    label: str,
) -> Optional[dict]:
    if "diagnostic/probe_config" in str(config).replace("\\", "/"):
        print(f"refuse nonproduction config {config}", file=sys.stderr)
        return None
    env = os.environ.copy()
    # Keep map free of burial/boom/wall knobs
    for k in list(env):
        if k.startswith("FLEXAIDDS_") and any(
            x in k for x in ("BOOM", "WAL_COERCIVE", "PB_CLASH", "MEMETIC", "WALL_PILOT")
        ):
            env.pop(k, None)
    env["FLEXAIDDS_WAL_COERCIVE"] = "0"
    cmd = [
        str(probe),
        "--receptor",
        str(receptor),
        "--pose",
        str(pose),
        "--ligand",
        str(ligand),
        "--config",
        str(config),
        "--binary",
        str(binary),
        "--data-dir",
        str(data_dir),
        "--pdb",
        label,
    ]
    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=600, env=env, check=False
        )
    except (OSError, subprocess.TimeoutExpired) as e:
        print(f"probe error {label}: {e}", file=sys.stderr)
        return None
    for line in reversed(proc.stdout.splitlines()):
        line = line.strip()
        if not line.startswith("{"):
            continue
        try:
            obj = json.loads(line)
            if "cf_total" in obj:
                return obj
        except (json.JSONDecodeError, TypeError, ValueError):
            continue
    if proc.returncode != 0:
        print(f"probe fail {label} rc={proc.returncode} {proc.stderr[-200:]}", file=sys.stderr)
    return None


def fnum(obj: Optional[dict], key: str = "cf_total") -> Optional[float]:
    if not obj or key not in obj:
        return None
    try:
        return float(obj[key])
    except (TypeError, ValueError):
        return None


def classify(cf_native: float, cf_elected: float, eps: float) -> str:
    # lower CF is better
    if cf_elected + eps < cf_native:
        return "SCORING-LOCKED"
    if cf_native + eps < cf_elected:
        return "SEARCH-MISS"
    return "TIED"


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    repo = repo_root()
    ap.add_argument("--repo", type=Path, default=repo)
    ap.add_argument("--pilot-dir", type=Path, required=True)
    ap.add_argument("--probe-cf", type=Path, default=None)
    ap.add_argument("--binary", type=Path, default=None)
    ap.add_argument("--data-dir", type=Path, default=None)
    ap.add_argument("--panel", nargs="*", default=None)
    ap.add_argument("--eps", type=float, default=EPS_DEFAULT)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(list(argv) if argv is not None else None)

    repo = args.repo.expanduser().resolve()
    pilot = args.pilot_dir.expanduser().resolve()
    out_dir = args.out_dir.expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    probe = args.probe_cf or (repo / "build" / "probe_cf")
    binary = args.binary or (repo / "build" / "FlexAIDdS")
    data_dir = args.data_dir or repo
    if not probe.is_file() or not binary.is_file():
        print("missing probe_cf or FlexAIDdS", file=sys.stderr)
        return 2

    panel = list(args.panel) if args.panel else list(PANEL)
    rows: List[Dict[str, Any]] = []
    for pdb in panel:
        cfg = resolve_config(repo, pdb, pilot)
        rec = resolve_receptor(repo, pdb)
        nat = resolve_native(repo, pdb)
        ele = resolve_elected(pilot, pdb)
        meta = load_result_metrics(pilot, pdb)
        if not all((cfg, rec, nat, ele)):
            rows.append(
                {
                    "pdb": pdb,
                    "status": "skip_missing",
                    "config": str(cfg) if cfg else "",
                    "elected": str(ele) if ele else "",
                }
            )
            continue
        print(f"[{pdb}] native + elected …")
        nobj = probe_cf(
            probe,
            receptor=rec,
            pose=nat,
            ligand=nat,
            config=cfg,
            binary=binary,
            data_dir=data_dir,
            label=f"{pdb}_native",
        )
        eobj = probe_cf(
            probe,
            receptor=rec,
            pose=ele,
            ligand=nat,
            config=cfg,
            binary=binary,
            data_dir=data_dir,
            label=f"{pdb}_elected",
        )
        cn, ce = fnum(nobj), fnum(eobj)
        if cn is None or ce is None:
            rows.append({"pdb": pdb, "status": "probe_fail"})
            continue
        cls = classify(cn, ce, args.eps)
        dcf = cn - ce  # negative => native better
        bcr = meta.get("best_cluster_rmsd")
        elect_rmsd = meta.get("rmsd_hungarian")
        sub = "n/a"
        if isinstance(bcr, (int, float)):
            if cls == "SEARCH-MISS" and bcr <= 2.0:
                sub = "election_pool"
            elif cls == "SEARCH-MISS" and bcr > 2.0:
                sub = "sampling_ceiling"
            elif cls == "SCORING-LOCKED":
                sub = "false_min_attractor"
        rows.append(
            {
                "pdb": pdb,
                "status": "ok",
                "class": cls,
                "subclass": sub,
                "cf_native": cn,
                "cf_elected": ce,
                "dCF_native_minus_elected": dcf,
                "clash_native": fnum(nobj, "cf_clash"),
                "clash_elected": fnum(eobj, "cf_clash"),
                "wal_native": fnum(nobj, "cf_wal"),
                "wal_elected": fnum(eobj, "cf_wal"),
                "rmsd_hungarian": elect_rmsd,
                "best_cluster_rmsd": bcr,
                "seed_echo": meta.get("seed_echo"),
                "result_elected_cf": meta.get("elected_cf"),
                "elected_path": str(ele),
                "config": str(cfg),
            }
        )
        print(f"  {pdb}: {cls} cn={cn:.3f} ce={ce:.3f} dCF={dcf:+.3f} sub={sub}")

    ok = [r for r in rows if r.get("status") == "ok"]
    n = len(ok)
    n_lock = sum(1 for r in ok if r["class"] == "SCORING-LOCKED")
    n_miss = sum(1 for r in ok if r["class"] == "SEARCH-MISS")
    n_tied = sum(1 for r in ok if r["class"] == "TIED")
    n_samp = sum(1 for r in ok if r.get("subclass") == "sampling_ceiling")
    n_elect = sum(1 for r in ok if r.get("subclass") == "election_pool")

    # ACCEPT completeness
    if n >= 6:
        verdict = "PASS"
        reason = f"classified {n}/8; LOCKED={n_lock} MISS={n_miss} TIED={n_tied}"
    else:
        verdict = "FAIL"
        reason = f"only {n}/8 scored (need ≥6)"

    # regime recommendation
    if n_miss > n_lock:
        next_lever = "SEARCH-MISS dominant → FLEXAIDDS_COARSE_ORIENTATIONS=256 W1 pilot (matrix 9dc9)"
    elif n_lock > n_miss:
        next_lever = "SCORING-LOCKED dominant → strong deep-interpenetration decoys / scoring; do not BOOM thrash"
    else:
        next_lever = "mixed → dock coarse-orient only on SEARCH-MISS codes"

    written = datetime.now(timezone.utc).isoformat()
    payload = {
        "schema": "native_elected_cf_inversion_map/v1",
        "written_utc": written,
        "verdict": verdict,
        "reason": reason,
        "eps": args.eps,
        "n_scored": n,
        "n_scoring_locked": n_lock,
        "n_search_miss": n_miss,
        "n_tied": n_tied,
        "n_sampling_ceiling": n_samp,
        "n_election_pool": n_elect,
        "next_lever": next_lever,
        "pilot_dir": str(pilot),
        "rows": rows,
        "one_variable": "pose_role under fixed production LOCCLF (native vs elected)",
        "blocks": [
            "no full-85",
            "no memetic",
            "no WAL_COERCIVE re-panel",
            "no interval-only BOOM",
        ],
    }
    (out_dir / "inversion_map.json").write_text(json.dumps(payload, indent=2) + "\n")
    fields = [
        "pdb",
        "status",
        "class",
        "subclass",
        "cf_native",
        "cf_elected",
        "dCF_native_minus_elected",
        "rmsd_hungarian",
        "best_cluster_rmsd",
        "seed_echo",
    ]
    with (out_dir / "inversion_map.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    lines = [
        f"# Native–Elected CF inversion map — {verdict}",
        "",
        f"**Written:** {written}  ",
        f"**Pilot poses:** `{pilot}`  ",
        f"**One variable:** pose role under fixed production LOCCLF (native vs elected)  ",
        f"**ε:** {args.eps} CF units  ",
        f"**OUT:** `{out_dir}`  ",
        "",
        f"**Verdict:** **{verdict}** — {reason}",
        "",
        f"**Next lever (guided):** {next_lever}",
        "",
        "## Counts",
        "",
        f"| Class | N |",
        f"|-------|--:|",
        f"| SCORING-LOCKED | {n_lock} |",
        f"| SEARCH-MISS | {n_miss} |",
        f"| TIED | {n_tied} |",
        f"| of MISS: sampling_ceiling (BCR>2) | {n_samp} |",
        f"| of MISS: election_pool (BCR≤2) | {n_elect} |",
        "",
        "## Per-target",
        "",
        "| PDB | class | subclass | CF native | CF elected | dCF(n−e) | elect RMSD | BCR |",
        "|-----|-------|----------|----------:|-----------:|---------:|-----------:|----:|",
    ]
    for r in ok:
        lines.append(
            f"| {r['pdb']} | {r['class']} | {r['subclass']} | "
            f"{r['cf_native']:.3f} | {r['cf_elected']:.3f} | "
            f"{r['dCF_native_minus_elected']:+.3f} | "
            f"{r.get('rmsd_hungarian', 'n/a')} | {r.get('best_cluster_rmsd', 'n/a')} |"
        )
    lines += [
        "",
        "## Cadence",
        "",
        "- Phase: 1.5 diagnosis (score-only) / roadmap Phase 1 extension  ",
        "- Genuine/BCR: **not claimed** this step  ",
        f"- **{verdict}**  ",
        "- No full-85; no memetic; no WAL re-panel  ",
        "",
        "Methodology: `NEXT_CAMPAIGN_STEP.md` · BENCHMARKING_ROADMAP Phase 1 then Phase 4 sampling.  ",
        "",
    ]
    (out_dir / "inversion_map.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
