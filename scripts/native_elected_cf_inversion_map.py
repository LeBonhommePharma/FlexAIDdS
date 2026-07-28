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
PILOT8_PANEL = (
    "1G9V",
    "1GPK",
    "1MEH",
    "1P62",
    "1Q4G",
    "1R9O",
    "1T40",
    "2BYS",
)
EPS_DEFAULT = 0.5


def repo_root() -> Path:
    p = Path(__file__).resolve().parents[1]
    return p if (p / "AGENTS.md").exists() else Path.cwd()


def _as_path(raw: str) -> Optional[Path]:
    """Resolve a path string; tolerate missing ``.pdb`` suffix on elected_path."""
    if not raw or not str(raw).strip():
        return None
    p = Path(str(raw).strip())
    if p.is_file():
        return p
    # pilot8 parse may drop extension
    if not p.suffix:
        for ext in (".pdb", ".PDB"):
            cand = Path(str(p) + ext)
            if cand.is_file():
                return cand
    # relative to cwd
    if p.name and Path(p.name).is_file():
        return Path(p.name).resolve()
    return None


def resolve_config(
    repo: Path,
    pdb: str,
    pilot: Optional[Path] = None,
    *,
    config_template: Optional[Path] = None,
) -> Optional[Path]:
    candidates = [
        repo / "ops" / "gates" / "configs" / f"{pdb}_dock_config.json",
    ]
    if pilot is not None:
        candidates.extend(
            [
                pilot / pdb / "dock_config.json",
                pilot / "configs" / f"{pdb}_dock_config.json",
                pilot / "dock_configs" / f"{pdb}_dock_config.json",
            ]
        )
    for c in candidates:
        if c.is_file():
            return c
    # Materialize from production template when per-PDB ops config is absent
    # (pilot8 dual-zero panel). Prefer leaf S4-style production knobs.
    if config_template is not None and config_template.is_file() and pilot is not None:
        dest_dir = pilot / "dock_configs"
        dest_dir.mkdir(parents=True, exist_ok=True)
        dest = dest_dir / f"{pdb}_dock_config.json"
        if not dest.is_file():
            dest.write_text(config_template.read_text(encoding="utf-8"), encoding="utf-8")
        return dest
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
    # Prefer pilot8 result.csv elected_path (rank-0 elect), then explicit names,
    # then r0_0 only as last-resort fallback.
    rc = pilot / pdb / "result.csv"
    if rc.is_file():
        with rc.open() as f:
            row = next(csv.DictReader(f), None)
        if row:
            for key in ("elected_path", "elected_pose_path", "cf_top1_pose_path"):
                raw = row.get(key)
                if not raw:
                    continue
                p = _as_path(str(raw))
                if p is not None:
                    return p
    for c in (
        pilot / pdb / "elected_pose.pdb",
        pilot / pdb / f"{pdb}_elected.pdb",
        pilot / pdb / f"{pdb}_r0_0.pdb",
    ):
        if c.is_file():
            return c
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
    # Normalize pilot8 + classic column names into a stable set.
    aliases = {
        "rmsd_hungarian": ("rmsd_hungarian", "rmsd_top1"),
        "best_cluster_rmsd": ("best_cluster_rmsd", "rmsd_bcr"),
        "seed_echo": ("seed_echo",),
        "elected_cf": ("elected_cf", "score_top1", "best_score"),
        "elected_pose_path": ("elected_pose_path", "elected_path"),
    }
    for dest, keys in aliases.items():
        for k in keys:
            if k not in row or row[k] in (None, ""):
                continue
            try:
                out[dest] = (
                    float(row[k]) if dest != "elected_pose_path" else row[k]
                )
            except ValueError:
                out[dest] = row[k]
            break
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
    ap.add_argument(
        "--pilot8",
        action="store_true",
        help="Use dual-zero pilot8 panel (1G9V…2BYS)",
    )
    ap.add_argument(
        "--config-template",
        type=Path,
        default=None,
        help="Production dock_config.json template when per-PDB ops config is missing",
    )
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
        print(f"missing probe_cf or FlexAIDdS: probe={probe} binary={binary}", file=sys.stderr)
        return 2

    cfg_template = args.config_template.expanduser().resolve() if args.config_template else None
    if cfg_template is None:
        # Prefer leaf production axes (S4-style LOCCLF) over sparse ops/gates set.
        for cand in (
            Path.home()
            / "flexaidds_results/s4_pheno_unique_near_miss_20260727_211213/arm_control/1L7F/dock_config.json",
            repo / "ops" / "gates" / "configs" / "1L7F_dock_config.json",
        ):
            if cand.is_file():
                cfg_template = cand
                break

    if args.panel:
        panel = list(args.panel)
    elif args.pilot8:
        panel = list(PILOT8_PANEL)
    else:
        panel = list(PANEL)
    rows: List[Dict[str, Any]] = []
    for pdb in panel:
        cfg = resolve_config(repo, pdb, pilot, config_template=cfg_template)
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
                    "receptor": str(rec) if rec else "",
                    "native": str(nat) if nat else "",
                    "elected": str(ele) if ele else "",
                }
            )
            print(
                f"[{pdb}] skip_missing cfg={bool(cfg)} rec={bool(rec)} "
                f"nat={bool(nat)} ele={bool(ele)}",
                file=sys.stderr,
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
            rows.append(
                {
                    "pdb": pdb,
                    "status": "probe_fail",
                    "config": str(cfg),
                    "elected_path": str(ele),
                }
            )
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
                "native_path": str(nat),
                "receptor_path": str(rec),
                "config": str(cfg),
                "config_template": str(cfg_template) if cfg_template else "",
                "binary": str(binary),
                "data_dir": str(data_dir),
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
    # Binary identity (reconstruction labels if present beside binary)
    bin_id = ""
    for idp in (
        Path(binary).parent / "IDENTITY.txt",
        Path(binary).with_name("IDENTITY.txt"),
    ):
        if idp.is_file():
            bin_id = idp.read_text(encoding="utf-8", errors="ignore").strip()
            break
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
        "probe_cf": str(probe),
        "binary": str(binary),
        "binary_identity": bin_id,
        "data_dir": str(data_dir),
        "config_template": str(cfg_template) if cfg_template else "",
        "dCF_sign": "CF_native - CF_elected (negative ⇒ native better on CF proxy)",
        "rows": rows,
        "one_variable": "pose_role under fixed production LOCCLF (native vs elected)",
        "claim_success_rates": False,
        "full85_authorized": False,
        "blocks": [
            "no full-85",
            "no memetic",
            "no WAL_COERCIVE re-panel",
            "no interval-only BOOM",
            "does not claim docking success rates",
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
        f"**dCF sign:** CF_native − CF_elected (negative ⇒ native better)  ",
        f"**Binary:** `{binary}`  ",
        f"**Config template:** `{cfg_template}`  ",
        f"**OUT:** `{out_dir}`  ",
        "",
        f"**Verdict:** **{verdict}** — {reason}",
        "",
        f"**Next lever (guided):** {next_lever}",
        "",
        "**Does not claim docking success rates; does not authorize full-85.**",
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
