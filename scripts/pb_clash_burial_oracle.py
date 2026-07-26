#!/usr/bin/env python3
"""pb_clash burial / scoring oracle (score-only).

ROADMAP_v2 (2026-07-25) voids the SEARCH-MISS "clean probe" PASS: ΔdCF 1e-4..0.02
is noise; those targets already have native CF-min. Revised Phase 2:

  Panel: SCORING-LOCKED 1OQ5 1SQ5 1YGC
  Decoy: each target's actual elected pose (false-min attractor)
  One variable: FLEXAIDDS_PB_CLASH_WEIGHT (OFF=0 vs ON weight; ladder = one arm each)
  ACCEPT (magnitude floor):
    dCF decreases by ≥1.0 kcal AND flips sign on ≥2/3 inverted targets,
    AND no SEARCH-MISS probe regresses (native stays CF-min on all 5).

Legacy --mode search-miss keeps the old panel for historical re-runs only; its
PASS must not unlock memetic (workorder VOID).

Usage:
  python3 scripts/pb_clash_burial_oracle.py \\
    --mode scoring-locked \\
    --elected-leaf ~/flexaidds_results/pilot_w1_boom_interval_20260725_134740 \\
    --weight 1.0 \\
    --out-dir ~/flexaidds_results/workorders/pb_clash_scoring_locked
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

# SEARCH-MISS: native already CF-better; sampling only (never judge scoring here).
SEARCH_MISS_PANEL = ("1J3J", "1K3U", "1L7F", "1N1M", "1M2Z")
# SCORING-LOCKED: elected beats crystal while BCR≤2; scoring only.
SCORING_LOCKED_PANEL = ("1OQ5", "1SQ5", "1YGC")
# Alias kept for older callers / tests
CLEAN_PANEL = SEARCH_MISS_PANEL

MAGNITUDE_FLOOR_KCAL = 1.0
SIGN_FLIP_MIN_FRAC = 2  # of 3 SCORING-LOCKED


def find_repo_root() -> Path:
    p = Path(__file__).resolve().parents[1]
    if (p / "AGENTS.md").exists():
        return p
    return Path.cwd()


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve_dock_config(
    repo: Path, pdb: str, elected_leaf: Optional[Path]
) -> Optional[Path]:
    c = repo / "ops" / "gates" / "configs" / f"{pdb}_dock_config.json"
    if c.is_file():
        return c
    if elected_leaf is not None:
        p = elected_leaf / pdb / "dock_config.json"
        if p.is_file():
            return p
    return None


def resolve_receptor(repo: Path, pdb: str) -> Optional[Path]:
    for c in (
        repo / f"benchmarks/astex_diverse/astex_diverse/{pdb}/{pdb}_apo.pdb",
        repo / f"benchmarks/astex_diverse/{pdb}/{pdb}_apo.pdb",
        Path.home() / f".flexaidds/benchmarks/astex_diverse/{pdb}/{pdb}_apo.pdb",
    ):
        if c.is_file():
            return c
    return None


def resolve_native(repo: Path, pdb: str) -> Optional[Path]:
    for c in (
        repo / f"benchmarks/astex_diverse/astex_diverse/{pdb}/{pdb}_ligand.sdf",
        repo / f"benchmarks/astex_diverse/{pdb}/{pdb}_ligand.sdf",
        Path.home() / f".flexaidds/benchmarks/astex_diverse/{pdb}/{pdb}_ligand.sdf",
        repo / "diagnostic" / "refs" / pdb / "native.sdf",
    ):
        if c.is_file():
            return c
    return None


def resolve_decoy(
    repo: Path,
    pdb: str,
    decoy_dir: Optional[Path],
    elected_leaf: Optional[Path],
) -> Optional[Path]:
    if elected_leaf is not None:
        for name in ("elected_pose.pdb", f"{pdb}_elected.pdb"):
            p = elected_leaf / pdb / name
            if p.is_file():
                return p
    if decoy_dir is not None:
        for name in (
            f"{pdb}_buried.pdb",
            f"{pdb}_falsemin.pdb",
            f"{pdb}.pdb",
            "elected_pose.pdb",
        ):
            p = decoy_dir / pdb / name if (decoy_dir / pdb).is_dir() else decoy_dir / name
            if p.is_file():
                return p
            p2 = decoy_dir / name
            if p2.is_file():
                return p2
    d = repo / "diagnostic" / "refs" / pdb / "falsemin_armA.pdb"
    return d if d.is_file() else None


def com_pdb(path: Path) -> Optional[Tuple[float, float, float]]:
    xs: List[float] = []
    ys: List[float] = []
    zs: List[float] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        try:
            xs.append(float(line[30:38]))
            ys.append(float(line[38:46]))
            zs.append(float(line[46:54]))
        except ValueError:
            continue
    if not xs:
        return None
    n = float(len(xs))
    return (sum(xs) / n, sum(ys) / n, sum(zs) / n)


def translate_pdb(src: Path, dst: Path, dx: float, dy: float, dz: float) -> None:
    lines_out: List[str] = []
    for line in src.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                x = float(line[30:38]) + dx
                y = float(line[38:46]) + dy
                z = float(line[46:54]) + dz
                line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
            except ValueError:
                pass
        lines_out.append(line)
    dst.write_text("\n".join(lines_out) + "\n", encoding="utf-8")


def probe_json(
    probe: Path,
    *,
    receptor: Path,
    pose: Path,
    ligand: Path,
    binary: Path,
    data_dir: Path,
    config: Path,
    label: str,
    pb_clash_weight: float,
) -> Optional[dict]:
    """Score one pose under production config; one env variable pb_clash weight."""
    if not config.is_file():
        print(f"  missing config {config}", file=sys.stderr)
        return None
    if "diagnostic/probe_config" in str(config).replace("\\", "/"):
        print(f"  refuse non-production config {config}", file=sys.stderr)
        return None
    env = os.environ.copy()
    env["FLEXAIDDS_PB_CLASH_WEIGHT"] = str(pb_clash_weight)
    env["FLEXAIDDS_WAL_COERCIVE"] = "0"
    env.pop("FLEXAIDDS_COM_BURIAL_CAP", None)
    cmd = [
        str(probe),
        "--receptor",
        str(receptor),
        "--pose",
        str(pose),
        "--ligand",
        str(ligand),
        "--binary",
        str(binary),
        "--data-dir",
        str(data_dir),
        "--config",
        str(config),
        "--pdb",
        label,
    ]
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600,
            env=env,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as e:
        print(f"  probe_cf error {label}: {e}", file=sys.stderr)
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
        print(
            f"  probe_cf fail {label} rc={proc.returncode}: {proc.stderr[-400:]}",
            file=sys.stderr,
        )
    return None


def fnum(obj: Optional[dict], key: str) -> Optional[float]:
    if not obj or key not in obj:
        return None
    try:
        return float(obj[key])
    except (TypeError, ValueError):
        return None


def score_pair(
    probe: Path,
    *,
    rec: Path,
    nat: Path,
    dec: Path,
    cfg: Path,
    binary: Path,
    data_dir: Path,
    pdb: str,
    weight: float,
) -> Optional[Dict[str, Any]]:
    nat_j = probe_json(
        probe,
        receptor=rec,
        pose=nat,
        ligand=nat,
        binary=binary,
        data_dir=data_dir,
        config=cfg,
        label=f"{pdb}_nat_w{weight}",
        pb_clash_weight=weight,
    )
    dec_j = probe_json(
        probe,
        receptor=rec,
        pose=dec,
        ligand=nat,
        binary=binary,
        data_dir=data_dir,
        config=cfg,
        label=f"{pdb}_dec_w{weight}",
        pb_clash_weight=weight,
    )
    if not nat_j or not dec_j:
        return None
    cn = fnum(nat_j, "cf_total")
    cd = fnum(dec_j, "cf_total")
    if cn is None or cd is None:
        return None
    return {
        "cf_native": cn,
        "cf_decoy": cd,
        "dCF": cn - cd,  # positive => decoy better (scoring-locked fingerprint)
        "clash_native": fnum(nat_j, "cf_clash"),
        "clash_decoy": fnum(dec_j, "cf_clash"),
        "native_min": cn <= cd,
    }


def build_buried_decoy(
    probe: Path,
    *,
    pdb: str,
    rec: Path,
    nat: Path,
    falsemin: Path,
    cfg: Path,
    binary: Path,
    data_dir: Path,
    out_pose: Path,
    min_clash: float = 5.0,
) -> Optional[Path]:
    """Translate falsemin toward receptor COM until cf_clash is elevated."""
    rcom = com_pdb(rec)
    lcom = com_pdb(falsemin)
    if rcom is None or lcom is None:
        return None
    dx0, dy0, dz0 = rcom[0] - lcom[0], rcom[1] - lcom[1], rcom[2] - lcom[2]
    norm = (dx0 * dx0 + dy0 * dy0 + dz0 * dz0) ** 0.5
    if norm < 1e-6:
        return None
    ux, uy, uz = dx0 / norm, dy0 / norm, dz0 / norm
    best_pose = None
    best_clash = -1.0
    tmp = out_pose.parent / f"{pdb}_try.pdb"
    for step in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0]:
        translate_pdb(falsemin, tmp, ux * step, uy * step, uz * step)
        obj = probe_json(
            probe,
            receptor=rec,
            pose=tmp,
            ligand=nat,
            binary=binary,
            data_dir=data_dir,
            config=cfg,
            label=f"{pdb}_bury_{step}",
            pb_clash_weight=1.0,
        )
        clash = fnum(obj, "cf_clash") or 0.0
        total = fnum(obj, "cf_total")
        print(f"  bury {pdb} step={step:.1f} clash={clash:.3f} total={total}")
        if clash > best_clash:
            best_clash = clash
            out_pose.write_text(tmp.read_text(encoding="utf-8", errors="replace"))
            best_pose = out_pose
        if clash >= min_clash:
            break
    if tmp.is_file():
        tmp.unlink(missing_ok=True)  # type: ignore[arg-type]
    return best_pose if best_clash > 0 else best_pose


def verdict_search_miss_legacy(ok: List[Dict[str, Any]]) -> Tuple[str, str]:
    """Historical sign-only ACCEPT — VOID for science (ROADMAP_v2)."""
    n = len(ok)
    n_moved = sum(1 for r in ok if r.get("moved_toward_native"))
    n_regressed = sum(
        1 for r in ok if r.get("native_regressed") or r.get("native_min_lost")
    )
    n_identical = sum(1 for r in ok if r.get("identical_off_on"))
    accept_move = n_moved >= 4 if n >= 5 else (n_moved >= max(1, n - 1) if n else False)
    if n == 0:
        return "FAIL", "no scored targets"
    if n_identical == n:
        return "FAIL", "OFF≡ON on all targets"
    if accept_move and n_regressed == 0:
        return (
            "VOID_LEGACY_PASS",
            f"sign-only move {n_moved}/{n} — VOID under ROADMAP_v2 (no magnitude floor; SEARCH-MISS panel)",
        )
    return (
        "FAIL",
        f"moved={n_moved}/{n}; regressed={n_regressed}",
    )


def verdict_scoring_locked(
    ok: List[Dict[str, Any]],
    floor: float,
    min_flips: int,
) -> Tuple[str, str, Dict[str, Any]]:
    """ROADMAP_v2 ACCEPT: ΔdCF decrease ≥floor AND sign flip on ≥min_flips/3."""
    n = len(ok)
    stats = {
        "n_scored": n,
        "n_decrease_ge_floor": 0,
        "n_sign_flip": 0,
        "n_both": 0,
        "magnitude_floor": floor,
        "min_sign_flips": min_flips,
    }
    if n == 0:
        return "FAIL", "no scored SCORING-LOCKED targets", stats

    for r in ok:
        d0 = float(r["dCF_off"])
        d1 = float(r["dCF_on"])
        decrease = d0 - d1  # positive if dCF fell (moved toward/past native-better)
        ge_floor = decrease >= floor - 1e-12
        # Sign flip on inverted targets: OFF has elected better (dCF>0), ON native better (dCF<0)
        flipped = d0 > 0.0 and d1 < 0.0
        r["decrease"] = decrease
        r["ge_floor"] = ge_floor
        r["sign_flip"] = flipped
        r["both_accept"] = ge_floor and flipped
        if ge_floor:
            stats["n_decrease_ge_floor"] += 1
        if flipped:
            stats["n_sign_flip"] += 1
        if ge_floor and flipped:
            stats["n_both"] += 1

    if stats["n_both"] >= min_flips:
        return (
            "PASS",
            f"dCF decrease≥{floor} and sign flip on {stats['n_both']}/{n} (need ≥{min_flips})",
            stats,
        )
    return (
        "FAIL",
        f"both floor+flip on {stats['n_both']}/{n} (need ≥{min_flips}); "
        f"ge_floor={stats['n_decrease_ge_floor']} flips={stats['n_sign_flip']}",
        stats,
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    repo = find_repo_root()
    ap.add_argument("--repo", type=Path, default=repo)
    ap.add_argument("--probe-cf", type=Path, default=None)
    ap.add_argument("--binary", type=Path, default=None)
    ap.add_argument("--data-dir", type=Path, default=None)
    ap.add_argument(
        "--mode",
        choices=("scoring-locked", "search-miss"),
        default="scoring-locked",
        help="scoring-locked = ROADMAP_v2 Phase 2; search-miss = legacy VOID panel",
    )
    ap.add_argument("--panel", nargs="*", default=None)
    ap.add_argument(
        "--weight",
        type=float,
        default=1.0,
        help="Single ON value for FLEXAIDDS_PB_CLASH_WEIGHT (OFF is always 0)",
    )
    ap.add_argument(
        "--magnitude-floor",
        type=float,
        default=MAGNITUDE_FLOOR_KCAL,
        help="Min dCF decrease (kcal) for SCORING-LOCKED ACCEPT (default 1.0)",
    )
    ap.add_argument("--decoy-dir", type=Path, default=None)
    ap.add_argument(
        "--elected-leaf",
        type=Path,
        default=None,
        help="Leaf with {PDB}/elected_pose.pdb and dock_config.json",
    )
    ap.add_argument(
        "--build-buried",
        action="store_true",
        help="(search-miss only) translate falsemin toward receptor COM",
    )
    ap.add_argument(
        "--skip-search-miss-check",
        action="store_true",
        help="Skip SEARCH-MISS no-regression recheck (not for claim runs)",
    )
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(list(argv) if argv is not None else None)

    repo = args.repo.expanduser().resolve()
    probe = args.probe_cf
    if probe is None:
        c = repo / "build" / "probe_cf"
        probe = c if c.is_file() else None
    if probe is None or not Path(probe).is_file():
        print("error: probe_cf not found", file=sys.stderr)
        return 2
    binary = args.binary
    if binary is None:
        for c in (repo / "build" / "FlexAIDdS",):
            if c.is_file():
                binary = c
                break
    if binary is None or not Path(binary).is_file():
        print("error: FlexAIDdS binary not found", file=sys.stderr)
        return 2

    data_dir = args.data_dir or repo
    out_dir = args.out_dir.expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    weight_on = float(args.weight)
    if weight_on <= 0.0:
        print("error: --weight must be > 0 (OFF arm is always 0)", file=sys.stderr)
        return 2

    elected_leaf = args.elected_leaf
    if elected_leaf is not None:
        elected_leaf = elected_leaf.expanduser().resolve()
    if args.mode == "scoring-locked" and elected_leaf is None:
        # Default pilot leaf used by inversion map
        cand = (
            Path.home()
            / "flexaidds_results"
            / "pilot_w1_boom_interval_20260725_134740"
        )
        if cand.is_dir():
            elected_leaf = cand

    if args.mode == "scoring-locked":
        panel = list(args.panel) if args.panel else list(SCORING_LOCKED_PANEL)
    else:
        panel = list(args.panel) if args.panel else list(SEARCH_MISS_PANEL)

    buried_dir = out_dir / "buried"
    if args.build_buried:
        buried_dir.mkdir(parents=True, exist_ok=True)

    rows: List[Dict[str, Any]] = []
    for pdb in panel:
        cfg = resolve_dock_config(repo, pdb, elected_leaf)
        rec = resolve_receptor(repo, pdb)
        nat = resolve_native(repo, pdb)
        dec = resolve_decoy(repo, pdb, args.decoy_dir, elected_leaf)
        if not cfg or not rec or not nat:
            rows.append(
                {
                    "pdb": pdb,
                    "status": "skip_missing_inputs",
                    "config": str(cfg) if cfg else "",
                }
            )
            continue
        if "diagnostic/probe_config" in str(cfg).replace("\\", "/"):
            rows.append({"pdb": pdb, "status": "fail_nonproduction_config"})
            continue

        if args.build_buried and args.mode == "search-miss":
            if not dec or not Path(dec).is_file():
                rows.append({"pdb": pdb, "status": "skip_missing_decoy_for_bury"})
                continue
            buried = buried_dir / f"{pdb}_buried.pdb"
            print(f"[{pdb}] building buried decoy…")
            built = build_buried_decoy(
                Path(probe),
                pdb=pdb,
                rec=Path(rec),
                nat=Path(nat),
                falsemin=Path(dec),
                cfg=Path(cfg),
                binary=Path(binary),
                data_dir=Path(data_dir),
                out_pose=buried,
            )
            if built is not None:
                dec = built
        if not dec or not Path(dec).is_file():
            rows.append({"pdb": pdb, "status": "skip_missing_decoy"})
            continue

        print(
            f"[{pdb}] mode={args.mode} config={cfg} decoy={dec} weight_on={weight_on}"
        )
        off = score_pair(
            Path(probe),
            rec=Path(rec),
            nat=Path(nat),
            dec=Path(dec),
            cfg=Path(cfg),
            binary=Path(binary),
            data_dir=Path(data_dir),
            pdb=pdb,
            weight=0.0,
        )
        on = score_pair(
            Path(probe),
            rec=Path(rec),
            nat=Path(nat),
            dec=Path(dec),
            cfg=Path(cfg),
            binary=Path(binary),
            data_dir=Path(data_dir),
            pdb=pdb,
            weight=weight_on,
        )
        if not off or not on:
            rows.append({"pdb": pdb, "status": "probe_fail"})
            continue

        dcf0 = off["dCF"]
        dcf1 = on["dCF"]
        native_min_off = off["native_min"]
        native_min_on = on["native_min"]
        moved_toward_native = dcf1 < dcf0 - 1e-9
        native_regressed = native_min_off and (
            on["cf_native"] > off["cf_native"] + 1e-6
        )
        native_min_lost = native_min_off and not native_min_on

        rows.append(
            {
                "pdb": pdb,
                "status": "ok",
                "config": str(cfg),
                "decoy": str(dec),
                "weight_on": weight_on,
                "cf_native_off": off["cf_native"],
                "cf_decoy_off": off["cf_decoy"],
                "cf_native_on": on["cf_native"],
                "cf_decoy_on": on["cf_decoy"],
                "dCF_off": dcf0,
                "dCF_on": dcf1,
                "delta_dCF": dcf1 - dcf0,
                "clash_native_off": off["clash_native"],
                "clash_decoy_off": off["clash_decoy"],
                "clash_native_on": on["clash_native"],
                "clash_decoy_on": on["clash_decoy"],
                "native_min_off": native_min_off,
                "native_min_on": native_min_on,
                "moved_toward_native": moved_toward_native,
                "native_regressed": native_regressed,
                "native_min_lost": native_min_lost,
                "identical_off_on": (
                    abs(off["cf_native"] - on["cf_native"]) < 1e-9
                    and abs(off["cf_decoy"] - on["cf_decoy"]) < 1e-9
                ),
            }
        )

    ok = [r for r in rows if r.get("status") == "ok"]
    accept_stats: Dict[str, Any] = {}
    if args.mode == "scoring-locked":
        verdict, reason, accept_stats = verdict_scoring_locked(
            ok, float(args.magnitude_floor), SIGN_FLIP_MIN_FRAC
        )
    else:
        verdict, reason = verdict_search_miss_legacy(ok)

    # SEARCH-MISS no-regression under ON weight (required for scoring-locked claim)
    sm_rows: List[Dict[str, Any]] = []
    sm_regress = 0
    if args.mode == "scoring-locked" and not args.skip_search_miss_check:
        for pdb in SEARCH_MISS_PANEL:
            cfg = resolve_dock_config(repo, pdb, elected_leaf)
            rec = resolve_receptor(repo, pdb)
            nat = resolve_native(repo, pdb)
            dec = resolve_decoy(repo, pdb, args.decoy_dir, elected_leaf)
            if not all((cfg, rec, nat, dec)):
                sm_rows.append({"pdb": pdb, "status": "skip_missing"})
                continue
            on = score_pair(
                Path(probe),
                rec=Path(rec),  # type: ignore[arg-type]
                nat=Path(nat),  # type: ignore[arg-type]
                dec=Path(dec),  # type: ignore[arg-type]
                cfg=Path(cfg),  # type: ignore[arg-type]
                binary=Path(binary),
                data_dir=Path(data_dir),
                pdb=pdb,
                weight=weight_on,
            )
            if not on:
                sm_rows.append({"pdb": pdb, "status": "probe_fail"})
                continue
            reg = not on["native_min"]
            if reg:
                sm_regress += 1
            sm_rows.append(
                {
                    "pdb": pdb,
                    "status": "ok",
                    "dCF_on": on["dCF"],
                    "native_min_on": on["native_min"],
                    "regressed": reg,
                    "cf_native_on": on["cf_native"],
                    "cf_decoy_on": on["cf_decoy"],
                }
            )
        if sm_regress > 0 and verdict == "PASS":
            verdict = "FAIL"
            reason = f"{reason}; SEARCH-MISS native CF-min lost on {sm_regress}/5"
        elif sm_regress > 0:
            reason = f"{reason}; SEARCH-MISS regress={sm_regress}/5"

    written = datetime.now(timezone.utc).isoformat()
    bin_sha = sha256_file(Path(binary))
    try:
        git = subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "HEAD"],
            text=True,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        git = "unknown"

    summary = {
        "schema": "pb_clash_burial_oracle/v2",
        "written_utc": written,
        "mode": args.mode,
        "roadmap": "ROADMAP_v2_PANEL_CORRECTION",
        "one_variable": f"FLEXAIDDS_PB_CLASH_WEIGHT={weight_on} (OFF=0)",
        "panel": panel,
        "class": "SCORING-LOCKED" if args.mode == "scoring-locked" else "SEARCH-MISS",
        "magnitude_floor_kcal": float(args.magnitude_floor)
        if args.mode == "scoring-locked"
        else None,
        "accept_stats": accept_stats,
        "n_scored": len(ok),
        "verdict": verdict,
        "reason": reason,
        "rows": rows,
        "search_miss_no_regress": sm_rows,
        "search_miss_regress_count": sm_regress,
        "binary": str(binary),
        "binary_sha256": bin_sha,
        "probe_cf": str(probe),
        "git": git,
        "elected_leaf": str(elected_leaf) if elected_leaf else None,
        "memetic_unlock": verdict == "PASS" and args.mode == "scoring-locked",
        "memetic_env_if_pass": "FLEXAIDDS_PB_CLASH_PHASE2_PASS=1",
        "prior_search_miss_oracle": "VOID (wrong panel + no magnitude floor)",
    }

    csv_path = out_dir / "pb_clash_oracle.csv"
    fields = [
        "pdb",
        "status",
        "cf_native_off",
        "cf_decoy_off",
        "cf_native_on",
        "cf_decoy_on",
        "dCF_off",
        "dCF_on",
        "delta_dCF",
        "decrease",
        "ge_floor",
        "sign_flip",
        "both_accept",
        "clash_decoy_off",
        "clash_decoy_on",
        "native_min_off",
        "native_min_on",
        "moved_toward_native",
        "native_regressed",
        "identical_off_on",
    ]
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    (out_dir / "pb_clash_oracle.json").write_text(json.dumps(summary, indent=2) + "\n")
    if sm_rows:
        (out_dir / "search_miss_no_regress.json").write_text(
            json.dumps(sm_rows, indent=2) + "\n"
        )

    lines = [
        f"# pb_clash oracle ({args.mode}) — {verdict}",
        "",
        f"**Written:** {written}  ",
        f"**Mode:** `{args.mode}` (ROADMAP_v2)  ",
        f"**One variable:** `FLEXAIDDS_PB_CLASH_WEIGHT={weight_on}` (OFF arm = 0)  ",
        f"**Panel:** {', '.join(panel)}  ",
        f"**OUT:** `{out_dir}`  ",
        f"**Binary sha256:** `{bin_sha}`  ",
        f"**Git:** `{git}`  ",
        "",
        f"**Verdict:** **{verdict}** — {reason}",
        "",
    ]
    if args.mode == "scoring-locked":
        lines += [
            "## ACCEPT (ROADMAP_v2 magnitude floor)",
            "",
            f"- dCF decreases by **≥ {args.magnitude_floor} kcal** AND **sign flip** "
            f"(OFF dCF>0 → ON dCF<0) on **≥{SIGN_FLIP_MIN_FRAC}/3** SCORING-LOCKED targets",
            "- No SEARCH-MISS native CF-min regression under ON weight",
            f"- Memetic unlock env if PASS: `FLEXAIDDS_PB_CLASH_PHASE2_PASS=1` "
            f"(still requires `FLEXAIDDS_MEMETIC=1`; claim default stays OFF)",
            "",
            f"| Metric | Value |",
            f"|--------|------:|",
            f"| Scored | {len(ok)} |",
            f"| decrease ≥ floor | {accept_stats.get('n_decrease_ge_floor', 0)} |",
            f"| sign flip | {accept_stats.get('n_sign_flip', 0)} |",
            f"| both | {accept_stats.get('n_both', 0)} |",
            f"| SEARCH-MISS regress | {sm_regress} |",
            "",
        ]
    else:
        lines += [
            "## VOID note",
            "",
            "Legacy SEARCH-MISS panel + sign-only ACCEPT is **VOID** under ROADMAP_v2.",
            "Do not unlock memetic from this path.",
            "",
        ]

    lines += [
        "## Per-target (primary panel)",
        "",
        "| PDB | dCF_off | dCF_on | ΔdCF | decrease | ge_floor | sign_flip | clash_dec OFF/ON |",
        "|-----|--------:|-------:|-----:|---------:|:--------:|:---------:|-----------------:|",
    ]
    for r in ok:
        decr = r.get("decrease", r["dCF_off"] - r["dCF_on"])
        lines.append(
            f"| {r['pdb']} | {r['dCF_off']:+.4f} | {r['dCF_on']:+.4f} | "
            f"{r['delta_dCF']:+.4f} | {decr:+.4f} | "
            f"{'Y' if r.get('ge_floor') else 'N'} | "
            f"{'Y' if r.get('sign_flip') else 'N'} | "
            f"{r.get('clash_decoy_off')}/{r.get('clash_decoy_on')} |"
        )

    if sm_rows:
        lines += [
            "",
            "## SEARCH-MISS no-regression (ON weight)",
            "",
            "| PDB | dCF_on | native_min | regressed |",
            "|-----|-------:|:----------:|:---------:|",
        ]
        for r in sm_rows:
            if r.get("status") != "ok":
                lines.append(f"| {r['pdb']} | — | — | {r.get('status')} |")
            else:
                lines.append(
                    f"| {r['pdb']} | {r['dCF_on']:+.4f} | "
                    f"{'Y' if r['native_min_on'] else 'N'} | "
                    f"{'Y' if r['regressed'] else 'N'} |"
                )

    lines += [
        "",
        "## Cadence",
        "",
        f"- Phase: 2b′ SCORING-LOCKED pb_clash (ROADMAP_v2)"
        if args.mode == "scoring-locked"
        else "- Phase: legacy SEARCH-MISS (VOID)",
        f"- One variable: PB_CLASH_WEIGHT={weight_on}",
        f"- **{verdict}**",
        "- Class-matched: scoring levers only on SCORING-LOCKED; sampling only on SEARCH-MISS",
        "",
    ]
    md_path = out_dir / "pb_clash_oracle.md"
    md_path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nWrote {csv_path}, {out_dir / 'pb_clash_oracle.json'}, {md_path}")
    # Exit 0 for PASS only; VOID_LEGACY_PASS and FAIL are non-zero for CI honesty
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
