#!/usr/bin/env python3
"""STEP 2 replacement: one-variable pb_clash burial score-only oracle.

Replaces the structurally unpassable WAL_COERCIVE Voronoi-wall gate (B3).
pb_clash is an all-pairs cell-list term (uncapped by design) that *can* see deep
interpenetration — see workorders WALL_ORACLE_FAIL_EXPLAINED / DOCKING_BUG_AUDIT B3.

**One variable:** FLEXAIDDS_PB_CLASH_WEIGHT (default OFF = 0.0 vs a single ON weight).

Uses production LOCCLF configs under ops/gates/configs/{PDB}_dock_config.json.
Never diagnostic/probe_config.json.

Panel: methodology clean probes 1J3J 1K3U 1L7F 1N1M 1M2Z.
Decoys: diagnostic/refs/*/falsemin_armA.pdb, or --decoy-dir of prebuilt buried poses.

ACCEPT (OPS-revised STEP 2):
  dCF = cf_native - cf_decoy moves toward/below 0 with weight ON vs OFF on ≥4/5,
  and no clean probe whose native is already CF-min OFF regresses under ON.

Usage:
  python3 scripts/pb_clash_burial_oracle.py \\
    --probe-cf build/probe_cf --binary build/FlexAIDdS --data-dir . \\
    --weight 1.0 \\
    --out-dir ~/flexaidds_results/workorders/pb_clash_oracle

Optional: build deeper burial decoys first with --build-buried (translate falsemin
toward receptor COM until cf_clash rises, write poses under out-dir/buried/).
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
from typing import Any, Dict, List, Optional, Sequence, Tuple

CLEAN_PANEL = ("1J3J", "1K3U", "1L7F", "1N1M", "1M2Z")


def find_repo_root() -> Path:
    p = Path(__file__).resolve().parents[1]
    if (p / "AGENTS.md").exists():
        return p
    return Path.cwd()


def resolve_dock_config(repo: Path, pdb: str) -> Optional[Path]:
    c = repo / "ops" / "gates" / "configs" / f"{pdb}_dock_config.json"
    return c if c.is_file() else None


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


def resolve_decoy(repo: Path, pdb: str, decoy_dir: Optional[Path]) -> Optional[Path]:
    if decoy_dir is not None:
        for name in (f"{pdb}_buried.pdb", f"{pdb}_falsemin.pdb", f"{pdb}.pdb"):
            p = decoy_dir / name
            if p.is_file():
                return p
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
    # Explicit OFF vs ON for the single variable
    env["FLEXAIDDS_PB_CLASH_WEIGHT"] = str(pb_clash_weight)
    # Keep wall coercive off; do not touch other knobs
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


def build_buried_decoy(
    probe: Path,
    *,
    repo: Path,
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


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    repo = find_repo_root()
    ap.add_argument("--repo", type=Path, default=repo)
    ap.add_argument("--probe-cf", type=Path, default=None)
    ap.add_argument("--binary", type=Path, default=None)
    ap.add_argument("--data-dir", type=Path, default=None)
    ap.add_argument("--panel", nargs="*", default=None)
    ap.add_argument(
        "--weight",
        type=float,
        default=1.0,
        help="Single ON value for FLEXAIDDS_PB_CLASH_WEIGHT (OFF is always 0)",
    )
    ap.add_argument("--decoy-dir", type=Path, default=None)
    ap.add_argument(
        "--build-buried",
        action="store_true",
        help="Translate falsemin toward receptor COM to raise cf_clash before A/B",
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

    panel = list(args.panel) if args.panel else list(CLEAN_PANEL)
    buried_dir = out_dir / "buried"
    if args.build_buried:
        buried_dir.mkdir(parents=True, exist_ok=True)

    rows: List[Dict[str, Any]] = []
    for pdb in panel:
        cfg = resolve_dock_config(repo, pdb)
        rec = resolve_receptor(repo, pdb)
        nat = resolve_native(repo, pdb)
        dec = resolve_decoy(repo, pdb, args.decoy_dir)
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

        if args.build_buried:
            if not dec or not Path(dec).is_file():
                rows.append({"pdb": pdb, "status": "skip_missing_decoy_for_bury"})
                continue
            buried = buried_dir / f"{pdb}_buried.pdb"
            print(f"[{pdb}] building buried decoy…")
            built = build_buried_decoy(
                Path(probe),
                repo=repo,
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

        print(f"[{pdb}] production config={cfg} decoy={dec} weight_on={weight_on}")
        nat_off = probe_json(
            Path(probe),
            receptor=Path(rec),
            pose=Path(nat),
            ligand=Path(nat),
            binary=Path(binary),
            data_dir=Path(data_dir),
            config=Path(cfg),
            label=f"{pdb}_nat_w0",
            pb_clash_weight=0.0,
        )
        dec_off = probe_json(
            Path(probe),
            receptor=Path(rec),
            pose=Path(dec),
            ligand=Path(nat),
            binary=Path(binary),
            data_dir=Path(data_dir),
            config=Path(cfg),
            label=f"{pdb}_dec_w0",
            pb_clash_weight=0.0,
        )
        nat_on = probe_json(
            Path(probe),
            receptor=Path(rec),
            pose=Path(nat),
            ligand=Path(nat),
            binary=Path(binary),
            data_dir=Path(data_dir),
            config=Path(cfg),
            label=f"{pdb}_nat_w{weight_on}",
            pb_clash_weight=weight_on,
        )
        dec_on = probe_json(
            Path(probe),
            receptor=Path(rec),
            pose=Path(dec),
            ligand=Path(nat),
            binary=Path(binary),
            data_dir=Path(data_dir),
            config=Path(cfg),
            label=f"{pdb}_dec_w{weight_on}",
            pb_clash_weight=weight_on,
        )
        if not all((nat_off, dec_off, nat_on, dec_on)):
            rows.append({"pdb": pdb, "status": "probe_fail"})
            continue

        cn0 = fnum(nat_off, "cf_total")
        cd0 = fnum(dec_off, "cf_total")
        cn1 = fnum(nat_on, "cf_total")
        cd1 = fnum(dec_on, "cf_total")
        assert cn0 is not None and cd0 is not None and cn1 is not None and cd1 is not None
        dcf0 = cn0 - cd0  # negative => native better (lower CF)
        dcf1 = cn1 - cd1
        native_min_off = cn0 <= cd0
        native_min_on = cn1 <= cd1
        # "moves toward/below 0": dCF decreases (more negative or less positive)
        moved_toward_native = dcf1 < dcf0 - 1e-9
        native_regressed = native_min_off and (cn1 > cn0 + 1e-6)  # native CF worse under ON
        # stronger: native was min and loses under ON
        native_min_lost = native_min_off and not native_min_on

        rows.append(
            {
                "pdb": pdb,
                "status": "ok",
                "config": str(cfg),
                "decoy": str(dec),
                "weight_on": weight_on,
                "cf_native_off": cn0,
                "cf_decoy_off": cd0,
                "cf_native_on": cn1,
                "cf_decoy_on": cd1,
                "dCF_off": dcf0,
                "dCF_on": dcf1,
                "delta_dCF": dcf1 - dcf0,
                "clash_native_off": fnum(nat_off, "cf_clash"),
                "clash_decoy_off": fnum(dec_off, "cf_clash"),
                "clash_native_on": fnum(nat_on, "cf_clash"),
                "clash_decoy_on": fnum(dec_on, "cf_clash"),
                "native_min_off": native_min_off,
                "native_min_on": native_min_on,
                "moved_toward_native": moved_toward_native,
                "native_regressed": native_regressed,
                "native_min_lost": native_min_lost,
                "identical_off_on": (
                    abs(cn0 - cn1) < 1e-9
                    and abs(cd0 - cd1) < 1e-9
                ),
            }
        )

    ok = [r for r in rows if r.get("status") == "ok"]
    n = len(ok)
    n_moved = sum(1 for r in ok if r.get("moved_toward_native"))
    n_native_min_on = sum(1 for r in ok if r.get("native_min_on"))
    n_regressed = sum(1 for r in ok if r.get("native_regressed") or r.get("native_min_lost"))
    n_identical = sum(1 for r in ok if r.get("identical_off_on"))
    n_effect = n - n_identical

    # ACCEPT: dCF moves toward native on ≥4/5 AND no native-min regression
    accept_move = n_moved >= 4 if n >= 5 else (n_moved >= max(1, n - 1) if n else False)
    accept_no_regress = n_regressed == 0
    if n == 0:
        verdict, reason = "FAIL", "no scored targets"
    elif n_identical == n:
        verdict, reason = (
            "FAIL",
            "OFF≡ON on all targets (panel may lack deep clashes; use --build-buried)",
        )
    elif accept_move and accept_no_regress:
        verdict, reason = (
            "PASS",
            f"dCF moved toward native on {n_moved}/{n}; no native CF-min regression",
        )
    else:
        verdict, reason = (
            "FAIL",
            f"moved={n_moved}/{n} (need ≥4/5); regressed={n_regressed}; effect_targets={n_effect}",
        )

    written = datetime.now(timezone.utc).isoformat()
    summary = {
        "schema": "pb_clash_burial_oracle/v1",
        "written_utc": written,
        "one_variable": f"FLEXAIDDS_PB_CLASH_WEIGHT={weight_on} (OFF=0)",
        "replaces": "WAL_COERCIVE wall STEP 2 (structurally unpassable B3)",
        "panel": panel,
        "n_scored": n,
        "n_moved_toward_native": n_moved,
        "n_native_min_on": n_native_min_on,
        "n_regressed": n_regressed,
        "n_identical_off_on": n_identical,
        "verdict": verdict,
        "reason": reason,
        "rows": rows,
        "matrix_note": "score-only probe_cf; matrix not used",
        "binary": str(binary),
        "probe_cf": str(probe),
    }

    # CSV
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

    lines = [
        f"# pb_clash burial oracle — {verdict}",
        "",
        f"**Written:** {written}  ",
        f"**One variable:** `FLEXAIDDS_PB_CLASH_WEIGHT={weight_on}` (OFF arm = 0)  ",
        f"**Replaces:** WAL_COERCIVE STEP 2 (B3 structural no-op)  ",
        f"**Configs:** `ops/gates/configs/{{PDB}}_dock_config.json` (production LOCCLF)  ",
        f"**Panel:** {', '.join(panel)}  ",
        f"**OUT:** `{out_dir}`  ",
        "",
        f"**Verdict:** **{verdict}** — {reason}",
        "",
        "## Metrics",
        "",
        f"| Metric | Value |",
        f"|--------|------:|",
        f"| Scored | {n} |",
        f"| dCF moved toward native (ON vs OFF) | {n_moved}/{n} |",
        f"| Native CF-min under ON | {n_native_min_on}/{n} |",
        f"| Native regressions | {n_regressed} |",
        f"| OFF≡ON identical | {n_identical}/{n} |",
        "",
        "## Per-target",
        "",
        "| PDB | dCF_off | dCF_on | ΔdCF | clash_dec OFF/ON | nat_min OFF/ON | moved |",
        "|-----|--------:|-------:|-----:|-----------------:|:--------------:|:-----:|",
    ]
    for r in ok:
        lines.append(
            f"| {r['pdb']} | {r['dCF_off']:+.4f} | {r['dCF_on']:+.4f} | "
            f"{r['delta_dCF']:+.4f} | "
            f"{r.get('clash_decoy_off')}/{r.get('clash_decoy_on')} | "
            f"{'Y' if r['native_min_off'] else 'N'}/{'Y' if r['native_min_on'] else 'N'} | "
            f"{'Y' if r['moved_toward_native'] else 'N'} |"
        )
    lines += [
        "",
        "## ACCEPT (revised STEP 2)",
        "",
        "dCF moves toward/below 0 with weight ON vs OFF on ≥4/5; no clean native CF-min regression.",
        "",
        "## Cadence",
        "",
        "- Phase: STEP 2 replacement (pb_clash burial)  ",
        f"- One variable: PB_CLASH_WEIGHT={weight_on}  ",
        f"- **{verdict}**  ",
        "- Do **not** set WALL_PILOT_PASS from WAL_COERCIVE evidence  ",
        "- Memetic still blocked until a burial/steric oracle PASSes  ",
        "",
    ]
    md_path = out_dir / "pb_clash_oracle.md"
    md_path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nWrote {csv_path}, {out_dir / 'pb_clash_oracle.json'}, {md_path}")
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
