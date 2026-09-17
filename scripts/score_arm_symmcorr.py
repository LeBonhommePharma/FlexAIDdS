#!/usr/bin/env python3
"""Flat-arm symmetry-corrected RMSD sidecar driver (JOB B / claim-metric layer).

Scores three endpoints against an explicit crystal-reference directory using
``scripts/rmsd_symmcorr.py`` (imported, never edited):

  S1       elected rank-0 head
  S_top10  min symmcorr over ranks 0..9 of the elected restart
  BCR      min symmcorr over every emitted cluster head

Arm layout is flat: ``<arm_dir>/<TARGET>/result.csv``. Crystal refs:
``<crystal_ref_dir>/<TARGET>/<TARGET>_ligand.sdf``.

Resumable and atomic: one CSV per target under ``--out-dir``, skip if present,
write ``*.tmp`` then ``os.replace``. Propagates ``SymmCorrUnavailable`` —
never substitutes a weaker metric.

Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import math
import os
import re
import sys
from pathlib import Path
from typing import Any, Optional

# Allow ``python scripts/score_arm_symmcorr.py`` and import-as-module.
_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

import rmsd_symmcorr as symm  # noqa: E402

RMSD_SUCCESS_A = 2.0
TOP_N_MODES = 10
PRODUCER = "scripts/score_arm_symmcorr.py"

# Filename patterns (stem without .pdb) — same contract as parse_flexaid_arm_results.
_RESTART_FO = re.compile(
    r"^(?P<pdb>[A-Za-z0-9]+)_r(?P<restart>\d+)_(?P<minpts>\d+)_(?P<rank>\d+)$",
    re.I,
)
_RESTART_CF = re.compile(
    r"^(?P<pdb>[A-Za-z0-9]+)_r(?P<restart>\d+)_(?P<rank>\d+)$",
    re.I,
)
_FO_DUAL = re.compile(
    r"^(?P<pdb>[A-Za-z0-9]+)_(?P<minpts>\d+)_(?P<rank>\d+)$",
    re.I,
)
_CF_RANK = re.compile(
    r"^(?P<pdb>[A-Za-z0-9]+)_(?P<rank>\d+)$",
    re.I,
)

TARGET_FIELDS = (
    "pdb_id",
    "rmsd_symmcorr_s1",
    "status_s1",
    "pose_s1",
    "pose_sha256_s1",
    "sha_verified_s1",
    "rmsd_symmcorr_s_top10",
    "status_s_top10",
    "n_top10_scored",
    "rmsd_symmcorr_bcr",
    "status_bcr",
    "n_bcr_scored",
    "csv_rmsd_top1",
    "csv_rmsd_bcr",
    "csv_mode_rmsd_min",
    "success_s1_csv",
    "success_s1_symm",
    "success_s_top10_csv",
    "success_s_top10_symm",
    "success_s3_csv",
    "success_bcr_symm",
    "flip_s1",
    "flip_s_top10",
    "flip_bcr",
    "metric",
    "producer",
    "spyrmsd_version",
    "crystal_sdf",
    "n_heavy",
)


def parse_pose_filename(path: Path, pdb: str) -> Optional[tuple[Optional[int], int]]:
    """Return (restart_or_None, emitted_rank) or None if not a cluster head."""
    stem = path.stem
    if stem.upper().endswith("_INI") or stem.upper() == f"{pdb.upper()}_INI":
        return None
    if not stem.upper().startswith(pdb.upper()):
        return None
    for pat in (_RESTART_FO, _RESTART_CF, _FO_DUAL, _CF_RANK):
        m = pat.match(stem)
        if not m or m.group("pdb").upper() != pdb.upper():
            continue
        restart = int(m.group("restart")) if "restart" in m.groupdict() and m.group("restart") is not None else None
        return restart, int(m.group("rank"))
    return None


def collect_heads(out_dir: Path, pdb: str) -> list[tuple[Optional[int], int, Path]]:
    rows: list[tuple[Optional[int], int, Path]] = []
    for path in sorted(out_dir.glob(f"{pdb}*.pdb")):
        parsed = parse_pose_filename(path, pdb)
        if parsed is None:
            continue
        restart, rank = parsed
        rows.append((restart, rank, path))
    # Also accept gzipped cluster heads next to the uncompressed name.
    for path in sorted(out_dir.glob(f"{pdb}*.pdb.gz")):
        stem_path = Path(str(path)[:-3])  # strip .gz → *.pdb path identity
        parsed = parse_pose_filename(stem_path, pdb)
        if parsed is None:
            continue
        restart, rank = parsed
        # Prefer uncompressed if both exist
        if any(p == stem_path for _, _, p in rows):
            continue
        rows.append((restart, rank, path))
    return rows


def _finite(text: str) -> Optional[float]:
    text = (text or "").strip()
    if not text or text.upper() in {"NA", "NAN", "NONE"}:
        return None
    try:
        value = float(text)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def _success(rmsd: Optional[float]) -> str:
    if rmsd is None:
        return ""
    return "1" if rmsd <= RMSD_SUCCESS_A else "0"


def _flip(csv_flag: str, sym_flag: str) -> str:
    if csv_flag not in {"0", "1"} or sym_flag not in {"0", "1"}:
        return ""
    return "1" if csv_flag != sym_flag else "0"


def load_result_row(csv_path: Path) -> dict[str, str]:
    with csv_path.open(newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    if len(rows) != 1:
        raise ValueError(f"{csv_path}: expected exactly one observation, got {len(rows)}")
    return {k: (v if v is not None else "") for k, v in rows[0].items()}


def crystal_sdf_for(crystal_ref_dir: Path, pdb_id: str) -> Path:
    return crystal_ref_dir / pdb_id / f"{pdb_id}_ligand.sdf"


def score_pose(crystal_sdf: Path, pose_path: Path) -> dict[str, Any]:
    """Score one pose artifact. Propagates SymmCorrUnavailable."""
    raw, actual = symm.read_pose_artifact(str(pose_path))
    out: dict[str, Any] = {
        "rmsd": None,
        "status": "",
        "pose": actual or str(pose_path),
        "pose_sha256": "",
        "sha_verified": "0",
        "n_heavy": 0,
    }
    if raw is None:
        out["status"] = "no_pose_artifact"
        return out
    out["pose_sha256"] = hashlib.sha256(raw).hexdigest()
    if not crystal_sdf.is_file():
        out["status"] = "no_crystal_sdf"
        return out
    res = symm.symmcorr_rmsd(str(crystal_sdf), raw)
    out["status"] = res["status"]
    out["n_heavy"] = res["n_heavy"]
    out["rmsd"] = res["rmsd"]
    return out


def elect_restart(heads: list[tuple[Optional[int], int, Path]], elected_path: str) -> Optional[int]:
    """Infer elected restart from elected_path filename, else first restart with rank-0."""
    if elected_path:
        base = Path(elected_path).name
        if base.endswith(".gz"):
            base = base[:-3]
        elected_stem = Path(base).stem
        for restart, _rank, path in heads:
            # Compare against the *elected* basename only. Building `names`
            # from `path.name` made `path.name in names` true for every head,
            # so the first glob hit (restart 0) always won and S_top10 scored
            # the wrong restart.
            if path.name == base or path.stem == elected_stem or path.name == base + ".gz":
                return restart
    restarts = sorted({r for r, rank, _ in heads if r is not None and rank == 0})
    if restarts:
        return restarts[0]
    return None


def select_s1_pose(
    heads: list[tuple[Optional[int], int, Path]],
    elected_path: str,
    target_dir: Path,
) -> Optional[Path]:
    if elected_path:
        cand = Path(elected_path)
        if not cand.is_absolute():
            cand = target_dir / cand.name
        raw, actual = symm.read_pose_artifact(str(cand))
        if raw is not None:
            return Path(actual or cand)
        # try basename in target dir
        raw, actual = symm.read_pose_artifact(str(target_dir / Path(elected_path).name))
        if raw is not None:
            return Path(actual or (target_dir / Path(elected_path).name))
    # Fall back: first rank-0 head (restart order)
    rank0 = [(r, p) for r, rank, p in heads if rank == 0]
    if not rank0:
        return None
    rank0.sort(key=lambda t: (999999 if t[0] is None else t[0], str(t[1])))
    return rank0[0][1]


def select_top10_poses(
    heads: list[tuple[Optional[int], int, Path]],
    elected_restart: Optional[int],
) -> list[Path]:
    restarts_present = {r for r, _, _ in heads if r is not None}
    if restarts_present:
        target = elected_restart if elected_restart in restarts_present else min(restarts_present)
        group = [(rank, path) for r, rank, path in heads if r == target]
    else:
        group = [(rank, path) for r, rank, path in heads if r is None]
    by_rank: dict[int, Path] = {}
    for rank, path in group:
        if 0 <= rank < TOP_N_MODES and rank not in by_rank:
            by_rank[rank] = path
    return [by_rank[i] for i in range(TOP_N_MODES) if i in by_rank]


def min_mode_rmsd_csv(row: dict[str, str]) -> Optional[float]:
    vals = []
    for i in range(TOP_N_MODES):
        v = _finite(row.get(f"mode_rmsd_{i}", ""))
        if v is not None:
            vals.append(v)
    return min(vals) if vals else None


def score_target(
    target_dir: Path,
    crystal_ref_dir: Path,
    spy_version: str,
) -> dict[str, str]:
    pdb_id = target_dir.name.strip().upper()
    result_csv = target_dir / "result.csv"
    if not result_csv.is_file():
        raise FileNotFoundError(f"missing result.csv: {result_csv}")
    row = load_result_row(result_csv)
    csv_pdb = (row.get("pdb_id") or row.get("pdb") or row.get("target") or "").strip().upper()
    if csv_pdb and csv_pdb != pdb_id:
        # Flat layout uses directory name; warn via status later if mismatch.
        pdb_id = pdb_id  # keep directory identity for pose globs
    if not csv_pdb:
        csv_pdb = pdb_id

    crystal = crystal_sdf_for(crystal_ref_dir, pdb_id)
    heads = collect_heads(target_dir, pdb_id)
    elected_path = (row.get("elected_path") or row.get("elected_pose_path") or "").strip()
    elected_restart = elect_restart(heads, elected_path)

    s1_path = select_s1_pose(heads, elected_path, target_dir)
    top10_paths = select_top10_poses(heads, elected_restart)
    bcr_paths = [p for _, _, p in heads]

    out: dict[str, str] = {k: "" for k in TARGET_FIELDS}
    out.update(
        {
            "pdb_id": csv_pdb,
            "metric": symm.METRIC_NAME,
            "producer": PRODUCER,
            "spyrmsd_version": spy_version,
            "crystal_sdf": str(crystal),
            "csv_rmsd_top1": (row.get("rmsd_top1") or "").strip(),
            "csv_rmsd_bcr": (row.get("rmsd_bcr") or "").strip(),
            "success_s1_csv": (row.get("success_s1") or "").strip(),
            "success_s_top10_csv": (row.get("success_s_top10") or "").strip(),
            "success_s3_csv": (row.get("success_s3") or "").strip(),
        }
    )
    mode_min = min_mode_rmsd_csv(row)
    out["csv_mode_rmsd_min"] = "" if mode_min is None else f"{mode_min:.4f}"

    # --- S1 ---
    if s1_path is None:
        out["status_s1"] = "no_s1_pose"
    else:
        s1 = score_pose(crystal, s1_path)
        out["status_s1"] = s1["status"]
        out["pose_s1"] = s1["pose"]
        out["pose_sha256_s1"] = s1["pose_sha256"]
        out["n_heavy"] = str(s1["n_heavy"] or "")
        declared = (row.get("pose_sha256") or "").strip().lower()
        if declared and s1["pose_sha256"] and declared == s1["pose_sha256"]:
            out["sha_verified_s1"] = "1"
        else:
            out["sha_verified_s1"] = "0"
        if s1["rmsd"] is not None:
            out["rmsd_symmcorr_s1"] = f"{s1['rmsd']:.5f}"

    # --- S_top10 ---
    top10_vals: list[float] = []
    top10_status = "ok"
    if not top10_paths:
        top10_status = "no_top10_poses"
    for path in top10_paths:
        res = score_pose(crystal, path)
        if res["rmsd"] is None:
            if top10_status == "ok":
                top10_status = res["status"] or "top10_pose_failed"
            continue
        top10_vals.append(float(res["rmsd"]))
    out["n_top10_scored"] = str(len(top10_vals))
    out["status_s_top10"] = top10_status if top10_vals or top10_status != "ok" else "ok"
    if top10_vals:
        out["rmsd_symmcorr_s_top10"] = f"{min(top10_vals):.5f}"
        out["status_s_top10"] = "ok"

    # --- BCR ---
    bcr_vals: list[float] = []
    bcr_status = "ok"
    if not bcr_paths:
        bcr_status = "no_cluster_heads"
    for path in bcr_paths:
        res = score_pose(crystal, path)
        if res["rmsd"] is None:
            if bcr_status == "ok":
                bcr_status = res["status"] or "bcr_pose_failed"
            continue
        bcr_vals.append(float(res["rmsd"]))
    out["n_bcr_scored"] = str(len(bcr_vals))
    out["status_bcr"] = bcr_status if bcr_vals or bcr_status != "ok" else "ok"
    if bcr_vals:
        out["rmsd_symmcorr_bcr"] = f"{min(bcr_vals):.5f}"
        out["status_bcr"] = "ok"

    s1_v = _finite(out["rmsd_symmcorr_s1"])
    top_v = _finite(out["rmsd_symmcorr_s_top10"])
    bcr_v = _finite(out["rmsd_symmcorr_bcr"])
    out["success_s1_symm"] = _success(s1_v)
    out["success_s_top10_symm"] = _success(top_v)
    out["success_bcr_symm"] = _success(bcr_v)
    out["flip_s1"] = _flip(out["success_s1_csv"], out["success_s1_symm"])
    out["flip_s_top10"] = _flip(out["success_s_top10_csv"], out["success_s_top10_symm"])
    out["flip_bcr"] = _flip(out["success_s3_csv"], out["success_bcr_symm"])
    return out


def atomic_write_csv(path: Path, row: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(TARGET_FIELDS))
        writer.writeheader()
        writer.writerow({k: row.get(k, "") for k in TARGET_FIELDS})
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def discover_targets(arm_dir: Path) -> list[Path]:
    targets = []
    for child in sorted(arm_dir.iterdir()):
        if child.is_dir() and (child / "result.csv").is_file():
            targets.append(child)
    return targets


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arm-dir", type=Path, required=True, help="Flat arm directory containing <TARGET>/result.csv")
    ap.add_argument(
        "--crystal-ref-dir",
        type=Path,
        required=True,
        help="Directory with <TARGET>/<TARGET>_ligand.sdf (explicit; no hard-coded cache)",
    )
    ap.add_argument("--out-dir", type=Path, required=True, help="Sidecar output directory (one CSV per target)")
    ap.add_argument("--targets", help="Optional comma-separated target allow-list")
    ap.add_argument("--force", action="store_true", help="Rescore targets even if sidecar CSV exists")
    args = ap.parse_args(argv)

    if not args.arm_dir.is_dir():
        print(f"error: arm dir not found: {args.arm_dir}", file=sys.stderr)
        return 2
    if not args.crystal_ref_dir.is_dir():
        print(f"error: crystal ref dir not found: {args.crystal_ref_dir}", file=sys.stderr)
        return 2

    try:
        _, spy_version = symm._spyrmsd()
    except symm.SymmCorrUnavailable as exc:
        print(f"[SYMMCORR-FAIL] {exc}", file=sys.stderr)
        raise

    allow = None
    if args.targets:
        allow = {t.strip().upper() for t in args.targets.split(",") if t.strip()}

    targets = discover_targets(args.arm_dir)
    if allow is not None:
        targets = [t for t in targets if t.name.upper() in allow]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    scored = skipped = failed = 0
    for target_dir in targets:
        out_csv = args.out_dir / f"{target_dir.name.upper()}.csv"
        if out_csv.is_file() and not args.force:
            skipped += 1
            print(f"skip {target_dir.name} (exists)")
            continue
        try:
            row = score_target(target_dir, args.crystal_ref_dir, spy_version)
            atomic_write_csv(out_csv, row)
            scored += 1
            print(
                f"ok {row['pdb_id']}: S1={row['rmsd_symmcorr_s1'] or row['status_s1']} "
                f"top10={row['rmsd_symmcorr_s_top10'] or row['status_s_top10']} "
                f"bcr={row['rmsd_symmcorr_bcr'] or row['status_bcr']}"
            )
        except symm.SymmCorrUnavailable:
            raise
        except Exception as exc:
            failed += 1
            print(f"FAIL {target_dir.name}: {exc}", file=sys.stderr)

    print(f"done scored={scored} skipped={skipped} failed={failed} denom={len(targets)}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except symm.SymmCorrUnavailable as exc:
        print(f"[SYMMCORR-FAIL] {exc}", file=sys.stderr)
        raise SystemExit(2) from exc
