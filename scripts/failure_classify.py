#!/usr/bin/env python3
"""
failure_classify.py — classify FAILED benchmark targets into causal failure modes.

After each run, assigns every FAILED target (rmsd >= lo_rmsd threshold) a failure
mode tag so the next run's priority queue can be populated based on which failures
a given fix is causally linked to.

Usage:
    python3 scripts/failure_classify.py <result_dir> [--lo-rmsd 2.0] [--hi-near-miss 2.5]

Outputs:
    <result_dir>/failure_modes.json    — per-target classification
    <result_dir>/failure_summary.txt   — human-readable summary table
    stdout: summary table

## Failure Mode Taxonomy (priority order — first match wins)

1. seed_echo       — GA echoed the oracle seed unmodified (no real search occurred).
                     Detected via explicit seed_echo column (0/1) or reconstructed
                     from |best_score − cf_native| ≤ 0.01 (DatasetRunner.cpp logic).
2. selection_miss  — oracle min RMSD < lo_rmsd BUT top-1 reported RMSD >= lo_rmsd.
                     The near-native pose EXISTED in the ensemble but cluster
                     selection chose the wrong one.  This is what Z+H composite fixes.
3. harness_artifact — data/runner red flags (v46 cache, cf_native mismatch, etc.)
4. CF_scoring_failure_deep — RMSD >= 6 Å AND best_score < cf_native - cf_delta
5. CF_scoring_failure_near — 2 <= RMSD < 6 Å with scoring/selection pathology
6. CF_false_minimum — legacy alias; deep/near modes subsume this when RMSD known
7. selection_miss  — (moved after harness) oracle had sub-2 Å pose, selector did not
8. timeout         — wall_time_s >= 900 s, or RMSD sentinel < 0 (−1.0 = not computed/failed).
9. search_failure  — everything else; sub-tagged by DoF:
                       :rigid    num_genes ≤ 4   (rigid-body only — scoring suspect)
                       :flexible num_genes ≥ 8   (high-flex — search space explosion)
                       :medium   5 ≤ num_genes ≤ 7 (default)
                       (no tag)  num_genes not determinable

## ⚠  Column-semantics correction (spec vs. DatasetRunner.h ground truth)
The original spec documentation had `rmsd_hungarian` and `best_cluster_rmsd`
semantics SWAPPED.  DatasetRunner.h lines 129 and 147 are authoritative:

    rmsd_hungarian    → top-1 SELECTED pose RMSD (Hungarian symmetry-corrected)
                        = THE benchmark / reported metric;
                        failure is defined as rmsd_hungarian >= lo_rmsd.
    best_cluster_rmsd → oracle minimum RMSD across ALL emitted cluster poses (0–19)
                        = what a PERFECT selector could have achieved.

Therefore selection_miss is:
    best_cluster_rmsd < lo_rmsd  AND  rmsd_hungarian >= lo_rmsd
NOT the inverted condition written in the spec.

Copyright 2026 Le Bonhomme Pharma.  Apache-2.0.
"""

import argparse
import csv
import glob
import json
import math
import os
import re
import sys

# ── Tunables ──────────────────────────────────────────────────────────────────
DEFAULT_LO_RMSD        = 2.0   # Å — failure threshold (rmsd_hungarian >= this)
DEFAULT_HI_NEAR_MISS   = 2.5   # Å — near-miss band upper bound (for summary)
DEFAULT_CF_DELTA       = 5.0   # CF units — false-minimum gap threshold
TIMEOUT_WALL_S         = 900.0 # seconds — hard timeout
SEED_ECHO_CF_TOL       = 0.01  # CF units — tolerance for seed_echo reconstruction
RIGID_GENES_MAX        = 4     # num_genes threshold for :rigid tag
FLEX_GENES_MIN         = 8     # num_genes threshold for :flexible tag


# ── Low-level helpers ─────────────────────────────────────────────────────────

def _fget(row, keys, default=None):
    """Return the first non-empty value from `row` that matches any of `keys`."""
    for k in keys:
        v = row.get(k)
        if v is not None and str(v).strip() not in ("", "N/A", "nan"):
            return str(v).strip()
    return default


def _float(row, keys, default=None):
    raw = _fget(row, keys)
    if raw is None:
        return default
    try:
        v = float(raw)
        return v if math.isfinite(v) else default
    except (ValueError, TypeError):
        return default


def _int(row, keys, default=None):
    raw = _fget(row, keys)
    if raw is None:
        return default
    try:
        return int(float(raw))
    except (ValueError, TypeError):
        return default


SEED_ECHO_RMSD_TOL = 0.5  # Å — max pose RMSD for seed_echo reconstruction to fire


def _seed_echo(row):
    """Detect seed_echo from explicit column, or reconstruct from CF proximity.

    DatasetRunner.cpp logic (line 5142):
        seed_echo = (cf_native != 0) && isfinite(best_cf)
                    && (fabs(best_cf - cf_native) <= 0.01f)

    ⚠ False-positive risk (CF degeneracy): when the scoring landscape has many
    near-degenerate minima, unrelated poses can satisfy |best_cf − cf_native| ≤ 0.01
    without the GA having literally echoed the oracle seed.

    Guard (reconstruction path only):
        Require ALSO rmsd_to_crystal < SEED_ECHO_RMSD_TOL (0.5 Å).
        A true seed echo of the oracle placement means the reported pose IS the seed,
        which should be near-native.  High RMSD + CF match → CF degeneracy, not echo.
    """
    se_raw = _fget(row, ["seed_echo"])
    if se_raw is not None:
        sl = se_raw.lower()
        if sl in ("1", "true", "yes"):
            return True
        if sl in ("0", "false", "no"):
            return False
        # Unknown value — fall through to reconstruction

    # Reconstruct: |best_score − cf_native| ≤ SEED_ECHO_CF_TOL
    #              AND rmsd_to_crystal < SEED_ECHO_RMSD_TOL
    # The RMSD gate eliminates CF-degeneracy false positives.
    cf_native   = _float(row, ["cf_native"])
    best_score  = _float(row, ["best_score", "best_cf"])
    rmsd        = _float(row, ["rmsd_hungarian", "rmsd_to_crystal", "best_cluster_rmsd"])
    if cf_native is not None and best_score is not None and cf_native != 0.0:
        if abs(best_score - cf_native) <= SEED_ECHO_CF_TOL:
            # Only fire if the pose is near-native; otherwise it's CF degeneracy.
            if rmsd is not None and rmsd < SEED_ECHO_RMSD_TOL:
                return True
    return False


def _harness_flags(row, target_dir, rmsd):
    """Runner/data red flags — scoring re-tune must not proceed on these."""
    flags = []
    count_delta = _int(row, ["count_delta"])
    if count_delta is not None and count_delta > 2:
        flags.append("v46_cache_count_delta")
    if rmsd is not None and rmsd >= 999.0:
        flags.append("rmsd_sentinel_999")
    cf_native = _float(row, ["cf_native"])
    best_score = _float(row, ["best_score", "best_cf"])
    if cf_native is not None and best_score is not None:
        if abs(best_score - cf_native) <= SEED_ECHO_CF_TOL and rmsd is not None and rmsd >= 6.0:
            flags.append("cf_match_high_rmsd_native_ic")
    if cf_native == 0.0 and _int(row, ["num_poses"], 0) > 0:
        flags.append("v51_cf_native_zero")
    if cf_native is not None and cf_native > 1000.0:
        flags.append("cf_native_wall_clash_overflow")
    if cf_native is not None and cf_native < -50000.0:
        flags.append("cf_native_sanity_underflow")
    se_raw = _fget(row, ["seed_echo"])
    if se_raw in ("1", "true", "yes") and cf_native is not None and cf_native > 0:
        flags.append("seed_echo_clash_attractor")
    return flags


def _parse_num_genes_from_dir(target_dir):
    """Infer num_genes from per-target files:
    1. dock_config.json: force_rigid=true → 4
    2. stdout.log:  'N rotatable bonds' → 7 + N
       (7 rigid-body genes: tx, ty, tz + 4-component rotation)

    Returns int or None if undetermined.
    """
    force_rigid = False
    cfg_path = os.path.join(target_dir, "dock_config.json")
    if os.path.isfile(cfg_path):
        try:
            with open(cfg_path) as fh:
                cfg = json.load(fh)
            force_rigid = cfg.get("flexibility", {}).get("force_rigid", False)
        except (OSError, json.JSONDecodeError):
            pass

    if force_rigid:
        return 4  # only tx, ty, tz, rot_angle (hip_eval.hip semantics)

    log_path = os.path.join(target_dir, "stdout.log")
    if not os.path.isfile(log_path):
        return None
    try:
        with open(log_path) as fh:
            for line in fh:
                m = re.search(r"(\d+)\s+rotatable\s+bond", line)
                if m:
                    return 7 + int(m.group(1))
    except OSError:
        pass
    return None


# ── Per-target classification ─────────────────────────────────────────────────

def classify_one(pdb_id, row, target_dir, lo_rmsd, cf_delta):
    """Classify a single target row.  Returns a result dict."""

    # ------------------------------------------------------------------
    # Benchmark RMSD: rmsd_hungarian = top-1 selected pose (the metric).
    # Oracle   RMSD: best_cluster_rmsd = best possible across all poses.
    # See module docstring for column-semantics correction vs. spec.
    # ------------------------------------------------------------------
    rmsd   = _float(row, ["rmsd_hungarian", "best_cluster_rmsd"], default=-1.0)
    oracle = _float(row, ["best_cluster_rmsd", "rmsd_hungarian"],  default=-1.0)

    # If both columns are equal (e.g. older CSVs where oracle wasn't tracked
    # separately), oracle_rmsd == rmsd — selection_miss cannot be detected.
    # We still try, but note it in missing_cols later.

    success = (rmsd >= 0.0) and (rmsd < lo_rmsd)

    out = {
        "rmsd":           round(rmsd,   4),
        "oracle_rmsd":    round(oracle, 4),
        "success":        success,
        "failure_mode":   None,
    }

    if success:
        return out   # no further classification needed

    # ── 1. seed_echo ────────────────────────────────────────────────────
    if _seed_echo(row):
        out["failure_mode"] = "seed_echo"
        out["cf_best"]   = _float(row, ["best_score", "best_cf"])
        out["cf_native"] = _float(row, ["cf_native"])
        return out

    # ── 2. harness_artifact ─────────────────────────────────────────────
    harness = _harness_flags(row, target_dir, rmsd) if target_dir else []
    if harness:
        out["failure_mode"] = "harness_artifact"
        out["harness_flags"] = harness
        return out

    # ── 3. selection_miss ───────────────────────────────────────────────
    if oracle >= 0.0 and oracle < lo_rmsd and rmsd >= lo_rmsd:
        out["failure_mode"]     = "selection_miss"
        out["rmsd_in_ensemble"] = round(oracle, 4)
        return out

    # ── 4–6. CF scoring failures (deep vs near) ─────────────────────────
    cf_best   = _float(row, ["best_score", "best_cf"])
    cf_native = _float(row, ["cf_native"])
    if cf_best is not None and cf_native is not None:
        if cf_best < cf_native - cf_delta:
            out["cf_best"]   = round(cf_best,   4)
            out["cf_native"] = round(cf_native, 4)
            out["cf_gap"]    = round(cf_best - cf_native, 4)
            if rmsd >= 6.0:
                out["failure_mode"] = "CF_scoring_failure_deep"
            elif rmsd >= lo_rmsd:
                out["failure_mode"] = "CF_scoring_failure_near"
            else:
                out["failure_mode"] = "CF_false_minimum"
            return out
    if rmsd >= lo_rmsd and 2.0 <= rmsd < 6.0 and oracle < lo_rmsd:
        out["failure_mode"] = "CF_scoring_failure_near"
        out["rmsd_in_ensemble"] = round(oracle, 4)
        return out

    # ── 4. timeout ──────────────────────────────────────────────────────
    wall = _float(row, ["wall_time_s", "wall_time"], default=0.0)
    status_str = (row.get("status") or row.get("success") or "")
    if wall >= TIMEOUT_WALL_S or rmsd >= 998.0 or rmsd < 0.0 or "timeout" in str(status_str).lower():
        out["failure_mode"] = "timeout"
        out["wall_time_s"]  = round(wall, 1)
        return out

    # ── 5. search_failure (sub-tagged by DoF) ───────────────────────────
    num_genes = _int(row, ["num_genes"])
    if num_genes is None and target_dir:
        num_genes = _parse_num_genes_from_dir(target_dir)

    if num_genes is not None:
        if num_genes <= RIGID_GENES_MAX:
            mode = "search_failure:rigid"
        elif num_genes >= FLEX_GENES_MIN:
            mode = "search_failure:flexible"
        else:
            mode = "search_failure:medium"
        out["num_genes"] = num_genes
    else:
        mode = "search_failure"
    out["failure_mode"] = mode
    return out


# ── Data loading ──────────────────────────────────────────────────────────────

def _load_per_target_csvs(result_dir):
    """Load per-target result.csv files (one row per target).
    Returns (rows_dict, missing_cols_set) where rows_dict maps pdb_id → row dict.
    """
    rows = {}
    for csv_path in sorted(glob.glob(os.path.join(result_dir, "*/result.csv"))):
        pdb_id = os.path.basename(os.path.dirname(csv_path))
        try:
            with open(csv_path, newline="") as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    rows[pdb_id] = row
                    break  # one data row per file
        except OSError:
            pass
    return rows


def _load_summary_csv(result_dir):
    """Fallback: load from a top-level summary CSV (astex_diverse_results.csv, etc.)."""
    candidates = sorted(
        glob.glob(os.path.join(result_dir, "*results*.csv"))
        + glob.glob(os.path.join(result_dir, "*summary*.csv"))
    )
    for csv_path in candidates:
        try:
            with open(csv_path, newline="") as fh:
                reader = csv.DictReader(fh)
                rows = {}
                for row in reader:
                    pid = (row.get("pdb_id") or "").strip()
                    if pid:
                        rows[pid] = row
                if rows:
                    return rows
        except OSError:
            pass
    return {}


def _detect_missing_columns(rows):
    """Report which spec-defined columns were absent across the dataset."""
    spec_cols = {
        "rmsd_hungarian", "best_cluster_rmsd", "seed_echo",
        "cf_native", "best_score", "best_cf",
        "wall_time_s", "wall_time", "num_genes",
    }
    seen = set()
    for row in rows.values():
        seen.update(row.keys())
    missing = spec_cols - seen
    return missing


# ── Main classify pipeline ────────────────────────────────────────────────────

def classify_run(result_dir, lo_rmsd=DEFAULT_LO_RMSD, cf_delta=DEFAULT_CF_DELTA):
    """Classify all targets in result_dir.
    Returns (results_dict, missing_cols_set).
    """
    # 1. Load rows: per-target CSVs preferred, summary CSV as fallback
    rows = _load_per_target_csvs(result_dir)
    if not rows:
        rows = _load_summary_csv(result_dir)
    if not rows:
        raise FileNotFoundError(
            f"No per-target result.csv files or summary CSV found in {result_dir!r}"
        )

    missing_cols = _detect_missing_columns(rows)

    # 2. Classify each target
    results = {}
    for pdb_id, row in sorted(rows.items()):
        target_dir = os.path.join(result_dir, pdb_id)
        if not os.path.isdir(target_dir):
            target_dir = None
        results[pdb_id] = classify_one(pdb_id, row, target_dir, lo_rmsd, cf_delta)

    return results, missing_cols


# ── Reporting ─────────────────────────────────────────────────────────────────

_MODE_ORDER = [
    "seed_echo",
    "harness_artifact",
    "selection_miss",
    "CF_scoring_failure_deep",
    "CF_scoring_failure_near",
    "CF_false_minimum",
    "timeout",
    "search_failure:rigid",
    "search_failure:medium",
    "search_failure:flexible",
    "search_failure",
]


def build_summary(results, lo_rmsd, hi_near_miss, missing_cols):
    """Return a list of strings forming the summary table."""
    lines = []
    n_total   = len(results)
    n_success = sum(1 for r in results.values() if r["success"])
    n_fail    = n_total - n_success
    near_miss = sum(
        1 for r in results.values()
        if not r["success"] and r["rmsd"] < hi_near_miss
    )

    lines.append("=" * 72)
    lines.append("FlexAIDdS Failure-Mode Classification")
    lines.append("=" * 72)
    lines.append(f"  Result dir : {results and '(see JSON)' or '—'}")
    lines.append(f"  Targets    : {n_total}  |  Success: {n_success}  "
                 f"({100*n_success/n_total:.1f}%)  |  Failed: {n_fail}")
    lines.append(f"  Near-misses: {near_miss}  (RMSD in [{lo_rmsd:.1f}, {hi_near_miss:.1f}))")
    if missing_cols:
        lines.append(f"  Missing CSV cols (fell back to defaults): {', '.join(sorted(missing_cols))}")
    lines.append("")

    # ── Per-mode counts ──────────────────────────────────────────────────
    mode_counts = {}
    for r in results.values():
        if not r["success"]:
            m = r["failure_mode"] or "unknown"
            mode_counts[m] = mode_counts.get(m, 0) + 1

    lines.append(f"  {'Failure mode':<28s}  {'Count':>5s}  {'% of failures':>13s}")
    lines.append("  " + "-" * 50)
    for mode in _MODE_ORDER:
        c = mode_counts.pop(mode, 0)
        if c:
            lines.append(f"  {mode:<28s}  {c:>5d}  {100*c/n_fail:>12.1f}%")
    for mode, c in sorted(mode_counts.items()):  # catch any unlisted modes
        lines.append(f"  {mode:<28s}  {c:>5d}  {100*c/n_fail:>12.1f}%")
    lines.append("")

    # ── Per-target table (failures only) ────────────────────────────────
    lines.append(
        f"  {'PDB':6s}  {'RMSD':>6s}  {'Oracle':>6s}  "
        f"{'Mode':<28s}  {'Extra'}"
    )
    lines.append("  " + "-" * 72)

    fail_rows = [
        (pid, r) for pid, r in sorted(results.items()) if not r["success"]
    ]
    fail_rows.sort(key=lambda x: x[1]["rmsd"])

    for pid, r in fail_rows:
        mode  = r.get("failure_mode") or "?"
        extra = ""
        if mode == "selection_miss":
            extra = f"ensemble_best={r.get('rmsd_in_ensemble', '?'):.4f}"
        elif mode in ("CF_false_minimum", "CF_scoring_failure_deep", "CF_scoring_failure_near"):
            extra = (f"cf_best={r.get('cf_best', '?'):.2f}  "
                     f"cf_native={r.get('cf_native', '?'):.2f}  "
                     f"gap={r.get('cf_gap', '?'):.2f}")
        elif mode == "harness_artifact":
            extra = ",".join(r.get("harness_flags", [])) or "?"
        elif mode == "seed_echo":
            extra = (f"cf_best={r.get('cf_best', '?')}  "
                     f"cf_native={r.get('cf_native', '?')}")
        elif mode == "timeout":
            extra = f"wall={r.get('wall_time_s', '?')}s"
        elif "search_failure" in mode:
            ng = r.get("num_genes")
            extra = f"num_genes={ng}" if ng is not None else "num_genes=?"
        lines.append(
            f"  {pid:<6s}  {r['rmsd']:>6.4f}  {r.get('oracle_rmsd', 999):>6.4f}  "
            f"{mode:<28s}  {extra}"
        )

    lines.append("=" * 72)
    return lines


# ── Entry point ───────────────────────────────────────────────────────────────

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Classify FlexAIDdS benchmark failures by causal mode."
    )
    parser.add_argument("result_dir", help="Root output directory of a benchmark run")
    parser.add_argument(
        "--lo-rmsd", type=float, default=DEFAULT_LO_RMSD,
        help=f"RMSD failure threshold in Å (default: {DEFAULT_LO_RMSD})"
    )
    parser.add_argument(
        "--hi-near-miss", type=float, default=DEFAULT_HI_NEAR_MISS,
        help=f"Near-miss band upper bound in Å (default: {DEFAULT_HI_NEAR_MISS})"
    )
    parser.add_argument(
        "--cf-delta", type=float, default=DEFAULT_CF_DELTA,
        help=f"CF gap threshold for false-minimum detection (default: {DEFAULT_CF_DELTA})"
    )
    args = parser.parse_args(argv)

    result_dir = os.path.expanduser(args.result_dir)
    if not os.path.isdir(result_dir):
        sys.exit(f"ERROR: result_dir not found: {result_dir!r}")

    results, missing_cols = classify_run(
        result_dir, lo_rmsd=args.lo_rmsd, cf_delta=args.cf_delta
    )

    # ── Write JSON ───────────────────────────────────────────────────────
    json_path = os.path.join(result_dir, "failure_modes.json")
    with open(json_path, "w") as fh:
        json.dump(results, fh, indent=2)
        fh.write("\n")
    print(f"  → {json_path}", file=sys.stderr)

    # ── Write + print summary ────────────────────────────────────────────
    summary_lines = build_summary(results, args.lo_rmsd, args.hi_near_miss, missing_cols)
    summary_lines_with_dir = [
        l.replace("(see JSON)", result_dir) if "(see JSON)" in l else l
        for l in summary_lines
    ]

    txt_path = os.path.join(result_dir, "failure_summary.txt")
    with open(txt_path, "w") as fh:
        fh.write("\n".join(summary_lines_with_dir) + "\n")
    print(f"  → {txt_path}", file=sys.stderr)

    print("\n".join(summary_lines_with_dir))


if __name__ == "__main__":
    main()
