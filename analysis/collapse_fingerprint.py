#!/usr/bin/env python3
"""collapse_fingerprint.py — deterministic CMA-ES entropy-trace fingerprint.

CLI
---
  python3 analysis/collapse_fingerprint.py TRACE.csv [--out fingerprint.json] [--tol 1e-3]
  python3 analysis/collapse_fingerprint.py --compare A.json B.json [--tol 1e-3]

Pure Python 3.9+ stdlib only (csv, json, math, hashlib, argparse).
Exit 0 on successful parse / successful compare (INVARIANT or DIVERGENT).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

# Column aliases (case-insensitive flexible header match).
# First alias in each list is the canonical key used internally.
_COLUMN_ALIASES: Dict[str, Tuple[str, ...]] = {
    "gen": ("gen", "generation", "g", "iter", "iteration", "step"),
    "H_search": ("h_search", "hsearch", "search_entropy", "hs", "shannon_search"),
    "H_energy": ("h_energy", "henergy", "energy_entropy", "he", "shannon_energy"),
    "F": ("f", "free_energy", "helmholtz", "f_free"),
    "best_cf": ("best_cf", "bestcf", "cf_best", "best_score", "cf", "score"),
    "n_evals": ("n_evals", "nevals", "evals", "n_eval", "evaluations"),
}

_REQUIRED = ("gen", "H_search", "H_energy", "F", "best_cf")
_OPTIONAL = ("n_evals",)

# Numeric fingerprint keys compared under --tol (excluding sha256 / string meta).
_NUMERIC_KEYS = (
    "n_rows",
    "n_gens",
    "H_search_start",
    "H_search_end",
    "H_search_min",
    "H_search_max",
    "H_search_delta",
    "H_energy_start",
    "H_energy_end",
    "H_energy_mean",
    "H_energy_std",
    "F_end",
    "best_cf_end",
    "best_cf_min",
    "collapse_ratio",
)

_ROUND_DECIMALS = 6
_EPS = 1e-12


def _normalize_header(name: str) -> str:
    return name.strip().lower().replace("-", "_").replace(" ", "_")


def _resolve_columns(fieldnames: Sequence[str]) -> Dict[str, str]:
    """Map canonical keys → actual CSV header names (flexible match)."""
    norm_to_actual = {_normalize_header(h): h for h in fieldnames if h is not None}
    resolved: Dict[str, str] = {}
    for canonical, aliases in _COLUMN_ALIASES.items():
        for alias in aliases:
            if alias in norm_to_actual:
                resolved[canonical] = norm_to_actual[alias]
                break
        # Also accept exact canonical (normalized) if present under different casing.
        if canonical not in resolved:
            cn = _normalize_header(canonical)
            if cn in norm_to_actual:
                resolved[canonical] = norm_to_actual[cn]
    missing = [k for k in _REQUIRED if k not in resolved]
    if missing:
        raise ValueError(
            "CSV missing required columns (flexible match failed): "
            + ", ".join(missing)
            + f"; found headers={list(fieldnames)}"
        )
    return resolved


def _to_float(value: Any, field: str, row_idx: int) -> float:
    if value is None or (isinstance(value, str) and value.strip() == ""):
        raise ValueError(f"row {row_idx}: empty value for {field}")
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"row {row_idx}: cannot parse {field}={value!r}") from exc


def load_trace(path: Path) -> List[Dict[str, float]]:
    """Load entropy trace CSV into list of row dicts with canonical keys."""
    with path.open(newline="", encoding="utf-8") as fh:
        # Skip blank / comment lines that might precede the header.
        sample_lines: List[str] = []
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            sample_lines.append(line)
        if not sample_lines:
            raise ValueError(f"empty or comment-only CSV: {path}")
        reader = csv.DictReader(sample_lines)
        if not reader.fieldnames:
            raise ValueError(f"no header row in {path}")
        colmap = _resolve_columns(list(reader.fieldnames))
        rows: List[Dict[str, float]] = []
        for i, raw in enumerate(reader, start=1):
            if raw is None:
                continue
            # Skip fully empty rows.
            if all((v is None or str(v).strip() == "") for v in raw.values()):
                continue
            row: Dict[str, float] = {}
            for key in _REQUIRED:
                row[key] = _to_float(raw.get(colmap[key]), key, i)
            if "n_evals" in colmap:
                val = raw.get(colmap["n_evals"])
                if val is not None and str(val).strip() != "":
                    row["n_evals"] = _to_float(val, "n_evals", i)
            rows.append(row)
    if not rows:
        raise ValueError(f"no data rows in {path}")
    return rows


def _mean(xs: Sequence[float]) -> float:
    return sum(xs) / len(xs) if xs else float("nan")


def _std_sample(xs: Sequence[float]) -> float:
    """Population std (ddof=0) for determinism / invariance across runs."""
    n = len(xs)
    if n == 0:
        return float("nan")
    if n == 1:
        return 0.0
    m = _mean(xs)
    return math.sqrt(sum((x - m) ** 2 for x in xs) / n)


def _round6(x: float) -> float:
    return round(float(x), _ROUND_DECIMALS)


def series_sha256(rows: Sequence[Mapping[str, float]]) -> str:
    """SHA256 of rounded (gen, H_search, H_energy, F, best_cf[, n_evals]) series.

    Values rounded to 6 decimals; fields joined with commas, rows with newlines.
    Deterministic UTF-8 digest, hex-encoded.
    """
    lines: List[str] = []
    for r in rows:
        parts = [
            f"{_round6(r['gen']):.6f}",
            f"{_round6(r['H_search']):.6f}",
            f"{_round6(r['H_energy']):.6f}",
            f"{_round6(r['F']):.6f}",
            f"{_round6(r['best_cf']):.6f}",
        ]
        if "n_evals" in r:
            parts.append(f"{_round6(r['n_evals']):.6f}")
        lines.append(",".join(parts))
    payload = "\n".join(lines).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def series_vectors(
    rows: Sequence[Mapping[str, float]],
) -> Dict[str, List[float]]:
    """Rounded series vectors used for optional L2 sha-fallback compare."""
    keys = ["gen", "H_search", "H_energy", "F", "best_cf"]
    out: Dict[str, List[float]] = {k: [_round6(r[k]) for r in rows] for k in keys}
    if any("n_evals" in r for r in rows):
        out["n_evals"] = [_round6(r.get("n_evals", 0.0)) for r in rows]
    return out


def compute_fingerprint(rows: Sequence[Mapping[str, float]]) -> Dict[str, Any]:
    """Build deterministic JSON-serializable fingerprint from ordered rows."""
    h_search = [float(r["H_search"]) for r in rows]
    h_energy = [float(r["H_energy"]) for r in rows]
    f_vals = [float(r["F"]) for r in rows]
    best_cf = [float(r["best_cf"]) for r in rows]
    gens = [float(r["gen"]) for r in rows]

    h0 = h_search[0]
    hN = h_search[-1]
    collapse = hN / max(h0, _EPS)

    # Unique generation count (supports non-integer gens).
    n_gens = len({_round6(g) for g in gens})

    fp: Dict[str, Any] = {
        "n_rows": len(rows),
        "n_gens": n_gens,
        "H_search_start": _round6(h0),
        "H_search_end": _round6(hN),
        "H_search_min": _round6(min(h_search)),
        "H_search_max": _round6(max(h_search)),
        "H_search_delta": _round6(hN - h0),
        "H_energy_start": _round6(h_energy[0]),
        "H_energy_end": _round6(h_energy[-1]),
        "H_energy_mean": _round6(_mean(h_energy)),
        "H_energy_std": _round6(_std_sample(h_energy)),
        "F_end": _round6(f_vals[-1]),
        "best_cf_end": _round6(best_cf[-1]),
        "best_cf_min": _round6(min(best_cf)),
        "collapse_ratio": _round6(collapse),
        "sha256": series_sha256(rows),
        # Embedded rounded series for offline L2 fallback when comparing dumps
        # that predate sha (optional; compare tolerates absence).
        "series_rounded": series_vectors(rows),
    }
    return fp


def _finite_close(a: float, b: float, tol: float) -> bool:
    if math.isnan(a) and math.isnan(b):
        return True
    if math.isinf(a) or math.isinf(b):
        return a == b
    return abs(a - b) <= tol


def _series_l2(a: Mapping[str, List[float]], b: Mapping[str, List[float]]) -> float:
    """RMS L2 over shared keys with equal length; inf if shapes mismatch."""
    keys = sorted(set(a) & set(b))
    if not keys:
        return float("inf")
    total = 0.0
    count = 0
    for k in keys:
        va, vb = a[k], b[k]
        if len(va) != len(vb):
            return float("inf")
        for x, y in zip(va, vb):
            total += (x - y) ** 2
            count += 1
    if count == 0:
        return float("inf")
    return math.sqrt(total / count)


def compare_fingerprints(
    fa: Mapping[str, Any],
    fb: Mapping[str, Any],
    tol: float,
) -> Tuple[str, List[str]]:
    """Return (status, detail_lines). status is INVARIANT or DIVERGENT."""
    details: List[str] = []
    numeric_ok = True
    for key in _NUMERIC_KEYS:
        if key not in fa or key not in fb:
            numeric_ok = False
            details.append(f"missing key: {key} in A={key in fa} B={key in fb}")
            continue
        va, vb = float(fa[key]), float(fb[key])
        if not _finite_close(va, vb, tol):
            numeric_ok = False
            details.append(f"{key}: {va} vs {vb} (Δ={abs(va - vb):.6g} > tol={tol})")

    sha_a = str(fa.get("sha256", ""))
    sha_b = str(fb.get("sha256", ""))
    sha_match = bool(sha_a) and sha_a == sha_b

    series_ok = False
    if not sha_match:
        sa = fa.get("series_rounded")
        sb = fb.get("series_rounded")
        if isinstance(sa, dict) and isinstance(sb, dict):
            l2 = _series_l2(sa, sb)  # type: ignore[arg-type]
            series_ok = l2 <= tol
            details.append(f"sha mismatch; series L2={l2:.6g} (tol={tol})")
        else:
            details.append("sha mismatch and no series_rounded for L2 fallback")
    else:
        details.append(f"sha256 match: {sha_a}")

    if numeric_ok and (sha_match or series_ok):
        return "INVARIANT", details
    return "DIVERGENT", details


def load_fingerprint_json(path: Path) -> Dict[str, Any]:
    with path.open(encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        raise ValueError(f"fingerprint JSON must be an object: {path}")
    return data


def dump_fingerprint(fp: Mapping[str, Any], path: Optional[Path]) -> str:
    text = json.dumps(fp, indent=2, sort_keys=True) + "\n"
    if path is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
    return text


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="collapse_fingerprint.py",
        description="Compute or compare CMA-ES entropy-trace collapse fingerprints.",
    )
    p.add_argument(
        "trace",
        nargs="?",
        default=None,
        help="Entropy trace CSV (gen, H_search, H_energy, F, best_cf [, n_evals])",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Write fingerprint JSON to this path (trace mode)",
    )
    p.add_argument(
        "--compare",
        nargs=2,
        metavar=("A.json", "B.json"),
        default=None,
        help="Compare two fingerprint JSON files",
    )
    p.add_argument(
        "--tol",
        type=float,
        default=1e-3,
        help="Absolute tolerance for numeric field compare (default: 1e-3)",
    )
    p.add_argument(
        "--no-series",
        action="store_true",
        help="Omit series_rounded from fingerprint JSON (smaller dumps)",
    )
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _build_parser().parse_args(argv)

    if args.compare is not None:
        path_a = Path(args.compare[0])
        path_b = Path(args.compare[1])
        fa = load_fingerprint_json(path_a)
        fb = load_fingerprint_json(path_b)
        status, details = compare_fingerprints(fa, fb, tol=args.tol)
        print(status)
        for line in details:
            print(f"  {line}")
        # Exit 0 always for successful parse (INVARIANT or DIVERGENT).
        return 0

    if args.trace is None:
        _build_parser().error("TRACE.csv is required unless --compare is used")

    rows = load_trace(Path(args.trace))
    fp = compute_fingerprint(rows)
    if args.no_series:
        fp = {k: v for k, v in fp.items() if k != "series_rounded"}

    text = dump_fingerprint(fp, args.out)
    # Always print fingerprint to stdout for piping / harness capture.
    sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, json.JSONDecodeError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
