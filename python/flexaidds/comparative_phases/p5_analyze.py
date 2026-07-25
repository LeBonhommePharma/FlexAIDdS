# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
"""Phase P5: S_top10 bootstrap scaffolding + comparative table + thin sync.

Implements docs/implementation/COMPARATIVE_GOAL_METHODOLOGY.md Phase 5:

  - For each arm with result.csv under
    ``$FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine/<ARM>/<campaign>/``
    invoke the **shipped** ``scripts/bootstrap_3dsig_s_top10.py`` via subprocess
    (never reimplement the bootstrap median here).
  - Emit ``COMPARATIVE_TABLE.md`` under
    ``campaigns/three_engine/analysis/<campaign>/``.
  - Optionally invoke ``scripts/sync_three_engine_local_to_icloud.sh``
    (supports ``--dry-run``).

Dry-run always validates the bootstrap script is invokable (``--help``)
and still writes a table skeleton (N/A metrics when no arm data).
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from .p1_binaries import DEFAULT_LOCAL_ROOT, load_arm_pins, resolve_repo_root

# Claim contract RMSD threshold (same as bootstrap_3dsig_s_top10.DEFAULT_THRESH)
RMSD_THRESH = 2.0
DEFAULT_BOOTSTRAPS = 10_000
DEFAULT_ARMS: Sequence[str] = ("A", "B", "C")

TABLE_COLUMNS = (
    "arm",
    "binary_sha",
    "commit",
    "temper",
    "clusta",
    "N",
    "S_top10_median",
    "S1",
    "BCR",
    "matrix_md5",
    "reconstruction",
)


def local_root(override: Optional[str] = None) -> Path:
    if override:
        return Path(override).expanduser().resolve()
    return Path(
        os.environ.get("FLEXAIDDS_LOCAL_ROOT", str(DEFAULT_LOCAL_ROOT))
    ).expanduser().resolve()


def bootstrap_script(root: Path) -> Path:
    return root / "scripts" / "bootstrap_3dsig_s_top10.py"


def sync_script(root: Path) -> Path:
    return root / "scripts" / "sync_three_engine_local_to_icloud.sh"


def arm_campaign_dir(local: Path, arm: str, campaign: str) -> Path:
    return local / "campaigns" / "three_engine" / arm / campaign


def analysis_dir_for(local: Path, campaign: str) -> Path:
    return local / "campaigns" / "three_engine" / "analysis" / campaign


def find_result_csvs(arm_dir: Path) -> List[Path]:
    """Non-recursive-on-CloudDocs: only walk local arm_dir for result.csv."""
    if not arm_dir.is_dir():
        return []
    return sorted(p for p in arm_dir.rglob("result.csv") if p.is_file())


def validate_bootstrap_invokable(root: Path) -> Dict[str, Any]:
    """Prove shipped bootstrap is invokable via --help and import s_top10."""
    script = bootstrap_script(root)
    out: Dict[str, Any] = {
        "script": str(script),
        "exists": script.is_file(),
        "help_ok": False,
        "s_top10_import_ok": False,
        "error": None,
    }
    if not script.is_file():
        out["error"] = f"missing {script}"
        return out

    help_proc = subprocess.run(
        [sys.executable, str(script), "--help"],
        capture_output=True,
        text=True,
    )
    out["help_ok"] = help_proc.returncode == 0
    out["help_returncode"] = help_proc.returncode
    if help_proc.returncode != 0:
        out["error"] = (help_proc.stderr or help_proc.stdout or "")[:500]
        return out

    # Import the real module and touch s_top10 (no reimplementation).
    import importlib.util

    try:
        spec = importlib.util.spec_from_file_location(
            "bootstrap_3dsig_s_top10_p5", script
        )
        if spec is None or spec.loader is None:
            out["error"] = "importlib failed to load bootstrap module"
            return out
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        # Smoke: inclusive ≤2.0 claim contract
        ok = (
            callable(getattr(mod, "s_top10", None))
            and mod.s_top10([2.0]) is True
            and mod.s_top10([2.0001]) is False
            and getattr(mod, "DEFAULT_BOOTSTRAPS", 0) == DEFAULT_BOOTSTRAPS
        )
        out["s_top10_import_ok"] = bool(ok)
        if not ok:
            out["error"] = "s_top10 contract smoke failed"
    except Exception as exc:  # noqa: BLE001 — surface import errors to CLI
        out["error"] = f"import bootstrap failed: {exc}"
    return out


def run_bootstrap(
    root: Path,
    arm_dir: Path,
    json_out: Path,
    *,
    bootstraps: int = DEFAULT_BOOTSTRAPS,
) -> Dict[str, Any]:
    """Invoke shipped bootstrap_3dsig_s_top10.py on arm-dir (subprocess only)."""
    script = bootstrap_script(root)
    if not script.is_file():
        raise FileNotFoundError(f"missing {script}")
    json_out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(script),
        "--arm-dir",
        str(arm_dir),
        "--bootstraps",
        str(bootstraps),
        "--json-out",
        str(json_out),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    result: Dict[str, Any] = {
        "cmd": cmd,
        "returncode": proc.returncode,
        "stdout": (proc.stdout or "")[-2000:],
        "stderr": (proc.stderr or "")[-2000:],
        "json_out": str(json_out),
    }
    if proc.returncode == 0 and json_out.is_file():
        try:
            result["bootstrap"] = json.loads(json_out.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            result["bootstrap"] = None
    return result


def _finite(val: object) -> Optional[float]:
    if val is None or val == "" or val == "NA":
        return None
    try:
        v = float(val)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    if v != v or v < 0.0:  # NaN or negative
        return None
    return v


def scan_s1_bcr(arm_dir: Path, *, thresh: float = RMSD_THRESH) -> Dict[str, Any]:
    """Scan result.csv for S1 (rmsd_top1) and BCR (rmsd_bcr) rates.

    Does **not** compute S_top10 — that is exclusively the bootstrap script's job.
    """
    csvs = find_result_csvs(arm_dir)
    n = 0
    s1_ok = 0
    bcr_ok = 0
    for path in csvs:
        try:
            with path.open(newline="", encoding="utf-8") as fh:
                rows = list(csv.DictReader(fh))
        except OSError:
            continue
        if not rows:
            continue
        row = rows[0]
        keys = {k.lower(): k for k in row}
        n += 1
        # Prefer explicit success flags when present
        if "success_s1" in keys:
            try:
                s1_ok += int(float(row[keys["success_s1"]]))
            except (TypeError, ValueError):
                pass
        else:
            top1_key = keys.get("rmsd_top1")
            if top1_key is not None:
                r = _finite(row[top1_key])
                if r is not None and r <= thresh:
                    s1_ok += 1
        if "success_s3" in keys:
            try:
                bcr_ok += int(float(row[keys["success_s3"]]))
            except (TypeError, ValueError):
                pass
        else:
            bcr_key = keys.get("rmsd_bcr") or keys.get("rmsd_s3")
            if bcr_key is not None:
                r = _finite(row[bcr_key])
                if r is not None and r <= thresh:
                    bcr_ok += 1
    return {
        "n_csv": n,
        "s1_success": s1_ok,
        "bcr_success": bcr_ok,
        "S1": f"{s1_ok}/{n}" if n else "N/A",
        "BCR": f"{bcr_ok}/{n}" if n else "N/A",
    }


def _fmt_median_ci(bootstrap: Dict[str, Any]) -> str:
    med = bootstrap.get("median")
    p05 = bootstrap.get("p05")
    p95 = bootstrap.get("p95")
    if med is None:
        return "N/A"
    try:
        med_f = float(med)
    except (TypeError, ValueError):
        return str(med)
    if p05 is not None and p95 is not None:
        try:
            return f"{med_f:.4f} [{float(p05):.4f},{float(p95):.4f}]"
        except (TypeError, ValueError):
            pass
    return f"{med_f:.4f}"


def write_comparative_table(
    out_path: Path,
    rows: List[Dict[str, Any]],
    *,
    campaign: str,
    matrix_md5: str,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"# Comparative table — {campaign}",
        "",
        "Source: Phase P5 (`scripts/comparative_p5_analyze.py`).",
        "Headline metric: **S_top10** bootstrap median via "
        "`scripts/bootstrap_3dsig_s_top10.py` (not S1/BCR).",
        "",
        f"matrix_md5: `{matrix_md5}`",
        "",
        "| arm | binary_sha | commit | temper | clusta | N | S_top10_median | S1 | BCR | matrix_md5 | reconstruction |",
        "|-----|------------|--------|--------|--------|---|----------------|----|-----|------------|----------------|",
    ]
    for r in rows:
        lines.append(
            "| {arm} | {binary_sha} | {commit} | {temper} | {clusta} | {N} | "
            "{S_top10_median} | {S1} | {BCR} | {matrix_md5} | {reconstruction} |".format(
                arm=r.get("arm", ""),
                binary_sha=(str(r.get("binary_sha") or "") or "")[:12] or "—",
                commit=(str(r.get("commit") or "") or "")[:8] or "—",
                temper=r.get("temper", "—"),
                clusta=r.get("clusta", "—"),
                N=r.get("N", ""),
                S_top10_median=r.get("S_top10_median", "N/A"),
                S1=r.get("S1", "N/A"),
                BCR=r.get("BCR", "N/A"),
                matrix_md5=(str(r.get("matrix_md5") or matrix_md5) or "")[:12] or "—",
                reconstruction=r.get("reconstruction", "") or "—",
            )
        )
    lines.append("")
    lines.append(
        f"Columns: {' | '.join(TABLE_COLUMNS)}"
    )
    lines.append("")
    out_path.write_text("\n".join(lines), encoding="utf-8")


def _load_receipt(local: Path, arm: str) -> Dict[str, Any]:
    path = local / "campaigns" / "three_engine" / "receipts" / f"arm_{arm}_binary.json"
    if not path.is_file():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except (OSError, json.JSONDecodeError):
        return {}


def _reconstruction_cell(receipt: Dict[str, Any]) -> str:
    if receipt.get("reconstruction") is True:
        return receipt.get("reconstruction_label") or "true"
    if receipt.get("reconstruction") is False:
        return "false"
    status = receipt.get("status") or ""
    if "MISSING" in str(status) or "RECONSTRUCTION" in str(status) or "BUILD_FROM" in str(status):
        return receipt.get("reconstruction_label") or str(status)
    return str(receipt.get("reconstruction") or "")


def maybe_run_sync(
    root: Path,
    campaign: str,
    *,
    dry_run: bool,
) -> Dict[str, Any]:
    """Call sync_three_engine_local_to_icloud.sh; always pass --dry-run when dry_run."""
    script = sync_script(root)
    info: Dict[str, Any] = {
        "script": str(script),
        "exists": script.is_file(),
        "invoked": False,
        "returncode": None,
        "note": "",
    }
    if not script.is_file():
        info["note"] = "sync script missing — skipped"
        return info

    cmd = ["bash", str(script), "--campaign", campaign]
    if dry_run:
        cmd.append("--dry-run")
    # Bound hang risk on CloudDocs (sync uses rsync; dry-run only logs)
    env = os.environ.copy()
    env["FLEXAIDDS_ROOT"] = str(root)
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            env=env,
            timeout=180 if dry_run else 300,
        )
        info["invoked"] = True
        info["returncode"] = proc.returncode
        info["stdout"] = (proc.stdout or "")[-1500:]
        info["stderr"] = (proc.stderr or "")[-1500:]
        info["note"] = (
            f"sync --dry-run exit={proc.returncode}"
            if dry_run
            else f"sync exit={proc.returncode}"
        )
    except subprocess.TimeoutExpired:
        info["invoked"] = True
        info["note"] = "sync timed out (CloudDocs hang risk) — skipped remainder"
    except OSError as exc:
        info["note"] = f"sync invoke failed: {exc}"
    return info


def run_p5(
    campaign: str,
    *,
    local_root_path: Optional[str] = None,
    arms: Optional[Sequence[str]] = None,
    dry_run: bool = False,
    bootstraps: int = DEFAULT_BOOTSTRAPS,
    run_sync: bool = True,
) -> Dict[str, Any]:
    """Execute Phase P5 analysis scaffolding.

    Returns a status dict with ``phase``, ``status``, ``table_path``, ``rows``, etc.
    """
    root = resolve_repo_root()
    local = local_root(local_root_path)
    want = list(arms) if arms else list(DEFAULT_ARMS)

    try:
        pins = load_arm_pins(root)
    except (FileNotFoundError, ValueError, OSError, json.JSONDecodeError) as exc:
        return {
            "phase": "P5",
            "status": "fail",
            "reason": f"arm_pins load failed: {exc}",
            "campaign": campaign,
        }

    matrix_md5 = (pins.get("matrix") or {}).get("md5", "")
    pin_arms = pins.get("arms") or {}

    # Always prove bootstrap is invokable (dry-run and live).
    boot_check = validate_bootstrap_invokable(root)
    if not boot_check.get("help_ok") or not boot_check.get("s_top10_import_ok"):
        return {
            "phase": "P5",
            "status": "fail",
            "reason": (
                f"bootstrap script not invokable: {boot_check.get('error')}"
            ),
            "campaign": campaign,
            "bootstrap_check": boot_check,
        }

    analysis_dir = analysis_dir_for(local, campaign)
    analysis_dir.mkdir(parents=True, exist_ok=True)

    rows: List[Dict[str, Any]] = []
    boot_results: Dict[str, Any] = {}

    for arm in want:
        arm_spec = pin_arms.get(arm) or {}
        receipt = _load_receipt(local, arm)
        arm_dir = arm_campaign_dir(local, arm, campaign)
        csvs = find_result_csvs(arm_dir) if arm_dir.is_dir() else []

        row: Dict[str, Any] = {
            "arm": arm,
            "binary_sha": receipt.get("binary_sha256"),
            "commit": (
                receipt.get("git_commit")
                or arm_spec.get("source_commit")
                or arm_spec.get("source_commit_at_doc")
            ),
            "temper": arm_spec.get("temper"),
            "clusta": arm_spec.get("clusta"),
            "N": 0,
            "S_top10_median": "N/A",
            "S1": "N/A",
            "BCR": "N/A",
            "matrix_md5": matrix_md5,
            "reconstruction": _reconstruction_cell(receipt),
            "arm_dir": str(arm_dir),
            "n_result_csv": len(csvs),
        }

        if dry_run:
            # Skeleton only; do not run 10k bootstrap (even if data present).
            row["N"] = "dry-run"
            if csvs:
                row["S1"] = "dry-run"
                row["BCR"] = "dry-run"
                row["S_top10_median"] = "dry-run (data present; not bootstrapped)"
            rows.append(row)
            continue

        if not arm_dir.is_dir() or not csvs:
            row["N"] = 0
            row["S_top10_median"] = "N/A (no result.csv)"
            rows.append(row)
            continue

        # Auxiliary S1 / BCR rates from result.csv (not S_top10).
        aux = scan_s1_bcr(arm_dir)
        row["S1"] = aux["S1"]
        row["BCR"] = aux["BCR"]
        row["N"] = aux["n_csv"]

        json_out = analysis_dir / f"{arm}_s_top10.json"
        try:
            br = run_bootstrap(root, arm_dir, json_out, bootstraps=bootstraps)
            boot_results[arm] = {
                "returncode": br["returncode"],
                "json_out": br["json_out"],
                "stdout_tail": br.get("stdout"),
                "stderr_tail": br.get("stderr"),
            }
            if br.get("bootstrap"):
                b = br["bootstrap"]
                row["S_top10_median"] = _fmt_median_ci(b)
                if b.get("n_cases") is not None:
                    row["N"] = b["n_cases"]
                boot_results[arm]["bootstrap_keys"] = sorted(b.keys())
            elif br["returncode"] != 0:
                row["S_top10_median"] = f"error:{br['returncode']}"
                boot_results[arm]["error"] = br.get("stderr") or br.get("stdout")
        except FileNotFoundError as exc:
            row["S_top10_median"] = f"error:{exc}"

        rows.append(row)

    table_path = analysis_dir / "COMPARATIVE_TABLE.md"
    write_comparative_table(
        table_path, rows, campaign=campaign, matrix_md5=matrix_md5
    )

    sync_info: Dict[str, Any] = {"note": "sync disabled"}
    if run_sync:
        sync_info = maybe_run_sync(root, campaign, dry_run=dry_run)
    else:
        sync_info = {
            "note": (
                "skipped (--no-sync); manual: "
                f"bash scripts/sync_three_engine_local_to_icloud.sh --campaign {campaign}"
            )
        }

    status = "pass"
    reason = f"table written dry_run={dry_run} path={table_path}"

    return {
        "phase": "P5",
        "status": status,
        "reason": reason,
        "campaign": campaign,
        "local_root": str(local),
        "table_path": str(table_path),
        "rows": rows,
        "bootstrap_check": boot_check,
        "bootstrap_help_ok": True,
        "bootstrap_results": boot_results,
        "sync": sync_info,
        "dry_run": dry_run,
        "arms": want,
        "matrix_md5": matrix_md5,
    }


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Phase P5: bootstrap S_top10 per arm (subprocess to "
            "bootstrap_3dsig_s_top10.py), write COMPARATIVE_TABLE.md, "
            "optional iCloud thin sync."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python3 scripts/comparative_p5_analyze.py "
            "--campaign comparative_pilot8 --dry-run\n"
            "  python3 scripts/comparative_p5_analyze.py "
            "--campaign comparative_full85 --arms A,B,C\n"
        ),
    )
    ap.add_argument(
        "--campaign",
        required=True,
        help="Campaign id under campaigns/three_engine/<ARM>/<campaign>/",
    )
    ap.add_argument(
        "--local-root",
        default=None,
        help="Override FLEXAIDDS_LOCAL_ROOT (default ~/flexaidds_results)",
    )
    ap.add_argument(
        "--arms",
        default="A,B,C",
        help="Comma-separated arms (default A,B,C; B0 allowed)",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="Write table skeleton; validate bootstrap --help; sync --dry-run",
    )
    ap.add_argument(
        "--bootstraps",
        type=int,
        default=DEFAULT_BOOTSTRAPS,
        help=f"Bootstrap resamples passed to bootstrap script (default {DEFAULT_BOOTSTRAPS})",
    )
    ap.add_argument(
        "--no-sync",
        action="store_true",
        help="Do not invoke sync_three_engine_local_to_icloud.sh",
    )
    ap.add_argument(
        "--json-out",
        type=Path,
        default=None,
        help="Optional path for full P5 status JSON",
    )
    return ap


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(list(argv) if argv is not None else None)
    arms = [a.strip() for a in str(args.arms).split(",") if a.strip()]
    if not arms:
        print("error: --arms empty", file=sys.stderr)
        return 2

    result = run_p5(
        args.campaign,
        local_root_path=args.local_root,
        arms=arms,
        dry_run=bool(args.dry_run),
        bootstraps=int(args.bootstraps),
        run_sync=not bool(args.no_sync),
    )

    status = result.get("status", "fail")
    print(f"PHASE=P5 status={status}")
    print(f"campaign={result.get('campaign')}")
    if result.get("table_path"):
        print(f"table={result['table_path']}")
    print(f"reason={result.get('reason', '')}")
    sync = result.get("sync") or {}
    if isinstance(sync, dict):
        print(f"sync={sync.get('note') or sync}")
    else:
        print(f"sync={sync}")
    bc = result.get("bootstrap_check") or {}
    if bc:
        print(
            f"bootstrap_help_ok={bc.get('help_ok')} "
            f"s_top10_import_ok={bc.get('s_top10_import_ok')}"
        )
    for row in result.get("rows") or []:
        print(
            f"  arm={row.get('arm')} N={row.get('N')} "
            f"S_top10={row.get('S_top10_median')} "
            f"S1={row.get('S1')} BCR={row.get('BCR')}"
        )

    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(
            json.dumps(result, indent=2, default=str) + "\n", encoding="utf-8"
        )
        print(f"json_out={args.json_out}")

    return 0 if status == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
