#!/usr/bin/env python3
"""Enforce a priori / a posteriori benchmark self-eval contract.

Defers floors to PHASE4 / METHODOLOGY narrative; validates structure and
on-disk metrics without reimplementing docking.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

REQUIRED_APRIORI = [
    "one_variable",
    "panel_class",
    "codes",
    "matrix_pin",
    "no_sec",
    "sol9",
    "matched_control",
    "magnitude_floor",
    "report_tiers_separately",
]

VALID_STATUS = {
    "PASS",
    "FAIL",
    "PASS_LIVENESS",
    "VOID",
    "INVALID",
    "MISSING_OUT",
    "IN_PROGRESS",
    "NOT_RUN",
}

MATRIX_9DC9 = "9dc93717dfed0698006d88dd6a9627bc"

# S2 closed-gate pin pack (publication / audit trail)
# evidence/accept.txt + per-arm binary SHA256 (arm_pins.json or resolvable stamp).
SHA256_HEX_RE = re.compile(r"\b([0-9a-fA-F]{64})\b")


def load_apriori(path: Path) -> dict:
    return json.loads(path.read_text())


def _first_sha256_hex(text: str) -> str | None:
    m = SHA256_HEX_RE.search(text or "")
    return m.group(1).lower() if m else None


def discover_arm_names(out: Path) -> list[str]:
    """Arm directory basenames under OUT (arm_control → control)."""
    names: list[str] = []
    if not out.is_dir():
        return names
    for p in sorted(out.glob("arm_*")):
        if p.is_dir():
            names.append(p.name.removeprefix("arm_"))
    return names


def resolve_arm_binary_sha256(out: Path, arm: str) -> str | None:
    """Resolve per-arm binary SHA256 from pin files or stamped binary hash.

    Precedence:
    1. evidence/arm_pins.json → arms[arm].binary_sha256
    2. evidence/binary_<arm>.sha256 or evidence/<arm>.sha256
    3. arm_<arm>/binary.sha256 or arm_<arm>/bin/*.sha256 text
    4. sha256 of arm_<arm>/bin/FlexAIDdS.stamped (or OUT/bin if shared)
    5. OUT/binary.sha256 / evidence/binary_after_dup_fix.sha256 (shared stamp;
       only accepted when a single arm is being validated, or arm_pins maps it)
    """
    import hashlib

    pins = out / "evidence" / "arm_pins.json"
    if pins.is_file():
        try:
            data = json.loads(pins.read_text())
            arms = data.get("arms") or {}
            ent = arms.get(arm) or arms.get(f"arm_{arm}")
            if isinstance(ent, dict):
                h = ent.get("binary_sha256") or ent.get("sha256")
                if h and _first_sha256_hex(str(h)):
                    return _first_sha256_hex(str(h))
            elif isinstance(ent, str) and _first_sha256_hex(ent):
                return _first_sha256_hex(ent)
        except (json.JSONDecodeError, OSError):
            pass

    for cand in (
        out / "evidence" / f"binary_{arm}.sha256",
        out / "evidence" / f"{arm}.sha256",
        out / f"arm_{arm}" / "binary.sha256",
        out / f"arm_{arm}" / "bin" / "binary.sha256",
    ):
        if cand.is_file():
            h = _first_sha256_hex(cand.read_text(errors="replace"))
            if h:
                return h

    for bin_path in (
        out / f"arm_{arm}" / "bin" / "FlexAIDdS.stamped",
        out / f"arm_{arm}" / "FlexAIDdS.stamped",
        out / "bin" / "FlexAIDdS.stamped",
    ):
        if bin_path.is_file():
            h = hashlib.sha256(bin_path.read_bytes()).hexdigest()
            return h

    # Shared OUT-level sha files: use only as last resort when pins file
    # explicitly sets "shared_binary": true, or when discovering single arm.
    for cand in (
        out / "binary.sha256",
        out / "evidence" / "binary_after_dup_fix.sha256",
        out / "evidence" / "binary.sha256",
    ):
        if cand.is_file():
            h = _first_sha256_hex(cand.read_text(errors="replace"))
            if h:
                # Only valid if arm_pins declares shared, or only one arm.
                if pins.is_file():
                    try:
                        data = json.loads(pins.read_text())
                        if data.get("shared_binary") is True:
                            return h
                    except (json.JSONDecodeError, OSError):
                        pass
                if len(discover_arm_names(out)) <= 1:
                    return h
    return None


def validate_closed_gate_pins(
    out: Path,
    *,
    arms: list[str] | None = None,
    require_accept: bool = True,
    require_arm_sha: bool = True,
) -> list[str]:
    """Validate closed-gate publication pin pack under OUT.

    Required (S2):
    - evidence/accept.txt exists and is non-empty
    - each arm has a resolvable binary SHA256 (see resolve_arm_binary_sha256)

    Returns a list of human-readable error strings (empty ⇒ pass).
    """
    errs: list[str] = []
    if not out.is_dir():
        return [f"OUT is not a directory: {out}"]

    if require_accept:
        accept = out / "evidence" / "accept.txt"
        if not accept.is_file():
            errs.append("missing evidence/accept.txt")
        else:
            body = accept.read_text(errors="replace").strip()
            if not body:
                errs.append("evidence/accept.txt is empty")

    if require_arm_sha:
        arm_list = list(arms) if arms is not None else discover_arm_names(out)
        if not arm_list:
            errs.append("no arm_* directories found and --arms not provided")
        for arm in arm_list:
            sha = resolve_arm_binary_sha256(out, arm)
            if not sha:
                errs.append(
                    f"missing per-arm binary SHA256 for arm={arm!r} "
                    f"(need evidence/arm_pins.json arms.{arm}.binary_sha256 "
                    f"or arm_{arm}/bin/FlexAIDdS.stamped)"
                )
    return errs


def validate_pins_report(
    out: Path,
    *,
    arms: list[str] | None = None,
) -> dict:
    """Machine-readable pin validation result for CLI / tests."""
    arm_list = list(arms) if arms is not None else discover_arm_names(out)
    errs = validate_closed_gate_pins(out, arms=arm_list)
    resolved = {a: resolve_arm_binary_sha256(out, a) for a in arm_list}
    accept_path = out / "evidence" / "accept.txt"
    return {
        "ok": len(errs) == 0,
        "out": str(out),
        "errors": errs,
        "accept_present": accept_path.is_file()
        and bool(accept_path.read_text(errors="replace").strip())
        if accept_path.is_file()
        else False,
        "arms": arm_list,
        "arm_sha256": resolved,
    }


def validate_apriori(ap: dict) -> list[str]:
    errs = []
    for k in REQUIRED_APRIORI:
        if k not in ap:
            errs.append(f"missing a priori field: {k}")
    if ap.get("matrix_pin") not in (MATRIX_9DC9, "9dc9", "9dc93717dfed0698006d88dd6a9627bc"):
        errs.append(f"matrix_pin must be 9dc9, got {ap.get('matrix_pin')}")
    if ap.get("panel_class") in ("SEARCH_MISS", "NEAR_MISS", "GROSS_MISS"):
        if not ap.get("no_sec"):
            errs.append("Phase-4 SEARCH sampling docks require no_sec=true")
    codes = ap.get("codes") or []
    if isinstance(codes, str):
        codes = [c.strip() for c in codes.split(",") if c.strip()]
    if ap.get("panel_class") == "NEAR_MISS":
        bad = set(codes) - {"1N1M", "1L7F"}
        if bad:
            errs.append(f"NEAR_MISS codes must be subset of 1N1M,1L7F; extra {bad}")
    if ap.get("report_tiers_separately") is not True:
        errs.append("report_tiers_separately must be true")
    return errs


def result_metrics(out: Path, codes: list[str]) -> dict[str, dict]:
    m = {}
    for c in codes:
        p = out / c / "result.csv"
        if not p.is_file():
            m[c] = {"missing": True}
            continue
        row = next(csv.DictReader(p.open()))
        m[c] = {
            "bcr": float(row["best_cluster_rmsd"]),
            "s3": float(row["conditional_scanned_pool_ceiling"]),
            "elect": float(row["rmsd_to_crystal"]),
            "missing": False,
        }
    return m


def iter_engine_logs(out: Path):
    """Yield engine log paths under an arm OUT.

    Production FlexAIDdS writes [BOOM] and other L4 markers to **stderr.log**
    (often only there). Restarts land under ``r0/``, ``r1/``, …  Scanning
    only top-level stdout.log false-negatives real G4.1 OUTs.
    """
    if not out.is_dir():
        return
    # Prefer explicit names, then any *.log (covers r*/stderr.log etc.)
    seen: set[Path] = set()
    for name in ("stderr.log", "stdout.log", "driver.log"):
        for log in out.rglob(name):
            if log.is_file() and log not in seen:
                seen.add(log)
                yield log
    for log in out.rglob("*.log"):
        if log.is_file() and log not in seen:
            seen.add(log)
            yield log


def count_marker(out: Path, marker: str) -> int:
    n = 0
    for log in iter_engine_logs(out):
        n += log.read_text(errors="replace").count(marker)
    return n


def wipeout_signature(out: Path) -> bool:
    for log in iter_engine_logs(out):
        t = log.read_text(errors="replace")
        if re.search(r"wipeout|CF\s*[=~]\s*0\.0+", t, re.I):
            return True
    return False


def posteriori_g4_1_style(
    control: Path,
    treatments: dict[str, Path],
    codes: list[str],
    *,
    marker: str = "[BOOM]",
    magnitude_floor: float = -0.5,
) -> dict:
    ctrl = result_metrics(control, codes)
    if any(ctrl[c].get("missing") for c in codes):
        return {"status": "IN_PROGRESS", "reason": "control result.csv missing"}
    tx_out = {}
    best_name = None
    best_mean = 999.0
    for name, path in treatments.items():
        m = result_metrics(path, codes)
        if any(m[c].get("missing") for c in codes):
            tx_out[name] = {"status": "IN_PROGRESS", "metrics": m}
            continue
        deltas = [m[c]["bcr"] - ctrl[c]["bcr"] for c in codes]
        mean_d = sum(deltas) / len(deltas)
        n_sub2 = sum(1 for c in codes if m[c]["bcr"] <= 2.0)
        n_boom = count_marker(path, marker)
        wo = wipeout_signature(path)
        mag = mean_d <= magnitude_floor or n_sub2 >= 1
        l4_ok = n_boom > 0
        st = "PASS" if mag and l4_ok and not wo else (
            "PASS_LIVENESS" if l4_ok and not mag and not wo else "FAIL"
        )
        if not l4_ok:
            st = "FAIL"
        tx_out[name] = {
            "status": st,
            "mean_dBCR": mean_d,
            "n_bcr_sub2": n_sub2,
            "n_markers": n_boom,
            "wipeout": wo,
            "metrics": m,
            "per_target_dBCR": {c: m[c]["bcr"] - ctrl[c]["bcr"] for c in codes},
        }
        if mean_d < best_mean and st != "IN_PROGRESS":
            best_mean = mean_d
            best_name = name
    ctrl_boom = count_marker(control, marker)
    accept = False
    if best_name and tx_out[best_name]["status"] == "PASS" and ctrl_boom == 0:
        accept = True
    overall = "PASS" if accept else "FAIL"
    if any(tx_out[n].get("status") == "IN_PROGRESS" for n in tx_out):
        overall = "IN_PROGRESS"
    return {
        "status": overall,
        "accept_magnitude": accept,
        "best_treatment": best_name,
        "control_markers": ctrl_boom,
        "control": ctrl,
        "treatments": tx_out,
        "codes": codes,
        "magnitude_floor": magnitude_floor,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("preflight")
    p.add_argument("--apriori", type=Path, required=True)
    p.add_argument("--write-out", type=Path, help="copy apriori into OUT/APRIORI.json")

    q = sub.add_parser("posteriori")
    q.add_argument("--control", type=Path, required=True)
    q.add_argument("--treatment", action="append", nargs=2, metavar=("NAME", "PATH"), required=True)
    q.add_argument("--codes", default="1N1M,1L7F")
    q.add_argument("--marker", default="[BOOM]")
    q.add_argument("--json", action="store_true")

    v = sub.add_parser("validate-contract-doc")
    v.add_argument("--path", type=Path, default=Path("workorders/BENCHMARK_SELF_EVAL_CONTRACT.md"))

    pins = sub.add_parser(
        "validate-pins",
        help="S2: require evidence/accept.txt + per-arm binary SHA256 under OUT",
    )
    pins.add_argument("--out", type=Path, required=True, help="Campaign OUT root")
    pins.add_argument(
        "--arms",
        default="",
        help="Comma-separated arm names (default: discover arm_* dirs)",
    )
    pins.add_argument("--json", action="store_true")

    args = ap.parse_args(argv)

    if args.cmd == "preflight":
        ap_data = load_apriori(args.apriori)
        errs = validate_apriori(ap_data)
        if errs:
            print("APRIORI_FAIL")
            for e in errs:
                print(" ", e)
            return 2
        print("APRIORI_OK")
        if args.write_out:
            args.write_out.mkdir(parents=True, exist_ok=True)
            (args.write_out / "APRIORI.json").write_text(json.dumps(ap_data, indent=2) + "\n")
        return 0

    if args.cmd == "posteriori":
        codes = [c.strip() for c in args.codes.split(",") if c.strip()]
        treatments = {n: Path(p) for n, p in args.treatment}
        rec = posteriori_g4_1_style(
            args.control, treatments, codes, marker=args.marker
        )
        human = (
            f"status={rec['status']} accept={rec['accept_magnitude']} "
            f"best={rec['best_treatment']} ctrl_markers={rec['control_markers']}"
        )
        if args.json:
            print(human, file=sys.stderr)
            print(json.dumps(rec, indent=2))
        else:
            print(human)
        return 0 if rec["status"] == "PASS" else 1

    if args.cmd == "validate-pins":
        arm_arg = [a.strip() for a in args.arms.split(",") if a.strip()] or None
        rep = validate_pins_report(args.out, arms=arm_arg)
        if args.json:
            print(json.dumps(rep, indent=2))
        else:
            if rep["ok"]:
                print("PINS_OK", args.out)
                for a, h in (rep.get("arm_sha256") or {}).items():
                    print(f"  arm={a} sha256={h}")
            else:
                print("PINS_FAIL", args.out)
                for e in rep["errors"]:
                    print(" ", e)
        return 0 if rep["ok"] else 2

    if args.cmd == "validate-contract-doc":
        text = args.path.read_text()
        need = [
            "a priori",
            "a posteriori",
            "PHASE4_GATES_ACTUALIZED",
            "METHODOLOGY.md",
            "magnitude",
            "NEAR_MISS",
            "Sol #9",
            "9dc9",
            "accept.txt",
            "binary_sha256",
            "arm_pins",
        ]
        missing = []
        for n in need:
            if n not in text and n.lower() not in text.lower():
                missing.append(n)
        if missing:
            print("CONTRACT_DOC_FAIL", missing)
            return 2
        print("CONTRACT_DOC_OK", args.path)
        return 0

    return 2


if __name__ == "__main__":
    raise SystemExit(main())
