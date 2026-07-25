"""Phase 1 — pin / receipt source-pinned binaries for comparative arms A, B, C.

See ``docs/implementation/COMPARATIVE_GOAL_METHODOLOGY.md`` § Phase 1 and
``docs/implementation/arm_pins.json``.

Does **not** compile FlexAID. Only inspects existing Mach-Os, SHA256s them,
writes receipts under ``$FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine/receipts/``.
Missing A/B require ``--allow-reconstruction`` for a labeled stub receipt
(no invented binary digests).

Copyright 2026 Le Bonhomme Pharma
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

# Science arms audited in Phase 1 (B0 is binary-control of B, not pinned here).
SCIENCE_ARMS: Tuple[str, ...] = ("A", "B", "C")

DEFAULT_LOCAL_ROOT = Path.home() / "flexaidds_results"
PINS_REL = Path("docs/implementation/arm_pins.json")
RECEIPTS_REL = Path("campaigns/three_engine/receipts")
BIN_REL = Path("three_engine_entropy_q1/bin")

# Candidate executable basenames per arm (order = preference).
_ARM_BINARY_NAMES: Mapping[str, Sequence[str]] = {
    "A": ("FlexAID", "FlexAIDdS"),
    "B": ("FlexAID", "FlexAIDdS"),
    "C": ("FlexAIDdS", "FlexAID"),
}


def resolve_repo_root(start: Optional[Path] = None) -> Path:
    """Resolve FlexAIDdS repo root via git or package layout."""
    if start is None:
        # p1_binaries.py → comparative_phases → flexaidds → python → repo
        start = Path(__file__).resolve()
        for parent in start.parents:
            if (parent / "docs" / "implementation" / "arm_pins.json").is_file():
                return parent
            if (parent / ".git").exists() and (parent / "CMakeLists.txt").is_file():
                return parent

    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"],
            cwd=str(start if start.is_dir() else start.parent),
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
        if out:
            return Path(out)
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        pass

    # Fallback: walk up from this file
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / PINS_REL).is_file():
            return parent
    raise FileNotFoundError(
        "Could not resolve FlexAIDdS repo root (need docs/implementation/arm_pins.json)"
    )


def _expand_env_path(raw: str, *, local_root: Path) -> Path:
    """Expand $FLEXAIDDS_LOCAL_ROOT / $HOME / ~ in pin paths."""
    s = raw.strip()
    s = s.replace("$FLEXAIDDS_LOCAL_ROOT", str(local_root))
    s = s.replace("${FLEXAIDDS_LOCAL_ROOT}", str(local_root))
    s = os.path.expandvars(s)
    s = os.path.expanduser(s)
    return Path(s)


def load_arm_pins(repo_root: Optional[Path] = None) -> Dict[str, Any]:
    """Load ``docs/implementation/arm_pins.json`` from the repo."""
    root = repo_root or resolve_repo_root()
    pins_path = root / PINS_REL
    if not pins_path.is_file():
        raise FileNotFoundError(f"arm_pins.json not found: {pins_path}")
    with pins_path.open("r", encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict) or "arms" not in data:
        raise ValueError(f"Invalid arm_pins schema at {pins_path}")
    return data


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _git_commit_for_arm(arm_spec: Mapping[str, Any]) -> Optional[str]:
    """Pinned source commit for the arm (full SHA when available)."""
    for key in ("source_commit", "source_commit_at_doc", "git_commit"):
        val = arm_spec.get(key)
        if isinstance(val, str) and val.strip():
            return val.strip()
    return None


def _expected_binary_path(
    arm: str,
    arm_spec: Mapping[str, Any],
    *,
    local_root: Path,
) -> Path:
    """Resolve expected binary path from pin ``binary_install_path`` or layout default."""
    install = arm_spec.get("binary_install_path")
    if isinstance(install, str) and install.strip():
        return _expand_env_path(install, local_root=local_root)

    names = _ARM_BINARY_NAMES.get(arm, ("FlexAID", "FlexAIDdS"))
    return local_root / BIN_REL / arm / names[0]


def _find_present_binary(
    arm: str,
    arm_spec: Mapping[str, Any],
    *,
    local_root: Path,
) -> Optional[Path]:
    """Return first existing executable path for the arm, else None."""
    candidates: List[Path] = []
    install = arm_spec.get("binary_install_path")
    if isinstance(install, str) and install.strip():
        candidates.append(_expand_env_path(install, local_root=local_root))
    for name in _ARM_BINARY_NAMES.get(arm, ("FlexAID", "FlexAIDdS")):
        p = local_root / BIN_REL / arm / name
        if p not in candidates:
            candidates.append(p)
    for p in candidates:
        if p.is_file():
            return p
    return None


def _reconstruction_note(arm: str, arm_spec: Mapping[str, Any]) -> str:
    label = arm_spec.get("reconstruction_label_if_missing")
    if isinstance(label, str) and label.strip():
        return label.strip()
    if arm == "A":
        return (
            "CF reconstruction on current FlexAIDdS --legacy TEMPER0 CLUSTA CF "
            "— not historical A SHA"
        )
    if arm == "B":
        return (
            "Entropy reconstruction — modern engine TEMPER21+FO, same axes "
            "— not first-entropy SHA / 3Dsig binary replica"
        )
    return "BUILD_FROM_CURRENT_TREE / reconstruction stub — binary not yet staged"


def write_receipt(path: Path, receipt: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(receipt, fh, indent=2, sort_keys=True)
        fh.write("\n")


def pin_arm_binaries(
    *,
    local_root: Optional[Path] = None,
    repo_root: Optional[Path] = None,
    allow_reconstruction: bool = False,
    arms: Sequence[str] = SCIENCE_ARMS,
) -> Dict[str, Any]:
    """Inspect arm binaries, write receipts, return structured result.

    Returns a dict with keys:
      - phase: "P1"
      - receipts: {arm: path}
      - arms: {arm: receipt_dict}
      - exit_code: int (0 or 2, or 1 for identical A/B)
      - messages: list[str]
    """
    root = repo_root or resolve_repo_root()
    lr = Path(
        local_root
        if local_root is not None
        else os.environ.get("FLEXAIDDS_LOCAL_ROOT", str(DEFAULT_LOCAL_ROOT))
    ).expanduser().resolve()
    pins = load_arm_pins(root)
    pin_arms: Mapping[str, Any] = pins.get("arms") or {}

    receipts_dir = lr / RECEIPTS_REL
    receipts_dir.mkdir(parents=True, exist_ok=True)

    ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    arm_results: Dict[str, Dict[str, Any]] = {}
    receipt_paths: Dict[str, str] = {}
    messages: List[str] = []
    missing_required: List[str] = []  # A/B missing without reconstruction path

    for arm in arms:
        if arm not in SCIENCE_ARMS:
            messages.append(f"skip unknown arm {arm}")
            continue
        spec = pin_arms.get(arm)
        if not isinstance(spec, Mapping):
            messages.append(f"WARN: arm {arm} not in arm_pins.json")
            spec = {}

        expected = _expected_binary_path(arm, spec, local_root=lr)
        present = _find_present_binary(arm, spec, local_root=lr)
        git_commit = _git_commit_for_arm(spec)
        science_id = spec.get("science_identity")
        receipt_path = receipts_dir / f"arm_{arm}_binary.json"

        if present is not None:
            digest = sha256_file(present)
            receipt: Dict[str, Any] = {
                "phase": "P1",
                "arm": arm,
                "binary_path": str(present.resolve()),
                "binary_sha256": digest,
                "git_commit": git_commit,
                "reconstruction": False,
                "status": "present",
                "science_identity": science_id,
                "expected_path": str(expected),
                "pins_file": str((root / PINS_REL).resolve()),
                "local_root": str(lr),
                "ts_utc": ts,
            }
            write_receipt(receipt_path, receipt)
            arm_results[arm] = receipt
            receipt_paths[arm] = str(receipt_path)
            messages.append(f"OK arm {arm}: present sha256={digest[:12]}… path={present}")
        else:
            # Missing binary — fail closed unless reconstruction allowed
            pin_status = spec.get("status") or "SOURCE_PINNED_BINARY_MISSING"
            if arm in ("A", "B") and pin_status == "BUILD_FROM_CURRENT_TREE":
                pin_status = "SOURCE_PINNED_BINARY_MISSING"

            receipt = {
                "phase": "P1",
                "arm": arm,
                "binary_path": str(expected),
                "binary_sha256": None,  # never invent digests for missing files
                "git_commit": git_commit,
                "reconstruction": bool(allow_reconstruction),
                "status": "SOURCE_PINNED_BINARY_MISSING"
                if arm in ("A", "B")
                else (pin_status if isinstance(pin_status, str) else "SOURCE_PINNED_BINARY_MISSING"),
                "science_identity": science_id,
                "expected_path": str(expected),
                "reconstruction_label": _reconstruction_note(arm, spec)
                if allow_reconstruction
                else None,
                "pins_file": str((root / PINS_REL).resolve()),
                "local_root": str(lr),
                "ts_utc": ts,
            }

            if allow_reconstruction:
                write_receipt(receipt_path, receipt)
                arm_results[arm] = receipt
                receipt_paths[arm] = str(receipt_path)
                messages.append(
                    f"RECONSTRUCTION arm {arm}: binary missing; labeled stub receipt → {receipt_path}"
                )
            else:
                # Still write nothing for missing without flag (fail closed).
                # Do not leave a present-looking receipt.
                arm_results[arm] = receipt
                if arm in ("A", "B"):
                    missing_required.append(arm)
                messages.append(
                    f"MISSING arm {arm}: expected {expected} "
                    f"(use --allow-reconstruction for labeled stub only)"
                )

    # Fail closed: A and B SHAs must differ when both present
    a_rec = arm_results.get("A") or {}
    b_rec = arm_results.get("B") or {}
    a_sha = a_rec.get("binary_sha256")
    b_sha = b_rec.get("binary_sha256")
    a_present = a_rec.get("status") == "present" and a_sha
    b_present = b_rec.get("status") == "present" and b_sha
    identical_ab = bool(a_present and b_present and a_sha == b_sha)

    if identical_ab:
        messages.append(
            "FAIL: arm A and B binary_sha256 are identical — not a claim split "
            "(B0 would be a deterministic twin of A)"
        )

    # Exit code policy (task contract):
    #   0 — A and B both present and SHAs differ
    #   0 — --allow-reconstruction and status documents missing (reconstruction receipts)
    #   2 — missing A/B without reconstruction flag
    #   1 — A/B both present but same SHA (fail closed claim split)
    if identical_ab:
        exit_code = 1
    elif missing_required:
        exit_code = 2
    elif a_present and b_present:
        exit_code = 0
    elif allow_reconstruction:
        # Reconstruction path: missing documented; exit 0 only when receipts label missing
        exit_code = 0
    else:
        # C-only or partial without A+B present and without reconstruction
        exit_code = 2

    summary = {
        "phase": "P1",
        "phase_status_line": "PHASE=P1",
        "local_root": str(lr),
        "repo_root": str(root),
        "receipts_dir": str(receipts_dir),
        "receipts": receipt_paths,
        "arms": arm_results,
        "a_present": bool(a_present),
        "b_present": bool(b_present),
        "claim_binary_split_ok": bool(a_present and b_present and not identical_ab),
        "identical_ab": identical_ab,
        "allow_reconstruction": allow_reconstruction,
        "exit_code": exit_code,
        "messages": messages,
        "ts_utc": ts,
    }
    return summary


def main(argv: Optional[Sequence[str]] = None) -> int:
    """CLI entry for Phase 1 binary pin/receipt helper."""
    import argparse

    ap = argparse.ArgumentParser(
        description=(
            "Phase P1: pin comparative arm binaries (A/B/C). "
            "SHA256 existing Mach-Os and write receipts; never invent digests."
        )
    )
    ap.add_argument(
        "--local-root",
        type=Path,
        default=None,
        help="Override FLEXAIDDS_LOCAL_ROOT (default: env or ~/flexaidds_results)",
    )
    ap.add_argument(
        "--repo-root",
        type=Path,
        default=None,
        help="FlexAIDdS repo root (default: git rev-parse / package parents)",
    )
    ap.add_argument(
        "--allow-reconstruction",
        action="store_true",
        help=(
            "If A/B binary missing, write labeled reconstruction stub receipts "
            "(binary_sha256=null; status=SOURCE_PINNED_BINARY_MISSING)"
        ),
    )
    ap.add_argument(
        "--json",
        action="store_true",
        help="Print full summary JSON to stdout",
    )
    args = ap.parse_args(list(argv) if argv is not None else None)

    # Always emit phase status line first for log grepping.
    print("PHASE=P1")

    try:
        result = pin_arm_binaries(
            local_root=args.local_root,
            repo_root=args.repo_root,
            allow_reconstruction=args.allow_reconstruction,
        )
    except (FileNotFoundError, ValueError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    for msg in result.get("messages") or []:
        print(msg)

    if args.json:
        # Compact serializable dump (arms already plain dicts)
        print(json.dumps(result, indent=2, sort_keys=True))

    print(
        f"claim_binary_split_ok={result.get('claim_binary_split_ok')} "
        f"exit_code={result.get('exit_code')}"
    )
    return int(result.get("exit_code", 2))


if __name__ == "__main__":
    sys.exit(main())
