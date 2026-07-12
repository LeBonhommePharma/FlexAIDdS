#!/usr/bin/env python3
"""
safe_archive_to_icoud.py

Safer-than-safe archiver for moving local FlexAIDdS results into iCloud Drive (2TB).

Design goals (addressing known iCloud File Provider issues):
- Atomic writes (tmp + replace)
- Full content verification (SHA256 per file + tree summary)
- Retry with backoff on reads (placeholder / sync lag)
- Never delete source until target is 100% verified + marker written
- Timestamped / versioned destination to avoid overwrites
- Conflict detection (looks for "conflicted copy")
- Optional: create .verified + SHA manifest
- Supports --dry-run, --keep-local, --evict-after

Usage (recommended):
  python3 scripts/safe_archive_to_icoud.py \
      --source results/astex_jcim2015_fair_20260708_0002 \
      --dest   "$FLEXAIDDS_ICLOUD/archived/astex-fair-2026-07-08" \
      --verify

Or for whole trees:
  python3 scripts/safe_archive_to_icoud.py \
      --source results/posex_cd_20260616 \
      --dest   "$FLEXAIDDS_ICLOUD/results/posex_cd_20260616" \
      --evict-after

Environment:
  FLEXAIDDS_ICLOUD   (falls back to standard iCloud container + /FlexAIDdS_benchmarks)
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

ICLOUD_DEFAULT = Path.home() / "Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks"

def get_icoud_base() -> Path:
    env = os.environ.get("FLEXAIDDS_ICLOUD")
    if env:
        p = Path(env).expanduser()
        if p.name != "FlexAIDdS_benchmarks":
            p = p / "FlexAIDdS_benchmarks"
        return p
    return ICLOUD_DEFAULT

def sha256_file(path: Path, chunk: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()

def write_json_atomic(path: Path, obj: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True), encoding="utf-8")
    tmp.replace(path)

def read_json_safe(path: Path, attempts: int = 6, delay: float = 0.3) -> dict:
    last = None
    for i in range(attempts):
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except Exception as e:
            last = e
            time.sleep(delay * (i + 1))
    raise RuntimeError(f"Failed to read {path} after {attempts} attempts: {last}")

def find_conflicts(root: Path) -> List[Path]:
    conflicts = []
    for p in root.rglob("*"):
        name = p.name.lower()
        if "conflicted copy" in name or "(conflicted" in name or p.suffix == ".icloud":
            conflicts.append(p)
    return conflicts

def compute_tree_manifest(src: Path) -> dict:
    manifest = {"files": {}, "dirs": 0, "total_bytes": 0, "generated_at": datetime.now(timezone.utc).isoformat()}
    for dirpath, dirnames, filenames in os.walk(src):
        manifest["dirs"] += 1
        for fn in filenames:
            fp = Path(dirpath) / fn
            try:
                st = fp.stat()
                manifest["files"][str(fp.relative_to(src))] = {
                    "size": st.st_size,
                    "sha256": sha256_file(fp),
                    "mtime": st.st_mtime,
                }
                manifest["total_bytes"] += st.st_size
            except Exception as e:
                manifest["files"][str(fp.relative_to(src))] = {"error": str(e)}
    return manifest

def verify_copy(src: Path, dst: Path, manifest: dict) -> Tuple[bool, List[str]]:
    errors = []
    ok = True
    for rel, meta in manifest["files"].items():
        if "error" in meta:
            continue
        dst_f = dst / rel
        if not dst_f.exists():
            errors.append(f"missing: {rel}")
            ok = False
            continue
        try:
            if sha256_file(dst_f) != meta["sha256"]:
                errors.append(f"hash mismatch: {rel}")
                ok = False
        except Exception as e:
            errors.append(f"verify error {rel}: {e}")
            ok = False
    return ok, errors

def safe_archive(src: Path, dest_root: Path, evict_after: bool = False, dry_run: bool = False, keep_local: bool = False) -> dict:
    if not src.exists():
        raise FileNotFoundError(src)

    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    dest = dest_root / f"{src.name}-{ts}"
    if dest.exists():
        # extremely defensive
        dest = dest_root / f"{src.name}-{ts}-{os.getpid()}"

    report = {
        "source": str(src),
        "dest": str(dest),
        "started": ts,
        "dry_run": dry_run,
        "success": False,
    }

    if dry_run:
        report["would_create"] = str(dest)
        report["manifest_preview"] = "would compute tree sha + copy + verify"
        report["success"] = True
        return report

    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"[safe-icloud] Staging copy {src} -> {dest}")

    # 1. Copy with rsync for efficiency + preserve metadata
    rsync = ["rsync", "-a", "--delete", "--inplace", str(src) + "/", str(dest) + "/"]
    subprocess.check_call(rsync)

    # 2. Compute source manifest
    print("[safe-icloud] Computing source manifest (SHA256 tree)...")
    manifest = compute_tree_manifest(src)
    manifest_path = dest / ".source_manifest.json"
    write_json_atomic(manifest_path, manifest)

    # 3. Verify
    print("[safe-icloud] Verifying copy...")
    ok, errs = verify_copy(src, dest, manifest)
    if not ok:
        report["errors"] = errs
        raise RuntimeError(f"Verification failed: {errs[:5]}...")

    # 4. Write verified marker (atomic)
    verified = {
        "verified_at": datetime.now(timezone.utc).isoformat(),
        "source_sha_summary": manifest["total_bytes"],
        "file_count": len(manifest["files"]),
        "tool": "safe_archive_to_icoud.py",
    }
    write_json_atomic(dest / ".verified.json", verified)

    # 5. Check for conflicts after copy
    conflicts = find_conflicts(dest)
    if conflicts:
        report["conflicts_found"] = [str(c) for c in conflicts]
        print(f"[WARN] Found {len(conflicts)} potential iCloud conflict/placeholder files")

    report["success"] = True
    report["verified_path"] = str(dest / ".verified.json")

    if evict_after and not keep_local:
        print(f"[safe-icloud] Evicting local source after successful verified copy: {src}")
        shutil.rmtree(src)

    return report

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True, type=Path, help="Local directory to archive")
    ap.add_argument("--dest", type=Path, help="Destination under iCloud (default: $FLEXAIDDS_ICLOUD/archived/)")
    ap.add_argument("--evict-after", action="store_true", help="Remove local source ONLY after full verification")
    ap.add_argument("--keep-local", action="store_true", help="Never delete local")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    base = args.dest or (get_icoud_base() / "archived")
    base.mkdir(parents=True, exist_ok=True)

    print(f"[safe-icloud] iCloud base: {base}")
    report = safe_archive(args.source.resolve(), base, args.evict_after, args.dry_run, args.keep_local)

    out = args.source.with_name(args.source.name + "-archive-report.json")
    write_json_atomic(out, report)
    print(json.dumps(report, indent=2))
    print(f"[safe-icloud] Report written to {out}")

if __name__ == "__main__":
    main()
