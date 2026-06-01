#!/usr/bin/env python3
"""
iCloud FS Sanity Checker for RE-DOCK / FlexAIDdS benchmarks on M3 Pro.

Purpose:
  Verify that read / write / executable-bit / rapid round-trips work reliably
  when all campaign state lives on the 2 TB iCloud Drive (File Provider).

  This directly addresses latency, placeholder eviction, sync races,
  permission stripping on write, and delayed visibility that can bite
  Codex / Grok Build / Claude Code resume workflows.

Usage (recommended on this machine):
  python benchmarks/re-dock/icloud_fs_check.py \
      --path ~/Library/Mobile\ Documents/com~apple~CloudDocs/FlexAIDdS/re-dock-test-fs

  # Or using the standard env
  source ~/.flexaidds_env
  python benchmarks/re-dock/icloud_fs_check.py \
      --path "$FLEXAIDDS_ICLOUD/re-dock-test-fs-$(date +%s)"

It will:
- Create the directory (if needed) under iCloud
- Perform multiple write/read/compare cycles for JSON (checkpoint-like)
- Write a self-contained worker-style .py script, chmod +x, execute it
- Do rapid churn (write + immediate re-read + re-exec) to surface races
- Report timing, any errors, file provider state if detectable
- Clean up its own test files (unless --keep)

Exit code 0 = all tests passed on the real iCloud volume.
Non-zero = problems found that would affect real Codex/local resume campaigns.

This tool itself lives in the source tree. All test data it generates
must be (and is) written under the iCloud tree per the Storage Invariant.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def write_json_atomic(path: Path, obj: Dict[str, Any]) -> None:
    """Atomic write to reduce partial-file races with iCloud sync."""
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True), encoding="utf-8")
    tmp.replace(path)  # atomic on APFS


def read_json_safe(path: Path, max_attempts: int = 5, delay: float = 0.2) -> Dict[str, Any]:
    """Read with small retry loop — iCloud can briefly return placeholder or stale data."""
    last_err = None
    for attempt in range(max_attempts):
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except Exception as e:
            last_err = e
            time.sleep(delay * (attempt + 1))
    raise RuntimeError(f"Failed to read {path} after {max_attempts} attempts: {last_err}")


def make_test_script(path: Path, chunk_id: str) -> Path:
    """Create a realistic self-contained worker script (like to_worker_script output)."""
    script = f'''#!/usr/bin/env python3
"""iCloud FS test worker script — chunk {chunk_id}"""
import json, random, hashlib, time
print("iCloud FS test worker starting...")
seed = int(hashlib.md5(b"{chunk_id}").hexdigest()[:8], 16)
random.seed(seed)
energies = [random.gauss(-8.5, 1.2) for _ in range(50)]
result = {{
    "chunk_id": "{chunk_id}",
    "best_energy": min(energies),
    "mean_energy": sum(energies) / len(energies),
    "n_poses": len(energies),
    "timestamp": time.time(),
    "icloud_test": True,
}}
print(json.dumps(result))
'''
    path.write_text(script, encoding="utf-8")
    path.chmod(0o755)
    return path


def run_script(path: Path, timeout: float = 10.0) -> Dict[str, Any]:
    """Execute the script and parse its JSON stdout (like real Codex ingestion)."""
    proc = subprocess.run(
        [sys.executable, str(path)],
        capture_output=True,
        text=True,
        timeout=timeout,
        cwd=path.parent,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"Script failed (rc={proc.returncode}):\n{proc.stderr}")
    # Find the last JSON-looking line
    for line in reversed(proc.stdout.strip().splitlines()):
        line = line.strip()
        if line.startswith("{") and line.endswith("}"):
            return json.loads(line)
    raise RuntimeError(f"No JSON result found in script output:\n{proc.stdout}")


def test_icloud_path(base: Path, keep: bool = False) -> Dict[str, Any]:
    """Run the full battery of read/write/exec tests under the given iCloud path."""
    base.mkdir(parents=True, exist_ok=True)
    test_dir = base / f"fs-test-{datetime.now().strftime('%Y%m%d-%H%M%S')}"
    test_dir.mkdir()

    results: Dict[str, Any] = {
        "test_dir": str(test_dir),
        "started_at": time.time(),
        "tests": [],
        "success": True,
        "errors": [],
    }

    try:
        # 1. JSON checkpoint-style roundtrip (multiple times + rapid churn)
        for i in range(3):
            cp = test_dir / f"checkpoint-gen{i}.json"
            data = {
                "campaign_id": "icloud-fs-test",
                "generation": i,
                "timestamp": time.time(),
                "temperatures": [298.0 + x * 40 for x in range(8)],
                "note": "Testing iCloud File Provider behavior",
            }
            write_json_atomic(cp, data)
            time.sleep(0.05)  # tiny gap to let sync breathe
            read_back = read_json_safe(cp)
            if read_back != data:
                raise AssertionError(f"JSON roundtrip mismatch on iteration {i}")
            results["tests"].append(f"json_roundtrip_{i}: OK")

        # Rapid churn test (write + immediate read 5 times)
        churn_file = test_dir / "churn.json"
        for j in range(5):
            payload = {"iteration": j, "value": f"data-{j}", "t": time.time()}
            write_json_atomic(churn_file, payload)
            back = read_json_safe(churn_file, max_attempts=3, delay=0.1)
            if back.get("iteration") != j:
                results["errors"].append(f"churn iteration {j} saw stale data")
                results["success"] = False
        results["tests"].append("rapid_churn: OK" if results["success"] else "rapid_churn: PARTIAL")

        # 2. Executable script write + chmod + execution
        script_path = test_dir / "test_worker_chunk.py"
        script_result = run_script(make_test_script(script_path, "icloud-test-001"))
        if not script_result.get("icloud_test"):
            raise AssertionError("Script did not report icloud_test flag")
        results["tests"].append("executable_script: OK")

        # 3. Subdirectory + nested result (mimics gen_N / chunk results)
        gen_dir = test_dir / "gen_042"
        gen_dir.mkdir()
        nested_result = gen_dir / "1a30_T298_g42.json"
        write_json_atomic(nested_result, {"pdb_id": "1a30", "best_energy": -12.34, "from": "codex-sim"})
        reread = read_json_safe(nested_result)
        assert reread["pdb_id"] == "1a30"
        results["tests"].append("nested_subdir: OK")

        # 4. Permission / executable bit persistence after write
        mode_after = script_path.stat().st_mode & 0o777
        if mode_after != 0o755:
            results["errors"].append(f"Executable bit not preserved: got {oct(mode_after)}")
            results["success"] = False
        results["tests"].append("executable_bit_persisted: " + ("OK" if mode_after == 0o755 else "FAIL"))

    except Exception as e:
        results["success"] = False
        results["errors"].append(str(e))
    finally:
        results["finished_at"] = time.time()
        results["duration_s"] = results["finished_at"] - results["started_at"]

        if not keep:
            try:
                shutil.rmtree(test_dir)
            except Exception:
                pass  # best effort cleanup

    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="iCloud FS read/write/exec sanity checker for RE-DOCK")
    parser.add_argument(
        "--path",
        type=Path,
        default=None,
        help="iCloud directory to test inside (will be created). "
             "Defaults to $FLEXAIDDS_ICLOUD/re-dock-fs-tests or a safe fallback.",
    )
    parser.add_argument("--keep", action="store_true", help="Do not delete test files after run")
    args = parser.parse_args()

    if args.path is None:
        icloud_base = os.environ.get("FLEXAIDDS_ICLOUD")
        if icloud_base:
            base = Path(icloud_base) / "re-dock-fs-tests"
        else:
            base = Path.home() / "Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS/re-dock-fs-tests"
    else:
        base = args.path

    print(f"=== iCloud FS Sanity Check ===")
    print(f"Target base: {base}")
    print(f"Started: {datetime.now().isoformat()}")
    print()

    report = test_icloud_path(base, keep=args.keep)

    print(json.dumps(report, indent=2, default=str))

    if report["success"]:
        print("\n✅ All iCloud FS tests PASSED on this volume.")
        sys.exit(0)
    else:
        print("\n❌ iCloud FS problems detected. See errors above.")
        print("This would likely break Codex → Grok Build resume or local campaign runs.")
        sys.exit(1)


if __name__ == "__main__":
    main()
