#!/usr/bin/env python3
"""Bonhomme Fleet Status HTTP Server — Lightweight API for FlexAIDdS Dashboard.

Serves fleet status JSON files and deep-sanned progress information via a
REST API. Designed to be compatible with the BonhommeViewer TypeScript PWA.

Endpoints:
  GET /                 — API index with endpoint listing
  GET /api/status       — Merged fleet status from all runners
  GET /api/progress/<dataset>/<run> — Deep-scan progress for a specific run
  GET /api/fleet/<filename>         — Raw fleet status JSON file

CORS headers (Access-Control-Allow-Origin: *) are added to all responses.

Usage:
  # Start on default port 8787
  python3 benchmarks/m3pro/dashboard/fleet_status_server.py

  # Custom port
  python3 benchmarks/m3pro/dashboard/fleet_status_server.py --port 9000

  # Override paths
  python3 benchmarks/m3pro/dashboard/fleet_status_server.py \\
      --icloud ~/Library/Mobile\\ Documents/com~apple~CloudDocs/FlexAIDdS \\
      --results ~/.flexaidds_fast/results

Environment variables:
  FLEXAIDDS_FAST_BASE  Override fast results path
  FLEXAIDDS_ICLOUD     Override iCloud FlexAIDdS path

Apache-2.0 (c) 2026 NRGlab, Universite de Montreal
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from datetime import datetime, timezone
from glob import glob
from http.server import HTTPServer, BaseHTTPRequestHandler
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import unquote

__version__ = "1.0.0"

ICLOUD_DEFAULT = os.path.expanduser(
    "~/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS"
)
FAST_DEFAULT = os.path.expanduser("~/.flexaidds_fast")
RESULTS_SUBDIR = "results/tier2"


def safe_read_json(filepath: str) -> Optional[Dict[str, Any]]:
    """Read and parse a JSON file, handling partial writes and iCloud eviction."""
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError, UnicodeDecodeError):
        return None


def discover_fleet_status_files(icloud_dir: str) -> List[str]:
    """Find all fleet_status*.json files in the iCloud directory."""
    pattern = os.path.join(icloud_dir, "fleet_status*.json")
    return sorted(glob(pattern))


def merge_fleet_status(icloud_dir: str) -> Dict[str, Any]:
    """Merge all fleet status JSON files into a unified response.

    Args:
        icloud_dir: Path to iCloud FlexAIDdS directory.

    Returns:
        Merged status dictionary with per-runner breakdowns.
    """
    files = discover_fleet_status_files(icloud_dir)
    merged: Dict[str, Any] = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "runners": {},
        "total_workers": 0,
        "fleet_files": [os.path.basename(f) for f in files],
    }

    for fpath in files:
        data = safe_read_json(fpath)
        if data is None:
            continue
        runner_name = data.get("runner", os.path.basename(fpath))
        merged["runners"][runner_name] = data

    return merged


def deep_scan_run(
    results_dirs: List[str],
    dataset: str,
    run_name: str,
) -> Dict[str, Any]:
    """Deep-scan result directories for per-target progress.

    Scans all configured result directories for a specific dataset+run
    combination, reporting per-target status (done/in_progress/stuck/queued).

    Args:
        results_dirs: List of base result directories.
        dataset: Dataset name (e.g., 'astex').
        run_name: Run directory name (e.g., 'run31').

    Returns:
        Dictionary with per-target progress and summary stats.
    """
    result: Dict[str, Any] = {
        "dataset": dataset,
        "run": run_name,
        "targets": {},
        "summary": {
            "done": 0,
            "in_progress": 0,
            "stuck": 0,
            "queued": 0,
            "total": 0,
        },
        "scanned_dirs": [],
    }

    for base in results_dirs:
        tier2 = base if base.endswith(RESULTS_SUBDIR) else os.path.join(base, RESULTS_SUBDIR)
        run_path = os.path.join(tier2, dataset, run_name)

        if not os.path.isdir(run_path):
            continue

        result["scanned_dirs"].append(run_path)

        try:
            entries = sorted(os.listdir(run_path))
        except OSError:
            continue

        for entry in entries:
            entry_path = os.path.join(run_path, entry)
            if not os.path.isdir(entry_path):
                continue
            if entry.startswith("."):
                continue
            if entry in ("logs", "tmp"):
                continue

            target_status = _classify_target(entry_path)
            result["targets"][entry] = target_status
            result["summary"][target_status["state"]] += 1
            result["summary"]["total"] += 1

    return result


def _classify_target(target_path: str) -> Dict[str, Any]:
    """Classify a single target directory's status.

    Args:
        target_path: Absolute path to the target directory.

    Returns:
        Dictionary with 'state', 'has_result', 'is_stuck' keys.
    """
    status: Dict[str, Any] = {
        "state": "queued",
        "has_result": False,
        "is_stuck": False,
    }

    try:
        entries = set(os.listdir(target_path))
    except OSError:
        return status

    has_result_csv = "result.csv" in entries
    has_dock_config = "dock_config.json" in entries
    has_stdout = "stdout.log" in entries

    has_output_pdb = any(
        e.endswith(".pdb") and not e.endswith("_INI.pdb")
        for e in entries
    )

    status["has_result"] = has_result_csv or has_output_pdb

    if has_stdout:
        stdout_path = os.path.join(target_path, "stdout.log")
        try:
            with open(stdout_path, "r", encoding="utf-8", errors="replace") as f:
                content = f.read()
            if "[STUCK]" in content:
                status["is_stuck"] = True
                status["state"] = "stuck"
                return status
        except OSError:
            pass

    if has_result_csv or has_output_pdb:
        status["state"] = "done"
    elif has_dock_config:
        status["state"] = "in_progress"
    else:
        status["state"] = "queued"

    return status


class FleetStatusHandler(BaseHTTPRequestHandler):
    """HTTP request handler for the fleet status API."""

    icloud_dir: str = ICLOUD_DEFAULT
    results_dirs: List[str] = []

    def log_message(self, format: str, *args: Any) -> None:
        """Override to use a cleaner log format."""
        sys.stderr.write(
            f"[fleet_status_server] {self.address_string()} - {format % args}\n"
        )

    def _send_json(self, data: Any, status_code: int = 200) -> None:
        """Send a JSON response with CORS headers."""
        body = json.dumps(data, indent=2, default=str).encode("utf-8")
        self.send_response(status_code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.send_header("Cache-Control", "no-cache")
        self.end_headers()
        self.wfile.write(body)

    def _send_error_json(self, message: str, status_code: int = 404) -> None:
        """Send a JSON error response."""
        self._send_json({"error": message, "status": status_code}, status_code)

    def do_OPTIONS(self) -> None:
        """Handle CORS preflight requests."""
        self.send_response(204)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.end_headers()

    def do_GET(self) -> None:
        """Route GET requests to appropriate handlers."""
        path = unquote(self.path).rstrip("/")

        if path == "" or path == "/":
            self._handle_index()
        elif path == "/api/status":
            self._handle_status()
        elif path.startswith("/api/progress/"):
            self._handle_progress(path)
        elif path.startswith("/api/fleet/"):
            self._handle_fleet_file(path)
        else:
            self._send_error_json(f"Unknown endpoint: {path}", 404)

    def _handle_index(self) -> None:
        """Handle GET / — API index."""
        self._send_json({
            "service": "Bonhomme Fleet Status Server",
            "version": __version__,
            "endpoints": {
                "GET /": "This index",
                "GET /api/status": "Merged fleet status from all runners",
                "GET /api/progress/<dataset>/<run>": "Deep-scan progress for a run",
                "GET /api/fleet/<filename>": "Raw fleet status JSON file",
            },
            "icloud_dir": self.icloud_dir,
            "results_dirs": self.results_dirs,
        })

    def _handle_status(self) -> None:
        """Handle GET /api/status — merged fleet status."""
        merged = merge_fleet_status(self.icloud_dir)
        self._send_json(merged)

    def _handle_progress(self, path: str) -> None:
        """Handle GET /api/progress/<dataset>/<run> — deep scan."""
        parts = path.replace("/api/progress/", "").split("/")
        if len(parts) < 2:
            self._send_error_json(
                "Usage: /api/progress/<dataset>/<run>", 400
            )
            return

        dataset = parts[0]
        run_name = parts[1]

        progress = deep_scan_run(self.results_dirs, dataset, run_name)
        self._send_json(progress)

    def _handle_fleet_file(self, path: str) -> None:
        """Handle GET /api/fleet/<filename> — raw fleet status file."""
        filename = path.replace("/api/fleet/", "")
        if not re.match(r"^fleet_status[\w.-]*\.json$", filename):
            self._send_error_json("Invalid filename pattern", 400)
            return

        filepath = os.path.join(self.icloud_dir, filename)
        filepath = os.path.realpath(filepath)

        if not filepath.startswith(os.path.realpath(self.icloud_dir)):
            self._send_error_json("Path traversal denied", 403)
            return

        data = safe_read_json(filepath)
        if data is None:
            self._send_error_json(f"File not found or unreadable: {filename}", 404)
            return

        self._send_json(data)


def build_argparser() -> argparse.ArgumentParser:
    """Build the CLI argument parser."""
    parser = argparse.ArgumentParser(
        prog="fleet_status_server",
        description="Bonhomme Fleet Status HTTP Server",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=8787,
        help="HTTP server port (default: 8787)",
    )
    parser.add_argument(
        "--host",
        type=str,
        default="0.0.0.0",
        help="Bind address (default: 0.0.0.0)",
    )
    parser.add_argument(
        "--icloud",
        type=str,
        default=None,
        help="Override iCloud status directory",
    )
    parser.add_argument(
        "--results",
        type=str,
        nargs="*",
        default=None,
        help="Override results directories (can specify multiple)",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    """Main entry point for the fleet status server.

    Args:
        argv: Command-line arguments (defaults to sys.argv[1:]).

    Returns:
        Exit code (0 for success).
    """
    parser = build_argparser()
    args = parser.parse_args(argv)

    icloud_dir = args.icloud or os.environ.get(
        "FLEXAIDDS_ICLOUD", ICLOUD_DEFAULT
    )
    fast_base = os.environ.get("FLEXAIDDS_FAST_BASE", FAST_DEFAULT)

    results_dirs: List[str] = []
    if args.results:
        results_dirs = args.results
    else:
        fast_tier2 = os.path.join(fast_base, RESULTS_SUBDIR)
        if os.path.isdir(fast_tier2):
            results_dirs.append(fast_base)
        icloud_results = os.path.join(icloud_dir, "results")
        icloud_tier2 = os.path.join(icloud_results, "tier2")
        if os.path.isdir(icloud_tier2):
            results_dirs.append(icloud_results)
        elif os.path.isdir(os.path.join(icloud_results, RESULTS_SUBDIR)):
            results_dirs.append(icloud_results)

    FleetStatusHandler.icloud_dir = icloud_dir
    FleetStatusHandler.results_dirs = results_dirs

    server = HTTPServer((args.host, args.port), FleetStatusHandler)

    print(f"Bonhomme Fleet Status Server v{__version__}")
    print(f"  Listening: http://{args.host}:{args.port}")
    print(f"  iCloud dir: {icloud_dir}")
    print(f"  Results dirs: {results_dirs}")
    print(f"  Endpoints:")
    print(f"    GET /api/status       — Merged fleet status")
    print(f"    GET /api/progress/<dataset>/<run> — Deep scan")
    print(f"    GET /api/fleet/<file> — Raw fleet JSON")
    print(f"  Press Ctrl+C to stop")
    print()

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down...")
        server.shutdown()

    return 0


if __name__ == "__main__":
    sys.exit(main())
