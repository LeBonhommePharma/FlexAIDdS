#!/usr/bin/env python3
"""Assert FlexAIDdS hidden stat markers match a repo-stats.json snapshot.

Used by update-site.yml both before upload (local HTML) and after deploy-pages
(live https://thebonhomme.com/FlexAIDdS/). Exits 0 on match.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.request

MARKER_IDS = (
    "stat-commits",
    "stat-langs",
    "stat-stars",
    "last-updated",
    "latest-release",
)


def parse_markers(html: str) -> dict[str, str]:
    found: dict[str, str] = {}
    for marker_id in MARKER_IDS:
        match = re.search(
            rf'id="{re.escape(marker_id)}"[^>]*>([^<]*)<',
            html,
        )
        if match:
            found[marker_id] = match.group(1).strip()
    return found


def load_snapshot(path: str) -> dict:
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, dict) or "commits" not in data:
        raise SystemExit(f"{path} is not a valid repo-stats snapshot")
    return data


def expected_markers(snapshot: dict) -> dict[str, str]:
    expected = {
        "stat-commits": str(snapshot["commits"]),
        "stat-langs": str(snapshot.get("languageCount", "")),
    }
    if snapshot.get("stars") is not None:
        expected["stat-stars"] = str(snapshot["stars"])
    if snapshot.get("latestRelease"):
        expected["latest-release"] = str(snapshot["latestRelease"])
    if snapshot.get("lastUpdated"):
        expected["last-updated"] = str(snapshot["lastUpdated"])
    return {key: value for key, value in expected.items() if value != ""}


def compare(actual: dict[str, str], expected: dict[str, str]) -> list[str]:
    errors: list[str] = []
    for key, want in expected.items():
        got = actual.get(key)
        if got != want:
            errors.append(f"{key}: live={got!r} expected={want!r}")
    return errors


def fetch_html(url: str, timeout: int) -> str:
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "FlexAIDdS-site-stats-verify",
            "Cache-Control": "no-cache",
            "Pragma": "no-cache",
        },
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8", errors="replace")


def cache_busted_url(url: str) -> str:
    """Append a cache-buster so CDN/HTML shells are not served stale."""
    sep = "&" if "?" in url else "?"
    return f"{url}{sep}t={int(time.time())}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", required=True, help="Path to repo-stats.json")
    parser.add_argument("--html", help="Local HTML to check (no network)")
    parser.add_argument("--url", help="Live URL to fetch")
    parser.add_argument("--local-only", action="store_true")
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument("--interval", type=int, default=15)
    args = parser.parse_args()

    snapshot = load_snapshot(args.json)
    expected = expected_markers(snapshot)
    print("expected:", expected)

    if args.local_only or args.html:
        html_path = args.html
        if not html_path:
            print("Error: --html is required with --local-only", file=sys.stderr)
            return 1
        with open(html_path, encoding="utf-8") as f:
            html = f.read()
        errors = compare(parse_markers(html), expected)
        if errors:
            print("local HTML markers do not match snapshot:", file=sys.stderr)
            print("\n".join(errors), file=sys.stderr)
            return 1
        print(f"local HTML markers match {html_path}")
        if args.local_only:
            return 0

    if not args.url:
        return 0

    deadline = time.time() + args.timeout
    last_errors: list[str] = ["live fetch not attempted"]
    while time.time() < deadline:
        try:
            html = fetch_html(cache_busted_url(args.url), timeout=30)
            last_errors = compare(parse_markers(html), expected)
            if not last_errors:
                print(f"live markers match {args.url}")
                return 0
            print("not yet:", "; ".join(last_errors), flush=True)
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as exc:
            last_errors = [str(exc)]
            print(f"fetch failed: {exc}", flush=True)
        time.sleep(args.interval)

    print("live markers still stale after timeout:", file=sys.stderr)
    print("\n".join(last_errors), file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
