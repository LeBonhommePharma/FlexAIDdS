#!/usr/bin/env python3
"""Update hardcoded repository stats in site/index.html.

Reads commit count from git and language breakdown from GitHub API,
then patches site/index.html in-place. Uses only the standard library.

Usage:
    python scripts/update_site_stats.py [--repo OWNER/REPO] [--html PATH]
"""

from __future__ import annotations

import argparse
import datetime
import json
import os
import re
import subprocess
import sys
import urllib.error
import urllib.request

# GitHub language name → (CSS variable suffix, short display name)
LANG_MAP = {
    "C++": ("cpp", "C++"),
    "Python": ("python", "Python"),
    "Swift": ("swift", "Swift"),
    "Objective-C++": ("objcpp", "Obj-C++"),
    "CMake": ("cmake", "CMake"),
    "TypeScript": ("ts", "TypeScript"),
    "Cuda": ("cuda", "CUDA"),
    "CUDA": ("cuda", "CUDA"),
    "C": ("c", "C"),
    "Shell": ("other", "Shell"),
    "Metal": ("other", "Metal"),
    "JavaScript": ("other", "JavaScript"),
}

MIN_PERCENT = 1.0  # Languages below this are grouped into "Other"


def get_commit_count() -> int:
    """Return total commit count on the current branch."""
    result = subprocess.run(
        ["git", "rev-list", "--count", "HEAD"],
        capture_output=True,
        text=True,
        check=True,
    )
    return int(result.stdout.strip())


def fetch_languages(repo: str) -> dict[str, int]:
    """Fetch language byte counts from GitHub API."""
    data = _github_api_get(f"/repos/{repo}/languages")
    return data if isinstance(data, dict) else {}


def _github_api_get(path: str) -> dict | list | None:
    """Fetch a GitHub API endpoint. Returns parsed JSON or None on failure."""
    url = f"https://api.github.com{path}"
    req = urllib.request.Request(url)
    req.add_header("Accept", "application/vnd.github+json")
    req.add_header("User-Agent", "FlexAIDdS-site-updater")

    token = os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
    if token:
        req.add_header("Authorization", f"Bearer {token}")

    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            return json.loads(resp.read().decode())
    except (urllib.error.URLError, urllib.error.HTTPError) as e:
        print(f"Warning: GitHub API request failed ({path}): {e}", file=sys.stderr)
        return None


def fetch_stars(repo: str) -> int | None:
    """Fetch the current stargazers count."""
    data = _github_api_get(f"/repos/{repo}")
    if data and isinstance(data, dict):
        return data.get("stargazers_count")
    return None


def fetch_latest_release(repo: str) -> str | None:
    """Fetch the latest release tag name (e.g. 'v2.0.0')."""
    data = _github_api_get(f"/repos/{repo}/releases/latest")
    if data and isinstance(data, dict):
        return data.get("tag_name")
    return None


def compute_percentages(languages: dict[str, int]) -> list[tuple[str, str, float]]:
    """Compute (css_var_suffix, display_name, percentage) sorted by percentage desc.

    Languages below MIN_PERCENT are grouped into 'Other'.
    """
    total = sum(languages.values())
    if total == 0:
        return []

    entries: list[tuple[str, str, float]] = []
    other_pct = 0.0

    for lang, bytes_count in sorted(languages.items(), key=lambda x: -x[1]):
        pct = round(bytes_count / total * 100, 1)
        if lang in LANG_MAP:
            css_suffix, display = LANG_MAP[lang]
            if pct < MIN_PERCENT:
                other_pct += pct
            else:
                entries.append((css_suffix, display, pct))
        else:
            other_pct += pct

    if other_pct > 0:
        entries.append(("other", "Other", round(other_pct, 1)))

    return entries


def count_source_languages(languages: dict[str, int]) -> int:
    """Count source languages for the stats badge.

    GitHub's languages API can include generated dependency artifacts under
    "Makefile"; keep that out of the user-facing source-language total.
    """
    return sum(1 for lang in languages if lang != "Makefile")


def build_lang_bar(entries: list[tuple[str, str, float]]) -> str:
    """Build the lang-bar HTML block."""
    lines = ['        <div class="lang-bar" aria-label="Language breakdown">']
    for css_suffix, display, pct in entries:
        lines.append(
            f'          <div class="lang-segment" style="width:{pct}%;'
            f'background:var(--lang-{css_suffix})" title="{display} {pct}%"></div>'
        )
    lines.append("        </div>")
    return "\n".join(lines)


def build_lang_legend(entries: list[tuple[str, str, float]]) -> str:
    """Build the lang-legend HTML block."""
    lines = ['        <div class="lang-legend">']
    for css_suffix, display, pct in entries:
        lines.append(
            f'          <span><i style="background:var(--lang-{css_suffix})"></i>'
            f"{display} {pct}%</span>"
        )
    lines.append("        </div>")
    return "\n".join(lines)


def update_html(
    html_path: str,
    commit_count: int,
    lang_entries: list[tuple[str, str, float]],
    *,
    language_count: int | None = None,
    stars: int | None = None,
    release: str | None = None,
) -> bool:
    """Patch the HTML file in-place. Returns True if changes were made."""
    with open(html_path, "r", encoding="utf-8") as f:
        content = f.read()

    original = content

    # 1. Update commit count in the new semantic marker (id="stat-commits") — any tag
    content = re.sub(
        r'(<[^>]*id="stat-commits"[^>]*>)\d+(</[^>]+>)',
        rf"\g<1>{commit_count}\g<2>",
        content,
        count=1,
    )

    # 2. Update language count in the new semantic marker (id="stat-langs") — any tag
    if language_count is not None:
        content = re.sub(
            r'(<[^>]*id="stat-langs"[^>]*>)\d+(</[^>]+>)',
            rf"\g<1>{language_count}\g<2>",
            content,
            count=1,
        )

    # 3. (Optional) Stars - only if a matching span exists in future revisions
    if stars is not None:
        content = re.sub(
            r'(<span[^>]*id="stat-stars"[^>]*>)\d+(</span>)',
            rf"\g<1>{stars}\g<2>",
            content,
            count=1,
        )

    # 4. Update "last updated" date if a marker span is present (future-proof)
    today = datetime.date.today().isoformat()
    content = re.sub(
        r'(<span[^>]*id="last-updated"[^>]*>)[^<]*(</span>)',
        rf"\g<1>{today}\g<2>",
        content,
        count=1,
    )

    # 5. Update latest release version if a marker span is present
    if release:
        content = re.sub(
            r'(<span[^>]*id="latest-release"[^>]*>)[^<]*(</span>)',
            rf"\g<1>{release}\g<2>",
            content,
            count=1,
        )

    changed = content != original

    if changed:
        with open(html_path, "w", encoding="utf-8") as f:
            f.write(content)

    return changed


def main() -> int:
    parser = argparse.ArgumentParser(description="Update site stats in index.html")
    parser.add_argument(
        "--repo", default="LeBonhommePharma/FlexAIDdS", help="GitHub repo (owner/name)"
    )
    parser.add_argument("--html", default="site/index.html", help="Path to index.html")
    args = parser.parse_args()

    # Get commit count
    try:
        commit_count = get_commit_count()
        print(f"Commit count: {commit_count}")
    except (subprocess.CalledError, FileNotFoundError) as e:
        print(f"Error getting commit count: {e}", file=sys.stderr)
        return 1

    # Get language breakdown
    languages = fetch_languages(args.repo)
    language_count = count_source_languages(languages) if languages else None
    lang_entries = compute_percentages(languages) if languages else []
    if lang_entries:
        print("Language breakdown:")
        for _, display, pct in lang_entries:
            print(f"  {display}: {pct}%")
    else:
        print("Skipping language update (API unavailable)")

    # Get stars
    stars = fetch_stars(args.repo)
    if stars is not None:
        print(f"Stars: {stars}")

    # Get latest release
    release = fetch_latest_release(args.repo)
    if release:
        print(f"Latest release: {release}")

    # Update HTML
    if not os.path.isfile(args.html):
        print(f"Error: {args.html} not found", file=sys.stderr)
        return 1

    changed = update_html(
        args.html,
        commit_count,
        lang_entries,
        language_count=language_count,
        stars=stars,
        release=release,
    )
    if changed:
        print(f"Updated {args.html}")
    else:
        print("No changes needed")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
