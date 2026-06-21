#!/usr/bin/env python3
"""Update hardcoded repository stats across all site pages simultaneously.

Patches the apex homepage (site/index.html) and the React full site
(site/FlexAIDdS/index.html), and writes a shared JSON snapshot at
site/assets/repo-stats.json consumed by both pages at runtime.

Usage:
    python scripts/update_site_stats.py [--repo OWNER/REPO]
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

# GitHub language name → (CSS variable suffix, short display name, bar color)
LANG_MAP = {
    "C++": ("cpp", "C++", "#f34b7d"),
    "Python": ("python", "Python", "#3572A5"),
    "Swift": ("swift", "Swift", "#F05138"),
    "Objective-C++": ("objcpp", "Obj-C++", "#438eff"),
    "CMake": ("cmake", "CMake", "#8b949e"),
    "TypeScript": ("ts", "TypeScript", "#3178c6"),
    "Cuda": ("cuda", "CUDA", "#76B900"),
    "CUDA": ("cuda", "CUDA", "#76B900"),
    "C": ("c", "C", "#555555"),
    "Shell": ("other", "Shell", "#89e051"),
    "Metal": ("other", "Metal", "#c4c4c4"),
    "JavaScript": ("other", "JavaScript", "#f1e05a"),
    "HTML": ("other", "HTML", "#e34c26"),
}

MIN_PERCENT = 1.0  # Languages below this are grouped into "Other"
OTHER_COLOR = "#555555"

# All HTML files that receive the same semantic stat markers.
HTML_TARGETS = [
    "site/index.html",
    "site/FlexAIDdS/index.html",
]

STATS_JSON_PATH = "site/assets/repo-stats.json"


def fetch_commit_count(repo: str) -> int | None:
    """Return total commit count from GitHub pagination headers."""
    url = f"https://api.github.com/repos/{repo}/commits?per_page=1"
    req = urllib.request.Request(url)
    req.add_header("Accept", "application/vnd.github+json")
    req.add_header("User-Agent", "FlexAIDdS-site-updater")

    token = os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
    if token:
        req.add_header("Authorization", f"Bearer {token}")

    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            link = resp.headers.get("Link", "")
            match = re.search(r'page=(\d+)>;\s*rel="last"', link)
            if match:
                return int(match.group(1))
            data = json.loads(resp.read().decode())
            if isinstance(data, list):
                return len(data)
    except (urllib.error.URLError, urllib.error.HTTPError) as e:
        print(f"Warning: GitHub commit count request failed: {e}", file=sys.stderr)

    return None


def get_commit_count(repo: str) -> int:
    """Return total commit count from GitHub, falling back to local git."""
    remote_count = fetch_commit_count(repo)
    if remote_count is not None and remote_count > 0:
        return remote_count

    try:
        result = subprocess.run(
            ["git", "rev-list", "--count", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        )
        return int(result.stdout.strip())
    except (subprocess.CalledProcessError, FileNotFoundError):
        return 0


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


def fetch_languages(repo: str) -> dict[str, int]:
    """Fetch language byte counts from GitHub API."""
    data = _github_api_get(f"/repos/{repo}/languages")
    return data if isinstance(data, dict) else {}


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


def compute_percentages(
    languages: dict[str, int],
) -> list[tuple[str, str, str, float]]:
    """Build bar segments from GitHub bytes; each lang ≥ MIN_PERCENT gets its own slice."""
    total = sum(languages.values())
    if total == 0:
        return []

    entries: list[tuple[str, str, str, float]] = []
    other_pct = 0.0

    for lang, bytes_count in sorted(languages.items(), key=lambda x: -x[1]):
        pct = round(bytes_count / total * 100, 1)
        if lang in LANG_MAP:
            css_suffix, display, color = LANG_MAP[lang]
            if pct < MIN_PERCENT:
                other_pct += pct
            else:
                entries.append((css_suffix, display, color, pct))
        else:
            other_pct += pct

    if other_pct > 0:
        entries.append(("other", "Other", OTHER_COLOR, round(other_pct, 1)))

    return entries


def count_source_languages(languages: dict[str, int]) -> int:
    """Count languages reported by GitHub (badge matches repo stats)."""
    return sum(1 for lang in languages if lang != "Makefile")


def load_cached_stats() -> dict | None:
    """Load the last known good stats snapshot when GitHub API is unavailable."""
    if not os.path.isfile(STATS_JSON_PATH):
        return None
    try:
        with open(STATS_JSON_PATH, "r", encoding="utf-8") as f:
            data = json.load(f)
        return data if isinstance(data, dict) else None
    except (OSError, json.JSONDecodeError):
        return None


def build_lang_bar(entries: list[tuple[str, str, str, float]]) -> str:
    """Build the lang-bar HTML block for the apex homepage."""
    lines = ['        <div class="lang-bar" aria-label="Language breakdown">']
    for css_suffix, display, _color, pct in entries:
        lines.append(
            f'          <div class="lang-segment" style="width:{pct}%;'
            f'background:var(--lang-{css_suffix})" title="{display} {pct}%"></div>'
        )
    lines.append("        </div>")
    return "\n".join(lines)


def build_lang_legend(entries: list[tuple[str, str, str, float]]) -> str:
    """Build the lang-legend HTML block for the apex homepage."""
    lines = ['        <div class="lang-legend">']
    for css_suffix, display, _color, pct in entries:
        lines.append(
            f'          <span><i style="background:var(--lang-{css_suffix})"></i>'
            f"{display} {pct}%</span>"
        )
    lines.append("        </div>")
    return "\n".join(lines)


def patch_html_markers(
    content: str,
    commit_count: int,
    language_count: int | None,
    *,
    stars: int | None = None,
    release: str | None = None,
) -> str:
    """Patch shared semantic stat markers in an HTML file."""
    content = re.sub(
        r'(<[^>]*id="stat-commits"[^>]*>)\d+(</[^>]+>)',
        rf"\g<1>{commit_count}\g<2>",
        content,
        count=1,
    )

    if language_count is not None:
        content = re.sub(
            r'(<[^>]*id="stat-langs"[^>]*>)\d+(</[^>]+>)',
            rf"\g<1>{language_count}\g<2>",
            content,
            count=1,
        )

    if stars is not None:
        content = re.sub(
            r'(<span[^>]*id="stat-stars"[^>]*>)\d*(</span>)',
            rf"\g<1>{stars}\g<2>",
            content,
            count=1,
        )

    today = datetime.date.today().isoformat()
    content = re.sub(
        r'(<span[^>]*id="last-updated"[^>]*>)[^<]*(</span>)',
        rf"\g<1>{today}\g<2>",
        content,
        count=1,
    )

    if release:
        content = re.sub(
            r'(<span[^>]*id="latest-release"[^>]*>)[^<]*(</span>)',
            rf"\g<1>{release}\g<2>",
            content,
            count=1,
        )

    return content


def patch_apex_stats(
    content: str,
    commit_count: int,
    language_count: int,
    lang_entries: list[tuple[str, str, str, float]],
) -> str:
    """Patch apex-only visible stats (data-count, lang bar/legend)."""
    content = patch_html_markers(content, commit_count, language_count)

    content = re.sub(
        r'(<span class="stat-value[^"]*" data-count=")\d+(">)',
        rf"\g<1>{commit_count}\g<2>",
        content,
        count=1,
    )

    content = re.sub(
        r'(<span class="stat-value" id="stat-langs-display">)\d+(</span>)',
        rf"\g<1>{language_count}\g<2>",
        content,
        count=1,
    )

    if lang_entries:
        bar = build_lang_bar(lang_entries)
        bar_legend = bar + "\n" + build_lang_legend(lang_entries)
        content = re.sub(
            r'        <div class="lang-bar" aria-label="Language breakdown">.*?        </div>\n        <div class="lang-legend">.*?</div>',
            bar_legend,
            content,
            count=1,
            flags=re.DOTALL,
        )

    return content


def write_stats_json(
    path: str,
    commit_count: int,
    language_count: int,
    lang_entries: list[tuple[str, str, str, float]],
    *,
    stars: int | None = None,
    release: str | None = None,
) -> None:
    """Write the shared JSON snapshot consumed by both pages."""
    payload = {
        "commits": commit_count,
        "languageCount": language_count,
        "lastUpdated": datetime.date.today().isoformat(),
        "languages": [
            {
                "id": css_suffix,
                "name": display,
                "percent": pct,
                "color": color,
            }
            for css_suffix, display, color, pct in lang_entries
        ],
    }
    if stars is not None:
        payload["stars"] = stars
    if release:
        payload["latestRelease"] = release

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
        f.write("\n")


def update_all(
    html_targets: list[str],
    commit_count: int,
    lang_entries: list[tuple[str, str, str, float]],
    *,
    language_count: int | None = None,
    stars: int | None = None,
    release: str | None = None,
) -> list[str]:
    """Patch all HTML targets. Returns list of changed file paths."""
    changed: list[str] = []
    lang_count = language_count if language_count is not None else len(lang_entries)

    for html_path in html_targets:
        if not os.path.isfile(html_path):
            print(f"Warning: {html_path} not found, skipping", file=sys.stderr)
            continue

        with open(html_path, "r", encoding="utf-8") as f:
            original = f.read()

        if html_path.endswith("site/index.html"):
            updated = patch_apex_stats(
                original, commit_count, lang_count, lang_entries
            )
        else:
            updated = patch_html_markers(
                original,
                commit_count,
                lang_count,
                stars=stars,
                release=release,
            )

        if updated != original:
            with open(html_path, "w", encoding="utf-8") as f:
                f.write(updated)
            changed.append(html_path)

    json_original = None
    if os.path.isfile(STATS_JSON_PATH):
        with open(STATS_JSON_PATH, "r", encoding="utf-8") as f:
            json_original = f.read()

    write_stats_json(
        STATS_JSON_PATH,
        commit_count,
        lang_count,
        lang_entries,
        stars=stars,
        release=release,
    )

    with open(STATS_JSON_PATH, "r", encoding="utf-8") as f:
        json_updated = f.read()

    if json_updated != json_original:
        changed.append(STATS_JSON_PATH)

    return changed


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Update site stats across apex + FlexAIDdS pages"
    )
    parser.add_argument(
        "--repo", default="LeBonhommePharma/FlexAIDdS", help="GitHub repo (owner/name)"
    )
    args = parser.parse_args()

    cached = load_cached_stats()
    remote_commits = fetch_commit_count(args.repo)
    local_commits = 0
    try:
        result = subprocess.run(
            ["git", "rev-list", "--count", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        )
        local_commits = int(result.stdout.strip())
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass

    commit_count = remote_commits or (cached or {}).get("commits") or local_commits
    if commit_count <= 0:
        print("Error getting commit count", file=sys.stderr)
        return 1
    print(f"Commit count: {commit_count}")

    languages = fetch_languages(args.repo)
    language_count: int | None = None
    lang_entries: list[tuple[str, str, str, float]] = []

    if languages:
        language_count = count_source_languages(languages)
        lang_entries = compute_percentages(languages)
        print("Language breakdown:")
        for _, display, _color, pct in lang_entries:
            print(f"  {display}: {pct}%")
    elif cached and cached.get("languages"):
        language_count = cached.get("languageCount")
        lang_entries = [
            (item["id"], item["name"], item["color"], item["percent"])
            for item in cached["languages"]
        ]
        print("Using cached language breakdown (API unavailable)")
    else:
        print("Skipping language update (API unavailable, no cache)", file=sys.stderr)
        return 1

    stars = fetch_stars(args.repo)
    if stars is None and cached:
        stars = cached.get("stars")

    release = fetch_latest_release(args.repo)
    if not release and cached:
        release = cached.get("latestRelease")

    if stars is not None:
        print(f"Stars: {stars}")
    if release:
        print(f"Latest release: {release}")

    if not lang_entries:
        print("No language data to publish", file=sys.stderr)
        return 1

    changed = update_all(
        HTML_TARGETS,
        commit_count,
        lang_entries,
        language_count=language_count,
        stars=stars,
        release=release,
    )

    if changed:
        for path in changed:
            print(f"Updated {path}")
    else:
        print("No changes needed")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())