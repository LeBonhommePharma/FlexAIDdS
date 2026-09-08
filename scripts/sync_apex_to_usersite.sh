#!/usr/bin/env bash
# Surgically patch repo-stat markers on lebonhommepharma.github.io.
#
# The apex repo is a full brand site (rive, drugs, entropy, …). Never rsync
# --delete site/ onto it: that would wipe unrelated routes. A FlexAIDdS
# GITHUB_TOKEN also cannot push there (403), so the durable publisher is the
# apex workflow `.github/workflows/update-flexaidds-stats.yml`, which pulls
# GitHub API stats itself.
#
# This script is an optional same-day push when USER_SITE_TOKEN (a PAT / fine-
# grained token with contents:write on the user-site repo) is set.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
USER_REPO="${USER_SITE_REPO:-LeBonhommePharma/lebonhommepharma.github.io}"
WORKDIR="/tmp/usersite-sync"
TOKEN="${USER_SITE_TOKEN:-}"
JSON="$ROOT/site/assets/repo-stats.json"

if [ -z "$TOKEN" ]; then
  echo "usersite sync: skipped (no USER_SITE_TOKEN)."
  echo "Apex stats are pulled by ${USER_REPO} workflow update-flexaidds-stats.yml."
  exit 0
fi

if [ ! -f "$JSON" ]; then
  echo "usersite sync: missing $JSON (run update_site_stats.py first)" >&2
  exit 1
fi

rm -rf "$WORKDIR"
git clone --depth 1 "https://x-access-token:${TOKEN}@github.com/${USER_REPO}.git" "$WORKDIR"

python3 "$ROOT/scripts/update_site_stats.py" --from-json "$JSON" --patch-tree "$WORKDIR"

cd "$WORKDIR"
git add -- FlexAIDdS/index.html flexaid-ds/index.html assets/repo-stats.json

if git diff --staged --quiet; then
  echo "user-site sync: no changes"
else
  git -c user.name="github-actions[bot]" \
      -c user.email="github-actions[bot]@users.noreply.github.com" \
      commit -m "chore: refresh FlexAIDdS repo stats from engine snapshot"
  git push origin HEAD:main
  echo "user-site sync: pushed stats-only patch"
fi

rm -rf "$WORKDIR"
